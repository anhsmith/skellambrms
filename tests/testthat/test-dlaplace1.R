# tests/testthat/test-dlaplace1.R
#
# Validates dlaplace1(), the discrete Laplace family (location fixed at
# 0, free scale), with the same rigor as skellam1/skellam2's test suites.
#
# Reference note: extraDistr::ddlaplace() implements a *different*
# discrete Laplace (its `scale` argument is a decay probability p for the
# exact closed form (1-p)/(1+p)*p^|z|, not a continuous-Laplace b) --
# confirmed numerically to NOT match this family's CDF-differenced PMF,
# so it is not usable as a reference here. Per the spec's documented
# fallback, validation instead uses a hand-derived CDF-difference R
# reference, computed in log-space (not a naive direct CDF subtraction,
# which underflows to -Inf once 0.5*exp(-x/b) is too small to perturb
# 1.0 in double precision -- e.g. sigma=0.5, z=14 -- well before the
# Stan implementation itself, which works in log-space throughout via
# log_diff_exp(), has any such problem).

# Stable log-CDF of the continuous Laplace(0, b): avoids the 1-tiny~1.0
# rounding loss a naive log(1 - 0.5*exp(-x/b)) suffers in the tail.
log_laplace_cdf <- function(x, b) {
  ifelse(x < 0, log(0.5) + x / b, log1p(-0.5 * exp(-x / b)))
}

# Stable log-CCDF (log(1-F(x))) of the continuous Laplace(0, b). For
# x >= 0 this is the *exact* closed form log(0.5) - x/b (the upper-tail
# survival function), not "1 - F(x)" computed via exp()/log1p() --
# computing it that way (log1p(-exp(log_laplace_cdf(x,b)))) re-introduces
# the same 1-tiny~1.0 cancellation log_laplace_cdf's own log1p form was
# built to avoid, this time one level up: once log_laplace_cdf(x,b) is
# within ~1e-16 of 0 (i.e. F(x) is within machine epsilon of 1),
# exp(that) itself rounds to exactly 1.0, and log1p(-1.0) = -Inf.
# Confirmed: this broke at sigma=0.2, y=6 (only ~46/b SDs out) using the
# naive form; the direct closed form below does not, for any x >= 0.
log_laplace_lccdf <- function(x, b) {
  ifelse(x >= 0, log(0.5) - x / b, log1p(-0.5 * exp(x / b)))
}

# log(exp(a) - exp(bb)) for a > bb, computed stably (R-side mirror of
# Stan's log_diff_exp()).
log_diff_exp_r <- function(a, bb) a + log1p(-exp(bb - a))

# R-side log-PMF: same computation as log_lik_dlaplace1, but in
# numerically stable log-space (see note above) rather than the
# direct-CDF-subtraction form used in log_lik_dlaplace1 itself -- that
# simpler form is fine for log_lik_dlaplace1's realistic operating range
# (matches this reference to 1e-6 there), but a "trusted reference for a
# wide validation grid" needs the more careful version.
r_lpmf_dl1 <- function(z, sigma) {
  b <- sigma / sqrt(2)
  log_diff_exp_r(log_laplace_cdf(z + 0.5, b), log_laplace_cdf(z - 0.5, b))
}

sigma_vals <- c(0.2, 0.5, 1, 2, 5, 10, 25, 50, 100)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R log-PMF normalises to 1 across grid", {
  for (sigma in sigma_vals) {
    # Range scaled to sigma -- the Laplace distribution has exponential
    # (not Gaussian) tails, so a fixed +-500 window misses a
    # non-negligible amount of mass once sigma is itself comparable to
    # that window (confirmed: a fixed range left ~1e-3 of mass
    # uncounted at sigma=100).
    kmax <- ceiling(40 * sigma + 50)
    # suppressWarnings: R's ifelse() eagerly evaluates both branches of
    # log_laplace_cdf/log_diff_exp_r for every element, including the
    # not-selected branch -- e.g. the x>=0 branch's exp(-x/b) overflows
    # for very negative x, producing a spurious (and discarded) NaN. The
    # value ifelse() actually selects is unaffected.
    total <- suppressWarnings(sum(exp(r_lpmf_dl1(-kmax:kmax, sigma))))
    expect_equal(total, 1, tolerance = 1e-8, label = paste0("sigma = ", sigma))
  }
})

test_that("R log-PMF matches log_lik_dlaplace1's direct-subtraction form", {
  # Cross-check the package's actual R-side log_lik formula (simple
  # direct CDF subtraction) against the log-space reference above, in
  # the realistic range where the direct form doesn't underflow.
  k_vals <- -10L:10L
  for (sigma in c(0.5, 1, 5, 20)) {
    for (k in k_vals) {
      direct <- {
        b <- sigma / sqrt(2)
        laplace_cdf <- function(x) ifelse(x < 0, 0.5 * exp(x / b), 1 - 0.5 * exp(-x / b))
        log(laplace_cdf(k + 0.5) - laplace_cdf(k - 0.5))
      }
      expect_equal(direct, r_lpmf_dl1(k, sigma), tolerance = 1e-6,
                   label = paste0("sigma = ", sigma, ", k = ", k))
    }
  }
})

test_that("R log-PMF is numerically stable, realistic-but-stressed range", {
  # k scaled to sigma (a few SD out, not a fixed absolute value) -- a
  # fixed k paired with a small sigma is an arbitrary, disconnected
  # extreme (e.g. k=50 at sigma=0.2 is 250 SDs out, where even the
  # log-space form above underflows: both CDF terms round to the same
  # double once the true PMF is far below 1e-308, which is a
  # double-precision floor, not a bug -- not worth chasing, mirrors the
  # skellam2 mu=-500 case already treated as out of scope).
  for (sigma in c(0.2, 1, 10, 100)) {
    for (k in round(sigma * c(-10, 0, 10))) {
      val <- r_lpmf_dl1(k, sigma)
      expect_false(is.nan(val),      label = paste0("NaN at sigma=", sigma, " k=", k))
      expect_false(is.infinite(val), label = paste0("Inf at sigma=", sigma, " k=", k))
    }
  }
})

log_sum_exp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# Brute-force log P(Z > y): sums from y+1 (not -y -- that's the
# skellam-style "S(y) = P(>= -y)" truncation quantity used in
# test-lccdf.R/test-skellam2.R, a different thing from the lccdf itself
# unless offset by -y-1; dlaplace1_lccdf(y, sigma) below is defined as
# log P(Z>y) directly, so its reference must match that directly too).
logS_bruteforce_dl1 <- function(y, sigma, K = NULL) {
  if (is.null(K)) K <- y + ceiling(40 * sigma + 50)
  log_sum_exp(r_lpmf_dl1((y + 1):K, sigma))
}

# Closed-form: log P(Z > y) = log(1 - F(y+0.5)), via the stable lccdf
# (y+0.5 is always >= 0 here, so this always takes the exact closed-form
# branch -- see log_laplace_lccdf above).
logS_closed_dl1 <- function(y, sigma) {
  b <- sigma / sqrt(2)
  log_laplace_lccdf(y + 0.5, b)
}

# y as a multiple of sigma (SD units), not a fixed absolute value -- a
# fixed y crossed against the smallest sigma (0.2) would put some grid
# points at an extreme, disconnected ratio. Capped at 10 SDs: even the
# log-space R reference (r_lpmf_dl1) underflows to -Inf somewhat beyond
# that for small sigma (confirmed at sigma=0.2, ~30 SDs out) -- a true
# double-precision floor on the reference, not a family bug (Stan's own
# dlaplace1_lpmf stays finite further out than this R mimic of it does),
# and not worth chasing per the realistic-but-stressed scoping agreed
# for skellam2.
y_multiples <- c(0, 1, 2, 3, 5, 10)
trunc_grid  <- expand.grid(sigma = sigma_vals, y_mult = y_multiples)
trunc_grid$y <- round(trunc_grid$sigma * trunc_grid$y_mult)

test_that("R brute-force log-CCDF agrees with the closed form across grid", {
  diffs <- mapply(
    function(sigma, y) logS_closed_dl1(y, sigma) - logS_bruteforce_dl1(y, sigma),
    trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs)), 1e-8)
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

lpmf_stan_code <- paste0("functions {\n", dlaplace1_stan_funs, "}\nmodel {}\n")

stan_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lpmf_stan_code)
      rstan::expose_stan_functions(sm)
    })
    stan_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan log-PMF matches R log-space reference across grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  # k scaled to sigma (up to 10 SDs out), not a fixed absolute range --
  # see the y_multiples comment above for why.
  for (sigma in sigma_vals) {
    k_vals    <- unique(round(sigma * (-10:10)))
    stan_vals <- vapply(k_vals, function(k) dlaplace1_lpmf(as.integer(k), sigma), numeric(1))
    r_vals    <- vapply(k_vals, r_lpmf_dl1, numeric(1), sigma = sigma)
    expect_equal(stan_vals, r_vals, tolerance = 1e-6, label = paste0("sigma = ", sigma))
  }
})

lccdf_stan_code <- paste0("functions {\n", dlaplace1_lccdf_stan, "}\nmodel {}\n")

lccdf_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm_lccdf <- rstan::stan_model(model_code = lccdf_stan_code)
      rstan::expose_stan_functions(sm_lccdf)
    })
    lccdf_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan dlaplace1_lccdf matches R closed-form / brute-force references", {
  skip_if_not(lccdf_ready, "rstan unavailable or Stan compilation failed")
  # dlaplace1_lccdf(y, sigma) is defined directly as log P(Z > y) (no
  # skellam-style -y-1 offset needed -- see logS_bruteforce_dl1 above).
  diffs_closed <- mapply(
    function(sigma, y) dlaplace1_lccdf(as.integer(y), sigma) - logS_closed_dl1(y, sigma),
    trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs_closed)), 1e-6)

  diffs_brute <- mapply(
    function(sigma, y) dlaplace1_lccdf(as.integer(y), sigma) - logS_bruteforce_dl1(y, sigma),
    trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs_brute)), 1e-6)
})

test_that("Stan dlaplace1_lpmf and dlaplace1_lccdf are numerically stable, realistic-but-stressed range", {
  skip_if_not(stan_ready && lccdf_ready, "rstan unavailable or Stan compilation failed")
  # mirrors the realistic-but-stressed range agreed for skellam2: sigma
  # up to 100 is comfortably beyond this project's real fitted scale, k
  # scaled to sigma (a few SD out) rather than a fixed absolute value --
  # see the R-side version of this test for why a fixed k is the wrong
  # comparison here.
  for (sigma in c(0.2, 1, 10, 100)) {
    for (k in round(sigma * c(-10, 0, 10))) {
      lpmf_val  <- dlaplace1_lpmf(as.integer(k), sigma)
      lccdf_val <- dlaplace1_lccdf(as.integer(k), sigma)
      expect_false(is.nan(lpmf_val) || is.infinite(lpmf_val),
                   label = paste0("lpmf at sigma=", sigma, " k=", k))
      expect_false(is.nan(lccdf_val) || is.infinite(lccdf_val),
                   label = paste0("lccdf at sigma=", sigma, " k=", k))
    }
  }
})

# -----------------------------------------------------------------------
# End-to-end smoke test: resp_trunc() recovery
# -----------------------------------------------------------------------

test_that("dlaplace1_stanvars() + dlaplace1_lccdf_stanvars() recovers sigma_true under resp_trunc()", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  set.seed(3)
  n          <- 50
  sigma_true <- 6
  b_true     <- sigma_true / sqrt(2)
  y_lb       <- sample(0:8, n, replace = TRUE)

  draw_trunc <- function(y) {
    repeat {
      d <- round(rexp(1, rate = 1 / b_true) - rexp(1, rate = 1 / b_true))
      if (d >= -y) return(d)
    }
  }
  delta <- vapply(y_lb, draw_trunc, numeric(1))
  dat   <- data.frame(delta = delta, neg_bound = -y_lb)

  sane_prior <- brms::prior(normal(1, 1.5), class = "Intercept")

  suppressMessages({
    fit <- brms::brm(
      brms::bf(delta | trunc(lb = neg_bound) ~ 1),
      family   = dlaplace1(),
      stanvars = dlaplace1_stanvars() + dlaplace1_lccdf_stanvars(),
      data     = dat,
      prior    = sane_prior,
      backend  = "rstan",
      chains   = 2,
      iter     = 800,
      warmup   = 400,
      seed     = 3,
      refresh  = 0
    )
  })

  expect_s3_class(fit, "brmsfit")

  draws     <- as.data.frame(fit)
  sigma_hat <- exp(draws[["b_Intercept"]])
  sigma_q   <- quantile(sigma_hat, c(0.025, 0.975))
  expect_true(
    sigma_true >= sigma_q[[1]] && sigma_true <= sigma_q[[2]],
    label = paste0("sigma_true = ", sigma_true,
                   ", 95% CI: [", round(sigma_q[[1]], 3),
                   ", ", round(sigma_q[[2]], 3), "]")
  )
})
