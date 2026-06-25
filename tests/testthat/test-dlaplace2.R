# tests/testthat/test-dlaplace2.R
#
# Validates dlaplace2(), the free-location/free-scale discrete Laplace
# family, with the same rigor as dlaplace1's test suite. See
# test-dlaplace1.R's file header for why extraDistr::ddlaplace() is not
# a usable reference, and for the numerically stable log-space reference
# pattern (direct closed-form survival function for the lccdf, not
# "1 - exp(log CDF)") used here too -- both lessons carried straight over
# from getting dlaplace1's tests right.

log_laplace_cdf <- function(x, b) {
  ifelse(x < 0, log(0.5) + x / b, log1p(-0.5 * exp(-x / b)))
}

log_laplace_lccdf <- function(x, b) {
  ifelse(x >= 0, log(0.5) - x / b, log1p(-0.5 * exp(x / b)))
}

log_diff_exp_r <- function(a, bb) a + log1p(-exp(bb - a))

# R-side log-PMF: same computation as log_lik_dlaplace2, in stable
# log-space (see test-dlaplace1.R header).
r_lpmf_dl2 <- function(z, mu, sigma) {
  b <- sigma / sqrt(2)
  log_diff_exp_r(log_laplace_cdf(z - mu + 0.5, b), log_laplace_cdf(z - mu - 0.5, b))
}

mu_vals    <- c(-20, -5, 0, 5, 20)
sigma_vals <- c(0.5, 1, 5, 20, 50)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R log-PMF normalises to 1 across grid, including nonzero mu", {
  for (mu in mu_vals) {
    for (sigma in sigma_vals) {
      half  <- ceiling(40 * sigma + 50)
      total <- suppressWarnings(sum(exp(r_lpmf_dl2((mu - half):(mu + half), mu, sigma))))
      expect_equal(total, 1, tolerance = 1e-8,
                   label = paste0("mu = ", mu, ", sigma = ", sigma))
    }
  }
})

test_that("mu = 0 reduces exactly to dlaplace1's lpmf", {
  r_lpmf_dl1 <- function(z, sigma) {
    b <- sigma / sqrt(2)
    log_diff_exp_r(log_laplace_cdf(z + 0.5, b), log_laplace_cdf(z - 0.5, b))
  }
  for (sigma in c(0.5, 5, 20)) {
    for (z in -5:5) {
      expect_equal(r_lpmf_dl2(z, 0, sigma), r_lpmf_dl1(z, sigma), tolerance = 1e-10,
                   label = paste0("sigma = ", sigma, ", z = ", z))
    }
  }
})

test_that("mu shifts the distribution exactly: r_lpmf_dl2(z, mu, sigma) == r_lpmf_dl2(z - mu, 0, sigma)", {
  for (mu in mu_vals) {
    for (sigma in c(1, 10)) {
      for (z_offset in -5:5) {
        z <- round(mu) + z_offset
        expect_equal(r_lpmf_dl2(z, mu, sigma), r_lpmf_dl2(z - mu, 0, sigma), tolerance = 1e-10,
                     label = paste0("mu = ", mu, ", sigma = ", sigma, ", z = ", z))
      }
    }
  }
})

test_that("mu and sigma are not coupled: small sigma with large |mu| is valid (unlike skellam2)", {
  # No sigma >= |mu| floor for this family -- confirm a tiny sigma
  # combined with a large mu produces a well-defined, normalised PMF.
  mu <- 1000
  sigma <- 0.5
  total <- suppressWarnings(sum(exp(r_lpmf_dl2((mu - 50):(mu + 50), mu, sigma))))
  expect_equal(total, 1, tolerance = 1e-6)
})

test_that("R log-PMF is numerically stable, realistic-but-stressed range", {
  for (mu in c(-50, 0, 50)) {
    for (sigma in c(0.2, 1, 10, 100)) {
      for (k in round(mu + sigma * c(-10, 0, 10))) {
        val <- r_lpmf_dl2(k, mu, sigma)
        expect_false(is.nan(val),      label = paste0("NaN at mu=", mu, " sigma=", sigma, " k=", k))
        expect_false(is.infinite(val), label = paste0("Inf at mu=", mu, " sigma=", sigma, " k=", k))
      }
    }
  }
})

log_sum_exp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# Brute-force log P(Z > y), summing from y+1 (dlaplace2_lccdf is defined
# directly as log P(Z>y), no skellam-style offset -- see test-dlaplace1.R).
logS_bruteforce_dl2 <- function(y, mu, sigma, K = NULL) {
  if (is.null(K)) K <- y + ceiling(40 * sigma + 50)
  log_sum_exp(r_lpmf_dl2((y + 1):K, mu, sigma))
}

# Closed-form log P(Z > y) via the exact upper-tail survival function
# (y - mu + 0.5 is not guaranteed >= 0 here since mu can be large and
# positive -- log_laplace_lccdf handles both signs correctly).
logS_closed_dl2 <- function(y, mu, sigma) {
  b <- sigma / sqrt(2)
  log_laplace_lccdf(y - mu + 0.5, b)
}

# y as mu + a multiple of sigma (SD units), capped at 10 SDs -- see
# test-dlaplace1.R for why (a true double-precision floor on the
# CDF-differencing reference beyond that, not a family bug).
y_multiples <- c(-10, -5, -3, -1, 0, 1, 3, 5, 10)
trunc_grid  <- expand.grid(mu = mu_vals, sigma = sigma_vals, y_mult = y_multiples)
trunc_grid$y <- round(trunc_grid$mu + trunc_grid$sigma * trunc_grid$y_mult)

test_that("R brute-force log-CCDF agrees with the closed form across grid", {
  # suppressWarnings: same benign ifelse()-eager-evaluation noise noted
  # in test-dlaplace1.R's normalisation test.
  diffs <- suppressWarnings(mapply(
    function(mu, sigma, y) logS_closed_dl2(y, mu, sigma) - logS_bruteforce_dl2(y, mu, sigma),
    trunc_grid$mu, trunc_grid$sigma, trunc_grid$y
  ))
  expect_lt(max(abs(diffs)), 1e-8)
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

lpmf_stan_code <- paste0("functions {\n", dlaplace2_stan_funs, "}\nmodel {}\n")

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
  for (mu in mu_vals) {
    for (sigma in sigma_vals) {
      k_vals    <- unique(round(mu + sigma * (-10:10)))
      stan_vals <- vapply(k_vals, function(k) dlaplace2_lpmf(as.integer(k), mu, sigma), numeric(1))
      r_vals    <- vapply(k_vals, r_lpmf_dl2, numeric(1), mu = mu, sigma = sigma)
      expect_equal(stan_vals, r_vals, tolerance = 1e-6,
                   label = paste0("mu = ", mu, ", sigma = ", sigma))
    }
  }
})

lccdf_stan_code <- paste0("functions {\n", dlaplace2_lccdf_stan, "}\nmodel {}\n")

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

test_that("Stan dlaplace2_lccdf matches R closed-form / brute-force references", {
  skip_if_not(lccdf_ready, "rstan unavailable or Stan compilation failed")
  diffs_closed <- mapply(
    function(mu, sigma, y) dlaplace2_lccdf(as.integer(y), mu, sigma) - logS_closed_dl2(y, mu, sigma),
    trunc_grid$mu, trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs_closed)), 1e-6)

  diffs_brute <- suppressWarnings(mapply(
    function(mu, sigma, y) dlaplace2_lccdf(as.integer(y), mu, sigma) - logS_bruteforce_dl2(y, mu, sigma),
    trunc_grid$mu, trunc_grid$sigma, trunc_grid$y
  ))
  expect_lt(max(abs(diffs_brute)), 1e-6)
})

test_that("Stan dlaplace2_lpmf/lccdf are numerically stable, realistic-but-stressed range", {
  skip_if_not(stan_ready && lccdf_ready, "rstan unavailable or Stan compilation failed")
  for (mu in c(-50, 0, 50)) {
    for (sigma in c(0.2, 1, 10, 100)) {
      for (k in round(mu + sigma * c(-10, 0, 10))) {
        lpmf_val  <- dlaplace2_lpmf(as.integer(k), mu, sigma)
        lccdf_val <- dlaplace2_lccdf(as.integer(k), mu, sigma)
        expect_false(is.nan(lpmf_val) || is.infinite(lpmf_val),
                     label = paste0("lpmf at mu=", mu, " sigma=", sigma, " k=", k))
        expect_false(is.nan(lccdf_val) || is.infinite(lccdf_val),
                     label = paste0("lccdf at mu=", mu, " sigma=", sigma, " k=", k))
      }
    }
  }
})

# -----------------------------------------------------------------------
# Structural check: no constraint coupling mu and sigma
# -----------------------------------------------------------------------

test_that("make_stancode() shows no constraint coupling mu and sigma", {
  skip_if_not_installed("brms")

  fake_dat <- data.frame(delta = sample(-5:5, 20, replace = TRUE), x = rnorm(20))
  code <- brms::make_stancode(
    brms::bf(delta ~ x, sigma ~ x),
    family   = dlaplace2(),
    stanvars = dlaplace2_stanvars(),
    data     = fake_dat
  )
  # No reject() anywhere, and sigma's only bound is its own plain lb=0
  # (from the log link), not anything depending on mu.
  expect_false(grepl("reject\\(", code))
})

# -----------------------------------------------------------------------
# End-to-end smoke tests: resp_trunc() recovery, zero and nonzero mu_true
# -----------------------------------------------------------------------

fit_and_check_dlaplace2 <- function(mu_true, sigma_true, seed) {
  set.seed(seed)
  n    <- 60
  b_true <- sigma_true / sqrt(2)
  y_lb <- sample(0:10, n, replace = TRUE)

  draw_trunc <- function(y) {
    repeat {
      d <- round(mu_true + rexp(1, rate = 1 / b_true) - rexp(1, rate = 1 / b_true))
      if (d >= -y) return(d)
    }
  }
  delta <- vapply(y_lb, draw_trunc, numeric(1))
  dat   <- data.frame(delta = delta, neg_bound = -y_lb)

  sane_prior <- brms::prior(normal(0, 10), class = "Intercept") +
    brms::prior(normal(1, 1.5), class = "Intercept", dpar = "sigma")

  fit <- brms::brm(
    brms::bf(delta | trunc(lb = neg_bound) ~ 1, sigma ~ 1),
    family   = dlaplace2(),
    stanvars = dlaplace2_stanvars() + dlaplace2_lccdf_stanvars(),
    data     = dat,
    prior    = sane_prior,
    backend  = "rstan",
    chains   = 2,
    iter     = 1000,
    warmup   = 500,
    seed     = seed,
    refresh  = 0
  )

  draws     <- as.data.frame(fit)
  mu_hat    <- draws[["b_Intercept"]]
  mu_q      <- quantile(mu_hat, c(0.025, 0.975))
  sigma_hat <- exp(draws[["b_sigma_Intercept"]])
  sigma_q   <- quantile(sigma_hat, c(0.025, 0.975))

  list(
    fit  = fit,
    mu_ok = mu_true >= mu_q[[1]] && mu_true <= mu_q[[2]],
    mu_label = paste0("mu_true = ", mu_true, ", 95% CI: [", round(mu_q[[1]], 3), ", ", round(mu_q[[2]], 3), "]"),
    sigma_ok = sigma_true >= sigma_q[[1]] && sigma_true <= sigma_q[[2]],
    sigma_label = paste0("sigma_true = ", sigma_true, ", 95% CI: [", round(sigma_q[[1]], 3), ", ", round(sigma_q[[2]], 3), "]")
  )
}

test_that("dlaplace2 recovers a genuinely nonzero mu_true under resp_trunc()", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  suppressMessages({
    res <- fit_and_check_dlaplace2(mu_true = 4, sigma_true = 5, seed = 4)
  })

  expect_s3_class(res$fit, "brmsfit")
  expect_true(res$mu_ok, label = res$mu_label)
  expect_true(res$sigma_ok, label = res$sigma_label)
})

test_that("dlaplace2 does not spuriously detect bias when mu_true = 0", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  suppressMessages({
    res <- fit_and_check_dlaplace2(mu_true = 0, sigma_true = 5, seed = 5)
  })

  expect_s3_class(res$fit, "brmsfit")
  expect_true(res$mu_ok, label = res$mu_label)
  expect_true(res$sigma_ok, label = res$sigma_label)
})
