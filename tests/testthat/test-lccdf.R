# tests/testthat/test-lccdf.R
#
# Ports the R-side and Stan-side validation of skellam1_lccdf from the
# source investigation (05-04-skellam-truncation-investigation.qmd in
# tnc001-belize-em): agreement against skellam::pskellam() and a
# brute-force log_sum_exp tail-sum, both on the R side and, once exposed
# via rstan::expose_stan_functions(), on the Stan side too.

# R-side log-PMF (same computation as log_lik_skellam1 / r_lpmf in
# test-lpmf.R; duplicated locally so this file has no cross-file
# dependency on test execution order).
lccdf_r_lpmf <- function(k, mu) log(besselI(2 * mu, abs(k), expon.scaled = TRUE))

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# log S(y, mu) = log P(delta >= -y), the truncated normalising constant.
# lccdf(k, mu) = log P(delta > k), so logS(y, mu) = lccdf(-y - 1, mu).
logS_bruteforce <- function(y, mu, K = NULL) {
  if (is.null(K)) K <- ceiling(mu + 40 * sqrt(mu) + 50)
  log_sum_exp(lccdf_r_lpmf((-y):K, mu))
}

logS_pkg <- function(y, mu) {
  log1p(-skellam::pskellam(-y - 1, lambda1 = mu, lambda2 = mu))
}

mu_grid <- c(0.2, 0.5, 1, 2, 5, 10, 25, 50, 100)
y_grid  <- c(0, 1, 2, 3, 5, 10, 20, 30)
grid    <- expand.grid(mu = mu_grid, y = y_grid)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R brute-force log-CCDF agrees with skellam::pskellam across grid", {
  diffs <- mapply(
    function(mu, y) logS_pkg(y, mu) - logS_bruteforce(y, mu),
    grid$mu, grid$y
  )
  expect_lt(max(abs(diffs)), 1e-8)
})

test_that("leakage at y = 0 matches the documented structural finding", {
  # At y = 0, leakage = 1 - S(0, mu) should be at least ~15% for any mu --
  # the key structural finding motivating truncation in the first place.
  leakage <- vapply(mu_grid, function(mu) 1 - exp(logS_bruteforce(0, mu)), numeric(1))
  expect_true(all(leakage > 0.14),
              label = paste0("min leakage = ", round(min(leakage), 3)))
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

lccdf_stan_block <- function(threshold = 100) {
  paste0("functions {\n", skellam1_lccdf_stan(threshold), "}\nmodel {}\n")
}

stan_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lccdf_stan_block(100))
      rstan::expose_stan_functions(sm)
    })
    stan_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan skellam1_lccdf (exact branch) matches R brute-force reference", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  # Grid is entirely mu <= 100, so every point uses the exact (Bessel-sum)
  # branch under the default normal_approx_threshold.
  diffs <- mapply(
    function(mu, y) skellam1_lccdf(as.integer(-y - 1), mu) - logS_bruteforce(y, mu),
    grid$mu, grid$y
  )
  expect_lt(max(abs(diffs)), 1e-6)
})

test_that("Stan skellam1_lccdf matches package pskellam reference", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  diffs <- mapply(
    function(mu, y) skellam1_lccdf(as.integer(-y - 1), mu) - logS_pkg(y, mu),
    grid$mu, grid$y
  )
  expect_lt(max(abs(diffs)), 1e-6)
})

test_that("Stan skellam1_lccdf is numerically stable at large mu (normal-approx branch)", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  for (mu in c(500, 5000)) {
    for (y in c(-5L, 0L, 5L)) {
      val <- skellam1_lccdf(y, mu)
      expect_false(is.nan(val),      label = paste0("NaN at mu=", mu, " y=", y))
      expect_false(is.infinite(val), label = paste0("Inf at mu=", mu, " y=", y))
    }
  }
})

test_that("normal_approx_threshold parameter actually moves the cutover point", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  low_thresh_env <- new.env()
  ok <- tryCatch({
    suppressMessages({
      sm_low <- rstan::stan_model(model_code = lccdf_stan_block(5))
      rstan::expose_stan_functions(sm_low, env = low_thresh_env)
    })
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "Stan compilation of custom-threshold variant failed")

  mu <- 10  # > 5 (custom threshold) but <= 100 (default threshold)
  y  <- 3L

  exact_val  <- skellam1_lccdf(y, mu)                  # default threshold: exact branch
  approx_val <- low_thresh_env$skellam1_lccdf(y, mu)   # threshold = 5: normal-approx branch

  expect_false(is.nan(approx_val))
  expect_false(is.infinite(approx_val))
  # Different branches at the same (y, mu) should give visibly different
  # results -- confirms the parameter actually moves the cutover, not
  # just that the generated code still parses.
  expect_gt(abs(exact_val - approx_val), 1e-4)
})

# -----------------------------------------------------------------------
# End-to-end smoke test: the resp_trunc() usage example from the docs
# -----------------------------------------------------------------------

test_that("skellam1_stanvars() + skellam1_lccdf_stanvars() works with resp_trunc()", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  set.seed(1)
  n       <- 40
  mu_true <- 5
  y_lb    <- sample(0:8, n, replace = TRUE)

  # Rejection-sample from the truncated Skellam(mu_true, mu_true) subject
  # to delta >= -y_lb, to build a dataset whose generative process matches
  # what resp_trunc() should correct for.
  draw_trunc <- function(y) {
    repeat {
      d <- skellam::rskellam(1, lambda1 = mu_true, lambda2 = mu_true)
      if (d >= -y) return(d)
    }
  }
  delta <- vapply(y_lb, draw_trunc, numeric(1))
  dat   <- data.frame(delta = delta, neg_bound = -y_lb)

  suppressMessages({
    fit <- brms::brm(
      brms::bf(delta | trunc(lb = neg_bound) ~ 1),
      family   = skellam1(),
      stanvars = skellam1_stanvars() + skellam1_lccdf_stanvars(),
      data     = dat,
      backend  = "rstan",
      chains   = 2,
      iter     = 800,
      warmup   = 400,
      seed     = 1,
      refresh  = 0
    )
  })

  expect_s3_class(fit, "brmsfit")
  expect_true(is.finite(brms::fixef(fit)["Intercept", "Estimate"]))
})
