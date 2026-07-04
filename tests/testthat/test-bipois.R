# tests/testthat/test-bipois.R

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("R brute-force lpmf normalises to 1 across parameter sets", {
  norm_check <- function(mu, lambdaem, lambdalb, K = 60) {
    ys <- 0:K
    yg <- expand.grid(y_em = ys, y_lb = ys)
    lp <- bipois_lpmf_r(yg$y_em, yg$y_lb, mu, lambdaem, lambdalb)
    sum(exp(lp))
  }
  for (p in list(c(1, 1, 1), c(5, 2, 3), c(0.5, 10, 10))) {
    expect_equal(norm_check(p[1], p[2], p[3]), 1, tolerance = 1e-8,
                 label = paste(p, collapse = ","))
  }
})

test_that("R brute-force lpmf reproduces Poisson(mu+lambdaem) marginal for y_em", {
  mu <- 3; lambdaem <- 2; lambdalb <- 4
  K <- 80
  marg_em <- vapply(0:K, function(r) {
    s_vals <- 0:K
    lp <- bipois_lpmf_r(rep(r, length(s_vals)), s_vals, mu, lambdaem, lambdalb)
    mx <- max(lp)
    mx + log(sum(exp(lp - mx)))
  }, numeric(1))
  ref_em <- dpois(0:K, mu + lambdaem, log = TRUE)
  expect_equal(marg_em, ref_em, tolerance = 1e-6)
})

test_that("N_shared cancels: d = y_em - y_lb marginal matches Skellam(lambdaem, lambdalb)", {
  mu <- 4; lambdaem <- 2; lambdalb <- 3
  K <- 60
  ys <- 0:K
  yg <- expand.grid(y_em = ys, y_lb = ys)
  lp <- bipois_lpmf_r(yg$y_em, yg$y_lb, mu, lambdaem, lambdalb)
  d  <- yg$y_em - yg$y_lb
  d_vals <- -10:10
  p_d <- vapply(d_vals, function(dd) sum(exp(lp[d == dd])), numeric(1))
  ref <- skellam::dskellam(d_vals, lambda1 = lambdaem, lambda2 = lambdalb)
  expect_equal(p_d, ref, tolerance = 1e-4)
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

stan_code <- paste0("functions {\n", bipois_stan_funs, "}\nmodel {}\n")

stan_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = stan_code)
      rstan::expose_stan_functions(sm)
    })
    stan_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan bipois_lpmf matches R brute-force reference across a grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu       = c(1e-6, 0.1, 1, 5, 20, 100),
    lambdaem = c(1e-6, 0.1, 1, 5, 20, 100),
    lambdalb = c(1e-6, 0.1, 1, 5, 20, 100)
  )

  check_one <- function(mu, lambdaem, lambdalb) {
    means <- c(mu + lambdaem, mu + lambdalb)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 5L), 0L))
    yg <- expand.grid(y_em = ys, y_lb = ys)
    stan_vals <- mapply(function(r, s) bipois_lpmf(r, mu, lambdaem, lambdalb, s),
                         yg$y_em, yg$y_lb)
    r_vals <- bipois_lpmf_r(yg$y_em, yg$y_lb, mu, lambdaem, lambdalb)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lambdaem, grid$lambdalb)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

test_that("Stan bipois_lpmf is numerically stable at near-zero and large rates", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  edge_cases <- list(
    c(mu = 1e-8, lambdaem = 1e-8, lambdalb = 1e-8),
    c(mu = 1e-8, lambdaem = 50,   lambdalb = 50),
    c(mu = 200,  lambdaem = 1e-8, lambdalb = 1e-8),
    c(mu = 200,  lambdaem = 200,  lambdalb = 200)
  )

  for (ec in edge_cases) {
    ys <- 0:5
    yg <- expand.grid(y_em = ys, y_lb = ys)
    stan_vals <- mapply(function(r, s) bipois_lpmf(r, ec[["mu"]], ec[["lambdaem"]], ec[["lambdalb"]], s),
                         yg$y_em, yg$y_lb)
    r_vals <- bipois_lpmf_r(yg$y_em, yg$y_lb, ec[["mu"]], ec[["lambdaem"]], ec[["lambdalb"]])
    expect_false(any(!is.finite(stan_vals)), label = paste(ec, collapse = ","))
    expect_equal(stan_vals, r_vals, tolerance = 1e-8, label = paste(ec, collapse = ","))
  }
})

# -----------------------------------------------------------------------
# Parameter recovery from simulated hierarchical data (brms end-to-end)
# -----------------------------------------------------------------------

test_that("bipois parameter recovery from simulated vessel-level data", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(20260704)

  n_vessel     <- 8L
  n_per_vessel <- 25L
  n            <- n_vessel * n_per_vessel

  true_log_mu_int       <- log(3)
  true_log_lambdaem_int <- log(1.5)
  true_log_lambdalb_int <- log(4)
  true_sd_vessel        <- 0.4

  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)

  z_mu       <- rnorm(n_vessel)
  z_lambdaem <- rnorm(n_vessel)
  z_lambdalb <- rnorm(n_vessel)

  mu_i       <- exp(true_log_mu_int       + true_sd_vessel * z_mu[vessel])
  lambdaem_i <- exp(true_log_lambdaem_int + true_sd_vessel * z_lambdaem[vessel])
  lambdalb_i <- exp(true_log_lambdalb_int + true_sd_vessel * z_lambdalb[vessel])

  n_shared <- rpois(n, mu_i)
  n10      <- rpois(n, lambdaem_i)
  n01      <- rpois(n, lambdalb_i)

  dat <- data.frame(
    y_em   = n_shared + n10,
    y_lb   = n_shared + n01,
    vessel = factor(vessel)
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y_em | vint(y_lb) ~ 1,
        mu       ~ 1 + (1 | vessel),
        lambdaem ~ 1 + (1 | vessel),
        lambdalb ~ 1 + (1 | vessel)
      ),
      family   = bipois(),
      stanvars = bipois_stanvars(),
      data     = dat,
      backend  = "rstan",
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      seed     = 20260704,
      refresh  = 0,
      control  = list(adapt_delta = 0.95)
    )
  })

  draws <- as.data.frame(fit)

  check_recovery <- function(true_val, draws_col) {
    q <- quantile(draws[[draws_col]], c(0.05, 0.95))
    true_val >= q[[1]] && true_val <= q[[2]]
  }

  # brms treats "mu" as the family's canonical/default dpar and drops its
  # infix from generated column names (b_Intercept, sd_vessel__Intercept),
  # unlike the other two, plainly-named dpars (b_lambdaem_Intercept, etc.)
  # -- same convention already relied on in test-recovery.R's skellam1
  # check, where "mu" likewise stands in for a non-mean quantity.
  expect_true(check_recovery(true_log_mu_int,       "b_Intercept"))
  expect_true(check_recovery(true_log_lambdaem_int, "b_lambdaem_Intercept"))
  expect_true(check_recovery(true_log_lambdalb_int, "b_lambdalb_Intercept"))
  expect_true(check_recovery(true_sd_vessel, "sd_vessel__Intercept"))
  expect_true(check_recovery(true_sd_vessel, "sd_vessel__lambdaem_Intercept"))
  expect_true(check_recovery(true_sd_vessel, "sd_vessel__lambdalb_Intercept"))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.01, label = paste0("max Rhat = ", round(max_rhat, 4)))
})
