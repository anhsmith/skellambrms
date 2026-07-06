# tests/testthat/test-binegbin.R

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("R brute-force lpmf normalises to 1 across parameter sets", {
  norm_check <- function(mu, lem, llb, ss, sx, K = 120) {
    ys <- 0:K
    yg <- expand.grid(y_em = ys, y_lb = ys)
    lp <- binegbin_lpmf_r(yg$y_em, yg$y_lb, mu, lem, llb, ss, sx)
    sum(exp(lp))
  }
  params <- list(
    c(mu = 5,  lem = 3, llb = 3, ss = 2,   sx = 1),
    c(mu = 10, lem = 2, llb = 2, ss = 1.5, sx = 2),
    c(mu = 2,  lem = 4, llb = 4, ss = 5,   sx = 0.8)
  )
  for (p in params) {
    expect_equal(norm_check(p[["mu"]], p[["lem"]], p[["llb"]], p[["ss"]], p[["sx"]]),
                 1, tolerance = 1e-6, label = paste(p, collapse = ","))
  }
})

test_that("moment identities hold (mean, Var(y_em), Var(d), Cov)", {
  # For y_em = N_shared + N10, y_lb = N_shared + N01, NB2(m, phi) with
  # var = m + m^2/phi: E[y_em] = mu + lem; Var(y_em) = (mu+mu^2/ss)+(lem+lem^2/sx);
  # Var(d) = 2(lem+lem^2/sx) when lem=llb, sx shared; Cov = mu + mu^2/ss.
  mu <- 8; lem <- 3; llb <- 3; ss <- 2; sx <- 1.5
  set.seed(1); n <- 3e6
  Ns  <- rnbinom(n, size = ss, mu = mu)
  N10 <- rnbinom(n, size = sx, mu = lem)
  N01 <- rnbinom(n, size = sx, mu = llb)
  ye <- Ns + N10; yl <- Ns + N01; d <- ye - yl
  expect_equal(mean(ye), mu + lem,                         tolerance = 0.02)
  expect_equal(var(ye),  (mu + mu^2/ss) + (lem + lem^2/sx), tolerance = 0.05)
  expect_equal(var(d),   2 * (lem + lem^2/sx),              tolerance = 0.05)
  expect_equal(cov(ye, yl), mu + mu^2/ss,                   tolerance = 0.1)
})

test_that("binegbin reduces to bipois as shapes, shapex -> Inf (Poisson limit)", {
  # NB2 -> Poisson as phi -> Inf with O(1/phi) residual; a large-but-finite phi
  # leaves a small tail difference, so use a generous phi and a modest tolerance.
  mu <- 4; lem <- 2; llb <- 3
  ys <- 0:40
  yg <- expand.grid(y_em = ys, y_lb = ys)
  lp_nb   <- binegbin_lpmf_r(yg$y_em, yg$y_lb, mu, lem, llb, 1e7, 1e7)
  lp_pois <- bipois_lpmf_r(yg$y_em, yg$y_lb, mu, lem, llb)
  expect_equal(lp_nb, lp_pois, tolerance = 1e-3)
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

stan_code <- paste0("functions {\n", binegbin_stan_funs, "}\nmodel {}\n")

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

test_that("Stan binegbin_lpmf matches R brute-force reference across a grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu     = c(0.2, 1, 5, 20),
    lem    = c(0.2, 1, 5),
    llb    = c(0.2, 1, 5),
    ss     = c(0.5, 2, 50),
    sx     = c(0.5, 2, 50)
  )

  check_one <- function(mu, lem, llb, ss, sx) {
    means <- c(mu + lem, mu + llb)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y_em = ys, y_lb = ys)
    stan_vals <- mapply(function(r, s) binegbin_lpmf(r, mu, lem, llb, ss, sx, s),
                        yg$y_em, yg$y_lb)
    r_vals <- binegbin_lpmf_r(yg$y_em, yg$y_lb, mu, lem, llb, ss, sx)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lem, grid$llb, grid$ss, grid$sx)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

test_that("Stan binegbin_lpmf is numerically stable at extreme rates and shapes", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  edge_cases <- list(
    c(mu = 1e-6, lem = 1e-6, llb = 1e-6, ss = 0.2, sx = 0.2),
    c(mu = 1e-6, lem = 50,   llb = 50,   ss = 100, sx = 100),
    c(mu = 200,  lem = 1e-6, llb = 1e-6, ss = 0.2, sx = 0.2),
    c(mu = 200,  lem = 200,  llb = 200,  ss = 50,  sx = 50)
  )

  for (ec in edge_cases) {
    ys <- 0:5
    yg <- expand.grid(y_em = ys, y_lb = ys)
    stan_vals <- mapply(
      function(r, s) binegbin_lpmf(r, ec[["mu"]], ec[["lem"]], ec[["llb"]], ec[["ss"]], ec[["sx"]], s),
      yg$y_em, yg$y_lb)
    r_vals <- binegbin_lpmf_r(yg$y_em, yg$y_lb, ec[["mu"]], ec[["lem"]], ec[["llb"]], ec[["ss"]], ec[["sx"]])
    expect_false(any(!is.finite(stan_vals)), label = paste(ec, collapse = ","))
    expect_equal(stan_vals, r_vals, tolerance = 1e-8, label = paste(ec, collapse = ","))
  }
})

# -----------------------------------------------------------------------
# Parameter recovery from simulated hierarchical data (brms end-to-end)
# -----------------------------------------------------------------------

test_that("binegbin parameter recovery from simulated vessel-level data", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(20260705)

  n_vessel     <- 8L
  n_per_vessel <- 25L
  n            <- n_vessel * n_per_vessel

  true_log_mu_int <- log(8)
  true_sd_vessel  <- 0.3
  true_lem        <- 3           # shared private rate (lambdaem = lambdalb)
  true_shapes     <- 2
  true_shapex     <- 1.5

  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)
  z_mu   <- rnorm(n_vessel)
  mu_i   <- exp(true_log_mu_int + true_sd_vessel * z_mu[vessel])

  n_shared <- rnbinom(n, size = true_shapes, mu = mu_i)
  n10      <- rnbinom(n, size = true_shapex, mu = true_lem)
  n01      <- rnbinom(n, size = true_shapex, mu = true_lem)

  dat <- data.frame(
    y_em   = n_shared + n10,
    y_lb   = n_shared + n01,
    vessel = factor(vessel)
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y_em | vint(y_lb) ~ 1,
        mu ~ 1 + (1 | vessel),
        brms::nlf(lambdaem ~ lamx),
        brms::nlf(lambdalb ~ lamx),
        lamx ~ 1, shapes ~ 1, shapex ~ 1, nl = TRUE
      ),
      family   = binegbin(),
      stanvars = binegbin_stanvars(),
      data     = dat,
      backend  = "rstan",
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      seed     = 20260705,
      refresh  = 0,
      control  = list(adapt_delta = 0.95)
    )
  })

  draws <- as.data.frame(fit)
  check_recovery <- function(true_val, draws_col) {
    q <- quantile(draws[[draws_col]], c(0.05, 0.95))
    true_val >= q[[1]] && true_val <= q[[2]]
  }

  expect_true(check_recovery(true_log_mu_int,   "b_Intercept"))
  expect_true(check_recovery(log(true_lem),     "b_lamx_Intercept"))
  expect_true(check_recovery(true_sd_vessel,    "sd_vessel__Intercept"))
  # shapes/shapex are on the log link; brms reports them as b_<dpar>_Intercept
  expect_true(check_recovery(log(true_shapes),  "b_shapes_Intercept"))
  expect_true(check_recovery(log(true_shapex),  "b_shapex_Intercept"))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.01, label = paste0("max Rhat = ", round(max_rhat, 4)))
})
