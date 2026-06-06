# tests/testthat/test-recovery.R

test_that("parameter recovery from simulated hierarchical data", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(42)

  n_groups          <- 10L
  n_obs_per_group   <- 20L
  true_mu_intercept <- log(3)
  true_sigma_group  <- 0.5

  group_effects <- rnorm(n_groups, mean = 0, sd = true_sigma_group)
  group         <- rep(seq_len(n_groups), each = n_obs_per_group)
  mu_i          <- exp(true_mu_intercept + group_effects[group])
  y             <- skellam::rskellam(length(mu_i), lambda1 = mu_i, lambda2 = mu_i)
  dat           <- data.frame(y = y, group = factor(group))

  suppressMessages({
    fit <- brms::brm(
      y ~ 1 + (1 | group),
      data     = dat,
      family   = skellam1(),
      stanvars = skellam1_stanvars(),
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      seed     = 42,
      refresh  = 0
    )
  })

  draws <- as.data.frame(fit)

  # 1. True mu_intercept within 90% posterior CI for intercept
  intercept_q <- quantile(draws[["b_Intercept"]], c(0.05, 0.95))
  expect_true(
    true_mu_intercept >= intercept_q[[1]] && true_mu_intercept <= intercept_q[[2]],
    label = paste0("true intercept = ", round(true_mu_intercept, 3),
                   ", 90% CI: [", round(intercept_q[[1]], 3),
                   ", ",          round(intercept_q[[2]], 3), "]")
  )

  # 2. True sigma_group within 90% posterior CI
  sd_q <- quantile(draws[["sd_group__Intercept"]], c(0.05, 0.95))
  expect_true(
    true_sigma_group >= sd_q[[1]] && true_sigma_group <= sd_q[[2]],
    label = paste0("sigma_group = ", true_sigma_group,
                   ", 90% CI: [", round(sd_q[[1]], 3),
                   ", ", round(sd_q[[2]], 3), "]")
  )

  # 3. No divergences
  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0,
               label = paste0(n_div, " divergent transitions"))

  # 4. All Rhat < 1.01
  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(
    max_rhat < 1.01,
    label = paste0("max Rhat = ", round(max_rhat, 4))
  )
})
