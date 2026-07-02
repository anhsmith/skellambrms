# tests/testthat/test-skellam1.R
#
# Unit tests for posterior_predict_skellam1()/posterior_epred_skellam1()
# under truncation -- the R-side sampling/expectation layer, as opposed to
# the Stan-side lpmf/lccdf already covered by test-lpmf.R/test-lccdf.R (this
# file intentionally does not duplicate those). See R/truncation.R for the
# shared inverse-CDF search and truncated-mean-by-summation machinery, and
# family.R's posterior_predict_skellam1/posterior_epred_skellam1 for how
# skellam1_lccdf_r/skellam1_lpmf_r are wired in.

test_that("posterior_predict_skellam1 respects lb (repro of confirmed bug)", {
  # sigma=3 -> mu_skellam=4.5, well below the normal-approx threshold, so
  # this exercises the exact Bessel-tail-sum branch of skellam1_lccdf_r.
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(3, nrow = 1500, ncol = 1)),
    Y = 0, lb = -2
  )
  set.seed(123)
  draws <- posterior_predict_skellam1(1, prep)
  expect_true(all(draws >= -2), label = paste0("min draw = ", min(draws)))
})

test_that("posterior_predict_skellam1 without lb/ub matches untruncated behaviour (fast path)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(3, nrow = 3000, ncol = 1)),
    Y = 0
  )
  set.seed(1)
  draws <- posterior_predict_skellam1(1, prep)
  expect_true(is.numeric(draws) && length(draws) == 3000)
  # Symmetric about 0 with mu_skellam=4.5 (SD = 3): untruncated spread
  # should extend well below -2 to confirm the fast path was taken.
  expect_true(min(draws) < -3)
})

test_that("posterior_predict_skellam1 draws match the truncated PMF distributionally", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(3, nrow = 4000, ncol = 1)),
    Y = 0, lb = -2, ub = 10
  )
  set.seed(42)
  draws <- posterior_predict_skellam1(1, prep)
  expect_true(all(draws >= -2 & draws <= 10))

  support <- -2:10
  probs <- exp(skellam1_lpmf_r(support, 3))
  probs <- probs / sum(probs)
  emp <- as.numeric(table(factor(draws, levels = support))) / length(draws)
  expect_lt(max(abs(emp - probs)), 0.03)
})

test_that("posterior_epred_skellam1 differs from untruncated mean (0) when lb is tight", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(3, nrow = 5, ncol = 1)),
    Y = 0, lb = -2
  )
  epred <- posterior_epred_skellam1(prep)
  expect_equal(dim(epred), c(5, 1))
  expect_true(all(epred[, 1] > 0), label = paste0("epred = ", epred[1, 1]))
})

test_that("posterior_epred_skellam1 matches brute-force truncated mean", {
  sigma_val <- 3; lb_val <- -2
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(sigma_val, nrow = 1, ncol = 1)),
    Y = 0, lb = lb_val
  )
  epred <- posterior_epred_skellam1(prep)

  support <- lb_val:60
  probs <- exp(skellam1_lpmf_r(support, sigma_val))
  probs <- probs / sum(probs)
  brute_force_mean <- sum(support * probs)

  expect_equal(epred[1, 1], brute_force_mean, tolerance = 1e-6)
})

test_that("posterior_epred_skellam1 leaves untruncated observations exactly at 0 (no regression)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(c(3, 5), nrow = 2, ncol = 2)),
    Y = c(0, 0)
  )
  epred <- posterior_epred_skellam1(prep)
  expect_equal(epred, 0 * prep$dpars$mu)
})
