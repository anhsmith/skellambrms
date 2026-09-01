# ==========================================================================
# Recovery checks: smoke gates, not calibration statements
# ==========================================================================
#
# WHY THIS EXISTS. The recovery tests used to assert that the true value falls
# inside a single fit's 90% credible interval. For a CORRECT model with a
# well-calibrated posterior, that assertion is a Bernoulli(0.9) draw, so it
# fails 10% of the time by construction -- and across the ~18 such assertions
# in the suite this package was split from, roughly 1.8 spurious failures per
# run were expected. Two were observed on every run, and went unnoticed from the
# day they were written because skip_on_cran() keeps the fitting tests out of
# R CMD check and nothing in CI ever set NOT_CRAN.
#
# So a single-fit interval check is a SMOKE gate and is written as one: a wide
# interval, answering "is this fit grossly wrong?" rather than making a
# calibration claim. Widen the interval rather than tightening it if a correct
# model starts failing here.
#
# The matching CALIBRATION instrument -- refit on fresh simulated data R times
# and measure how often the truth lands inside the nominal interval, with a pass
# threshold derived from the Binomial(R, 0.9) null -- was built for the joint
# bivariate-count families and went with them to
# https://github.com/anhsmith/bicountbrms. It is not reproduced here because
# nothing in this suite currently calls it. If a difference family's
# parameterisation changes, port it back rather than tightening the gate below:
# a tighter single-fit interval measures luck, not calibration.

recovery_ok <- function(draws, true_val, draws_col, level = 0.99) {
  stopifnot(draws_col %in% names(draws))
  a <- (1 - level) / 2
  q <- stats::quantile(draws[[draws_col]], c(a, 1 - a), names = FALSE)
  true_val >= q[[1]] && true_val <= q[[2]]
}
