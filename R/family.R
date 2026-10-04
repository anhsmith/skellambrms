#' Symmetric Skellam custom family for brms
#'
#' @description
#' Returns a brms custom family for the symmetric Skellam distribution,
#' Skellam(mu_skellam, mu_skellam): the distribution of the difference of
#' two independent Poisson(mu_skellam) random variables. The mean is zero.
#' The single parameter is `sigma` (log link), the SD of the difference. The
#' Stan code derives `mu_skellam = sigma^2 / 2` and evaluates the Skellam
#' PMF, which is written in terms of a modified Bessel function, at that
#' value.
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = skellam1(), stanvars = skellam1_stanvars(), data = ...)
#'
#' @section Parameter named mu:
#' `brms::custom_family()` requires one distributional parameter to be named
#' `"mu"`, whatever that parameter represents. In `skellam1()`, the parameter
#' named `mu` is `sigma`, the SD of the difference, on the log link; the mean
#' of `skellam1()` is zero throughout. In formulas and priors, refer to
#' `sigma` as `mu`: for example, `prior(normal(1, 1.5), class = "Intercept")`
#' is a prior on the intercept of log(sigma). The post-processing functions
#' in this package read the parameter with `brms::get_dpar(prep, "mu")` and
#' assign it at once to a variable named `sigma`.
#'
#' @return
#' `skellam1()` returns a brms `custom_family` object. `skellam1_stanvars()`
#' returns a `stanvars` object holding the Stan code for `skellam1_lpmf`.
#' `log_lik_skellam1()` returns a numeric vector of log-densities, one per
#' posterior draw, for observation `i`. `posterior_predict_skellam1()`
#' returns a vector of simulated differences, one per posterior draw, for
#' observation `i`, drawn within the truncation bounds of that observation
#' if it has any. `posterior_epred_skellam1()` returns a draws x observations
#' matrix of means, taken over the truncated distribution on any row that is
#' bounded.
#' @seealso [skellam1_lccdf_stanvars()] for truncation; [skellam2()] for the
#'   free-mean Skellam; [dlaplace1()] and [dnorm1()] for the same fixed
#'   mean with heavier and lighter tails.
#' @export
skellam1 <- function() {
  brms::custom_family(
    name  = "skellam1",
    dpars = "mu",   # forced by brms; represents sigma here -- see Details
    links = "log",
    lb    = 0,
    type  = "int"
  )
}

#' @rdname skellam1
#' @export
skellam1_stanvars <- function() {
  brms::stanvar(block = "functions", scode = skellam1_stan_funs)
}

#' Log-CCDF of the symmetric Skellam distribution, for truncated fits
#'
#' @description
#' Returns a `brms::stanvar()` defining `skellam1_lccdf`, the log
#' complementary CDF of the symmetric Skellam(mu_skellam, mu_skellam)
#' distribution: `skellam1_lccdf(y, sigma)` = log P(delta > y), with
#' `mu_skellam = sigma^2 / 2` derived inside the function. For truncation
#' with `resp_trunc()`, brms finds the log-CCDF of a custom family by the
#' name `<family>_lccdf`. Adding this stanvar to `skellam1_stanvars()` is
#' therefore all that a truncated fit requires, including a fit whose lower
#' bound varies by row.
#'
#' @details
#' For `mu_skellam` above `normal_approx_threshold`, `skellam1_lccdf` uses a
#' normal approximation with variance `2 * mu_skellam`. At or below the
#' threshold, it sums the PMF exactly, over whichever tail lies away from
#' the mean: for `y >= 0`, the upper tail upward from `y + 1`; for `y < 0`,
#' the lower tail downward from `y`, returning the log of one minus that
#' sum. Each term of the sum evaluates a modified Bessel function, and the
#' log link places no upper bound on `sigma`, so warmup can propose values
#' of `mu_skellam` at which the exact sum is slow to evaluate. The sum stops
#' after 500 terms, or earlier once a term is more than 40 log-units below
#' the running sum; neither limit can be changed.
#'
#' A larger threshold applies the exact sum to more evaluations, and a
#' smaller one applies the normal approximation, which is less accurate in
#' the tails, to more. The threshold is on the `mu_skellam` scale, not the
#' `sigma` scale.
#'
#' @param normal_approx_threshold Numeric scalar. Values of `mu_skellam`
#'   above this threshold use the normal approximation instead of the exact
#'   tail sum. Default `100`.
#'
#' @return A `brms::stanvars` object defining the `skellam1_lccdf` Stan
#'   function, for combining with `skellam1_stanvars()` via `+`.
#'
#' @examples
#' \dontrun{
#' library(brms)
#'
#' brm(
#'   bf(y | trunc(lb = neg_bound) ~ x),
#'   family   = skellam1(),
#'   stanvars = skellam1_stanvars() + skellam1_lccdf_stanvars(),
#'   data     = dat
#' )
#' }
#'
#' @seealso [skellam1()] for the family itself; [skellam2_lccdf_stanvars()] for
#'   the free-mean counterpart.
#' @export
skellam1_lccdf_stanvars <- function(normal_approx_threshold = 100) {
  brms::stanvar(
    block = "functions",
    scode = skellam1_lccdf_stan(normal_approx_threshold)
  )
}

# --------------------------------------------------------------------------
# brms interface functions — found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname skellam1
#' @export
#' @keywords internal
log_lik_skellam1 <- function(i, prep) {
  sigma <- brms::get_dpar(prep, "mu", i = i)  # brms dpar name "mu" is sigma here -- see skellam1() Details
  mu    <- sigma^2 / 2
  y     <- prep$data$Y[i]
  # besselI(..., expon.scaled = TRUE) returns I_nu(x) * exp(-x), so
  # log(besselI(2*mu, |y|, expon.scaled=TRUE)) = log(I_|y|(2mu)) - 2mu
  # which equals the full log-PMF: -2mu + log(I_|y|(2mu))
  log(besselI(2 * mu, abs(y), expon.scaled = TRUE))
}

#' @rdname skellam1
#' @export
#' @keywords internal
posterior_predict_skellam1 <- function(i, prep, ...) {
  sigma <- brms::get_dpar(prep, "mu", i = i)  # brms dpar name "mu" is sigma here -- see skellam1() Details
  mu    <- sigma^2 / 2
  lb <- .get_bound(prep, "lb", i)
  ub <- .get_bound(prep, "ub", i)
  if (!is.finite(lb) && !is.finite(ub)) {
    return(skellam::rskellam(length(mu), lambda1 = mu, lambda2 = mu))
  }
  u <- stats::runif(length(sigma))
  .invert_truncated_cdf(function(y, idx) skellam1_lccdf_r(y, sigma[idx]), u, lb, ub)
}

#' @rdname skellam1
#' @export
#' @keywords internal
posterior_epred_skellam1 <- function(prep) {
  sigma <- brms::get_dpar(prep, "mu")  # brms dpar name "mu" is sigma here -- see skellam1() Details
  out <- 0 * sigma  # E[Skellam(mu, mu)] = 0; preserves draw x obs matrix dimensions
  lb_full <- prep$data$lb; ub_full <- prep$data$ub
  if (is.null(lb_full) && is.null(ub_full)) return(out)
  nobs <- ncol(sigma)
  lb_obs <- .get_bound(prep, "lb", seq_len(nobs))
  ub_obs <- .get_bound(prep, "ub", seq_len(nobs))
  trunc_obs <- which(is.finite(lb_obs) | is.finite(ub_obs))
  if (length(trunc_obs) == 0) return(out)
  for (j in trunc_obs) {
    lpmf <- function(y, idx) skellam1_lpmf_r(y, sigma[idx, j])
    out[, j] <- .truncated_mean_by_sum(lpmf, lb_obs[j], ub_obs[j], center = rep(0, nrow(sigma)))
  }
  out
}

# ==========================================================================
# skellam2: asymmetric Skellam, free mean (Koopman-style parameterisation)
# ==========================================================================

#' Asymmetric Skellam custom family for brms
#'
#' @description
#' Returns a brms custom family for the asymmetric Skellam distribution,
#' Skellam(theta1, theta2): the distribution of the difference of two
#' independent Poisson random variables with rates `theta1` and `theta2`,
#' which may differ. The two parameters are `mu` (identity link), the mean
#' of the difference, and `sigmaexcess` (log link, so non-negative). The SD
#' of the difference, `sigma`, and the rates `theta1` and `theta2` are
#' derived from `mu` and `sigmaexcess` (see Details).
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = skellam2(), stanvars = skellam2_stanvars(), data = ...)
#'
#' @section Parameter named sigmaexcess:
#' `brms::custom_family()` does not allow underscores or dots in the names
#' of distributional parameters, so the second parameter is named
#' `sigmaexcess` rather than `sigma_excess`.
#'
#' @section Constraint between mu and sigma:
#' The Skellam rates satisfy `theta1 + theta2 = sigma^2` and
#' `theta1 - theta2 = mu`, so `theta1 = (sigma^2 + mu) / 2` and
#' `theta2 = (sigma^2 - mu) / 2`. Both rates are non-negative only if
#' `sigma^2 >= |mu|`: the variance of a Skellam variable is at least the
#' absolute value of its mean. `skellam2()` satisfies this constraint by
#' construction, through
#'   sigma^2 = |mu| + sigmaexcess^2,
#' for every `mu` and every `sigmaexcess >= 0`, with equality at
#' `sigmaexcess = 0`. For `mu >= 0`, `theta1 = mu + sigmaexcess^2 / 2` and
#' `theta2 = sigmaexcess^2 / 2`. For `mu < 0`, `theta1 = sigmaexcess^2 / 2`
#' and `theta2 = |mu| + sigmaexcess^2 / 2`. Both rates are strictly positive
#' whenever `sigmaexcess > 0`, which the log link guarantees for any finite
#' linear predictor. At `mu = 0`, `skellam2()` reduces to `skellam1()`, with
#' `sigma = sigmaexcess` and `theta1 = theta2 = sigmaexcess^2 / 2`.
#'
#' The alternative construction `sigma = sqrt(mu^2 + sigmaexcess^2)`
#' guarantees only `sigma >= |mu|`, which does not imply `sigma^2 >= |mu|`
#' when `|mu| < 1`. It can then give a negative rate: `mu = 0.5` and
#' `sigmaexcess = 0` give `sigma = 0.5` and `theta2 = -0.125`.
#'
#' @section Derived quantities:
#' The Stan program that brms generates declares the per-observation vectors
#' of distributional parameters, here `mu` and `sigmaexcess`, as local
#' variables in the `model` block rather than as transformed parameters. A
#' `generated quantities` block therefore cannot refer to them, with or
#' without `loop = TRUE`. [skellam2_dpars()] instead computes `mu`, `sigma`, `sigma^2`, `theta1` and
#' `theta2` in R with `brms::get_dpar()`, which works for any formula.
#'
#' @return
#' `skellam2()` returns a brms `custom_family` object. `skellam2_stanvars()`
#' returns a `stanvars` object holding the Stan code for `skellam2_lpmf`.
#' `log_lik_skellam2()` returns a numeric vector of log-densities, one per
#' posterior draw, for observation `i`. `posterior_predict_skellam2()`
#' returns a vector of simulated differences, one per posterior draw, for
#' observation `i`, drawn within the truncation bounds of that observation
#' if it has any. `posterior_epred_skellam2()` returns a draws x observations
#' matrix of means, taken over the truncated distribution on any row that is
#' bounded.
#' @seealso [skellam2_lccdf_stanvars()] for truncation; [skellam2_dpars()] to
#'   recover `sigma`, `theta1` and `theta2` from a fit; [skellam1()] for
#'   the fixed-mean family this one reduces to at `mu = 0`; [dlaplace2()]
#'   and [dnorm2()] for free-mean families that leave mean and spread
#'   uncoupled.
#' @export
skellam2 <- function() {
  brms::custom_family(
    name  = "skellam2",
    dpars = c("mu", "sigmaexcess"),
    links = c("identity", "log"),
    lb    = c(NA, 0),
    type  = "int"
  )
}

#' @rdname skellam2
#' @export
skellam2_stanvars <- function() {
  brms::stanvar(block = "functions", scode = skellam2_stan_funs)
}

#' Log-CCDF of the asymmetric Skellam distribution, for truncated fits
#'
#' @description
#' Returns a `brms::stanvar()` defining `skellam2_lccdf`, the log
#' complementary CDF of the asymmetric Skellam(theta1, theta2)
#' distribution: `skellam2_lccdf(y, mu, sigmaexcess)` = log P(delta > y).
#' `skellam2_lccdf_stanvars()` is used in the same way as
#' [skellam1_lccdf_stanvars()], and the documentation of
#' `skellam1_lccdf_stanvars()` describes how `resp_trunc()` finds the
#' function, the normal approximation and the exact sum. In
#' `skellam2_lccdf`, the threshold is compared with
#' `mu_skellam = (theta1 + theta2) / 2`, which equals the `mu_skellam` of
#' `skellam1()` when `theta1 = theta2`, and the exact sum changes from the
#' upper to the lower tail at `y = mu` rather than at `y = 0`.
#'
#' @param normal_approx_threshold Numeric scalar, compared with
#'   `mu_skellam` as described above. Default `100`.
#'
#' @return A `brms::stanvars` object defining the `skellam2_lccdf` Stan
#'   function, for combining with `skellam2_stanvars()` via `+`.
#' @seealso [skellam2()] for the family itself; [skellam1_lccdf_stanvars()] for
#'   the normal approximation and the exact sum.
#' @export
skellam2_lccdf_stanvars <- function(normal_approx_threshold = 100) {
  brms::stanvar(
    block = "functions",
    scode = skellam2_lccdf_stan(normal_approx_threshold)
  )
}

#' Derived quantities of the asymmetric Skellam family from a fitted model
#'
#' @description
#' Returns `mu`, `sigma`, `sigma^2`, `theta1` and `theta2`, each as a
#' draws x observations matrix, from a `brmsfit` fitted with `skellam2()`.
#' The quantities are computed in R with `brms::get_dpar()`, because the
#' Stan program generated by brms cannot report them from a
#' `generated quantities` block. The section "Derived quantities" in
#' [skellam2()] gives the reason.
#'
#' @param fit A `brmsfit` fitted with `family = skellam2()`.
#' @param newdata Optional new data, passed to `brms::prepare_predictions()`.
#'
#' @return A named list of draws x observations matrices: `mu`, `sigma`,
#'   `sigmasq`, `theta1`, `theta2`.
#' @seealso [skellam2()] for the family and the algebra these quantities are
#'   derived from.
#' @export
skellam2_dpars <- function(fit, newdata = NULL) {
  prep        <- brms::prepare_predictions(fit, newdata = newdata)
  mu          <- brms::get_dpar(prep, "mu")
  sigmaexcess <- brms::get_dpar(prep, "sigmaexcess")
  sigmasq     <- abs(mu) + sigmaexcess^2
  list(
    mu      = mu,
    sigma   = sqrt(sigmasq),
    sigmasq = sigmasq,
    theta1  = (sigmasq + mu) / 2,
    theta2  = (sigmasq - mu) / 2
  )
}

# --------------------------------------------------------------------------
# brms interface functions — found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname skellam2
#' @export
#' @keywords internal
log_lik_skellam2 <- function(i, prep) {
  mu          <- brms::get_dpar(prep, "mu", i = i)
  sigmaexcess <- brms::get_dpar(prep, "sigmaexcess", i = i)
  sigmasq <- abs(mu) + sigmaexcess^2
  theta1  <- (sigmasq + mu) / 2
  theta2  <- (sigmasq - mu) / 2
  y       <- prep$data$Y[i]
  z       <- 2 * sqrt(theta1 * theta2)
  # Matches skellam::dskellam's own internal formula exactly (besselI's
  # expon.scaled=TRUE bakes in -z, added back here, same trick as skellam1).
  log(besselI(z, abs(y), expon.scaled = TRUE)) + z - theta1 - theta2 + (y / 2) * log(theta1 / theta2)
}

#' @rdname skellam2
#' @export
#' @keywords internal
posterior_predict_skellam2 <- function(i, prep, ...) {
  mu          <- brms::get_dpar(prep, "mu", i = i)
  sigmaexcess <- brms::get_dpar(prep, "sigmaexcess", i = i)
  lb <- .get_bound(prep, "lb", i)
  ub <- .get_bound(prep, "ub", i)
  if (!is.finite(lb) && !is.finite(ub)) {
    sigmasq <- abs(mu) + sigmaexcess^2
    theta1  <- (sigmasq + mu) / 2
    theta2  <- (sigmasq - mu) / 2
    return(skellam::rskellam(length(mu), lambda1 = theta1, lambda2 = theta2))
  }
  u <- stats::runif(length(mu))
  .invert_truncated_cdf(
    function(y, idx) skellam2_lccdf_r(y, mu[idx], sigmaexcess[idx]), u, lb, ub
  )
}

#' @rdname skellam2
#' @export
#' @keywords internal
posterior_epred_skellam2 <- function(prep) {
  mu  <- brms::get_dpar(prep, "mu")  # E[Skellam(theta1, theta2)] = theta1 - theta2 = mu
  out <- mu
  lb_full <- prep$data$lb; ub_full <- prep$data$ub
  if (is.null(lb_full) && is.null(ub_full)) return(out)
  sigmaexcess <- brms::get_dpar(prep, "sigmaexcess")
  nobs <- ncol(mu)
  lb_obs <- .get_bound(prep, "lb", seq_len(nobs))
  ub_obs <- .get_bound(prep, "ub", seq_len(nobs))
  trunc_obs <- which(is.finite(lb_obs) | is.finite(ub_obs))
  if (length(trunc_obs) == 0) return(out)
  for (j in trunc_obs) {
    lpmf <- function(y, idx) skellam2_lpmf_r(y, mu[idx, j], sigmaexcess[idx, j])
    out[, j] <- .truncated_mean_by_sum(lpmf, lb_obs[j], ub_obs[j], center = mu[, j])
  }
  out
}

# ==========================================================================
# dlaplace1: discrete Laplace, location fixed at 0, free scale
# ==========================================================================

#' Discrete-Laplace custom family for brms (location 0, free scale)
#'
#' @description
#' Returns a brms custom family for the discrete Laplace distribution with
#' location fixed at 0, obtained from the continuous Laplace(0, b) by CDF
#' differencing: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`. The single parameter
#' is `sigma` (log link), the SD of the continuous distribution before
#' discretisation (see "Conversion from sigma to b"); the mean is zero.
#' Unlike the PMF and CCDF
#' of `skellam1()` and `skellam2()`, those of `dlaplace1()` are closed-form,
#' built on the Stan function `double_exponential_lcdf` (Stan names the
#' Laplace distribution "double exponential"). They require no Bessel
#' function, no large-argument branch and no iteration cap.
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = dlaplace1(), stanvars = dlaplace1_stanvars(), data = ...)
#'
#' @section Parameter named mu:
#' As in `skellam1()`, `brms::custom_family()` requires a parameter named
#' `"mu"`, and in `dlaplace1()` that parameter is `sigma`, the SD, not a
#' mean. See [skellam1()] for the consequences for formulas and priors.
#'
#' @section Conversion from sigma to b:
#' `double_exponential_lcdf` takes the scale `b` of the continuous Laplace
#' distribution. The variance of Laplace(0, b) is `2 * b^2`, so its SD is
#' `b * sqrt(2)`, and `dlaplace1_lpmf` and `dlaplace1_lccdf` both begin by
#' computing `b = sigma / sqrt(2)`. The conversion treats `sigma` as the SD
#' of the continuous distribution. For `sigma >= 0.5`, the SD of the
#' discretised distribution is larger: by 7.9% at `sigma = 0.5`, 3.4% at
#' `sigma = 1`, 1.0% at `sigma = 2` and 0.2% at `sigma = 5`. At
#' `sigma = 0.25`, it is 2.2% smaller, and below that it falls far below
#' `sigma`, because nearly all of the probability falls on 0: at
#' `sigma = 0.1`, the discretised SD is 0.029. In `skellam1()` and
#' `skellam2()`, by contrast, `sigma` is the exact SD of the integer-valued
#' difference.
#'
#' @section Comparison with extraDistr::ddlaplace:
#' `extraDistr::ddlaplace()` implements a different discrete Laplace
#' distribution, with PMF `P(z) = (1 - p) / (1 + p) * p^|z|`, in which its
#' argument `scale` is the decay probability `p` rather than the continuous
#' scale `b`. The two PMFs differ: at `b = 3` and `p = exp(-1/3)`,
#' `P(0) = 0.1535` under `dlaplace1()` and `0.1651` under
#' `extraDistr::ddlaplace()`. The tests in this package therefore compare
#' `dlaplace1()` with an R implementation of the CDF difference rather than
#' with `extraDistr`.
#'
#' @return
#' `dlaplace1()` returns a brms `custom_family` object.
#' `dlaplace1_stanvars()` returns a `stanvars` object holding the Stan code
#' for `dlaplace1_lpmf`. `log_lik_dlaplace1()` returns a numeric vector of
#' log-densities, one per posterior draw, for observation `i`.
#' `posterior_predict_dlaplace1()` returns a vector of simulated
#' differences, one per posterior draw, for observation `i`, drawn within
#' the truncation bounds of that observation if it has any.
#' `posterior_epred_dlaplace1()` returns a draws x observations matrix of
#' means, taken over the truncated distribution on any row that is bounded.
#' @seealso [dlaplace1_lccdf_stanvars()] for truncation; [dlaplace2()] for the
#'   free-mean version; [skellam1()] and [dnorm1()] for the same fixed
#'   mean under a Skellam and a lighter-tailed alternative.
#' @export
dlaplace1 <- function() {
  brms::custom_family(
    name  = "dlaplace1",
    dpars = "mu",   # forced by brms; represents sigma here -- see Details
    links = "log",
    lb    = 0,
    type  = "int"
  )
}

#' @rdname dlaplace1
#' @export
dlaplace1_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dlaplace1_stan_funs)
}

#' Log-CCDF of the discrete Laplace distribution with location 0, for
#' truncated fits
#'
#' @description
#' Returns a `brms::stanvar()` defining `dlaplace1_lccdf`, the log
#' complementary CDF of the discrete Laplace(0, sigma) distribution:
#' `dlaplace1_lccdf(y, sigma)` = log P(Z > y). `dlaplace1_lccdf_stanvars()`
#' is used in the same way as [skellam1_lccdf_stanvars()], but takes no
#' threshold argument, because the closed form
#' `log1m_exp(double_exponential_lcdf(...))` involves no Bessel function and
#' no iterative sum.
#'
#' @return A `brms::stanvars` object defining the `dlaplace1_lccdf` Stan
#'   function, for combining with `dlaplace1_stanvars()` via `+`.
#' @seealso [dlaplace1()] for the family itself; [skellam1_lccdf_stanvars()] for
#'   how `resp_trunc()` locates these functions by name.
#' @export
dlaplace1_lccdf_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dlaplace1_lccdf_stan)
}

# --------------------------------------------------------------------------
# brms interface functions — found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname dlaplace1
#' @export
#' @keywords internal
log_lik_dlaplace1 <- function(i, prep) {
  sigma <- brms::get_dpar(prep, "mu", i = i)  # brms dpar name "mu" is sigma here -- see Details
  b     <- sigma / sqrt(2)
  z     <- prep$data$Y[i]
  # ifelse()'s output length follows its *test* argument, not the yes/no
  # branches -- since z is a single observation (scalar) while b varies by
  # draw (vector), the test x < 0 would be evaluated at length 1 and the
  # whole ifelse() would silently collapse to length 1 instead of ndraws.
  # Broadcasting x against b before the test fixes this.
  laplace_cdf <- function(x) {
    x <- x + 0 * b
    ifelse(x < 0, 0.5 * exp(x / b), 1 - 0.5 * exp(-x / b))
  }
  log(laplace_cdf(z + 0.5) - laplace_cdf(z - 0.5))
}

#' @rdname dlaplace1
#' @export
#' @keywords internal
posterior_predict_dlaplace1 <- function(i, prep, ...) {
  sigma <- brms::get_dpar(prep, "mu", i = i)  # brms dpar name "mu" is sigma here -- see Details
  b     <- sigma / sqrt(2)
  lb <- .get_bound(prep, "lb", i)
  ub <- .get_bound(prep, "ub", i)
  if (!is.finite(lb) && !is.finite(ub)) {
    # Difference of two iid Exponential(rate = 1/b) draws is Laplace(0, b)
    # exactly (the discrete-Laplace analogue of skellam's "difference of
    # two iid Poissons"); rounding a continuous Laplace(0,b) draw to the
    # nearest integer reproduces this family's CDF-differenced PMF exactly,
    # since P(round(X)=z) = P(z-0.5 <= X < z+0.5) = F(z+0.5) - F(z-0.5).
    n <- length(b)
    return(round(stats::rexp(n, rate = 1 / b) - stats::rexp(n, rate = 1 / b)))
  }
  u <- stats::runif(length(sigma))
  .invert_truncated_cdf(function(y, idx) dlaplace1_lccdf_r(y, sigma[idx]), u, lb, ub)
}

#' @rdname dlaplace1
#' @export
#' @keywords internal
posterior_epred_dlaplace1 <- function(prep) {
  sigma <- brms::get_dpar(prep, "mu")  # brms dpar name "mu" is sigma here -- see Details
  out <- 0 * sigma  # E[discrete Laplace(0, sigma)] = 0; preserves draw x obs matrix dimensions
  lb_full <- prep$data$lb; ub_full <- prep$data$ub
  if (is.null(lb_full) && is.null(ub_full)) return(out)
  nobs <- ncol(sigma)
  lb_obs <- .get_bound(prep, "lb", seq_len(nobs))
  ub_obs <- .get_bound(prep, "ub", seq_len(nobs))
  trunc_obs <- which(is.finite(lb_obs) | is.finite(ub_obs))
  if (length(trunc_obs) == 0) return(out)
  for (j in trunc_obs) {
    lpmf <- function(y, idx) dlaplace1_lpmf_r(y, sigma[idx, j])
    out[, j] <- .truncated_mean_by_sum(lpmf, lb_obs[j], ub_obs[j], center = rep(0, nrow(sigma)))
  }
  out
}

# ==========================================================================
# dlaplace2: discrete Laplace, free location AND free scale, uncoupled
# ==========================================================================

#' Discrete-Laplace custom family for brms (free location and scale)
#'
#' @description
#' Returns a brms custom family for the discrete Laplace distribution with
#' free location (`mu`, identity link) and free scale (`sigma`, log link),
#' obtained by CDF differencing as in `dlaplace1()` but centred at `mu`
#' rather than 0: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`, where `F` is the CDF
#' of the continuous Laplace(mu, b).
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = dlaplace2(), stanvars = dlaplace2_stanvars(), data = ...)
#'
#' @section Parameter named mu:
#' In `dlaplace2()`, the parameter named `mu` is the location, so the
#' requirement of `brms::custom_family()` for a parameter named `"mu"` (see
#' [skellam1()]) is met without reinterpreting it.
#'
#' @section Independence of mu and sigma:
#' `skellam2()` requires `sigma^2 >= |mu|`, the relation between the
#' variance and the mean of a Skellam variable (see [skellam2()]). The
#' discrete Laplace distribution has no such relation, and in `dlaplace2()`
#' the parameters `mu` and `sigma` vary independently. Fitting both
#' `skellam2()` and `dlaplace2()` to the same data therefore compares a
#' model in which the bias and the spread of the difference are coupled
#' with a model in which they are not. A constraint between `mu` and `sigma`
#' in `dlaplace2()` would remove that contrast.
#'
#' @section Conversion from sigma to b:
#' As in `dlaplace1()`, `b = sigma / sqrt(2)`, and `sigma` is the SD of the
#' continuous distribution before discretisation. `double_exponential_lcdf`
#' takes the location as an argument, as `normal_lcdf` does, so `mu` is
#' passed to it directly and the Stan code does not shift `z`.
#'
#' @section Mean of the discretised distribution:
#' The mean of the discretised distribution equals `mu` when `mu` is an
#' integer or a half-integer. For other values of `mu`, the mean differs
#' from `mu` by up to 0.015 at `sigma = 1`, 0.054 at `sigma = 0.5` and 0.14
#' at `sigma = 0.25`; the largest differences occur near `mu = 0.3` and
#' `mu = 0.7` modulo 1. `mu` is therefore the location of the continuous
#' distribution before discretisation, not the mean of the discretised one.
#' `posterior_epred_dlaplace2()` returns the mean.
#'
#' @return
#' `dlaplace2()` returns a brms `custom_family` object.
#' `dlaplace2_stanvars()` returns a `stanvars` object holding the Stan code
#' for `dlaplace2_lpmf`. `log_lik_dlaplace2()` returns a numeric vector of
#' log-densities, one per posterior draw, for observation `i`.
#' `posterior_predict_dlaplace2()` returns a vector of simulated
#' differences, one per posterior draw, for observation `i`, drawn within
#' the truncation bounds of that observation if it has any.
#' `posterior_epred_dlaplace2()` returns a draws x observations matrix of
#' means of the discretised distribution: computed in closed form on rows
#' without truncation bounds, and by summing the truncated PMF on bounded
#' rows.
#' @seealso [dlaplace2_lccdf_stanvars()] for truncation; [dlaplace1()] for the
#'   fixed-mean version; [skellam2()] for the coupled comparison;
#'   [dnorm2()] for the light-tailed alternative.
#' @export
dlaplace2 <- function() {
  brms::custom_family(
    name  = "dlaplace2",
    dpars = c("mu", "sigma"),
    links = c("identity", "log"),
    lb    = c(NA, 0),
    type  = "int"
  )
}

#' @rdname dlaplace2
#' @export
dlaplace2_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dlaplace2_stan_funs)
}

#' Log-CCDF of the discrete Laplace distribution with free location, for
#' truncated fits
#'
#' @description
#' Returns a `brms::stanvar()` defining `dlaplace2_lccdf`, the log
#' complementary CDF of the discrete Laplace(mu, sigma) distribution:
#' `dlaplace2_lccdf(y, mu, sigma)` = log P(Z > y).
#' `dlaplace2_lccdf_stanvars()` is used in the same way as
#' [dlaplace1_lccdf_stanvars()], and likewise takes no threshold argument.
#'
#' @return A `brms::stanvars` object defining the `dlaplace2_lccdf` Stan
#'   function, for combining with `dlaplace2_stanvars()` via `+`.
#' @seealso [dlaplace2()] for the family itself; [dlaplace1_lccdf_stanvars()] for
#'   the fixed-mean version.
#' @export
dlaplace2_lccdf_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dlaplace2_lccdf_stan)
}

# --------------------------------------------------------------------------
# brms interface functions — found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname dlaplace2
#' @export
#' @keywords internal
log_lik_dlaplace2 <- function(i, prep) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  b     <- sigma / sqrt(2)
  z     <- prep$data$Y[i]
  laplace_cdf <- function(x) ifelse(x < 0, 0.5 * exp(x / b), 1 - 0.5 * exp(-x / b))
  log(laplace_cdf(z - mu + 0.5) - laplace_cdf(z - mu - 0.5))
}

#' @rdname dlaplace2
#' @export
#' @keywords internal
posterior_predict_dlaplace2 <- function(i, prep, ...) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  b     <- sigma / sqrt(2)
  lb <- .get_bound(prep, "lb", i)
  ub <- .get_bound(prep, "ub", i)
  if (!is.finite(lb) && !is.finite(ub)) {
    n <- length(mu)
    return(round(mu + stats::rexp(n, rate = 1 / b) - stats::rexp(n, rate = 1 / b)))
  }
  u <- stats::runif(length(mu))
  .invert_truncated_cdf(function(y, idx) dlaplace2_lccdf_r(y, mu[idx], sigma[idx]), u, lb, ub)
}

#' @rdname dlaplace2
#' @export
#' @keywords internal
posterior_epred_dlaplace2 <- function(prep) {
  # E[Z] equals mu only for integer or half-integer mu; see .dlaplace_mean()
  # in truncation.R for the closed form used here.
  mu    <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  out   <- .dlaplace_mean(mu, sigma)
  lb_full <- prep$data$lb; ub_full <- prep$data$ub
  if (is.null(lb_full) && is.null(ub_full)) return(out)
  nobs <- ncol(mu)
  lb_obs <- .get_bound(prep, "lb", seq_len(nobs))
  ub_obs <- .get_bound(prep, "ub", seq_len(nobs))
  trunc_obs <- which(is.finite(lb_obs) | is.finite(ub_obs))
  if (length(trunc_obs) == 0) return(out)
  for (j in trunc_obs) {
    lpmf <- function(y, idx) dlaplace2_lpmf_r(y, mu[idx, j], sigma[idx, j])
    out[, j] <- .truncated_mean_by_sum(lpmf, lb_obs[j], ub_obs[j], center = mu[, j])
  }
  out
}

# ==========================================================================
# dnorm1: discrete normal, location fixed at 0, free scale
# ==========================================================================

#' Discrete-normal custom family for brms (location 0, free scale)
#'
#' @description
#' Returns a brms custom family for the discrete normal distribution with
#' location fixed at 0, obtained from the continuous Normal(0, sigma) by CDF
#' differencing: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`. The single parameter
#' is `sigma` (log link), the SD of the continuous distribution before
#' discretisation; the mean is zero. The PMF and CCDF are closed-form: the
#' Stan code computes the PMF with `normal_lcdf` for `z < 0` and with
#' `erfc()` for `z >= 0`, and the CCDF with `erfc()` (see "Evaluation in the
#' upper tail"). As in `dlaplace1()`, they require no Bessel function and
#' no iteration cap.
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = dnorm1(), stanvars = dnorm1_stanvars(), data = ...)
#'
#' @section Parameter named mu:
#' As in `skellam1()` and `dlaplace1()`, `brms::custom_family()` requires a
#' parameter named `"mu"`, and in `dnorm1()` that parameter is `sigma`, the
#' SD, not a mean. See [skellam1()] for the consequences for formulas and
#' priors.
#'
#' @section Scale parameter:
#' `sigma` is the SD of the continuous normal distribution and enters the
#' normal CDF unchanged. `dlaplace1()`, by contrast, converts `sigma` to the
#' Laplace scale `b` first. The SD of the discretised distribution is close
#' to `sqrt(sigma^2 + 1/12)` for `sigma >= 0.5`: larger than `sigma` by 4.1%
#' at `sigma = 1` and by 1.0% at `sigma = 2`.
#'
#' @section Evaluation in the upper tail:
#' The Stan function `normal_lccdf` returns negative infinity for
#' standardised arguments above about 8.25 (stan-dev/math#1985).
#' `dnorm1_lccdf` therefore computes the upper-tail probability as
#' `0.5 * erfc(x / (sigma * sqrt(2)))`, and `dnorm1_lpmf` uses the same form
#' for `z >= 0`. In R, `log_lik_dnorm1()` uses `pnorm(lower.tail = FALSE)`
#' for `z >= 0`, because `log(pnorm(z + 0.5) - pnorm(z - 0.5))` evaluates
#' to `log(0)` once `z` exceeds about 8.5 SDs. Both the Stan and the R
#' log-PMF are accurate until they underflow below about -708, near 38 SDs
#' from 0.
#'
#' @return
#' `dnorm1()` returns a brms `custom_family` object. `dnorm1_stanvars()`
#' returns a `stanvars` object holding the Stan code for `dnorm1_lpmf`.
#' `log_lik_dnorm1()` returns a numeric vector of log-densities, one per
#' posterior draw, for observation `i`. `posterior_predict_dnorm1()` returns
#' a vector of simulated differences, one per posterior draw, for
#' observation `i`, drawn within the truncation bounds of that observation
#' if it has any. `posterior_epred_dnorm1()` returns a draws x observations
#' matrix of means, taken over the truncated distribution on any row that is
#' bounded.
#' @seealso [dnorm1_lccdf_stanvars()] for truncation; [dnorm2()] for the
#'   free-mean version; [skellam1()] and [dlaplace1()] for the same fixed
#'   mean under a Skellam and a heavier-tailed alternative.
#' @export
dnorm1 <- function() {
  brms::custom_family(
    name  = "dnorm1",
    dpars = "mu",   # forced by brms; represents sigma here -- see Details
    links = "log",
    lb    = 0,
    type  = "int"
  )
}

#' @rdname dnorm1
#' @export
dnorm1_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dnorm1_stan_funs)
}

#' Log-CCDF of the discrete normal distribution with location 0, for
#' truncated fits
#'
#' @description
#' Returns a `brms::stanvar()` defining `dnorm1_lccdf`, the log
#' complementary CDF of the discrete Normal(0, sigma) distribution:
#' `dnorm1_lccdf(y, sigma)` = log P(Z > y). `dnorm1_lccdf_stanvars()` is
#' used in the same way as [dlaplace1_lccdf_stanvars()]. `dnorm1_lccdf`
#' computes the upper tail with `erfc()` rather than with the Stan function
#' `normal_lccdf`, which returns negative infinity for standardised
#' arguments above about 8.25 (see [dnorm1()]). `dnorm1_lccdf_stanvars()`
#' takes no threshold argument.
#'
#' @return A `brms::stanvars` object defining the `dnorm1_lccdf` Stan
#'   function, for combining with `dnorm1_stanvars()` via `+`.
#' @seealso [dnorm1()] for the family itself; [skellam1_lccdf_stanvars()] for how
#'   `resp_trunc()` locates these functions by name.
#' @export
dnorm1_lccdf_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dnorm1_lccdf_stan)
}

# --------------------------------------------------------------------------
# brms interface functions — found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname dnorm1
#' @export
#' @keywords internal
log_lik_dnorm1 <- function(i, prep) {
  sigma <- brms::get_dpar(prep, "mu", i = i)  # brms dpar name "mu" is sigma here -- see Details
  z     <- prep$data$Y[i]
  # Not simply log(pnorm(z+0.5) - pnorm(z-0.5)): both terms round to
  # exactly 1.0 from about 9 SDs into the positive tail, giving log(0) =
  # -Inf -- inside this package's "realistic but stressed" test range,
  # unlike the analogous direct-subtraction form in
  # log_lik_dlaplace1 (the Laplace's heavier tail keeps that one accurate
  # much further out). Same z >= 0 branch as dnorm1_lpmf in
  # stanfunctions.R: difference two survival probabilities (small, hence
  # distinguishable) instead of two CDF probabilities (both near 1) when
  # z is on the far side of the mean.
  if (z >= 0) {
    log(stats::pnorm(z - 0.5, sd = sigma, lower.tail = FALSE) -
        stats::pnorm(z + 0.5, sd = sigma, lower.tail = FALSE))
  } else {
    log(stats::pnorm(z + 0.5, sd = sigma) - stats::pnorm(z - 0.5, sd = sigma))
  }
}

#' @rdname dnorm1
#' @export
#' @keywords internal
posterior_predict_dnorm1 <- function(i, prep, ...) {
  sigma <- brms::get_dpar(prep, "mu", i = i)  # brms dpar name "mu" is sigma here -- see Details
  lb <- .get_bound(prep, "lb", i)
  ub <- .get_bound(prep, "ub", i)
  if (!is.finite(lb) && !is.finite(ub)) {
    # P(round(X)=z) = P(z-0.5 <= X < z+0.5) = F(z+0.5) - F(z-0.5) for
    # X ~ Normal(0, sigma) reproduces this family's CDF-differenced PMF
    # exactly, the same rounding identity used in posterior_predict_dlaplace1.
    return(round(stats::rnorm(length(sigma), mean = 0, sd = sigma)))
  }
  u <- stats::runif(length(sigma))
  .invert_truncated_cdf(function(y, idx) dnorm1_lccdf_r(y, sigma[idx]), u, lb, ub)
}

#' @rdname dnorm1
#' @export
#' @keywords internal
posterior_epred_dnorm1 <- function(prep) {
  sigma <- brms::get_dpar(prep, "mu")  # brms dpar name "mu" is sigma here -- see Details
  out <- 0 * sigma  # E[discrete Normal(0, sigma)] = 0; preserves draw x obs matrix dimensions
  lb_full <- prep$data$lb; ub_full <- prep$data$ub
  if (is.null(lb_full) && is.null(ub_full)) return(out)
  nobs <- ncol(sigma)
  lb_obs <- .get_bound(prep, "lb", seq_len(nobs))
  ub_obs <- .get_bound(prep, "ub", seq_len(nobs))
  trunc_obs <- which(is.finite(lb_obs) | is.finite(ub_obs))
  if (length(trunc_obs) == 0) return(out)
  for (j in trunc_obs) {
    lpmf <- function(y, idx) dnorm1_lpmf_r(y, sigma[idx, j])
    out[, j] <- .truncated_mean_by_sum(lpmf, lb_obs[j], ub_obs[j], center = rep(0, nrow(sigma)))
  }
  out
}

# ==========================================================================
# dnorm2: discrete normal, free location AND free scale, uncoupled
# ==========================================================================

#' Discrete-normal custom family for brms (free location and scale)
#'
#' @description
#' Returns a brms custom family for the discrete normal distribution with
#' free location (`mu`, identity link) and free scale (`sigma`, log link),
#' obtained by CDF differencing as in `dnorm1()` but centred at `mu` rather
#' than 0: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`, where `F` is the CDF of the
#' continuous Normal(mu, sigma).
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = dnorm2(), stanvars = dnorm2_stanvars(), data = ...)
#'
#' @section Parameter named mu:
#' In `dnorm2()`, the parameter named `mu` is the location, so the
#' requirement of `brms::custom_family()` for a parameter named `"mu"` (see
#' [skellam1()]) is met without reinterpreting it.
#'
#' @section Independence of mu and sigma:
#' As in [dlaplace2()], `mu` and `sigma` vary independently, with no
#' constraint between them. Fitting `skellam2()`, `dlaplace2()` and
#' `dnorm2()` to the same data compares a model in which the bias and the
#' spread of the difference are coupled with two models in which they are
#' not.
#'
#' @section Evaluation in the upper tail:
#' As in [dnorm1()], the PMF and CCDF compute the upper tail with `erfc()`,
#' with the PMF branching on whether `z` lies above or below `mu` rather
#' than 0.
#'
#' @section Mean of the discretised distribution:
#' The mean of the discretised distribution equals `mu` when `mu` is an
#' integer or a half-integer. For other values of `mu`, the mean differs
#' from `mu` by less than 1e-9 at `sigma = 1` and by less than 5e-6 at
#' `sigma = 0.75`, but by up to 0.0023 at `sigma = 0.5` and 0.093 at
#' `sigma = 0.25`. `mu` is therefore the location of the continuous
#' distribution before discretisation, not the mean of the discretised one.
#' `posterior_epred_dnorm2()` returns the mean.
#'
#' @return
#' `dnorm2()` returns a brms `custom_family` object. `dnorm2_stanvars()`
#' returns a `stanvars` object holding the Stan code for `dnorm2_lpmf`.
#' `log_lik_dnorm2()` returns a numeric vector of log-densities, one per
#' posterior draw, for observation `i`. `posterior_predict_dnorm2()` returns
#' a vector of simulated differences, one per posterior draw, for
#' observation `i`, drawn within the truncation bounds of that observation
#' if it has any. `posterior_epred_dnorm2()` returns a draws x observations
#' matrix of means of the discretised distribution: computed from a
#' rapidly converging series on rows without truncation bounds, and by
#' summing the truncated PMF on bounded rows.
#' @seealso [dnorm2_lccdf_stanvars()] for truncation; [dnorm1()] for the
#'   fixed-mean version; [skellam2()] for the coupled comparison;
#'   [dlaplace2()] for the heavy-tailed alternative.
#' @export
dnorm2 <- function() {
  brms::custom_family(
    name  = "dnorm2",
    dpars = c("mu", "sigma"),
    links = c("identity", "log"),
    lb    = c(NA, 0),
    type  = "int"
  )
}

#' @rdname dnorm2
#' @export
dnorm2_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dnorm2_stan_funs)
}

#' Log-CCDF of the discrete normal distribution with free location, for
#' truncated fits
#'
#' @description
#' Returns a `brms::stanvar()` defining `dnorm2_lccdf`, the log
#' complementary CDF of the discrete Normal(mu, sigma) distribution:
#' `dnorm2_lccdf(y, mu, sigma)` = log P(Z > y). `dnorm2_lccdf_stanvars()` is
#' used in the same way as [dnorm1_lccdf_stanvars()], and likewise takes no
#' threshold argument.
#'
#' @return A `brms::stanvars` object defining the `dnorm2_lccdf` Stan
#'   function, for combining with `dnorm2_stanvars()` via `+`.
#' @seealso [dnorm2()] for the family itself; [dnorm1_lccdf_stanvars()] for the
#'   fixed-mean version.
#' @export
dnorm2_lccdf_stanvars <- function() {
  brms::stanvar(block = "functions", scode = dnorm2_lccdf_stan)
}

# --------------------------------------------------------------------------
# brms interface functions — found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname dnorm2
#' @export
#' @keywords internal
log_lik_dnorm2 <- function(i, prep) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  z     <- prep$data$Y[i]
  # Same z-vs-mean branch as log_lik_dnorm1, but mu varies by draw here
  # (unlike dnorm1's fixed location 0), so the branch itself must be
  # vectorised over draws via ifelse(), not a scalar if().
  ifelse(
    z >= mu,
    log(stats::pnorm(z - mu - 0.5, sd = sigma, lower.tail = FALSE) -
        stats::pnorm(z - mu + 0.5, sd = sigma, lower.tail = FALSE)),
    log(stats::pnorm(z - mu + 0.5, sd = sigma) - stats::pnorm(z - mu - 0.5, sd = sigma))
  )
}

#' @rdname dnorm2
#' @export
#' @keywords internal
posterior_predict_dnorm2 <- function(i, prep, ...) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  lb <- .get_bound(prep, "lb", i)
  ub <- .get_bound(prep, "ub", i)
  if (!is.finite(lb) && !is.finite(ub)) {
    return(round(stats::rnorm(length(mu), mean = mu, sd = sigma)))
  }
  u <- stats::runif(length(mu))
  .invert_truncated_cdf(function(y, idx) dnorm2_lccdf_r(y, mu[idx], sigma[idx]), u, lb, ub)
}

#' @rdname dnorm2
#' @export
#' @keywords internal
posterior_epred_dnorm2 <- function(prep) {
  # E[Z] equals mu only for integer or half-integer mu; see .dnorm_mean()
  # in truncation.R for the series used here.
  mu    <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  out   <- .dnorm_mean(mu, sigma)
  lb_full <- prep$data$lb; ub_full <- prep$data$ub
  if (is.null(lb_full) && is.null(ub_full)) return(out)
  nobs <- ncol(mu)
  lb_obs <- .get_bound(prep, "lb", seq_len(nobs))
  ub_obs <- .get_bound(prep, "ub", seq_len(nobs))
  trunc_obs <- which(is.finite(lb_obs) | is.finite(ub_obs))
  if (length(trunc_obs) == 0) return(out)
  for (j in trunc_obs) {
    lpmf <- function(y, idx) dnorm2_lpmf_r(y, mu[idx, j], sigma[idx, j])
    out[, j] <- .truncated_mean_by_sum(lpmf, lb_obs[j], ub_obs[j], center = mu[, j])
  }
  out
}
