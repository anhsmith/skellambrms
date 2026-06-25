#' Symmetric Skellam custom family for brms
#'
#' @description
#' Returns a brms custom family for the symmetric Skellam distribution,
#' Skellam(mu_skellam, mu_skellam) — the distribution of the difference of
#' two independent Poisson(mu_skellam) random variables. The single
#' parameter is sigma (link = "log"), the SD of that difference; the mean
#' is always zero. Internally, mu_skellam = sigma^2 / 2 is derived as a
#' transformed parameter and fed to the underlying Bessel-function PMF,
#' which is otherwise unchanged.
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = skellam1(), stanvars = skellam1_stanvars(), data = ...)
#'
#' @details
#' This family was originally parameterised directly on mu_skellam
#' (link = "log"); it now samples on sigma instead, for a common
#' (mean, SD-scale) convention shared with skellam2(), dlaplace1(), and
#' dlaplace2(). Since sigma = sqrt(2 * mu_skellam), a prior previously
#' stated on log(mu_skellam) — e.g. normal(1, 1.5), used in
#' 05-04-skellam-laplace-truncation-validation.qmd — translates as:
#'   log(sigma) = 0.5 * log(2) + 0.5 * log(mu_skellam)
#' so an intercept of 1 on the old log(mu_skellam) scale corresponds to
#' an intercept of 0.5*log(2) + 0.5*1 ≈ 0.847 on the new log(sigma) scale,
#' and the old prior's SD of 1.5 becomes 0.75 on the new scale (a linear
#' transform of a normal is normal). This is a scale correspondence only —
#' slope-coefficient interpretations from the old parameterisation are
#' NOT carried forward; any offset-vs-free-slope diagnostic should be
#' redone fresh against this sigma-scale parameterisation.
#'
#' @details
#' **Naming note.** `brms::custom_family()` hard-requires one `dpars`
#' entry to be literally named `"mu"` (`stop2("All families must have a
#' 'mu' parameter.")`, unconditional, no override) — every family built
#' on it, including this one, must comply regardless of what that
#' parameter actually represents. For skellam1, the brms/Stan-level dpar
#' named `mu` IS sigma (the SD of the difference, log-linked); it is NOT
#' the distribution's mean, which is structurally zero throughout. This
#' is a forced naming collision with brms's API, not a reversion to the
#' pre-reparameterisation behaviour: internally, `mu_skellam = mu^2 / 2`
#' is still derived from it before reaching the Bessel-function PMF,
#' exactly as documented above for "sigma". All R-side helper functions
#' below immediately rebind this dpar to a variable called `sigma` so
#' that no code in this package, other than the literal `dpars`/
#' `get_dpar()` calls forced by brms, ever refers to it as `mu`.
#'
#' @return A brms custom_family object.
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

#' Truncated-Skellam log-CCDF for use with brms's resp_trunc()
#'
#' @description
#' Returns a `brms::stanvar()` defining `skellam1_lccdf`, the log
#' complementary CDF of the symmetric Skellam(mu_skellam, mu_skellam)
#' distribution — `skellam1_lccdf(y, sigma)` = log P(delta > y), where
#' `mu_skellam = sigma^2 / 2` is derived internally. brms's generic
#' truncation machinery (`resp_trunc()`) locates a custom family's
#' log-CCDF by name convention (`<family>_lccdf`), so adding this stanvar
#' alongside `skellam1_stanvars()` is sufficient to support truncated
#' fits, including a row-varying lower bound — no other wiring is
#' required.
#'
#' @details
#' For `mu_skellam` above `normal_approx_threshold`, the exact log-CCDF —
#' an iterative tail-sum of the Skellam PMF, each term a Bessel function
#' evaluation — is replaced by a normal approximation, using
#' Var(Skellam(mu_skellam, mu_skellam)) = 2 * mu_skellam. This guards
#' against two confirmed failure modes, both triggered by an unadapted
#' HMC proposal pushing `sigma` (and hence `mu_skellam`) to an extreme
#' value during warmup (the log link on `sigma` places no ceiling on it):
#'
#' - A crash (`std::bad_alloc`) from `log_modified_bessel_first_kind`
#'   being evaluated at an enormous Bessel order.
#' - A slow-motion version of the same problem: `mu_skellam` in the
#'   hundreds still triggers the expensive exact loop, and if many rows do
#'   this within a single deep NUTS tree the cost compounds
#'   multiplicatively rather than crashing outright — observed as 200+
#'   CPU-seconds and several GB of memory consumed without completing one
#'   iteration.
#'
#' The exact loop also carries a hard cap of 500 iterations past `y`,
#' with an early exit once the tail term becomes negligible (more than
#' ~40 log-units below the running sum). These guard the same two
#' failure modes as the threshold itself and are not configurable here.
#'
#' The default threshold of 100 is **not a universal constant** — it
#' was calibrated to one project's data, where real per-taxon
#' `mu_skellam` estimates topped out around 30 (this is on the
#' `mu_skellam` scale, unaffected by the sigma-reparameterisation). The 3x
#' margin above that (rather than setting the threshold at, say, 35)
#' exists because HMC warmup transiently proposes values well outside any
#' final posterior estimate, not because 30 itself needed padding. When
#' using this function with a different count scale, consider what
#' *implausible but reachable during warmup* looks like for your
#' `mu_skellam`, not just your expected posterior range, and set the
#' threshold a few-fold above that. Setting it too low pays for the
#' normal approximation's bias more often than necessary; setting it too
#' high re-exposes the crash/slow-blowup risk this exists to prevent.
#'
#' @param normal_approx_threshold Numeric scalar; `mu_skellam` values
#'   above this use the normal approximation instead of the exact
#'   Bessel-sum tail. Default `100`. See Details for how to choose this
#'   for your data.
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
  skellam::rskellam(length(mu), lambda1 = mu, lambda2 = mu)
}

#' @rdname skellam1
#' @export
#' @keywords internal
posterior_epred_skellam1 <- function(prep) {
  sigma <- brms::get_dpar(prep, "mu")  # brms dpar name "mu" is sigma here -- see skellam1() Details
  0 * sigma  # E[Skellam(mu, mu)] = 0; preserves draw x obs matrix dimensions
}

# ==========================================================================
# skellam2: asymmetric Skellam, free mean (Koopman-style parameterisation)
# ==========================================================================

#' Asymmetric Skellam custom family for brms
#'
#' @description
#' Returns a brms custom family for the general (asymmetric) Skellam
#' distribution, Skellam(theta1, theta2) — the distribution of the
#' difference of two independent Poisson(theta1), Poisson(theta2) random
#' variables with possibly unequal rates. Two parameters: `mu` (link =
#' "identity"), the mean of the difference, and `sigmaexcess` (link =
#' "log", so `>= 0`), from which `sigma`, the SD of the difference, and
#' the underlying rates `theta1`, `theta2` are derived as transformed
#' quantities (see Details).
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = skellam2(), stanvars = skellam2_stanvars(), data = ...)
#'
#' @details
#' **Naming note.** `brms::custom_family()` disallows underscores in
#' `dpars` (`stop2("Dots or underscores are not allowed in 'dpars'.")`),
#' so the second parameter is spelled `sigmaexcess`, not `sigma_excess`
#' as in the package's design notes and Stan code comments — the two
#' names refer to the same quantity.
#'
#' **Constraint algebra — corrected from the original design.** The
#' natural-seeming construction `sigma = sqrt(mu^2 + sigmaexcess^2)`
#' (Pythagorean in mu and sigmaexcess) only guarantees `sigma >= |mu|`.
#' That is NOT the condition Skellam validity actually needs. With
#' `theta1 = (sigma^2 + mu) / 2` and `theta2 = (sigma^2 - mu) / 2`
#' (from `theta1 + theta2 = sigma^2` and `theta1 - theta2 = mu`),
#' `theta1, theta2 >= 0` requires `sigma^2 >= |mu|` — i.e. Var >= |mean|,
#' the genuine Skellam constraint (sum of two nonnegative Poisson rates
#' is always >= their difference's absolute value). `sigma >= |mu|` and
#' `sigma^2 >= |mu|` coincide only when `|mu| >= 1`; for `|mu| < 1` they
#' diverge, and `sigma = sqrt(mu^2 + sigmaexcess^2)` can produce a
#' *negative* theta1 or theta2 — confirmed numerically, e.g.
#' `mu = 0.5, sigmaexcess = 0` gives `sigma = 0.5`, `theta2 = -0.125`.
#' This package instead uses:
#'   sigma^2 = |mu| + sigmaexcess^2
#' which guarantees `sigma^2 >= |mu|` directly (the right-hand side is
#' `|mu|` plus a nonnegative term), for every `mu` and every
#' `sigmaexcess >= 0`, with equality (the minimal-spread boundary) at
#' `sigmaexcess = 0`. `theta1` and `theta2` are then both sums of
#' nonnegative terms (verify: for `mu >= 0`, `theta1 = mu +
#' sigmaexcess^2/2 >= 0` and `theta2 = sigmaexcess^2/2 >= 0`; for `mu <
#' 0`, the roles swap) — strictly positive whenever `sigmaexcess > 0`,
#' which the log link guarantees for any finite linear predictor. At
#' `mu = 0` this reduces exactly to skellam1's symmetric family
#' (`sigma = sigmaexcess`, `theta1 = theta2 = sigmaexcess^2 / 2`).
#'
#' **Generated-quantities note.** This family does *not* expose `mu`,
#' `sigma`, `sigma^2`, `theta1`, `theta2` via a Stan `generated
#' quantities` block. Confirmed via `make_stancode()`: brms declares a
#' custom family's per-observation dpar vectors (`mu`, `sigmaexcess`
#' here) as local variables inside the generated model's `model` block,
#' not `transformed parameters` — out of Stan-scope for `generated
#' quantities`, regardless of `loop = TRUE/FALSE`. Reconstructing them
#' from brms's internal linear-predictor variable names (`Xc`, `b`,
#' `Intercept`, ...) would only work for simple fixed-effects-only
#' formulas, breaking silently for anything with random effects or
#' splines. `skellam2_dpars()` (below) reports the same five quantities
#' from R instead, via `brms::get_dpar()` — works for any formula.
#'
#' @return A brms custom_family object.
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

#' Truncated-asymmetric-Skellam log-CCDF for use with brms's resp_trunc()
#'
#' @description
#' Returns a `brms::stanvar()` defining `skellam2_lccdf`, the log
#' complementary CDF of the asymmetric Skellam(theta1, theta2)
#' distribution — `skellam2_lccdf(y, mu, sigmaexcess)` = log P(delta > y).
#' Same role and calling convention as `skellam1_lccdf_stanvars()`; see
#' that function's documentation for how `resp_trunc()` locates it and
#' for the rationale behind the normal-approximation threshold (here
#' checked against `mu_skellam = (theta1 + theta2) / 2`, the direct
#' generalisation of skellam1's threshold quantity to the asymmetric
#' case — see skellam2_lccdf_stan() in stanfunctions.R).
#'
#' @param normal_approx_threshold Numeric scalar; see
#'   `skellam1_lccdf_stanvars()` for how to choose this for your data.
#'   Default `100`.
#'
#' @return A `brms::stanvars` object defining the `skellam2_lccdf` Stan
#'   function, for combining with `skellam2_stanvars()` via `+`.
#' @export
skellam2_lccdf_stanvars <- function(normal_approx_threshold = 100) {
  brms::stanvar(
    block = "functions",
    scode = skellam2_lccdf_stan(normal_approx_threshold)
  )
}

#' Report skellam2's derived quantities from a fitted model
#'
#' @description
#' Returns `mu`, `sigma`, `sigma^2`, `theta1`, and `theta2` (each a
#' draws x observations matrix) from a `skellam2()` `brmsfit`, computed
#' in R via `brms::get_dpar()` rather than a Stan `generated quantities`
#' block — see "Generated-quantities note" in `?skellam2` for why the
#' latter isn't available for this family.
#'
#' @param fit A `brmsfit` fitted with `family = skellam2()`.
#' @param newdata Optional new data, passed to `brms::prepare_predictions()`.
#'
#' @return A named list of draws x observations matrices: `mu`, `sigma`,
#'   `sigmasq`, `theta1`, `theta2`.
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
  sigmasq <- abs(mu) + sigmaexcess^2
  theta1  <- (sigmasq + mu) / 2
  theta2  <- (sigmasq - mu) / 2
  skellam::rskellam(length(mu), lambda1 = theta1, lambda2 = theta2)
}

#' @rdname skellam2
#' @export
#' @keywords internal
posterior_epred_skellam2 <- function(prep) {
  brms::get_dpar(prep, "mu")  # E[Skellam(theta1, theta2)] = theta1 - theta2 = mu
}

# ==========================================================================
# dlaplace1: discrete Laplace, location fixed at 0, free scale
# ==========================================================================

#' Discrete-Laplace custom family for brms (location 0, free scale)
#'
#' @description
#' Returns a brms custom family for the discrete Laplace distribution,
#' location fixed at 0, discretised from the continuous Laplace(0, b) via
#' CDF differencing: `P(Z=z) = F(z+0.5) - F(z-0.5)`. One parameter, sigma
#' (link = "log"), the SD; the mean is always zero. Unlike skellam1/
#' skellam2, the PMF and CCDF are closed-form
#' (`double_exponential_lcdf`-based -- Stan's name for the Laplace
#' distribution -- no Bessel function, no large-argument branch or
#' iteration cap needed).
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = dlaplace1(), stanvars = dlaplace1_stanvars(), data = ...)
#'
#' @details
#' **Naming note.** Same forced naming as `skellam1()`: `brms::custom_family()`
#' requires a dpar literally named `"mu"`; here it represents sigma (the
#' SD), not a mean. See `?skellam1` Details for the full rationale.
#'
#' **sigma-to-b conversion.** Stan's `double_exponential_lcdf` expects
#' the continuous Laplace's own scale parameter, `b`. Var(Laplace(0,b)) = `2*b^2`, so SD
#' = `b*sqrt(2)`; treating sigma as exactly that SD (the discretisation
#' perturbs the true discrete variance only slightly, and this keeps
#' sigma on the same scale as the other three families) gives
#' `b = sigma / sqrt(2)`, computed first in both `dlaplace1_lpmf` and
#' `dlaplace1_lccdf`.
#'
#' **Validation note.** `extraDistr::ddlaplace()` implements a different
#' discrete Laplace — its `scale` argument is actually a decay
#' probability `p` for the exact closed form `P(z) = (1-p)/(1+p) *
#' p^|z|`, not a continuous-Laplace `b` — confirmed numerically to NOT
#' match this family's CDF-differenced PMF (e.g. at `b=3`,
#' `p=exp(-1/3)`: `P(0) = 0.1535` here vs `0.1651` there). This package's
#' tests validate against a hand-derived CDF-difference R reference
#' instead (the documented fallback for when a package reference isn't
#' applicable), matching the CDF-differencing already used in
#' `05-04-skellam-truncation-investigation.qmd`'s exploratory plots.
#'
#' @return A brms custom_family object.
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

#' Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc()
#'
#' @description
#' Returns a `brms::stanvar()` defining `dlaplace1_lccdf`, the log
#' complementary CDF of the discrete Laplace(0, sigma) family --
#' `dlaplace1_lccdf(y, sigma)` = log P(Z > y). Same role and calling
#' convention as `skellam1_lccdf_stanvars()`. Unlike the Skellam
#' families' lccdf stanvars, this takes no threshold argument: the
#' closed-form `log1m_exp(double_exponential_lcdf(...))` has no
#' large-argument failure mode to guard against.
#'
#' @return A `brms::stanvars` object defining the `dlaplace1_lccdf` Stan
#'   function, for combining with `dlaplace1_stanvars()` via `+`.
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
  laplace_cdf <- function(x) ifelse(x < 0, 0.5 * exp(x / b), 1 - 0.5 * exp(-x / b))
  log(laplace_cdf(z + 0.5) - laplace_cdf(z - 0.5))
}

#' @rdname dlaplace1
#' @export
#' @keywords internal
posterior_predict_dlaplace1 <- function(i, prep, ...) {
  sigma <- brms::get_dpar(prep, "mu", i = i)  # brms dpar name "mu" is sigma here -- see Details
  b     <- sigma / sqrt(2)
  # Difference of two iid Exponential(rate = 1/b) draws is Laplace(0, b)
  # exactly (the discrete-Laplace analogue of skellam's "difference of
  # two iid Poissons"); rounding a continuous Laplace(0,b) draw to the
  # nearest integer reproduces this family's CDF-differenced PMF exactly,
  # since P(round(X)=z) = P(z-0.5 <= X < z+0.5) = F(z+0.5) - F(z-0.5).
  n <- length(b)
  round(stats::rexp(n, rate = 1 / b) - stats::rexp(n, rate = 1 / b))
}

#' @rdname dlaplace1
#' @export
#' @keywords internal
posterior_epred_dlaplace1 <- function(prep) {
  sigma <- brms::get_dpar(prep, "mu")  # brms dpar name "mu" is sigma here -- see Details
  0 * sigma  # E[discrete Laplace(0, sigma)] = 0; preserves draw x obs matrix dimensions
}

# ==========================================================================
# dlaplace2: discrete Laplace, free location AND free scale, uncoupled
# ==========================================================================

#' Discrete-Laplace custom family for brms (free location and scale)
#'
#' @description
#' Returns a brms custom family for the discrete Laplace distribution
#' with both location (`mu`, link = "identity") and scale (`sigma`,
#' link = "log") free, discretised via CDF differencing exactly as
#' `dlaplace1()` but centred at `mu` instead of fixed at 0:
#' `P(Z=z) = F(z+0.5) - F(z-0.5)`, `F` the continuous Laplace(mu, b) CDF.
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = dlaplace2(), stanvars = dlaplace2_stanvars(), data = ...)
#'
#' @details
#' **No naming workaround needed.** Unlike `skellam1()`/`dlaplace1()`,
#' `mu` here genuinely is the family's mean, so brms's "must have a `mu`
#' parameter" requirement (see `?skellam1` Details) is satisfied
#' directly — no forced reinterpretation.
#'
#' **No constraint coupling mu and sigma.** This is a genuine structural
#' difference from `skellam2()`, which structurally requires `sigma >=
#' |mu|` (the Skellam family's actual mean/variance relationship — see
#' `?skellam2` Details). The discrete Laplace has no such relationship:
#' `mu` and `sigma` are free, independent parameters. The point of
#' having both an asymmetric-Skellam and a free-location discrete-Laplace
#' family in this package is to compare a model where bias and spread are
#' structurally coupled (skellam2) against one where they are not
#' (dlaplace2) — do not impose any artificial coupling here.
#'
#' **sigma-to-b conversion.** Same as `dlaplace1()`: `b = sigma /
#' sqrt(2)`. `mu` is passed straight through to
#' `double_exponential_lcdf`'s own location argument (it takes location
#' and scale directly, like `normal_lcdf`), so no manual shift of `z` is
#' needed in the Stan code.
#'
#' @return A brms custom_family object.
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

#' Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc()
#' (free location and scale)
#'
#' @description
#' Returns a `brms::stanvar()` defining `dlaplace2_lccdf`, the log
#' complementary CDF of the discrete Laplace(mu, sigma) family --
#' `dlaplace2_lccdf(y, mu, sigma)` = log P(Z > y). Same role, calling
#' convention, and no-threshold-argument rationale as
#' `dlaplace1_lccdf_stanvars()`.
#'
#' @return A `brms::stanvars` object defining the `dlaplace2_lccdf` Stan
#'   function, for combining with `dlaplace2_stanvars()` via `+`.
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
  n     <- length(mu)
  round(mu + stats::rexp(n, rate = 1 / b) - stats::rexp(n, rate = 1 / b))
}

#' @rdname dlaplace2
#' @export
#' @keywords internal
posterior_epred_dlaplace2 <- function(prep) {
  brms::get_dpar(prep, "mu")  # E[discrete Laplace(mu, sigma)] = mu
}
