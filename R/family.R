#' Symmetric Skellam custom family for brms
#'
#' @description
#' Returns a brms custom family for the symmetric Skellam distribution,
#' Skellam(mu, mu) — the distribution of the difference of two independent
#' Poisson(mu) random variables. The single parameter mu (link = "log")
#' controls dispersion; the mean is always zero.
#'
#' Use in a brm() call as:
#'   brm(y ~ ..., family = skellam1(), stanvars = skellam1_stanvars(), data = ...)
#'
#' @return A brms custom_family object.
#' @export
skellam1 <- function() {
  brms::custom_family(
    name  = "skellam1",
    dpars = "mu",
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
#' complementary CDF of the symmetric Skellam(mu, mu) distribution —
#' `skellam1_lccdf(y, mu)` = log P(delta > y). brms's generic truncation
#' machinery (`resp_trunc()`) locates a custom family's log-CCDF by name
#' convention (`<family>_lccdf`), so adding this stanvar alongside
#' `skellam1_stanvars()` is sufficient to support truncated fits,
#' including a row-varying lower bound — no other wiring is required.
#'
#' @details
#' For `mu` above `normal_approx_threshold`, the exact log-CCDF — an
#' iterative tail-sum of the Skellam PMF, each term a Bessel function
#' evaluation — is replaced by a normal approximation, using
#' Var(Skellam(mu, mu)) = 2 * mu. This guards against two confirmed
#' failure modes, both triggered by an unadapted HMC proposal pushing
#' `mu` to an extreme value during warmup (the log link places no
#' ceiling on it):
#'
#' - A crash (`std::bad_alloc`) from `log_modified_bessel_first_kind`
#'   being evaluated at an enormous Bessel order.
#' - A slow-motion version of the same problem: `mu` in the hundreds
#'   still triggers the expensive exact loop, and if many rows do this
#'   within a single deep NUTS tree the cost compounds multiplicatively
#'   rather than crashing outright — observed as 200+ CPU-seconds and
#'   several GB of memory consumed without completing one iteration.
#'
#' The exact loop also carries a hard cap of 500 iterations past `y`,
#' with an early exit once the tail term becomes negligible (more than
#' ~40 log-units below the running sum). These guard the same two
#' failure modes as the threshold itself and are not configurable here.
#'
#' The default threshold of 100 is **not a universal constant** — it
#' was calibrated to one project's data, where real per-taxon `mu`
#' estimates topped out around 30. The 3x margin above that (rather
#' than setting the threshold at, say, 35) exists because HMC warmup
#' transiently proposes values well outside any final posterior
#' estimate, not because 30 itself needed padding. When using this
#' function with a different count scale, consider what *implausible
#' but reachable during warmup* looks like for your `mu`, not just
#' your expected posterior range, and set the threshold a few-fold
#' above that. Setting it too low pays for the normal approximation's
#' bias more often than necessary; setting it too high re-exposes the
#' crash/slow-blowup risk this exists to prevent.
#'
#' @param normal_approx_threshold Numeric scalar; `mu` values above this
#'   use the normal approximation instead of the exact Bessel-sum tail.
#'   Default `100`. See Details for how to choose this for your data.
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
  mu <- brms::get_dpar(prep, "mu", i = i)
  y  <- prep$data$Y[i]
  # besselI(..., expon.scaled = TRUE) returns I_nu(x) * exp(-x), so
  # log(besselI(2*mu, |y|, expon.scaled=TRUE)) = log(I_|y|(2mu)) - 2mu
  # which equals the full log-PMF: -2mu + log(I_|y|(2mu))
  log(besselI(2 * mu, abs(y), expon.scaled = TRUE))
}

#' @rdname skellam1
#' @export
#' @keywords internal
posterior_predict_skellam1 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  skellam::rskellam(length(mu), lambda1 = mu, lambda2 = mu)
}

#' @rdname skellam1
#' @export
#' @keywords internal
posterior_epred_skellam1 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  0 * mu  # E[Skellam(mu, mu)] = 0; preserves draw x obs matrix dimensions
}
