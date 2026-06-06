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
