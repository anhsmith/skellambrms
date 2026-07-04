# ==========================================================================
# bipois: joint bivariate Poisson via trivariate reduction
# (Holgate 1964; recurrence: "A fast way to calculate the bivariate poisson
# in STAN", stan-users Google Group, March 2016 -- thread SUcp-ktkXn4,
# posted by Andre, refined by Bob Carpenter -- adapted below from
# log-rate/`poisson_log_lpmf` parameterisation to the natural-scale-rate
# parameterisation this package's dpars use post-link.)
#
# Unlike skellam1/2, dlaplace1/2, dnorm1/2 (all of which model a single
# *difference* response, y_em - y_lb), bipois models the *joint* pair
# (y_em, y_lb) directly -- the escape from the d ~ y_lb regression-to-the-
# mean problem documented in 05-07-audit-track-step2-redesign.qmd. Both
# y_em and y_lb are generated as y_em = N_shared + N10, y_lb = N_shared +
# N01, with N_shared ~ Poisson(mu), N10 ~ Poisson(lambdaem), N01 ~
# Poisson(lambdalb) mutually independent given their rates. N_shared is
# unobserved and marginalised out analytically (05-07 @sec-likelihood):
#
#   P(y_em=x, y_lb=y) = sum_{k=0}^{min(x,y)}
#     P(N_shared=k) P(N10=x-k) P(N01=y-k)
#
# rather than evaluated as a naive per-k sum of three separate Poisson
# lpmfs, the Stan implementation below uses the cited bipois2 incremental
# recurrence: the k=0 term is computed directly, and each subsequent term
# is obtained from the previous one by a ratio update in log space,
# accumulated via log_sum_exp -- see bipois_stan_funs below for the
# derivation of that ratio and its correspondence to the cited source.

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Joint EM/logbook bivariate-Poisson custom family for brms
#'
#' @description
#' Returns a brms custom family for the joint distribution of a matched pair
#' of counts, `(y_em, y_lb)`, constructed via trivariate reduction: `y_em =
#' N_shared + N10`, `y_lb = N_shared + N01`, with `N_shared ~
#' Poisson(mu)`, `N10 ~ Poisson(lambdaem)`, `N01 ~ Poisson(lambdalb)`
#' mutually independent given their rates. All three rates are link = "log".
#' See `05-07-audit-track-step2-redesign.qmd` in the `tnc001-belize-em`
#' project for the full generative-model rationale (why this replaces a
#' `d = y_em - y_lb ~ y_lb` design).
#'
#' `y_em` is the family's response; `y_lb` is passed in as supplementary
#' integer data via brms's `vint()` addition term, since brms's
#' `custom_family()` machinery is built around a single declared response
#' column -- see Details.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y_em | vint(y_lb) ~ ...),
#'     family   = bipois(),
#'     stanvars = bipois_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **Naming note.** Same forced naming as `skellam1()`/`dlaplace1()`/
#' `dnorm1()`: `brms::custom_family()` requires a dpar literally named
#' `"mu"` (`stop2("All families must have a 'mu' parameter.")`,
#' unconditional). Here it is bound to `lambda_shared`, the rate of the
#' component shared between `y_em` and `y_lb` -- not a mean of either
#' response individually. `lambdaem` (EM-only rate) and `lambdalb`
#' (LB-only rate) are the other two dpars, plainly named (no forced
#' reinterpretation needed for those two).
#'
#' **Why `y_lb` travels via `vint()`, not as a second response.** brms's
#' `custom_family()` API supports exactly one declared response column
#' (`Y`) plus optional supplementary integer/real data (`vint()`/`vreal()`
#' addition terms) -- the same mechanism used for, e.g., binomial trial
#' counts in the brms custom-families vignette. There is no
#' *undeclared-response* concept for a genuinely joint two-count
#' likelihood; `vint(y_lb)` is the correct fit for that gap, not a
#' workaround. This does mean `y_lb` is *not* itself treated as
#' brms-modelled response data (no missing-value handling, no
#' resp_*() addition terms apply to it) -- it is fixed, observed
#' per-row data, consistent with the fact that every row used here comes
#' from the matched (both-observed) subset.
#'
#' **Order of dpars matters for the generated Stan call.** brms generates
#' `target += bipois_lpmf(Y[n] | mu[n], lambdaem[n], lambdalb[n],
#' vint1[n])` -- dpars in the order declared here, then vint/vreal args in
#' the order declared in `vars`. `bipois_stan_funs` (stanfunctions via
#' `bipois_stanvars()`) declares `bipois_lpmf` with exactly this argument
#' order; changing the order here without changing the Stan signature (or
#' vice versa) silently swaps which rate governs which count.
#'
#' @return A brms custom_family object.
#' @export
bipois <- function() {
  brms::custom_family(
    name  = "bipois",
    dpars = c("mu", "lambdaem", "lambdalb"),  # mu = lambda_shared -- see Details
    links = c("log", "log", "log"),
    lb    = c(0, 0, 0),
    type  = "int",
    vars  = "vint1[n]"
  )
}

#' @rdname bipois
#' @export
bipois_stanvars <- function() {
  brms::stanvar(block = "functions", scode = bipois_stan_funs)
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Derivation of the recurrence, in this package's natural-scale-rate terms
# (the cited stan-users thread works in log-rates via poisson_log_lpmf;
# translated here since this package's dpars are already inv-link-
# transformed to natural scale by the time they reach the lpmf, the same
# convention as every other family in this package -- see skellam1_lpmf,
# which likewise takes `sigma`, not `log(sigma)`).
#
# Write lambda_shared = mu, and let term(k) = P(N_shared=k) P(N10=r-k)
# P(N01=s-k) for r = y_em, s = y_lb, m = min(r,s). Then:
#
#   term(0) = exp(-mu) * [lambdaem^r exp(-lambdaem) / r!]
#                       * [lambdalb^s exp(-lambdalb) / s!]
#
# so log(term(0)) = poisson_lpmf(r|lambdaem) + poisson_lpmf(s|lambdalb) - mu
# (the `+lambdaem +lambdalb` inside each poisson_lpmf cancels algebraically
# against the `-lambdaem -lambdalb` that would otherwise appear from
# factoring exp(-mu-lambdaem-lambdalb) out front -- this is exactly the
# same cancellation used in the cited thread's `ss <- poisson_log_log(r,
# mu1) + poisson_log_log(s, mu2) - exp(mu3)` starting term).
#
# The ratio between consecutive terms is:
#   term(k) / term(k-1) = [(r-k+1)/lambdaem] * [(s-k+1)/lambdalb] * [mu/k]
# (one fewer factor of lambdaem in the N10 count as k increases by one,
# one fewer of lambdalb, one more of mu -- and the corresponding
# factorial/combinatorial adjustment (r-k+1), (s-k+1), 1/k). In log space:
#   log(term(k)) = log(term(k-1)) + log(r-k+1) + log(s-k+1) - log(k)
#                  + log(mu) - log(lambdaem) - log(lambdalb)
# which is exactly the cited thread's per-step update (`log_s <- log_s +
# log(r-k+1) + mus + log(s-k+1) - log(k)`, with `mus = -mu1-mu2+mu3`
# translating directly to `log(mu) - log(lambdaem) - log(lambdalb)` once
# mu1, mu2, mu3 are read as log-rates). Each term is accumulated into the
# running total via log_sum_exp -- no separate normalising-constant term is
# needed (unlike skellam1_lccdf_stan's tail sum), since this sum is finite
# (m+1 terms) and exactly equals the marginalised joint log-likelihood, not
# a survival function.
#
# No large-argument blowup risk analogous to skellam1/skellam2's Bessel
# function: m = min(y_em, y_lb) is bounded by the data itself (expected
# ~50-60 terms at most for this project's data, per 05-07
# @sec-likelihood), and every term here is a sum/difference of ordinary
# Poisson log-densities, not a Bessel function of large order. No
# iteration cap or normal-approximation branch is needed for that reason.
bipois_stan_funs <- "
  real bipois_lpmf(int y_em, real mu, real lambdaem, real lambdalb, int y_lb) {
    int m = min(y_em, y_lb);
    real ss = poisson_lpmf(y_em | lambdaem) + poisson_lpmf(y_lb | lambdalb) - mu;
    if (m > 0) {
      real log_ratio = log(mu) - log(lambdaem) - log(lambdalb);
      real log_term = ss;
      for (k in 1:m) {
        log_term += log(y_em - k + 1) + log(y_lb - k + 1) - log(k) + log_ratio;
        ss = log_sum_exp(ss, log_term);
      }
    }
    return ss;
  }
"

# --------------------------------------------------------------------------
# R-side reference implementation
# --------------------------------------------------------------------------

# Direct brute-force evaluation of the same marginalised sum, computed
# term-by-term from the original P(N_shared=k)*P(N10=x-k)*P(N01=y-k)
# definition (05-07 @sec-likelihood) rather than via the recurrence --
# an independent route to the same quantity, used to validate the Stan
# recurrence above and to power log_lik_bipois()/posterior_epred_bipois()
# (evaluated post-hoc, not inside the sampler's hot loop, so there is no
# reason to use the recurrence's algebraic shortcuts here). Not exported:
# internal reference only, exactly the role skellam1_lpmf_r/dlaplace1_lpmf_r
# etc. play in truncation.R for their families.
bipois_lpmf_r <- function(y_em, y_lb, mu, lambdaem, lambdalb) {
  n <- max(length(y_em), length(y_lb), length(mu), length(lambdaem), length(lambdalb))
  y_em     <- rep_len(y_em, n)
  y_lb     <- rep_len(y_lb, n)
  mu       <- rep_len(mu, n)
  lambdaem <- rep_len(lambdaem, n)
  lambdalb <- rep_len(lambdalb, n)

  vapply(seq_len(n), function(i) {
    m <- min(y_em[i], y_lb[i])
    k <- 0:m
    log_terms <- stats::dpois(k, mu[i], log = TRUE) +
      stats::dpois(y_em[i] - k, lambdaem[i], log = TRUE) +
      stats::dpois(y_lb[i] - k, lambdalb[i], log = TRUE)
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }, numeric(1))
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname bipois
#' @export
#' @keywords internal
log_lik_bipois <- function(i, prep) {
  mu       <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared -- see Details
  lambdaem <- brms::get_dpar(prep, "lambdaem", i = i)
  lambdalb <- brms::get_dpar(prep, "lambdalb", i = i)
  y_em <- prep$data$Y[i]
  y_lb <- prep$data$vint1[i]
  bipois_lpmf_r(y_em, y_lb, mu, lambdaem, lambdalb)
}

#' @rdname bipois
#' @export
#' @keywords internal
posterior_predict_bipois <- function(i, prep, ...) {
  mu       <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared -- see Details
  lambdaem <- brms::get_dpar(prep, "lambdaem", i = i)
  lambdalb <- brms::get_dpar(prep, "lambdalb", i = i)
  y_lb <- prep$data$vint1[i]
  # y_lb is fixed, observed data (not itself re-simulated -- see "Why y_lb
  # travels via vint()" in ?bipois). Consistent with that, y_em is
  # predicted *conditional on the real y_lb*, via the same closed-form
  # conditional split documented in 05-07 @sec-full-data-fit:
  # N_shared | y_lb ~ Binomial(y_lb, mu / (mu + lambdalb)); N10 fresh from
  # its own marginal; y_em = N_shared + N10.
  p_shared <- mu / (mu + lambdalb)
  n_shared <- stats::rbinom(length(mu), size = y_lb, prob = p_shared)
  n10      <- stats::rpois(length(mu), lambdaem)
  n_shared + n10
}

#' @rdname bipois
#' @export
#' @keywords internal
posterior_epred_bipois <- function(prep) {
  mu       <- brms::get_dpar(prep, "mu")        # lambda_shared -- see Details
  lambdaem <- brms::get_dpar(prep, "lambdaem")
  lambdalb <- brms::get_dpar(prep, "lambdalb")
  y_lb <- prep$data$vint1
  # E[y_em | y_lb] = E[N_shared | y_lb] + E[N10] = y_lb * mu/(mu+lambdalb) + lambdaem,
  # the same conditional split as posterior_predict_bipois above, in
  # expectation rather than simulated.
  y_lb_mat <- matrix(y_lb, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdalb)
  y_lb_mat * p_shared + lambdaem
}
