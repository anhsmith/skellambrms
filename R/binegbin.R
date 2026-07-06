# ==========================================================================
# binegbin: joint bivariate Negative-Binomial via trivariate reduction
#
# The overdispersed sibling of bipois (see bipois.R). Same trivariate-
# reduction construction -- y_em = N_shared + N10, y_lb = N_shared + N01,
# the three latent counts mutually independent given their rates -- but each
# latent count is Negative-Binomial rather than Poisson:
#
#   N_shared ~ NB2(mu,       shapes)   (shared component; drives correlation)
#   N10      ~ NB2(lambdaem, shapex)   (EM-only excess)
#   N01      ~ NB2(lambdalb, shapex)   (LB-only excess)
#
# NB2(m, phi) is Stan's neg_binomial_2 / R's dnbinom(size = phi, mu = m):
# mean m, variance m + m^2/phi. shapes is the shared-component dispersion,
# shapex the (shared) excess dispersion -- 05-09 Decision 2's "2 kappa"
# structure (kappa_ns separate; kappa_excess shared for N10/N01).
#
# WHY NEGBIN AND NOT AN OLRE ON bipois. The plain-Poisson bipois cannot be
# overdispersed (Var == mean for each latent count), so it underfits the
# real marginal variances by ~10x (05-10: YFT Var(y_em) fitted 16.6 vs
# observed 179) and Var(d) by ~3.5x. The obvious fix -- add a per-set
# observation-level random effect (OLRE) on the excess components -- FAILS
# synthetic recovery: with one bivariate observation per set but three
# per-set latent deviates (mu-OLRE + two excess OLREs), the excess deviates
# act as residual-absorbers, their population SD collapses toward the prior
# mode, and drawing fresh deviates does NOT regenerate the observed spread
# (recovered excess SD 0.37 vs true 0.85; fresh-deviate Var(d) 2.9 vs true
# 19.2). A conditional posterior-predictive check hides this completely --
# only a marginal (fresh-deviate) check exposes it. binegbin carries the
# dispersion in SCALAR shapes/shapex instead, estimated from aggregate
# mean-variance mismatch across sets -- identifiable, no per-set overfitting,
# clean marginal PPC, and consistent with the review-track NegBin models.
#
# LIKELIHOOD. N_shared is unobserved and marginalised out analytically,
# exactly as in bipois -- the sum structure is identical, only the component
# pmfs change from Poisson to NegBin:
#
#   P(y_em=x, y_lb=y) = sum_{k=0}^{min(x,y)}
#     NB2(k | mu, shapes) NB2(x-k | lambdaem, shapex) NB2(y-k | lambdalb, shapex)
#
# This is NOT the "two stacked marginalisations" problem flagged in 05-07
# (Gamma-mixing a Poisson while keeping it Poisson): the components are
# DIRECTLY NegBin, so the marginalisation sum is the same finite sum as
# bipois with neg_binomial_2_lpmf swapped in for poisson_lpmf. The bipois
# incremental recurrence does not carry over cleanly (NegBin consecutive-term
# ratios are less simple than Poisson's), so binegbin_lpmf evaluates the sum
# directly via log_sum_exp. m = min(y_em, y_lb) is bounded by the data
# (~50-60 terms at most for this project), so the direct sum is not a
# performance concern -- the same reasoning bipois's own docs give for why no
# large-argument branch is needed.
#
# Validation (grid cross-check of the Stan lpmf against the independent R
# brute-force reference to ~1e-14, normalisation to 1, moment identities, and
# synthetic parameter recovery with a marginal PPC) is documented in the
# tnc001-belize-em project alongside the bipois validation (05-08 lineage).

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Joint EM/logbook bivariate-Negative-Binomial custom family for brms
#'
#' @description
#' Overdispersed sibling of [bipois()]. Returns a brms custom family for the
#' joint distribution of a matched count pair `(y_em, y_lb)` via trivariate
#' reduction with Negative-Binomial (rather than Poisson) latent components:
#' `y_em = N_shared + N10`, `y_lb = N_shared + N01`, with
#' `N_shared ~ NB2(mu, shapes)`, `N10 ~ NB2(lambdaem, shapex)`,
#' `N01 ~ NB2(lambdalb, shapex)` mutually independent given their rates.
#' `NB2(m, phi)` has mean `m` and variance `m + m^2/phi` (Stan
#' `neg_binomial_2`; R `dnbinom(size = phi, mu = m)`).
#'
#' Five dpars: the three rates (`mu` = shared rate, `lambdaem`/`lambdalb` =
#' EM-/LB-only rates) plus two dispersions -- `shapes` for the shared
#' component and `shapex` shared across the two excess components. All five
#' use `link = "log"`. Supply the excess rates through a non-linear formula
#' without an explicit `exp()` (the log link applies it): `nlf(lambdaem ~ lamx)`
#' gives `lambdaem = exp(lamx)`.
#'
#' See the `binegbin.R` file header and the `tnc001-belize-em` project docs
#' (05-07 generative rationale; the OLRE-failure / NegBin-resolution finding)
#' for why NegBin components are used instead of an observation-level random
#' effect on [bipois()].
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y_em | vint(y_lb) ~ 1,
#'        mu ~ 1 + (1 | vessel),
#'        nlf(lambdaem ~ lamx), nlf(lambdalb ~ lamx), lamx ~ 1,
#'        shapes ~ 1, shapex ~ 1, nl = TRUE),
#'     family   = binegbin(),
#'     stanvars = binegbin_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **Forced `mu` naming, and `y_lb` via `vint()`.** Identical conventions to
#' [bipois()] -- `mu` is brms's mandatory dpar name, here bound to the shared
#' component's rate (`lambda_shared`), not a mean of either response; `y_lb`
#' travels as supplementary integer data through `vint()` because
#' `custom_family()` declares a single response column. See [bipois()] for
#' the full explanation.
#'
#' **Order of dpars matters for the generated Stan call.** brms generates
#' `target += binegbin_lpmf(Y[n] | mu[n], lambdaem[n], lambdalb[n],
#' shapes[n], shapex[n], vint1[n])` -- dpars in the order declared here, then
#' vint args. `binegbin_stan_funs` declares `binegbin_lpmf` with exactly this
#' signature; reordering one without the other silently swaps which rate or
#' dispersion governs which component.
#'
#' @return A brms custom_family object.
#' @export
binegbin <- function() {
  brms::custom_family(
    name  = "binegbin",
    dpars = c("mu", "lambdaem", "lambdalb", "shapes", "shapex"),
    links = c("log", "log", "log", "log", "log"),
    lb    = c(0, 0, 0, 0, 0),
    type  = "int",
    vars  = "vint1[n]"
  )
}

#' @rdname binegbin
#' @export
binegbin_stanvars <- function() {
  brms::stanvar(block = "functions", scode = binegbin_stan_funs)
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Direct marginalisation sum (see file header for why not a recurrence).
# Every term is a sum of three neg_binomial_2 log-densities; accumulated via
# log_sum_exp over k = 0..min(y_em, y_lb). neg_binomial_2_lpmf(0 | m, phi) is
# well-defined, so the k = 0 term and zero-count responses need no special
# casing.
binegbin_stan_funs <- "
  real binegbin_lpmf(int y_em, real mu, real lambdaem, real lambdalb,
                     real shapes, real shapex, int y_lb) {
    int m = min(y_em, y_lb);
    vector[m + 1] lp;
    for (k in 0:m) {
      lp[k + 1] = neg_binomial_2_lpmf(k        | mu,       shapes)
                + neg_binomial_2_lpmf(y_em - k | lambdaem, shapex)
                + neg_binomial_2_lpmf(y_lb - k | lambdalb, shapex);
    }
    return log_sum_exp(lp);
  }
"

# --------------------------------------------------------------------------
# R-side reference implementation
# --------------------------------------------------------------------------

# Independent brute-force evaluation of the same sum via R's dnbinom (an
# independent route from Stan's neg_binomial_2), used to validate the Stan
# lpmf and to power log_lik_binegbin()/posterior_epred_binegbin() post-hoc.
# Internal reference only, mirroring bipois_lpmf_r's role.
binegbin_lpmf_r <- function(y_em, y_lb, mu, lambdaem, lambdalb, shapes, shapex) {
  n <- max(length(y_em), length(y_lb), length(mu), length(lambdaem),
           length(lambdalb), length(shapes), length(shapex))
  y_em     <- rep_len(y_em, n)
  y_lb     <- rep_len(y_lb, n)
  mu       <- rep_len(mu, n)
  lambdaem <- rep_len(lambdaem, n)
  lambdalb <- rep_len(lambdalb, n)
  shapes   <- rep_len(shapes, n)
  shapex   <- rep_len(shapex, n)

  vapply(seq_len(n), function(i) {
    m <- min(y_em[i], y_lb[i])
    k <- 0:m
    log_terms <- stats::dnbinom(k,           size = shapes[i], mu = mu[i],       log = TRUE) +
      stats::dnbinom(y_em[i] - k,            size = shapex[i], mu = lambdaem[i], log = TRUE) +
      stats::dnbinom(y_lb[i] - k,            size = shapex[i], mu = lambdalb[i], log = TRUE)
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }, numeric(1))
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname binegbin
#' @export
#' @keywords internal
log_lik_binegbin <- function(i, prep) {
  mu       <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaem <- brms::get_dpar(prep, "lambdaem", i = i)
  lambdalb <- brms::get_dpar(prep, "lambdalb", i = i)
  shapes   <- brms::get_dpar(prep, "shapes", i = i)
  shapex   <- brms::get_dpar(prep, "shapex", i = i)
  y_em <- prep$data$Y[i]
  y_lb <- prep$data$vint1[i]
  binegbin_lpmf_r(y_em, y_lb, mu, lambdaem, lambdalb, shapes, shapex)
}

#' @rdname binegbin
#' @export
#' @keywords internal
posterior_predict_binegbin <- function(i, prep, ...) {
  mu       <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaem <- brms::get_dpar(prep, "lambdaem", i = i)
  lambdalb <- brms::get_dpar(prep, "lambdalb", i = i)
  shapes   <- brms::get_dpar(prep, "shapes", i = i)
  shapex   <- brms::get_dpar(prep, "shapex", i = i)
  y_lb <- prep$data$vint1[i]
  # y_em predicted conditional on the real observed y_lb. Unlike bipois, the
  # conditional split N_shared | y_lb is NOT Binomial (a NegBin sum condition
  # is not Binomial); it is P(N_shared = k | y_lb) proportional to
  # NB2(k | mu, shapes) NB2(y_lb - k | lambdalb, shapex) over k = 0..y_lb.
  # Sample that discrete conditional, then add a fresh N10 ~ NB2(lambdaem,
  # shapex).
  S <- length(mu)
  out <- integer(S)
  for (s in seq_len(S)) {
    k <- 0:y_lb
    lw <- stats::dnbinom(k,        size = shapes[s], mu = mu[s],       log = TRUE) +
          stats::dnbinom(y_lb - k, size = shapex[s], mu = lambdalb[s], log = TRUE)
    w <- exp(lw - max(lw))
    n_shared <- if (y_lb == 0) 0L else sample(k, 1, prob = w)
    out[s] <- n_shared + stats::rnbinom(1, size = shapex[s], mu = lambdaem[s])
  }
  out
}

#' @rdname binegbin
#' @export
#' @keywords internal
posterior_epred_binegbin <- function(prep) {
  mu       <- brms::get_dpar(prep, "mu")        # lambda_shared
  lambdaem <- brms::get_dpar(prep, "lambdaem")
  lambdalb <- brms::get_dpar(prep, "lambdalb")
  y_lb <- prep$data$vint1
  # E[y_em | y_lb] = E[N_shared | y_lb] + lambdaem. For binegbin there is no
  # closed-form E[N_shared | y_lb] as clean as bipois's y_lb * mu/(mu+lambdalb);
  # this returns the analogous point approximation using the marginal shared
  # fraction mu/(mu+lambdalb), adequate for epred display (posterior_predict
  # uses the exact discrete conditional).
  y_lb_mat <- matrix(y_lb, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdalb)
  y_lb_mat * p_shared + lambdaem
}
