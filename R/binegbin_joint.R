# ==========================================================================
# binegbin_joint: censoring-aware bivariate Negative-Binomial
#
# The censoring-aware extension of binegbin (see binegbin.R). Same generative
# model -- y_em = N_shared + N10, y_lb = N_shared + N01 with
#
#   N_shared ~ NB2(mu,       shapes)   (shared component; drives correlation)
#   N10      ~ NB2(lambdaem, shapex)   (EM-only excess)
#   N01      ~ NB2(lambdalb, shapex)   (LB-only excess)
#
# and the same five dpars and links -- but each row now carries a second
# supplementary integer, an observation flag em_obs, alongside y_lb:
#
#   em_obs == 1 (matched set):  full joint binegbin lpmf on (y_em, y_lb) --
#       identical to binegbin_lpmf, byte-for-byte.
#   em_obs == 0 (LB-only set):  the LB MARGINAL of the SAME bivariate model --
#       P(y_lb) = sum_k NB2(k | mu, shapes) NB2(y_lb - k | lambdalb, shapex),
#       i.e. the joint with the y_em (N10) term integrated out over all y_em.
#
# WHY A DEDICATED FAMILY AND NOT TWO SEPARATE FITS. The unmatched (LB-only)
# rows never observe y_em, so a plain binegbin fit could only use the matched
# rows. But the LB-only rows still carry information about the SHARED
# structure (mu, shapes, lambdalb) and the vessel/trip random effects: their
# y_lb is a draw from the same bivariate model, merely with its EM margin
# unobserved. Integrating y_em out (rather than dropping those rows, or --
# worse -- giving them their own single-dispersion neg_binomial_2 on y_lb,
# which is a DIFFERENT model inconsistent with the matched decomposition)
# lets one brm() call pool all rows under one coherent likelihood. lambdaem
# (the EM-only excess rate) and the directional bias between the two excess
# rates are then identified ONLY by the matched rows; the LB-only rows sharpen
# mu / shapes / lambdalb and the shared vessel+trip structure. This is the
# standard partially-observed-margin construction (a censored/"missing at
# random on the EM side" likelihood), not a heuristic.
#
# EXACT RELATIONSHIP TO binegbin. On the matched branch this family's lpmf is
# the binegbin lpmf, exactly (same sum, same neg_binomial_2 terms). The
# LB-only branch is the y_em-integrated marginal of that same joint. Two
# consequences the package tests pin down:
#   * sum over y_em of the matched-branch lpmf == the LB-only-branch lpmf
#     (marginal identity), and
#   * binegbin_joint_lpmf(em_obs == 1) == binegbin_lpmf on identical inputs
#     (equivalence).
# The second identity is what licenses reading a binegbin fit's cached draws
# as the em_obs == 1 slice of a binegbin_joint fit and vice versa -- the two
# families cannot silently drift apart without a test going red.
#
# Validated by the marginal identity, the binegbin equivalence, an
# expose_functions grid cross-check of the Stan lpmf against an independent R
# brute-force reference (~1e-14), and the conditional-prediction identity
# (posterior_predict draws == joint / marginal). See the tnc001-belize-em
# project docs (05-13/05-14 lineage; the full-year censored fit) for the
# modelling context.

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Censoring-aware joint EM/logbook bivariate-Negative-Binomial family for brms
#'
#' @description
#' Censoring-aware extension of [binegbin()]. Models the same trivariate-
#' reduction bivariate Negative-Binomial pair `(y_em, y_lb)` -- `y_em =
#' N_shared + N10`, `y_lb = N_shared + N01`, with `N_shared ~ NB2(mu,
#' shapes)`, `N10 ~ NB2(lambdaem, shapex)`, `N01 ~ NB2(lambdalb, shapex)`
#' mutually independent given their rates -- but allows the EM margin (`y_em`)
#' to be UNOBSERVED on some rows. Each row carries two supplementary integers
#' via `vint()`: `y_lb` (the always-observed logbook count) and `em_obs` (a
#' 0/1 flag marking whether `y_em` was observed for that row).
#'
#' On `em_obs == 1` (matched) rows the likelihood is the full joint
#' [binegbin()] lpmf on `(y_em, y_lb)`. On `em_obs == 0` (LB-only) rows it is
#' the `y_em`-integrated marginal of that same joint,
#' `P(y_lb) = sum_k NB2(k | mu, shapes) NB2(y_lb - k | lambdalb, shapex)` --
#' NOT a separate single-dispersion `neg_binomial_2` on `y_lb`, which would be
#' a different model inconsistent with the matched decomposition. This lets
#' one `brm()` call pool matched and LB-only rows under one coherent
#' likelihood: `lambdaem` and the EM/LB bias are identified only by the
#' matched rows, while the LB-only rows sharpen `mu`, `shapes`, `lambdalb`,
#' and the shared vessel/trip random-effect structure.
#'
#' Five dpars, identical to [binegbin()]: the three rates (`mu` = shared rate,
#' `lambdaem`/`lambdalb` = EM-/LB-only rates) plus two dispersions -- `shapes`
#' for the shared component and `shapex` shared across the two excess
#' components. All five use `link = "log"` (see [binegbin()]). To share the
#' excess level across the two rates and split them by a directional bias,
#' supply them through non-linear formulas *without* an explicit `exp()` --
#' the log link applies it, so `nlf(lambdaem ~ lamx + methd)` gives
#' `lambdaem = exp(lamx + methd)`.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y_em | vint(y_lb, em_obs) ~ 1,
#'        mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
#'        nlf(lambdaem ~ lamx + methd),
#'        nlf(lambdalb ~ lamx - methd),
#'        lamx ~ 1, methd ~ 1,
#'        shapes ~ 1, shapex ~ 1, nl = TRUE),
#'     family   = binegbin_joint(),
#'     stanvars = binegbin_joint_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **Two `vint()` arguments, in declared order.** brms appends `vint()`
#' integers to the generated lpmf call in the order they are listed in the
#' formula's `vint()` term, matching the `vars` declared here
#' (`c("vint1[n]", "vint2[n]")`): so `vint(y_lb, em_obs)` binds `vint1 = y_lb`
#' and `vint2 = em_obs`. brms generates `target += binegbin_joint_lpmf(Y[n] |
#' mu[n], lambdaem[n], lambdalb[n], shapes[n], shapex[n], vint1[n], vint2[n])`
#' -- dpars in the order declared here, then the two vint args.
#' `binegbin_joint_stan_funs` declares `binegbin_joint_lpmf` with exactly this
#' signature; reordering the dpars or the two `vint()` terms without matching
#' the Stan signature silently swaps which quantity governs which component or
#' which integer is the branch flag.
#'
#' **Forced `mu` naming, and the second count via `vint()`.** Identical
#' conventions to [binegbin()]/[bipois()] -- `mu` is brms's mandatory dpar
#' name, here bound to the shared component's rate, not a mean of either
#' response; `y_lb` (and `em_obs`) travel as supplementary integer data
#' through `vint()` because `custom_family()` declares a single response
#' column. See [bipois()] for the full explanation.
#'
#' **Relationship to [binegbin()].** On `em_obs == 1` rows this family's lpmf
#' equals the [binegbin()] lpmf exactly (same marginalisation sum). The
#' `em_obs == 0` branch is the `y_em`-integrated marginal of that same
#' bivariate model. The package tests pin both identities (marginal identity;
#' binegbin equivalence).
#'
#' @return A brms custom_family object.
#' @export
binegbin_joint <- function() {
  brms::custom_family(
    name  = "binegbin_joint",
    dpars = c("mu", "lambdaem", "lambdalb", "shapes", "shapex"),
    links = c("log", "log", "log", "log", "log"),
    lb    = c(0, 0, 0, 0, 0),
    type  = "int",
    vars  = c("vint1[n]", "vint2[n]")  # vint1 = y_lb, vint2 = em_obs
  )
}

#' @rdname binegbin_joint
#' @export
binegbin_joint_stanvars <- function() {
  brms::stanvar(block = "functions", scode = binegbin_joint_stan_funs)
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Branching marginalisation sum. The em_obs == 1 branch is byte-for-byte the
# binegbin_lpmf body (direct log_sum_exp over k = 0..min(y_em, y_lb); see
# binegbin.R for why a direct sum rather than a recurrence). The em_obs == 0
# branch drops the y_em (N10) term and sums over k = 0..y_lb, i.e. the same
# joint with the EM margin integrated out. neg_binomial_2_lpmf(0 | m, phi) is
# well-defined, so zero counts and the k = 0 term need no special casing.
binegbin_joint_stan_funs <- "
  real binegbin_joint_lpmf(int y_em, real mu, real lambdaem, real lambdalb,
                           real shapes, real shapex, int y_lb, int em_obs) {
    if (em_obs == 1) {
      int m = min(y_em, y_lb);
      vector[m + 1] lp;
      for (k in 0:m) {
        lp[k + 1] = neg_binomial_2_lpmf(k        | mu,       shapes)
                  + neg_binomial_2_lpmf(y_em - k | lambdaem, shapex)
                  + neg_binomial_2_lpmf(y_lb - k | lambdalb, shapex);
      }
      return log_sum_exp(lp);
    } else {
      vector[y_lb + 1] lp;
      for (k in 0:y_lb) {
        lp[k + 1] = neg_binomial_2_lpmf(k        | mu,       shapes)
                  + neg_binomial_2_lpmf(y_lb - k | lambdalb, shapex);
      }
      return log_sum_exp(lp);
    }
  }
"

# --------------------------------------------------------------------------
# R-side reference implementation
# --------------------------------------------------------------------------

# Independent brute-force evaluation of the same branching sum via R's
# dnbinom, used to validate the Stan lpmf and to power
# log_lik_binegbin_joint()/posterior_predict_binegbin_joint() post-hoc.
# Internal reference only, mirroring binegbin_lpmf_r's role. Vectorised over
# all arguments (recycled to common length); em_obs selects the branch
# per-row.
binegbin_joint_lpmf_r <- function(y_em, y_lb, em_obs, mu, lambdaem, lambdalb,
                                  shapes, shapex) {
  n <- max(length(y_em), length(y_lb), length(em_obs), length(mu),
           length(lambdaem), length(lambdalb), length(shapes), length(shapex))
  y_em     <- rep_len(y_em, n)
  y_lb     <- rep_len(y_lb, n)
  em_obs   <- rep_len(em_obs, n)
  mu       <- rep_len(mu, n)
  lambdaem <- rep_len(lambdaem, n)
  lambdalb <- rep_len(lambdalb, n)
  shapes   <- rep_len(shapes, n)
  shapex   <- rep_len(shapex, n)

  vapply(seq_len(n), function(i) {
    if (em_obs[i] == 1) {
      m <- min(y_em[i], y_lb[i])
      k <- 0:m
      log_terms <- stats::dnbinom(k,          size = shapes[i], mu = mu[i],       log = TRUE) +
        stats::dnbinom(y_em[i] - k,           size = shapex[i], mu = lambdaem[i], log = TRUE) +
        stats::dnbinom(y_lb[i] - k,           size = shapex[i], mu = lambdalb[i], log = TRUE)
    } else {
      k <- 0:y_lb[i]
      log_terms <- stats::dnbinom(k,          size = shapes[i], mu = mu[i],       log = TRUE) +
        stats::dnbinom(y_lb[i] - k,           size = shapex[i], mu = lambdalb[i], log = TRUE)
    }
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }, numeric(1))
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname binegbin_joint
#' @export
#' @keywords internal
log_lik_binegbin_joint <- function(i, prep) {
  mu       <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaem <- brms::get_dpar(prep, "lambdaem", i = i)
  lambdalb <- brms::get_dpar(prep, "lambdalb", i = i)
  shapes   <- brms::get_dpar(prep, "shapes", i = i)
  shapex   <- brms::get_dpar(prep, "shapex", i = i)
  y_em   <- prep$data$Y[i]
  y_lb   <- prep$data$vint1[i]
  em_obs <- prep$data$vint2[i]
  binegbin_joint_lpmf_r(y_em, y_lb, em_obs, mu, lambdaem, lambdalb, shapes, shapex)
}

#' @rdname binegbin_joint
#' @export
#' @keywords internal
posterior_predict_binegbin_joint <- function(i, prep, ...) {
  mu       <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaem <- brms::get_dpar(prep, "lambdaem", i = i)
  lambdalb <- brms::get_dpar(prep, "lambdalb", i = i)
  shapes   <- brms::get_dpar(prep, "shapes", i = i)
  shapex   <- brms::get_dpar(prep, "shapex", i = i)
  y_lb <- prep$data$vint1[i]
  # em_obs is deliberately IGNORED here: every set gets a y_em draw conditional
  # on its observed y_lb, matched and LB-only alike -- the fleet-wide y_em
  # simulation the deliverable needs (05-14 sec-fleetwide). The conditional
  # split N_shared | y_lb is NOT Binomial (a NegBin sum condition is not
  # Binomial); it is P(N_shared = k | y_lb) proportional to
  # NB2(k | mu, shapes) NB2(y_lb - k | lambdalb, shapex) over k = 0..y_lb.
  # Sample that discrete conditional, then add a fresh N10 ~ NB2(lambdaem,
  # shapex). Identical construction to posterior_predict_binegbin.
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
