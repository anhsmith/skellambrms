# Shared truncation-sampling / truncated-epred machinery for all six custom
# families' posterior_predict_<family>()/posterior_epred_<family>() (see
# family.R). None of the functions here are exported -- they exist purely to
# support the R-side truncated prediction/expectation logic, which brms does
# not provide automatically for custom families (confirmed: brms's own
# rejection-sampling/ntrys machinery and posterior_epred_trunc() correction
# are only ever invoked for *built-in* families; a custom family's own
# posterior_predict_<name>()/posterior_epred_<name>() must implement
# truncation-aware behaviour itself).
#
# Per-family PMF/lccdf math below is NOT re-derived -- each function is a
# direct transcription of the corresponding, already-validated Stan
# lpmf/lccdf in stanfunctions.R (skellam1/2_lpmf/lccdf, dlaplace1/2_lpmf/
# lccdf, dnorm1/2_lpmf/lccdf), duplicated fresh rather than calling
# log_lik_<family>() -- see family.R's log_lik_* functions, which are left
# untouched.

# --------------------------------------------------------------------------
# Small numerical utilities
# --------------------------------------------------------------------------

# erfc(x) = 2 * Phi(-x*sqrt(2)), the standard identity relating the
# complementary error function to the normal CDF. Used by dnorm1/2_lccdf_r
# and the skellam1/2 normal-approximation branch, mirroring the erfc() calls
# already used in the corresponding Stan lccdf code (stanfunctions.R).
.erfc <- function(x) 2 * stats::pnorm(x * sqrt(2), lower.tail = FALSE)

# Vectorized, numerically stable log(exp(a) + exp(b)); handles -Inf inputs
# (log_sum_exp(-Inf, -Inf) = -Inf, not NaN).
.log_sum_exp_pair <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[is.nan(out)] <- -Inf
  out
}

# Extracts the lower ("lb") or upper ("ub") truncation bound for observation
# index/indices `i` from a brms `prep` object, defaulting to -Inf/Inf when
# the field is absent entirely (confirmed: prep$data$lb/ub are NULL, not
# -Inf/Inf-filled, whenever no trunc()/resp_trunc() is used anywhere in the
# formula -- verified directly against brms::make_standata()) or non-finite.
.get_bound <- function(prep, which = c("lb", "ub"), i) {
  which <- match.arg(which)
  default <- if (which == "lb") -Inf else Inf
  b <- prep$data[[which]]
  if (is.null(b)) return(rep(default, length(i)))
  v <- b[i]
  v[!is.finite(v)] <- default
  v
}

# --------------------------------------------------------------------------
# Generic inverse-CDF search (for posterior_predict_<family>)
# --------------------------------------------------------------------------

# Finds the smallest integer y in [lb, ub] (either bound may be -Inf/Inf,
# but not both -- callers only invoke this once at least one is finite) with
# g(y) TRUE, where g is monotone: FALSE for small y, TRUE for large y (as is
# the case for g(y) = logS(y) <= target, since logS is a non-increasing
# survival function). Uses an exponential ("galloping") search from whatever
# bound is finite (doubling the offset each round), then integer bisection
# within the resulting bracket -- this handles lower-only, upper-only, and
# two-sided truncation uniformly, with no special-casing needed: if the
# anchor itself already satisfies (or, symmetrically, never satisfies) g,
# the gallop's own bound-clamping naturally resolves to the anchor.
#
# max_iter bounds each phase (gallop, bisection) separately. This mirrors
# the existing Stan tail-sum's iteration cap philosophy (skellam1_lccdf_stan/
# skellam2_lccdf_stan's 500-iteration cap) -- a safety net against
# pathological inputs, not a tuned parameter: with a doubling offset, 60
# rounds covers a search radius of 2^60, far beyond any realistic parameter
# scale.
.search_monotone <- function(g, lb, ub, max_iter = 60L) {
  anchor <- if (is.finite(lb)) lb else if (is.finite(ub)) ub else 0L

  if (isTRUE(g(anchor))) {
    prev <- anchor; offset <- 1L; lo <- anchor; hi <- anchor
    for (iter in seq_len(max_iter)) {
      cand <- anchor - offset
      if (is.finite(lb)) cand <- max(cand, lb)
      if (!isTRUE(g(cand))) { lo <- cand; hi <- prev; break }
      prev <- cand
      if (is.finite(lb) && cand <= lb) { lo <- cand; hi <- cand; break }
      offset <- offset * 2L
      lo <- cand; hi <- cand   # best effort if max_iter is exhausted below
    }
  } else {
    prev <- anchor; offset <- 1L; lo <- anchor; hi <- anchor
    for (iter in seq_len(max_iter)) {
      cand <- anchor + offset
      if (is.finite(ub)) cand <- min(cand, ub)
      if (isTRUE(g(cand))) { lo <- prev; hi <- cand; break }
      prev <- cand
      if (is.finite(ub) && cand >= ub) { lo <- cand; hi <- cand; break }
      offset <- offset * 2L
      lo <- cand; hi <- cand
    }
  }

  while (hi - lo > 1L) {
    mid <- lo + (hi - lo) %/% 2L
    if (isTRUE(g(mid))) hi <- mid else lo <- mid
  }
  hi
}

# Draws one truncated sample per posterior draw via inverse-CDF sampling
# (not rejection sampling -- rejection sampling was confirmed empirically
# slow/low-acceptance for tight bounds, especially given skellam2's
# expensive per-evaluation Bessel-function cost).
#
# `logS` must be a function(y, idx) returning log P(Y > y) for the draws
# indexed by `idx` (a subset of seq_along(u)), using each draw's own
# parameters (captured in `logS`'s closure, e.g.
# `function(y, idx) dnorm2_lccdf_r(y, mu[idx], sigma[idx])`). `u` is one
# Uniform(0,1) draw per posterior draw. `lb`/`ub` are the (per-observation,
# draw-invariant) truncation bounds, each length 1 or length(u).
#
# Derivation: for Y truncated to [lb, ub], the truncated CDF is
#   G(y) = [S(lb-1) - S(y)] / [S(lb-1) - S(ub)]
# Solving G(y) >= u for the smallest such integer y gives
#   S(y) <= S(lb-1)*(1-u) + S(ub)*u
# i.e. the smallest y with logS(y) <= log_target, where log_target is the
# (numerically stable) log of the right-hand side.
.invert_truncated_cdf <- function(logS, u, lb, ub, max_iter = 60L) {
  n  <- length(u)
  lb <- rep_len(lb, n)
  ub <- rep_len(ub, n)

  log_S_lb1 <- rep(0, n)
  has_lb <- is.finite(lb)
  if (any(has_lb)) {
    idx <- which(has_lb)
    log_S_lb1[idx] <- logS(lb[idx] - 1L, idx)
  }
  log_S_ub <- rep(-Inf, n)
  has_ub <- is.finite(ub)
  if (any(has_ub)) {
    idx <- which(has_ub)
    log_S_ub[idx] <- logS(ub[idx], idx)
  }
  log_target <- .log_sum_exp_pair(log(1 - u) + log_S_lb1, log(u) + log_S_ub)

  out <- integer(n)
  for (d in seq_len(n)) {
    g <- function(y) logS(y, d) <= log_target[d]
    out[d] <- .search_monotone(g, lb[d], ub[d], max_iter = max_iter)
  }
  out
}

# --------------------------------------------------------------------------
# Generic truncated-mean-by-summation (for posterior_epred_<family>)
# --------------------------------------------------------------------------

# No closed form exists for any of the six families' truncated means
# (Bessel-sum and CDF-differenced discrete distributions don't reduce to
# elementary closed forms), so the truncated conditional expectation is
# computed by deterministic two-sided numerical summation of the PMF,
# weighted by y, walking outward from `center` (the untruncated mean,
# clamped into [lb, ub]) until the per-step probability contribution is
# negligible relative to the running total for several consecutive steps,
# capped at max_iter per direction. Deterministic (not Monte Carlo): this
# keeps epred exact to a documented tolerance and reproducible, avoiding
# the seeding/averaging MC would need to keep "differs from the untruncated
# mean when the bound is tight" test assertions reliable. The two
# accumulators (num, den) are summed over the same support, so this is
# self-normalising -- no separate call to a family's lccdf is needed.
.truncated_mean_scalar <- function(lpmf_scalar, lb, ub, center, max_iter, tol) {
  start <- round(center)
  if (is.finite(lb)) start <- max(start, lb)
  if (is.finite(ub)) start <- min(start, ub)

  p0  <- exp(lpmf_scalar(start))
  num <- start * p0
  den <- p0

  y <- start + 1L
  consec_neg <- 0L
  iter <- 0L
  while ((!is.finite(ub) || y <= ub) && iter < max_iter) {
    p <- exp(lpmf_scalar(y))
    num <- num + y * p
    den <- den + p
    consec_neg <- if (p < tol * den) consec_neg + 1L else 0L
    if (consec_neg >= 5L) break
    y <- y + 1L
    iter <- iter + 1L
  }

  y <- start - 1L
  consec_neg <- 0L
  iter <- 0L
  while ((!is.finite(lb) || y >= lb) && iter < max_iter) {
    p <- exp(lpmf_scalar(y))
    num <- num + y * p
    den <- den + p
    consec_neg <- if (p < tol * den) consec_neg + 1L else 0L
    if (consec_neg >= 5L) break
    y <- y - 1L
    iter <- iter + 1L
  }

  num / den
}

# Vectorized (over draws) truncated mean. `lpmf` is function(y, idx)
# returning log P(Y=y) for the draws indexed by `idx`, e.g.
# `function(y, idx) dnorm2_lpmf_r(y, mu[idx], sigma[idx])`. `center` is one
# starting point per draw (length(center) draws total); `lb`/`ub` are
# scalar (draw-invariant, since bounds are per-observation).
.truncated_mean_by_sum <- function(lpmf, lb, ub, center, max_iter = 10000L, tol = 1e-12) {
  n <- length(center)
  vapply(seq_len(n), function(d) {
    .truncated_mean_scalar(function(y) lpmf(y, d), lb, ub, center[d], max_iter, tol)
  }, numeric(1))
}

# --------------------------------------------------------------------------
# skellam1: symmetric Skellam(mu_skellam, mu_skellam), mu_skellam = sigma^2/2
# --------------------------------------------------------------------------

# log-PMF at integer k, mirroring log_lik_skellam1's own formula exactly
# (besselI(..., expon.scaled=TRUE) already bakes in the -2*mu_skellam term).
skellam1_lpmf_r <- function(k, sigma) {
  mu_skellam <- sigma^2 / 2
  log(besselI(2 * mu_skellam, abs(k), expon.scaled = TRUE))
}

# Tail-sum log-CCDF below the normal-approx threshold: log P(K > y) via the
# same iterative Bessel-sum as skellam1_lccdf_stan/skellam2_lccdf_stan (500-
# iteration hard cap, early exit once a term is 40 log-units below the
# running sum). `lpmf_at_k` is a closure over the fixed per-draw parameters,
# called once per candidate integer k >= y+1.
.skellam_tailsum_lccdf <- function(y, lpmf_at_k) {
  acc <- -Inf
  k <- y + 1
  hard_cap <- y + 1 + 500
  repeat {
    lp_k <- lpmf_at_k(k)
    new_acc <- .log_sum_exp_pair(acc, lp_k)
    if (lp_k < new_acc - 40 && k > y + 5) { acc <- new_acc; break }
    acc <- new_acc
    k <- k + 1
    if (k >= hard_cap) break
  }
  acc
}

# log-CCDF (survival function), mirroring skellam1_lccdf_stan exactly: same
# normal-approximation branch above `normal_approx_threshold` (on the
# mu_skellam scale), same exact Bessel tail-sum below it.
skellam1_lccdf_r <- function(y, sigma, normal_approx_threshold = 100) {
  mu_skellam <- sigma^2 / 2
  n <- max(length(y), length(mu_skellam))
  y <- rep_len(y, n)
  mu_skellam <- rep_len(mu_skellam, n)
  out <- numeric(n)

  approx <- mu_skellam > normal_approx_threshold
  if (any(approx)) {
    z <- (y[approx] + 0.5) / sqrt(2 * mu_skellam[approx])
    out[approx] <- log(0.5) + log(.erfc(z / sqrt(2)))
  }
  if (any(!approx)) {
    idx <- which(!approx)
    out[idx] <- mapply(
      function(yy, mm) .skellam_tailsum_lccdf(yy, function(k) log(besselI(2 * mm, abs(k), expon.scaled = TRUE))),
      y[idx], mu_skellam[idx]
    )
  }
  out
}

# --------------------------------------------------------------------------
# skellam2: asymmetric Skellam(theta1, theta2)
# --------------------------------------------------------------------------

# log-PMF at integer k, mirroring log_lik_skellam2's own formula exactly.
skellam2_lpmf_r <- function(k, mu, sigmaexcess) {
  sigmasq <- abs(mu) + sigmaexcess^2
  theta1  <- (sigmasq + mu) / 2
  theta2  <- (sigmasq - mu) / 2
  z <- 2 * sqrt(theta1 * theta2)
  log(besselI(z, abs(k), expon.scaled = TRUE)) + z - theta1 - theta2 + (k / 2) * log(theta1 / theta2)
}

# log-CCDF, mirroring skellam2_lccdf_stan exactly: same threshold check on
# mu_skellam = (theta1+theta2)/2, same normal-approx branch (mean mu, sd
# sqrt(sigmasq)), same exact Bessel tail-sum below threshold.
skellam2_lccdf_r <- function(y, mu, sigmaexcess, normal_approx_threshold = 100) {
  n <- max(length(y), length(mu), length(sigmaexcess))
  y <- rep_len(y, n)
  mu <- rep_len(mu, n)
  sigmaexcess <- rep_len(sigmaexcess, n)

  sigmasq <- abs(mu) + sigmaexcess^2
  theta1  <- (sigmasq + mu) / 2
  theta2  <- (sigmasq - mu) / 2
  mu_skellam <- (theta1 + theta2) / 2
  out <- numeric(n)

  approx <- mu_skellam > normal_approx_threshold
  if (any(approx)) {
    sigma_a <- sqrt(sigmasq[approx])
    z <- (y[approx] + 0.5 - mu[approx]) / sigma_a
    out[approx] <- log(0.5) + log(.erfc(z / sqrt(2)))
  }
  if (any(!approx)) {
    idx <- which(!approx)
    out[idx] <- mapply(
      function(yy, t1, t2) {
        zz <- 2 * sqrt(t1 * t2)
        .skellam_tailsum_lccdf(yy, function(k) {
          log(besselI(zz, abs(k), expon.scaled = TRUE)) + zz - t1 - t2 + (k / 2) * log(t1 / t2)
        })
      },
      y[idx], theta1[idx], theta2[idx]
    )
  }
  out
}

# --------------------------------------------------------------------------
# dlaplace1 / dlaplace2: discrete Laplace (location 0 / free location)
# --------------------------------------------------------------------------

# log-PMF, mirroring log_lik_dlaplace1/dlaplace2's own CDF-differencing
# formula exactly (same laplace_cdf two-branch closed form).
dlaplace1_lpmf_r <- function(z, sigma) {
  b <- sigma / sqrt(2)
  # Broadcast x against b before the ifelse() test -- see log_lik_dlaplace1
  # for why (test-length, not yes/no-branch-length, drives ifelse() output).
  laplace_cdf <- function(x) {
    x <- x + 0 * b
    ifelse(x < 0, 0.5 * exp(x / b), 1 - 0.5 * exp(-x / b))
  }
  log(laplace_cdf(z + 0.5) - laplace_cdf(z - 0.5))
}
dlaplace2_lpmf_r <- function(z, mu, sigma) {
  b <- sigma / sqrt(2)
  laplace_cdf <- function(x) ifelse(x < 0, 0.5 * exp(x / b), 1 - 0.5 * exp(-x / b))
  log(laplace_cdf(z - mu + 0.5) - laplace_cdf(z - mu - 0.5))
}

# log-CCDF: log(0.5) - x/b for x = y(-mu)+0.5 >= 0 (right tail, direct);
# log1p(-0.5*exp(x/b)) for x < 0 (left tail -- NOT the same single-branch
# formula as the right tail: evaluating log(0.5) - x/b at negative x gives
# log(0.5) + |x|/b, i.e. a log-survival greater than log(1) = 0, which is
# nonsensical -- this is exactly the range that matters for S(lb-1) under a
# lower-truncation bound, so the branch is not optional). Both branches
# mirror dlaplace2_lccdf_stan's existing branch-free Stan form
# (log1m_exp(double_exponential_lcdf(...))), expanded here into the
# equivalent explicit closed form for each sign of x.
dlaplace1_lccdf_r <- function(y, sigma) {
  b <- sigma / sqrt(2)
  x <- y + 0.5
  # Broadcast x against b before the ifelse() test -- see log_lik_dlaplace1
  # for why (test-length, not yes/no-branch-length, drives ifelse() output).
  x <- x + 0 * b
  ifelse(x >= 0, log(0.5) - x / b, log1p(-0.5 * exp(x / b)))
}
dlaplace2_lccdf_r <- function(y, mu, sigma) {
  b <- sigma / sqrt(2)
  x <- y - mu + 0.5
  x <- x + 0 * b
  ifelse(x >= 0, log(0.5) - x / b, log1p(-0.5 * exp(x / b)))
}

# --------------------------------------------------------------------------
# dnorm1 / dnorm2: discrete normal (location 0 / free location)
# --------------------------------------------------------------------------

# log-PMF, mirroring log_lik_dnorm1/dnorm2's own z>=mean branch exactly
# (survival-difference form to avoid the log-CDF cancellation risk in the
# positive tail).
dnorm1_lpmf_r <- function(z, sigma) {
  # Broadcast z against sigma before the ifelse() test -- see
  # log_lik_dlaplace1 for why (test-length, not yes/no-branch-length,
  # drives ifelse() output).
  z <- z + 0 * sigma
  ifelse(
    z >= 0,
    log(stats::pnorm(z - 0.5, sd = sigma, lower.tail = FALSE) -
        stats::pnorm(z + 0.5, sd = sigma, lower.tail = FALSE)),
    log(stats::pnorm(z + 0.5, sd = sigma) - stats::pnorm(z - 0.5, sd = sigma))
  )
}
dnorm2_lpmf_r <- function(z, mu, sigma) {
  z <- z + 0 * mu + 0 * sigma
  ifelse(
    z >= mu,
    log(stats::pnorm(z - mu - 0.5, sd = sigma, lower.tail = FALSE) -
        stats::pnorm(z - mu + 0.5, sd = sigma, lower.tail = FALSE)),
    log(stats::pnorm(z - mu + 0.5, sd = sigma) - stats::pnorm(z - mu - 0.5, sd = sigma))
  )
}

# log-CCDF via erfc(), mirroring dnorm1_lccdf_stan/dnorm2_lccdf_stan exactly
# -- single branch, valid at all y (erfc() is well-conditioned in both
# tails, unlike Stan's normal_lccdf; no branch needed here, unlike the
# discrete-Laplace case above).
dnorm1_lccdf_r <- function(y, sigma) {
  z <- (y + 0.5) / (sigma * sqrt(2))
  log(0.5) + log(.erfc(z))
}
dnorm2_lccdf_r <- function(y, mu, sigma) {
  z <- (y + 0.5 - mu) / (sigma * sqrt(2))
  log(0.5) + log(.erfc(z))
}
