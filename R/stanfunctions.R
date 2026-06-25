# Stan function block for the symmetric Skellam log-PMF.
# Injected into the Stan model via stanvars.
skellam1_stan_funs <- "
  real skellam1_lpmf(int k, real mu) {
    return -2 * mu + log_modified_bessel_first_kind(abs(k), 2 * mu);
  }
"

# Stan function block for the symmetric Skellam log-CCDF (truncation
# support via brms's resp_trunc()). The normal-approximation threshold
# is templated in at call time rather than fixed, since the right value
# depends on the data's plausible mu range -- see
# skellam1_lccdf_stanvars() for the rationale and how to choose it.
# The iteration cap (500) and early-exit tolerance are fixed: they guard
# a confirmed std::bad_alloc crash in log_modified_bessel_first_kind and
# a confirmed multi-GB memory blowup at extreme mu, independent of where
# the normal-approximation threshold is set.
skellam1_lccdf_stan <- function(normal_approx_threshold = 100) {
  sprintf("
  real skellam1_lccdf(int y, real mu) {
    // log P(delta > y). Beyond the threshold, skip the exact Bessel-sum
    // tail and use a normal approximation (Skellam(mu,mu) has variance
    // 2*mu, so CLT applies).
    if (mu > %s) {
      real z = (y + 0.5) / sqrt(2 * mu);
      return normal_lccdf(z | 0, 1);
    }
    real acc = negative_infinity();
    int k = y + 1;
    int hard_cap = y + 1 + 500;
    while (k < hard_cap) {
      real lp_k = -2 * mu + log_modified_bessel_first_kind(abs(k), 2 * mu);
      real new_acc = log_sum_exp(acc, lp_k);
      if (lp_k < new_acc - 40 && k > y + 5) {
        acc = new_acc;
        break;
      }
      acc = new_acc;
      k += 1;
    }
    return acc;
  }
", normal_approx_threshold)
}
