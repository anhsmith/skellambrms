# Stan function block for the symmetric Skellam log-PMF.
# Injected into the Stan model via stanvars.
skellam1_stan_funs <- "
  real skellam1_lpmf(int k, real mu) {
    return -2 * mu + log_modified_bessel_first_kind(abs(k), 2 * mu);
  }
"
