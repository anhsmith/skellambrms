# Package index

## Difference families

Model the difference `d = y1 - y2` directly. The `1` variants fix the
location at zero (do the two sources agree on average?); the `2`
variants estimate it (by how much do they disagree?).

- [`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  [`skellam1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  [`log_lik_skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  [`posterior_predict_skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  [`posterior_epred_skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  : Symmetric Skellam custom family for brms
- [`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  [`skellam2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  [`log_lik_skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  [`posterior_predict_skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  [`posterior_epred_skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  : Asymmetric Skellam custom family for brms
- [`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
  [`dnorm1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
  [`log_lik_dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
  [`posterior_predict_dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
  [`posterior_epred_dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
  : Discrete-normal custom family for brms (location 0, free scale)
- [`dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  [`dnorm2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  [`log_lik_dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  [`posterior_predict_dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  [`posterior_epred_dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  : Discrete-normal custom family for brms (free location and scale)
- [`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  [`dlaplace1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  [`log_lik_dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  [`posterior_predict_dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  [`posterior_epred_dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  : Discrete-Laplace custom family for brms (location 0, free scale)
- [`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
  [`dlaplace2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
  [`log_lik_dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
  [`posterior_predict_dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
  [`posterior_epred_dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
  : Discrete-Laplace custom family for brms (free location and scale)

## Truncation support

Stan `lccdf` functions injected as stanvars so each family normalises
correctly under
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html).

- [`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
  : Log-CCDF of the symmetric Skellam distribution, for truncated fits
- [`skellam2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_lccdf_stanvars.md)
  : Log-CCDF of the asymmetric Skellam distribution, for truncated fits
- [`dnorm1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1_lccdf_stanvars.md)
  : Log-CCDF of the discrete normal distribution with location 0, for
  truncated fits
- [`dnorm2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm2_lccdf_stanvars.md)
  : Log-CCDF of the discrete normal distribution with free location, for
  truncated fits
- [`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md)
  : Log-CCDF of the discrete Laplace distribution with location 0, for
  truncated fits
- [`dlaplace2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2_lccdf_stanvars.md)
  : Log-CCDF of the discrete Laplace distribution with free location,
  for truncated fits

## Parameterisation helpers

Read a fitted family’s dpars under the names its documentation uses.

- [`skellam2_dpars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_dpars.md)
  : Derived quantities of the asymmetric Skellam family from a fitted
  model
