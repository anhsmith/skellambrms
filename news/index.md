# Changelog

## pairedcountbrms 0.6.0

- **Renamed: `skellambrms` is now `pairedcountbrms`.** The old name
  named one family; the package’s subject is the comparison of two
  paired count sources, by two complementary routes — the difference
  families (`skellam1`/`skellam2`, `dnorm1`/`dnorm2`,
  `dlaplace1`/`dlaplace2`) and the joint bivariate families (`bipois`,
  `binegbin`, `binegbin_joint`). Skellam is one member of that set, and
  no longer the most used one.

  **No family names change**, so no fitted model needs refitting:
  `binegbin`, `binegbin_joint`, `bipois`, `skellam1`/`skellam2`,
  `dnorm1`/`dnorm2` and `dlaplace1`/`dlaplace2` all keep their names,
  and brms continues to resolve each fit’s `log_lik_*` /
  `posterior_predict_*` / `posterior_epred_*` methods off the attached
  search path exactly as before. The only change a user needs to make is
  [`library(skellambrms)`](https://rdrr.io/r/base/library.html) →
  [`library(pairedcountbrms)`](https://github.com/anhsmith/pairedcountbrms)
  (and any `skellambrms::` prefix).

  The GitHub repository moves to `anhsmith/pairedcountbrms`. GitHub
  serves a permanent redirect from the old path for both web and git, so
  existing clones and `pak::pak("anhsmith/skellambrms")` calls keep
  working.

- New coordinate helpers,
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)
  and
  [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_dpars_to_mfd.md),
  converting between the rate dpars the joint families take (`mu`,
  `lambdaem`, `lambdalb`, plus `shapes`/`shapex`) and the interpretable
  `(M, f, delta)` coordinates — overall level, congruence, and method
  bias — along with the SD-scale dispersions `kappas`/`kappax`. Pure
  transforms; they fit nothing, and fitting in these coordinates still
  goes through
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  (documented on
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)).
  The inverse reports `delta` as `NA` at `f = 1`, where there is no
  excess to be biased and the bias is genuinely unidentified, rather
  than silently returning `0`.

- New vignette, `Getting started with pairedcountbrms`: simulates paired
  counts from known `binegbin` parameters, fits them with
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) plus
  [`binegbin_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md),
  checks that the five dpars recover the truth, and exercises
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  and
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html).
