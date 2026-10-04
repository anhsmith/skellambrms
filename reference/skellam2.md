# Asymmetric Skellam custom family for brms

Returns a brms custom family for the asymmetric Skellam distribution,
Skellam(theta1, theta2): the distribution of the difference of two
independent Poisson random variables with rates `theta1` and `theta2`,
which may differ. The two parameters are `mu` (identity link), the mean
of the difference, and `sigmaexcess` (log link, so non-negative). The SD
of the difference, `sigma`, and the rates `theta1` and `theta2` are
derived from `mu` and `sigmaexcess` (see Details).

Use in a brm() call as: brm(y ~ ..., family = skellam2(), stanvars =
skellam2_stanvars(), data = ...)

## Usage

``` r
skellam2()

skellam2_stanvars()

log_lik_skellam2(i, prep)

posterior_predict_skellam2(i, prep, ...)

posterior_epred_skellam2(prep)
```

## Value

`skellam2()` returns a brms `custom_family` object.
`skellam2_stanvars()` returns a `stanvars` object holding the Stan code
for `skellam2_lpmf`. `log_lik_skellam2()` returns a numeric vector of
log-densities, one per posterior draw, for observation `i`.
`posterior_predict_skellam2()` returns a vector of simulated
differences, one per posterior draw, for observation `i`, drawn within
the truncation bounds of that observation if it has any.
`posterior_epred_skellam2()` returns a draws x observations matrix of
means, taken over the truncated distribution on any row that is bounded.

## Parameter named sigmaexcess

[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
does not allow underscores or dots in the names of distributional
parameters, so the second parameter is named `sigmaexcess` rather than
`sigma_excess`.

## Constraint between mu and sigma

The Skellam rates satisfy `theta1 + theta2 = sigma^2` and
`theta1 - theta2 = mu`, so `theta1 = (sigma^2 + mu) / 2` and
`theta2 = (sigma^2 - mu) / 2`. Both rates are non-negative only if
`sigma^2 >= |mu|`: the variance of a Skellam variable is at least the
absolute value of its mean. `skellam2()` satisfies this constraint by
construction, through sigma^2 = \|mu\| + sigmaexcess^2, for every `mu`
and every `sigmaexcess >= 0`, with equality at `sigmaexcess = 0`. For
`mu >= 0`, `theta1 = mu + sigmaexcess^2 / 2` and
`theta2 = sigmaexcess^2 / 2`. For `mu < 0`, `theta1 = sigmaexcess^2 / 2`
and `theta2 = |mu| + sigmaexcess^2 / 2`. Both rates are strictly
positive whenever `sigmaexcess > 0`, which the log link guarantees for
any finite linear predictor. At `mu = 0`, `skellam2()` reduces to
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md),
with `sigma = sigmaexcess` and `theta1 = theta2 = sigmaexcess^2 / 2`.

The alternative construction `sigma = sqrt(mu^2 + sigmaexcess^2)`
guarantees only `sigma >= |mu|`, which does not imply `sigma^2 >= |mu|`
when `|mu| < 1`. It can then give a negative rate: `mu = 0.5` and
`sigmaexcess = 0` give `sigma = 0.5` and `theta2 = -0.125`.

## Derived quantities

The Stan program that brms generates declares the per-observation
vectors of distributional parameters, here `mu` and `sigmaexcess`, as
local variables in the `model` block rather than as transformed
parameters. A `generated quantities` block therefore cannot refer to
them, with or without `loop = TRUE`.
[`skellam2_dpars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_dpars.md)
instead computes `mu`, `sigma`, `sigma^2`, `theta1` and `theta2` in R
with
[`brms::get_dpar()`](https://paulbuerkner.com/brms/reference/get_dpar.html),
which works for any formula.

## See also

[`skellam2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_lccdf_stanvars.md)
for truncation;
[`skellam2_dpars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_dpars.md)
to recover `sigma`, `theta1` and `theta2` from a fit;
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
for the fixed-mean family this one reduces to at `mu = 0`;
[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
and
[`dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
for free-mean families that leave mean and spread uncoupled.
