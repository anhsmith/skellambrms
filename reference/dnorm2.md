# Discrete-normal custom family for brms (free location and scale)

Returns a brms custom family for the discrete normal distribution with
free location (`mu`, identity link) and free scale (`sigma`, log link),
obtained by CDF differencing as in
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
but centred at `mu` rather than 0: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`,
where `F` is the CDF of the continuous Normal(mu, sigma).

Use in a brm() call as: brm(y ~ ..., family = dnorm2(), stanvars =
dnorm2_stanvars(), data = ...)

## Usage

``` r
dnorm2()

dnorm2_stanvars()

log_lik_dnorm2(i, prep)

posterior_predict_dnorm2(i, prep, ...)

posterior_epred_dnorm2(prep)
```

## Value

`dnorm2()` returns a brms `custom_family` object. `dnorm2_stanvars()`
returns a `stanvars` object holding the Stan code for `dnorm2_lpmf`.
`log_lik_dnorm2()` returns a numeric vector of log-densities, one per
posterior draw, for observation `i`. `posterior_predict_dnorm2()`
returns a vector of simulated differences, one per posterior draw, for
observation `i`, drawn within the truncation bounds of that observation
if it has any. `posterior_epred_dnorm2()` returns a draws x observations
matrix of means of the discretised distribution: computed from a rapidly
converging series on rows without truncation bounds, and by summing the
truncated PMF on bounded rows.

## Parameter named mu

In `dnorm2()`, the parameter named `mu` is the location, so the
requirement of
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
for a parameter named `"mu"` (see
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md))
is met without reinterpreting it.

## Independence of mu and sigma

As in
[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md),
`mu` and `sigma` vary independently, with no constraint between them.
Fitting
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md),
[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
and `dnorm2()` to the same data compares a model in which the bias and
the spread of the difference are coupled with two models in which they
are not.

## Evaluation in the upper tail

As in
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md),
the PMF and CCDF compute the upper tail with `erfc()`, with the PMF
branching on whether `z` lies above or below `mu` rather than 0.

## Mean of the discretised distribution

The mean of the discretised distribution equals `mu` when `mu` is an
integer or a half-integer. For other values of `mu`, the mean differs
from `mu` by less than 1e-9 at `sigma = 1` and by less than 5e-6 at
`sigma = 0.75`, but by up to 0.0023 at `sigma = 0.5` and 0.093 at
`sigma = 0.25`. `mu` is therefore the location of the continuous
distribution before discretisation, not the mean of the discretised one.
`posterior_epred_dnorm2()` returns the mean.

## See also

[`dnorm2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm2_lccdf_stanvars.md)
for truncation;
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
for the fixed-mean version;
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
for the coupled comparison;
[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
for the heavy-tailed alternative.
