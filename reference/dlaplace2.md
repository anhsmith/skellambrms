# Discrete-Laplace custom family for brms (free location and scale)

Returns a brms custom family for the discrete Laplace distribution with
free location (`mu`, identity link) and free scale (`sigma`, log link),
obtained by CDF differencing as in
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
but centred at `mu` rather than 0: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`,
where `F` is the CDF of the continuous Laplace(mu, b).

Use in a brm() call as: brm(y ~ ..., family = dlaplace2(), stanvars =
dlaplace2_stanvars(), data = ...)

## Usage

``` r
dlaplace2()

dlaplace2_stanvars()

log_lik_dlaplace2(i, prep)

posterior_predict_dlaplace2(i, prep, ...)

posterior_epred_dlaplace2(prep)
```

## Value

`dlaplace2()` returns a brms `custom_family` object.
`dlaplace2_stanvars()` returns a `stanvars` object holding the Stan code
for `dlaplace2_lpmf`. `log_lik_dlaplace2()` returns a numeric vector of
log-densities, one per posterior draw, for observation `i`.
`posterior_predict_dlaplace2()` returns a vector of simulated
differences, one per posterior draw, for observation `i`, drawn within
the truncation bounds of that observation if it has any.
`posterior_epred_dlaplace2()` returns a draws x observations matrix of
means of the discretised distribution: computed in closed form on rows
without truncation bounds, and by summing the truncated PMF on bounded
rows.

## Parameter named mu

In `dlaplace2()`, the parameter named `mu` is the location, so the
requirement of
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
for a parameter named `"mu"` (see
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md))
is met without reinterpreting it.

## Independence of mu and sigma

[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
requires `sigma^2 >= |mu|`, the relation between the variance and the
mean of a Skellam variable (see
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)).
The discrete Laplace distribution has no such relation, and in
`dlaplace2()` the parameters `mu` and `sigma` vary independently.
Fitting both
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
and `dlaplace2()` to the same data therefore compares a model in which
the bias and the spread of the difference are coupled with a model in
which they are not. A constraint between `mu` and `sigma` in
`dlaplace2()` would remove that contrast.

## Conversion from sigma to b

As in
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md),
`b = sigma / sqrt(2)`, and `sigma` is the SD of the continuous
distribution before discretisation. `double_exponential_lcdf` takes the
location as an argument, as `normal_lcdf` does, so `mu` is passed to it
directly and the Stan code does not shift `z`.

## Mean of the discretised distribution

The mean of the discretised distribution equals `mu` when `mu` is an
integer or a half-integer. For other values of `mu`, the mean differs
from `mu` by up to 0.015 at `sigma = 1`, 0.054 at `sigma = 0.5` and 0.14
at `sigma = 0.25`; the largest differences occur near `mu = 0.3` and
`mu = 0.7` modulo 1. `mu` is therefore the location of the continuous
distribution before discretisation, not the mean of the discretised one.
`posterior_epred_dlaplace2()` returns the mean.

## See also

[`dlaplace2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2_lccdf_stanvars.md)
for truncation;
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
for the fixed-mean version;
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
for the coupled comparison;
[`dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
for the light-tailed alternative.
