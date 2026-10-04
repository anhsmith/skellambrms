# Discrete-Laplace custom family for brms (location 0, free scale)

Returns a brms custom family for the discrete Laplace distribution with
location fixed at 0, obtained from the continuous Laplace(0, b) by CDF
differencing: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`. The single parameter
is `sigma` (log link), the SD of the continuous distribution before
discretisation (see "Conversion from sigma to b"); the mean is zero.
Unlike the PMF and CCDF of
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
and
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md),
those of `dlaplace1()` are closed-form, built on the Stan function
`double_exponential_lcdf` (Stan names the Laplace distribution "double
exponential"). They require no Bessel function, no large-argument branch
and no iteration cap.

Use in a brm() call as: brm(y ~ ..., family = dlaplace1(), stanvars =
dlaplace1_stanvars(), data = ...)

## Usage

``` r
dlaplace1()

dlaplace1_stanvars()

log_lik_dlaplace1(i, prep)

posterior_predict_dlaplace1(i, prep, ...)

posterior_epred_dlaplace1(prep)
```

## Value

`dlaplace1()` returns a brms `custom_family` object.
`dlaplace1_stanvars()` returns a `stanvars` object holding the Stan code
for `dlaplace1_lpmf`. `log_lik_dlaplace1()` returns a numeric vector of
log-densities, one per posterior draw, for observation `i`.
`posterior_predict_dlaplace1()` returns a vector of simulated
differences, one per posterior draw, for observation `i`, drawn within
the truncation bounds of that observation if it has any.
`posterior_epred_dlaplace1()` returns a draws x observations matrix of
means, taken over the truncated distribution on any row that is bounded.

## Parameter named mu

As in
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md),
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
requires a parameter named `"mu"`, and in `dlaplace1()` that parameter
is `sigma`, the SD, not a mean. See
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
for the consequences for formulas and priors.

## Conversion from sigma to b

`double_exponential_lcdf` takes the scale `b` of the continuous Laplace
distribution. The variance of Laplace(0, b) is `2 * b^2`, so its SD is
`b * sqrt(2)`, and `dlaplace1_lpmf` and `dlaplace1_lccdf` both begin by
computing `b = sigma / sqrt(2)`. The conversion treats `sigma` as the SD
of the continuous distribution. For `sigma >= 0.5`, the SD of the
discretised distribution is larger: by 7.9% at `sigma = 0.5`, 3.4% at
`sigma = 1`, 1.0% at `sigma = 2` and 0.2% at `sigma = 5`. At
`sigma = 0.25`, it is 2.2% smaller, and below that it falls far below
`sigma`, because nearly all of the probability falls on 0: at
`sigma = 0.1`, the discretised SD is 0.029. In
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
and
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md),
by contrast, `sigma` is the exact SD of the integer-valued difference.

## Comparison with extraDistr

:ddlaplace: `extraDistr::ddlaplace()` implements a different discrete
Laplace distribution, with PMF `P(z) = (1 - p) / (1 + p) * p^|z|`, in
which its argument `scale` is the decay probability `p` rather than the
continuous scale `b`. The two PMFs differ: at `b = 3` and
`p = exp(-1/3)`, `P(0) = 0.1535` under `dlaplace1()` and `0.1651` under
`extraDistr::ddlaplace()`. The tests in this package therefore compare
`dlaplace1()` with an R implementation of the CDF difference rather than
with `extraDistr`.

## See also

[`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md)
for truncation;
[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
for the free-mean version;
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
and
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
for the same fixed mean under a Skellam and a lighter-tailed
alternative.
