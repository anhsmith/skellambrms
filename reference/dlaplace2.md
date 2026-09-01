# Discrete-Laplace custom family for brms (free location and scale)

Returns a brms custom family for the discrete Laplace distribution with
both location (`mu`, link = "identity") and scale (`sigma`, link =
"log") free, discretised via CDF differencing exactly as
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
but centred at `mu` instead of fixed at 0:
`P(Z=z) = F(z+0.5) - F(z-0.5)`, `F` the continuous Laplace(mu, b) CDF.

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
differences, one per posterior draw, for observation `i`, drawn subject
to that row's
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
bounds where it has any. `posterior_epred_dlaplace2()` returns a draws x
observations matrix of means, taken over the truncated distribution on
any row that is bounded.

## Details

**No naming workaround needed.** Unlike
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)/[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md),
`mu` here genuinely is the family's mean, so brms's "must have a `mu`
parameter" requirement (see
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
Details) is satisfied directly — no forced reinterpretation.

**No constraint coupling mu and sigma.** This is a genuine structural
difference from
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md),
which structurally requires `sigma >= |mu|` (the Skellam family's actual
mean/variance relationship — see
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
Details). The discrete Laplace has no such relationship: `mu` and
`sigma` are free, independent parameters. Fitting
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
against `dlaplace2()` compares a model where bias and spread are
structurally coupled against one where they are not. This package
supplies both families for that comparison. Do not impose any artificial
coupling here.

**sigma-to-b conversion.** Same as
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md):
`b = sigma / sqrt(2)`. `mu` is passed straight through to
`double_exponential_lcdf`'s own location argument (it takes location and
scale directly, like `normal_lcdf`), so no manual shift of `z` is needed
in the Stan code.

## See also

[`dlaplace2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2_lccdf_stanvars.md)
for truncation;
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
for the fixed-mean version;
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
for the coupled comparison;
[`dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
for the light-tailed alternative.
