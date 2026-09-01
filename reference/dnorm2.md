# Discrete-normal custom family for brms (free location and scale)

Returns a brms custom family for the discrete normal distribution with
both location (`mu`, link = "identity") and scale (`sigma`, link =
"log") free, discretised via CDF differencing exactly as
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
but centred at `mu` instead of fixed at 0:
`P(Z=z) = F(z+0.5) - F(z-0.5)`, `F` the continuous Normal(mu, sigma)
CDF.

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
observation `i`, drawn subject to that row's
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
bounds where it has any. `posterior_epred_dnorm2()` returns a draws x
observations matrix of means, taken over the truncated distribution on
any row that is bounded.

## Details

**No naming workaround needed.** Unlike
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)/[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)/
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md),
`mu` here genuinely is the family's mean, so brms's "must have a `mu`
parameter" requirement (see
[`?skellam1`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
Details) is satisfied directly – no forced reinterpretation.

**No constraint coupling mu and sigma.** Same structural contrast with
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
already documented for
[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
(see
[`?dlaplace2`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
Details): `mu` and `sigma` are free, independent parameters here, by
design. Fitting
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
against
[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
and `dnorm2()` compares a model where bias and spread are structurally
coupled against ones where they are not; this package supplies all three
families for that comparison.

**Cancellation in the PMF.** Same issue and fix as
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
(see its Details), generalised to branch on whether `z` is on the far
side of `mu` rather than of 0.
