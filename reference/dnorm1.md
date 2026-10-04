# Discrete-normal custom family for brms (location 0, free scale)

Returns a brms custom family for the discrete normal distribution with
location fixed at 0, obtained from the continuous Normal(0, sigma) by
CDF differencing: `P(Z = z) = F(z + 0.5) - F(z - 0.5)`. The single
parameter is `sigma` (log link), the SD of the continuous distribution
before discretisation; the mean is zero. The PMF and CCDF are
closed-form: the Stan code computes the PMF with `normal_lcdf` for
`z < 0` and with `erfc()` for `z >= 0`, and the CCDF with `erfc()` (see
"Evaluation in the upper tail"). As in
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md),
they require no Bessel function and no iteration cap.

Use in a brm() call as: brm(y ~ ..., family = dnorm1(), stanvars =
dnorm1_stanvars(), data = ...)

## Usage

``` r
dnorm1()

dnorm1_stanvars()

log_lik_dnorm1(i, prep)

posterior_predict_dnorm1(i, prep, ...)

posterior_epred_dnorm1(prep)
```

## Value

`dnorm1()` returns a brms `custom_family` object. `dnorm1_stanvars()`
returns a `stanvars` object holding the Stan code for `dnorm1_lpmf`.
`log_lik_dnorm1()` returns a numeric vector of log-densities, one per
posterior draw, for observation `i`. `posterior_predict_dnorm1()`
returns a vector of simulated differences, one per posterior draw, for
observation `i`, drawn within the truncation bounds of that observation
if it has any. `posterior_epred_dnorm1()` returns a draws x observations
matrix of means, taken over the truncated distribution on any row that
is bounded.

## Parameter named mu

As in
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
and
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md),
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
requires a parameter named `"mu"`, and in `dnorm1()` that parameter is
`sigma`, the SD, not a mean. See
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
for the consequences for formulas and priors.

## Scale parameter

`sigma` is the SD of the continuous normal distribution and enters the
normal CDF unchanged.
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md),
by contrast, converts `sigma` to the Laplace scale `b` first. The SD of
the discretised distribution is close to `sqrt(sigma^2 + 1/12)` for
`sigma >= 0.5`: larger than `sigma` by 4.1% at `sigma = 1` and by 1.0%
at `sigma = 2`.

## Evaluation in the upper tail

The Stan function `normal_lccdf` returns negative infinity for
standardised arguments above about 8.25 (stan-dev/math#1985).
`dnorm1_lccdf` therefore computes the upper-tail probability as
`0.5 * erfc(x / (sigma * sqrt(2)))`, and `dnorm1_lpmf` uses the same
form for `z >= 0`. In R, `log_lik_dnorm1()` uses
`pnorm(lower.tail = FALSE)` for `z >= 0`, because
`log(pnorm(z + 0.5) - pnorm(z - 0.5))` evaluates to `log(0)` once `z`
exceeds about 8.5 SDs. Both the Stan and the R log-PMF are accurate
until they underflow below about -708, near 38 SDs from 0.

## See also

[`dnorm1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1_lccdf_stanvars.md)
for truncation;
[`dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
for the free-mean version;
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
and
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
for the same fixed mean under a Skellam and a heavier-tailed
alternative.
