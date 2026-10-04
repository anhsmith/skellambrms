# Symmetric Skellam custom family for brms

Returns a brms custom family for the symmetric Skellam distribution,
Skellam(mu_skellam, mu_skellam): the distribution of the difference of
two independent Poisson(mu_skellam) random variables. The mean is zero.
The single parameter is `sigma` (log link), the SD of the difference.
The Stan code derives `mu_skellam = sigma^2 / 2` and evaluates the
Skellam PMF, which is written in terms of a modified Bessel function, at
that value.

Use in a brm() call as: brm(y ~ ..., family = skellam1(), stanvars =
skellam1_stanvars(), data = ...)

## Usage

``` r
skellam1()

skellam1_stanvars()

log_lik_skellam1(i, prep)

posterior_predict_skellam1(i, prep, ...)

posterior_epred_skellam1(prep)
```

## Value

`skellam1()` returns a brms `custom_family` object.
`skellam1_stanvars()` returns a `stanvars` object holding the Stan code
for `skellam1_lpmf`. `log_lik_skellam1()` returns a numeric vector of
log-densities, one per posterior draw, for observation `i`.
`posterior_predict_skellam1()` returns a vector of simulated
differences, one per posterior draw, for observation `i`, drawn within
the truncation bounds of that observation if it has any.
`posterior_epred_skellam1()` returns a draws x observations matrix of
means, taken over the truncated distribution on any row that is bounded.

## Parameter named mu

[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
requires one distributional parameter to be named `"mu"`, whatever that
parameter represents. In `skellam1()`, the parameter named `mu` is
`sigma`, the SD of the difference, on the log link; the mean of
`skellam1()` is zero throughout. In formulas and priors, refer to
`sigma` as `mu`: for example,
`prior(normal(1, 1.5), class = "Intercept")` is a prior on the intercept
of log(sigma). The post-processing functions in this package read the
parameter with `brms::get_dpar(prep, "mu")` and assign it at once to a
variable named `sigma`.

## See also

[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
for truncation;
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
for the free-mean Skellam;
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
and
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
for the same fixed mean with heavier and lighter tails.
