# Discrete-normal custom family for brms (location 0, free scale)

Returns a brms custom family for the discrete normal distribution,
location fixed at 0, discretised from the continuous Normal(0, sigma)
via CDF differencing: `P(Z=z) = F(z+0.5) - F(z-0.5)`. One parameter,
sigma (link = "log"), the SD; the mean is always zero. Same
CDF-differencing pattern as
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md),
using Stan's built-in `normal_lcdf`/`normal_lccdf` directly – no Bessel
function and no iteration cap needed, but see the cancellation note
below for a branch this family's PMF does need.

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
observation `i`, drawn subject to that row's
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
bounds where it has any. `posterior_epred_dnorm1()` returns a draws x
observations matrix of means, taken over the truncated distribution on
any row that is bounded.

## Details

**Naming note.** Same forced naming as
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)/[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md):
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
requires a dpar literally named `"mu"`; here it represents sigma (the
SD), not a mean. See
[`?skellam1`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
Details for the full rationale.

**No scale conversion needed.** Unlike
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md),
where Stan's `double_exponential_lcdf` expects the continuous Laplace's
own scale `b` (requiring `b = sigma / sqrt(2)` first), the continuous
normal's own SD parameter *is* sigma directly – `sigma` is passed
straight to `normal_lcdf`/`normal_lccdf` with no intermediate
conversion.

**Cancellation in the PMF, fixed by branching on z's sign.** The naive
`log_diff_exp(normal_lcdf(z+0.5), normal_lcdf(z-0.5))` fails once `z` is
far enough into the positive tail that both `normal_lcdf` calls round to
the same double (both within machine epsilon of `log(1)=0`) – confirmed
to occur at only ~10 SDs out, well inside this package's
realistic-but-stressed test range for the other families, and far sooner
than the analogous direct-subtraction form in
[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
(the normal's thinner tail saturates near 1 much faster per SD than the
Laplace's). Fixed in `dnorm1_lpmf` (and the R-side `log_lik_dnorm1`) by
differencing two *survival* values (`normal_lccdf`, both small and hence
distinguishable) instead of two *CDF* values when `z >= 0` – the same
exact-survival-form idea `dlaplace1_lccdf`/`dlaplace2_lccdf` already
use, applied here to the PMF rather than the CCDF, since CDF
differencing is itself the operation that creates the cancellation risk
in the first place.
