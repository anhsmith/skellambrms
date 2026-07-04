# skellambrms

Custom [brms](https://paul-buerkner.github.io/brms/) families for comparing
two paired counts. Two complementary approaches: **difference families** that
model the integer-valued difference on ℤ (one paired count minus the other),
and **joint bivariate-count families** that model the matched pair itself —
its correlation, marginal overdispersion, and difference together.

## Background

Many count-comparison problems produce a response that is an integer, but
can be negative: the difference between two observers' counts of the same
event, the goal difference between two teams, or any other set-level
difference between two paired count sources. Standard `brms` count
families (Poisson, negative binomial, ...) don't support a `Z`-valued
response.

This package provides six **difference** families, built around three
distributional families, plus two **joint bivariate-count** families
(documented in their own section below):

- **Skellam** — the distribution of the difference of two independent
  Poisson random variables.
- **Discrete Laplace** — a Laplace distribution discretised onto the
  integers via CDF differencing.
- **Discrete normal** — a normal distribution discretised the same way.

The joint families take a complementary view: rather than collapsing the
pair to its difference (which discards the pair's overall level and
correlation), they model both counts at once. Use a difference family when
only the disagreement matters and you want truncation support; use a joint
family when the pair's correlation, marginal overdispersion, or the need to
simulate one count conditional on the other matters.

Each of the three is available in two parameterisations: one with the mean
fixed at zero (for testing whether two sources agree on average) and one
with a free mean (for quantifying systematic disagreement), giving six
families in total — `skellam1`/`skellam2`, `dlaplace1`/`dlaplace2`, and
`dnorm1`/`dnorm2`. All six are sampled on a common (mean, SD-scale) spread
convention, and all support truncation via `brms`'s `resp_trunc()`.

## Families

| Family | Mean | Spread parameter | Mean–spread coupling |
|---|---|---|---|
| `skellam1()` | fixed at 0 | `sigma` (log link) | none exposed to the user; internally `theta1 = theta2 = sigma^2/2` |
| `skellam2()` | free, `mu` (identity link) | `sigmaexcess` (log link) | `sigma^2 = |mu| + sigmaexcess^2` — the genuine Skellam validity constraint (`theta1, theta2 >= 0`); reduces exactly to `skellam1()` at `mu = 0` |
| `dlaplace1()` | fixed at 0 | `sigma` (log link) | none |
| `dlaplace2()` | free, `mu` (identity link) | `sigma` (log link) | none — deliberately decoupled from `mu` |
| `dnorm1()` | fixed at 0 | `sigma` (log link) | none |
| `dnorm2()` | free, `mu` (identity link) | `sigma` (log link) | none — deliberately decoupled from `mu` |

`skellam1`, `dlaplace1`, and `dnorm1` fix the mean at zero and estimate a
single spread parameter, `sigma` — the SD of the response — on the log
scale. `skellam2`, `dlaplace2`, and `dnorm2` additionally estimate a free
mean `mu` on the identity scale.

For `skellam2()`, `sigma` and `mu` are not free: `sigma^2 = |mu| +
sigmaexcess^2` guarantees the underlying Poisson rates `theta1 =
(sigma^2 + mu)/2` and `theta2 = (sigma^2 - mu)/2` are both non-negative for
every `mu` and every `sigmaexcess >= 0` — the actual Skellam constraint
(variance ≥ |mean|). This is a genuine structural coupling: a Skellam
difference with a large mean *must* also have large variance. The discrete
Laplace and discrete normal families have no such constraint — `mu` and
`sigma` are free, independent parameters throughout. Comparing `skellam2()`
against `dlaplace2()`/`dnorm2()` lets you test whether bias and spread are
structurally coupled in your data or vary independently.

## Joint bivariate-count families

`bipois()` and `binegbin()` model a matched count pair `(y_em, y_lb)`
jointly, rather than its difference. Both use the same trivariate-reduction
construction: a shared latent count plus two private latent counts,

```
y_em = N_shared + N10       y_lb = N_shared + N01
```

with `N_shared`, `N10`, `N01` mutually independent given their rates and
`N_shared` marginalised out of the joint likelihood analytically. The shared
component induces positive correlation between the two counts; the private
components govern their difference (`N_shared` cancels from `y_em - y_lb`).

| Family | Latent components | Dispersion | Captures |
|---|---|---|---|
| `bipois()` | three independent **Poisson** (`mu`, `lambdaem`, `lambdalb`) | none — each component has `Var == mean` | correlation and difference of the pair, but only when the margins are *not* overdispersed |
| `binegbin()` | three independent **Negative-Binomial** | scalar `shapes` (shared) and `shapex` (private) | as above, plus marginal overdispersion and the associated over-spread of the difference |

`binegbin()` is the one to reach for on real count data, which is almost
always overdispersed: `bipois()` forces `Var == mean` on every component and
so underfits the marginal variances (and hence the difference variance) of
overdispersed pairs. `binegbin()` carries the extra dispersion in two
identifiable *scalar* parameters — a Negative-Binomial `shapes` for the
shared component and a `shapex` shared across the two private components.
(An observation-level random effect on the private components was tried as an
alternative and rejected: with one pair per unit it overfits and the
dispersion SD collapses, so a fresh-draw posterior-predictive check fails to
reproduce the difference variance even though a conditional one looks fine.)

The second count travels via brms's `vint()` addition term, since
`custom_family()` declares a single response column. There is no forced-`mu`
naming quirk here — `mu` is genuinely the shared component's rate — but note
that `mu` is *not* the mean of either response individually.

```r
library(brms)
library(skellambrms)

# bipois(): joint bivariate Poisson (non-overdispersed margins)
fit_bp <- brm(
  bf(y_em | vint(y_lb) ~ 1,
     mu ~ 1 + (1 | vessel), lambdaem ~ 1, lambdalb ~ 1),
  data     = dat,
  family   = bipois(),
  stanvars = bipois_stanvars(),
  chains   = 4
)

# binegbin(): joint bivariate Negative-Binomial (overdispersed margins)
fit_nb <- brm(
  bf(y_em | vint(y_lb) ~ 1,
     mu ~ 1 + (1 | vessel),
     nlf(lambdaem ~ exp(lamx)), nlf(lambdalb ~ exp(lamx)), lamx ~ 1,
     shapes ~ 1, shapex ~ 1, nl = TRUE),
  data     = dat,
  family   = binegbin(),
  stanvars = binegbin_stanvars(),
  chains   = 4
)
```

The `nlf(... exp(lamx))` idiom above ties the two private rates to a shared
value (a "no systematic bias" assumption, `E[y_em] = E[y_lb]`); drop it and
give `lambdaem`/`lambdalb` free intercepts to allow a mean difference. For
either family, `posterior_predict()` simulates `y_em` conditional on the
observed `y_lb` — for `bipois()` via the closed-form
`Binomial(y_lb, mu/(mu + lambdalb))` split, for `binegbin()` via the discrete
`N_shared | y_lb` conditional. Truncation (`resp_trunc()`) does not apply to
these joint families.

### A naming quirk to be aware of

`brms::custom_family()` unconditionally requires one `dpars` entry to be
literally named `"mu"`. For `skellam1()`, `dlaplace1()`, and `dnorm1()`,
that forced `"mu"` dpar is actually `sigma` — the mean is structurally
zero and isn't represented as a dpar at all. If you inspect
`brms::make_stancode()` output or call `brms::get_dpar(prep, "mu")` for one
of these three families, you're looking at `sigma`, not a mean. All R-side
functions in this package immediately rebind that dpar to a variable
called `sigma`, so nothing else in the package (or in this README) ever
calls it `mu`. `skellam2()`, `dlaplace2()`, and `dnorm2()` don't have this
issue — their `mu` genuinely is the mean.

Separately, `skellam2()`'s excess-spread parameter is spelled
`sigmaexcess`, not `sigma_excess`: `brms::custom_family()` disallows dots
and underscores in `dpars` names.

## Installation

```r
# install.packages("pak")
pak::pak("anhsmith/skellambrms")
```

Stan and a C++ toolchain are required. On Windows, install
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
Works with either rstan or cmdstanr as the brms backend.

## Usage

Every family follows the same pattern: pass `family = <family>()` and
`stanvars = <family>_stanvars()` to `brm()`. Add `<family>_lccdf_stanvars()`
to `stanvars` (combined via `+`) to also support truncation through
`resp_trunc()`.

### Skellam

```r
library(brms)
library(skellambrms)

# skellam1(): mean fixed at 0 -- do two sources agree on average?
fit1 <- brm(
  bf(delta | trunc(lb = neg_bound) ~ 1 + (1 | group)),
  data     = dat,
  family   = skellam1(),
  stanvars = skellam1_stanvars() + skellam1_lccdf_stanvars(),
  chains   = 4
)

# skellam2(): free mean -- how large and how uncertain is the disagreement?
fit2 <- brm(
  bf(delta | trunc(lb = neg_bound) ~ 1 + x, sigmaexcess ~ 1),
  data     = dat,
  family   = skellam2(),
  stanvars = skellam2_stanvars() + skellam2_lccdf_stanvars(),
  chains   = 4
)
```

### Discrete Laplace

```r
# dlaplace1(): mean fixed at 0
fit3 <- brm(
  bf(delta | trunc(lb = neg_bound) ~ 1 + (1 | group)),
  data     = dat,
  family   = dlaplace1(),
  stanvars = dlaplace1_stanvars() + dlaplace1_lccdf_stanvars(),
  chains   = 4
)

# dlaplace2(): free mean and free scale, structurally uncoupled
fit4 <- brm(
  bf(delta | trunc(lb = neg_bound) ~ 1 + x, sigma ~ 1),
  data     = dat,
  family   = dlaplace2(),
  stanvars = dlaplace2_stanvars() + dlaplace2_lccdf_stanvars(),
  chains   = 4
)
```

### Discrete normal

```r
# dnorm1(): mean fixed at 0
fit5 <- brm(
  bf(delta | trunc(lb = neg_bound) ~ 1 + (1 | group)),
  data     = dat,
  family   = dnorm1(),
  stanvars = dnorm1_stanvars() + dnorm1_lccdf_stanvars(),
  chains   = 4
)

# dnorm2(): free mean and free scale, structurally uncoupled
fit6 <- brm(
  bf(delta | trunc(lb = neg_bound) ~ 1 + x, sigma ~ 1),
  data     = dat,
  family   = dnorm2(),
  stanvars = dnorm2_stanvars() + dnorm2_lccdf_stanvars(),
  chains   = 4
)
```

`neg_bound` in these examples is a column giving a (possibly row-varying)
lower truncation bound, e.g. `-y_lb` where `y_lb` is how far below zero the
response could plausibly have gone for that row. All families support
arbitrary `brms` formula syntax, including random effects, non-linear
predictors on the spread dpars, and (for the free-mean families) on `mu`
too.

## Truncation

Each family exports a `<family>_lccdf_stanvars()` function that adds a
`brms::stanvar()` defining a Stan function named `<family>_lccdf` — the log
complementary CDF, `log P(Z > y)`. `brms`'s `resp_trunc()` machinery finds
this function purely by name convention (`<family>_lccdf`, matching the
family name) and uses it to compute the log-normalisation constant for the
truncated likelihood, including a row-varying lower bound. No further
wiring is required beyond adding the stanvar to your `stanvars` argument,
as shown above.

For `skellam1()` and `skellam2()`, the exact log-CCDF is an iterative
tail-sum over the Bessel-function PMF. Above a configurable
`normal_approx_threshold` (default `100`, measured on the underlying
`mu_skellam` scale — `sigma^2/2` for `skellam1()`, `(theta1 + theta2)/2`
for `skellam2()`), this exact sum is replaced by a normal approximation.
This guards against two confirmed numerical failure modes that can occur
during HMC warmup, when an unadapted proposal pushes the spread parameter
to an extreme value (the log link places no ceiling on it): a crash from
evaluating the Bessel function at an enormous order, and a much slower
blow-up in cost and memory when many rows within a single deep NUTS tree
each hit the expensive exact loop. Read `?skellam1_lccdf_stanvars` and the
comments in `R/family.R` for the full rationale and guidance on choosing a
threshold appropriate to your own data's scale.

`dlaplace1()`/`dlaplace2()` and `dnorm1()`/`dnorm2()` have closed-form
log-CCDFs (built on `double_exponential_lcdf` for the Laplace families and
an exact `erfc()`-based survival function for the normal families), so
their `_lccdf_stanvars()` functions take no threshold argument — there is
no large-argument failure mode to guard against.

`posterior_predict_<family>()` and `posterior_epred_<family>()` correctly
account for `resp_trunc()` bounds too (as of `skellambrms` 0.3.1) — see
Limitations below for an important caveat about `brms::posterior_epred()`
specifically on truncated fits.

## Limitations

**`brms::posterior_epred()` — and anything built on it, including
`fitted()` and `conditional_effects()` — errors on any truncated fit, for
all six families here.** This is a `brms` limitation, not this package's:
`brms`'s internal `posterior_epred.brmsprep()` checks whether the model is
truncated *before* checking whether the family is a custom one, and
unconditionally routes truncated fits to `brms:::posterior_epred_trunc()`.
That function has no fallback to a custom family's own
`posterior_epred_<family>()` on the truncated branch — it looks for a
generic `posterior_epred_trunc_custom()` inside `brms`'s own namespace,
which doesn't exist for *any* custom family, truncated or not. Confirmed
directly against the installed `brms` source
(`brms:::posterior_epred_trunc`); the error raised is `"posterior_epred
values on the respone scale not yet implemented for truncated 'custom'
models"`.

**Workaround:** call the family's `posterior_epred_<family>()` directly on
a real `prepare_predictions()` object, bypassing `brms`'s generic dispatch:

```r
prep  <- brms::prepare_predictions(fit)
epred <- posterior_epred_dnorm2(prep)  # or posterior_epred_skellam1(), etc.
```

Each family's `posterior_epred_<family>()` correctly accounts for
`resp_trunc()` bounds when called this way — the underlying computation
isn't the problem, only `brms`'s own generic dispatch is.

`brms::posterior_predict()` is unaffected by this and works correctly for
truncated fits of every family — its dispatch calls the family's own
`posterior_predict_<family>()` unconditionally and only checks truncation
bounds afterwards, with no analogous gate.

## Testing

The test suite (`tests/testthat/`) validates, for every family:

- The R-side log-PMF and log-CCDF against a trusted external reference —
  `skellam::dskellam()`/`pskellam()` for the Skellam families, and a
  hand-derived, numerically stable log-space CDF-differencing reference for
  the discrete Laplace and discrete normal families (`extraDistr::ddlaplace()`
  was checked and found to implement a *different* discrete Laplace, so it
  isn't usable as a reference).
- That the PMF sums to 1 across a grid of parameter values.
- Numerical stability (no `NaN`/`Inf`) across a "realistic-but-stressed"
  range of parameter values, including deep into the tails.
- That the Stan implementations of the log-PMF and log-CCDF agree with the
  R-side references, once exposed via `rstan::expose_stan_functions()`.
- For `skellam1()`/`skellam2()`, that the exact and normal-approximation
  branches of the log-CCDF agree closely at the threshold seam, and that
  changing `normal_approx_threshold` actually moves the cutover point.
- Structural (not rejection-based) enforcement of validity constraints, by
  inspecting `brms::make_stancode()` output directly: `skellam2()`'s
  `theta1, theta2 >= 0` constraint and `dlaplace2()`/`dnorm2()`'s lack of
  any constraint coupling `mu` and `sigma`.
- That the free-mean families (`skellam2()`, `dlaplace2()`, `dnorm2()`)
  reduce exactly to their fixed-mean counterparts at `mu = 0`.
- End-to-end parameter recovery: fitting each family with `brm()` (using
  `cmdstanr`) to simulated data — including truncated data via
  `resp_trunc()` — and checking that the true generating parameters fall
  within the posterior credible interval, alongside divergence and Rhat
  checks.
- That `log_lik_<family>()` and the internal R-side `_lpmf_r()`/`_lccdf_r()`
  helpers return one value per posterior draw for a single observation, not
  just for a vector of observations — the direction that silently broke
  `log_lik_dlaplace1()` (and hence `loo()`) prior to 0.3.2, since R's
  `ifelse()` takes its output length from its test argument alone.

For the joint families (`bipois()`, `binegbin()`) the suite instead validates
the marginalised joint log-PMF against an independent R brute-force reference
across a rate/shape grid and at extreme edge cases, checks that it normalises
to 1, checks the analytic moment identities (mean, marginal variance,
difference variance, covariance), confirms `binegbin()` reduces to `bipois()`
in the Poisson limit (`shapes`, `shapex` → ∞), and runs end-to-end parameter
recovery from simulated hierarchical data with divergence/Rhat checks.

## Functions

| Function | Purpose |
|---|---|
| `skellam1()` | Custom family object for `skellam1` (mean fixed at 0) |
| `skellam1_stanvars()` | Stan code block for `skellam1` |
| `skellam1_lccdf_stanvars()` | Stan log-CCDF block enabling `resp_trunc()` support for `skellam1` |
| `skellam2()` | Custom family object for `skellam2` (free mean) |
| `skellam2_stanvars()` | Stan code block for `skellam2` |
| `skellam2_lccdf_stanvars()` | Stan log-CCDF block enabling `resp_trunc()` support for `skellam2` |
| `skellam2_dpars()` | Reports `mu`, `sigma`, `sigmasq`, `theta1`, `theta2` from a fitted `skellam2()` model, computed in R via `get_dpar()` |
| `dlaplace1()` | Custom family object for `dlaplace1` (mean fixed at 0) |
| `dlaplace1_stanvars()` | Stan code block for `dlaplace1` |
| `dlaplace1_lccdf_stanvars()` | Stan log-CCDF block enabling `resp_trunc()` support for `dlaplace1` |
| `dlaplace2()` | Custom family object for `dlaplace2` (free mean and scale) |
| `dlaplace2_stanvars()` | Stan code block for `dlaplace2` |
| `dlaplace2_lccdf_stanvars()` | Stan log-CCDF block enabling `resp_trunc()` support for `dlaplace2` |
| `dnorm1()` | Custom family object for `dnorm1` (mean fixed at 0) |
| `dnorm1_stanvars()` | Stan code block for `dnorm1` |
| `dnorm1_lccdf_stanvars()` | Stan log-CCDF block enabling `resp_trunc()` support for `dnorm1` |
| `dnorm2()` | Custom family object for `dnorm2` (free mean and scale) |
| `dnorm2_stanvars()` | Stan code block for `dnorm2` |
| `dnorm2_lccdf_stanvars()` | Stan log-CCDF block enabling `resp_trunc()` support for `dnorm2` |
| `bipois()` | Custom family object for the joint bivariate Poisson |
| `bipois_stanvars()` | Stan code block for `bipois` |
| `binegbin()` | Custom family object for the joint bivariate Negative-Binomial |
| `binegbin_stanvars()` | Stan code block for `binegbin` |

Each family also exports `log_lik_<family>()`, `posterior_predict_<family>()`,
and `posterior_epred_<family>()` — the standard `brms` custom-family
interface functions, located by `brms` via name convention and not normally
called directly.

## References

Skellam JG (1946) The Frequency Distribution of the Difference Between Two
Poisson Variates Belonging to Different Populations. *Journal of the Royal
Statistical Society* 109:296.

Holgate P (1964) Estimation for the Bivariate Poisson Distribution.
*Biometrika* 51:241–245. (The trivariate-reduction construction underlying
`bipois()` and `binegbin()`.)

Karlis D, Ntzoufras I (2003) Analysis of Sports Data by Using Bivariate
Poisson Models. *Journal of the Royal Statistical Society: Series D (The
Statistician)* 52:381–393.

Karlis D, Ntzoufras I (2006) Bayesian Analysis of the Differences of Count
Data. *Statistics in Medicine* 25:1885–1905.

Karlis D, Michels R, Ötting M (2026) Modelling Handball Outcomes Using
Univariate and Bivariate Approaches. *Statistical Methods & Applications*
35:263–284.
