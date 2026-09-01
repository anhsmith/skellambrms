# skellambrms <img src="man/figures/logo.png" align="right" height="139" alt="skellambrms logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/anhsmith/skellambrms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anhsmith/skellambrms/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/1322796194.svg)](https://doi.org/10.5281/zenodo.22231870)
<!-- badges: end -->

Custom [brms](https://paulbuerkner.com/brms/) families (Bürkner 2017) for the
**difference between two paired counts**, $d = y_1 - y_2 \in \mathbb{Z}$: two
observers, two instruments, two reporting channels, where the question is how
much and how systematically they disagree.

Six families, three distributions in two flavours each — the Skellam (the
difference of two Poisson variates), the discrete Laplace, and the discrete
normal. The `1` variants fix the location at zero (do the two sources agree on
average?); the `2` variants estimate it (by how much do they disagree?). All six
are parameterised on a common (mean, SD-scale) convention so that fits are
directly comparable, and all six support truncation via `resp_trunc()`.

Standard `brms` count families (Poisson, negative binomial, …) model a
non-negative response and cannot take a $\mathbb{Z}$-valued one. This package
supplies that.

## Information discarded by reducing a pair to its difference

Reducing a pair to its difference discards the pair's overall level and the
correlation between its two sources. Under the trivariate-reduction view — the
two counts share an unobserved latent component, plus a private component each
— the shared component **cancels** in $d$. So $d$ alone identifies neither the
level nor the correlation; only the disagreement survives.

Discarding the level and the correlation is the right trade when the
disagreement *is* the question. Modelling the difference is also the only
option when the response must be truncated, since `resp_trunc()` takes a
univariate response. Where the level and the congruence matter too, model the
pair jointly instead: the companion package
[`bicountbrms`](https://github.com/anhsmith/bicountbrms) supplies bivariate
Poisson and bivariate negative-binomial families that retain the level, the
correlation and the disagreement together, and admit rows on which only one of
the two counts was observed. The two packages are connected: for the bivariate
Poisson, the induced difference $y_1 - y_2$ is exactly Skellam-distributed, so
the difference model implied by `bipois()` is `skellam2()`.

One thing not to do, whether you model the difference or the pair: regressing
the difference on one of the two counts places that count on both sides of the
equation. The fitted slope is then biased towards $-1$ by shared sampling error
alone (Bland and Altman 1986).

## Difference families

A difference family models $d = y_1 - y_2$ with a distribution on $\mathbb{Z}$.
All three underlying distributions are parameterised on a common **(mean,
SD-scale)** convention so that fits are directly comparable: a location $\mu$
(the mean of $d$) and a spread $\sigma$ (its standard deviation, on the log
scale). Each comes in two flavours — mean fixed at $0$ (does the pair agree on
average?) and free mean (how large is the systematic bias?):

| Family | Mean | Spread | Mean–spread coupling |
|---|---|---|---|
| `skellam1()` | fixed at $0$ | $\sigma$ (log) | internally $\theta_1=\theta_2=\sigma^2/2$ |
| `skellam2()` | free $\mu$ (identity) | $\sigma_{\text{excess}}$ (log) | $\sigma^2 = \lvert\mu\rvert + \sigma_{\text{excess}}^2$ — the genuine Skellam constraint; $\to$ `skellam1()` at $\mu=0$ |
| `dlaplace1()` | fixed at $0$ | $\sigma$ (log) | none |
| `dlaplace2()` | free $\mu$ (identity) | $\sigma$ (log) | none — deliberately decoupled |
| `dnorm1()` | fixed at $0$ | $\sigma$ (log) | none |
| `dnorm2()` | free $\mu$ (identity) | $\sigma$ (log) | none — deliberately decoupled |

### Skellam

If $Y_1 \sim \text{Poisson}(\theta_1)$ and $Y_2 \sim \text{Poisson}(\theta_2)$
independently, then $D = Y_1 - Y_2$ has the **Skellam** distribution
(Skellam 1946)

$$
P(D = k) = e^{-(\theta_1+\theta_2)}\left(\frac{\theta_1}{\theta_2}\right)^{k/2}
           I_{\lvert k\rvert}\!\left(2\sqrt{\theta_1\theta_2}\right),
\qquad k \in \mathbb{Z},
$$

where $I_\nu$ is the modified Bessel function of the first kind, with mean
$\theta_1-\theta_2$ and variance $\theta_1+\theta_2$. The package samples on
$(\mu,\sigma)$ rather than $(\theta_1,\theta_2)$:

- **`skellam1()`** fixes $\mu=0$, so $\theta_1=\theta_2=\sigma^2/2$ and
  $\mathrm{Var}(D)=\sigma^2$. One parameter, $\sigma$.
- **`skellam2()`** frees the mean via $\theta_1=(\sigma^2+\mu)/2$,
  $\theta_2=(\sigma^2-\mu)/2$. Because $\theta_1,\theta_2\ge 0$ requires
  $\sigma^2\ge\lvert\mu\rvert$ (variance $\ge \lvert\text{mean}\rvert$ — a sum of
  two non-negative rates can never be smaller than the size of their
  difference), the family sets $\sigma^2 = \lvert\mu\rvert +
  \sigma_{\text{excess}}^2$ with $\sigma_{\text{excess}}\ge 0$ free. This makes
  the constraint hold structurally for every $\mu$, and reduces exactly to
  `skellam1()` at $\mu=0$. It is a *genuine* coupling: a Skellam difference
  with a large mean must also have large variance.

### Discrete Laplace and discrete normal

The discrete Laplace and the discrete normal are both obtained by discretising
a continuous distribution $F$ onto the integers by **CDF differencing**,

$$
P(Z = z) = F\!\left(z + \tfrac12\right) - F\!\left(z - \tfrac12\right),
$$

with $F$ the $\mathrm{Laplace}(\mu, b)$ or $\mathrm{Normal}(\mu, \sigma)$ CDF. The scale is put
on the same SD footing as the Skellam families: for the Laplace,
$\mathrm{Var}=2b^2$, so $b=\sigma/\sqrt2$; for the normal, $\sigma$ is
already the SD. The `*1` versions fix $\mu=0$; the `*2` versions free it.

Unlike `skellam2()`, `dlaplace2()` and `dnorm2()` impose **no** coupling
between $\mu$ and $\sigma$ — they are free, independent parameters. That
contrast is deliberate. Fitting `skellam2()` (bias and spread structurally
coupled) against `dnorm2()` / `dlaplace2()` (uncoupled) tests whether your
data's disagreement obeys the Skellam $\sigma^2 \ge \lvert\text{mean}\rvert$
relationship or not. The discrete normal is the light-tailed reference; the
discrete Laplace the heavy-tailed one.

### Usage

Every family follows the same pattern: pass `family = <family>()` and
`stanvars = <family>_stanvars()` to `brm()`; add `<family>_lccdf_stanvars()`
(combined with `+`) to enable truncation.

```r
library(brms)
library(skellambrms)

# skellam1(): mean fixed at 0 -- do the two sources agree on average?
fit1 <- brm(
  bf(d | trunc(lb = neg_bound) ~ 1 + (1 | group)),
  data     = dat,
  family   = skellam1(),
  stanvars = skellam1_stanvars() + skellam1_lccdf_stanvars(),
  chains   = 4
)

# skellam2(): free mean -- how large, and how uncertain, is the disagreement?
fit2 <- brm(
  bf(d | trunc(lb = neg_bound) ~ 1 + x, sigmaexcess ~ 1),
  data     = dat,
  family   = skellam2(),
  stanvars = skellam2_stanvars() + skellam2_lccdf_stanvars(),
  chains   = 4
)
```

`dlaplace1/2()` and `dnorm1/2()` are drop-in replacements with the same
call shape (their free-scale dpar is `sigma`, e.g. `bf(d ~ 1 + x, sigma ~ 1)`).
`neg_bound` is a column giving a (possibly row-varying) lower truncation
bound — e.g. `-y_2`, if $d$ could not have fallen more than $y_2$ below zero
for that row. All families accept arbitrary `brms` formula syntax: random
effects, and non-linear or covariate-dependent predictors on the spread dpar
(and, for the free-mean families, on `mu`).

### Truncation

Each difference family exports `<family>_lccdf_stanvars()`, which defines a
Stan function `<family>_lccdf` — the log complementary CDF $\log P(Z > y)$.
`brms`'s `resp_trunc()` finds it purely by name convention and uses it for the
truncated likelihood's normalising constant, including a row-varying bound. No
wiring beyond adding the stanvar is needed.

For `skellam1()`/`skellam2()`, the exact log-CCDF is an iterative tail-sum over
the Bessel-function PMF. Above a configurable `normal_approx_threshold`
(default `100`, on the underlying $\mu_{\text{skellam}}$ scale — $\sigma^2/2$
for `skellam1()`, $(\theta_1+\theta_2)/2$ for `skellam2()`) the exact sum is
replaced by a normal approximation. This guards against two confirmed failure
modes that occur when HMC warmup pushes the (log-linked, hence unbounded)
spread to an extreme: a crash from a huge-order Bessel evaluation, and a slower
blow-up in cost/memory when many rows hit the exact loop inside one deep NUTS
tree. See `?skellam1_lccdf_stanvars` for how to pick a threshold for your
data's scale. The discrete Laplace and discrete normal families have
closed-form log-CCDFs (via `double_exponential_lcdf` and an `erfc`-based
survival function), so their `_lccdf_stanvars()` take no threshold — there is
no large-argument mode to guard.

## Parameterisation and naming notes

**The forced `"mu"` dpar.** `brms::custom_family()` unconditionally requires one
dpar to be named literally `"mu"`. For `skellam1()`, `dlaplace1()` and
`dnorm1()` — the fixed-mean families, whose mean is structurally $0$ and *not* a
parameter — that forced `"mu"` slot actually holds `sigma`. If you read
`make_stancode()` output or call `get_dpar(prep, "mu")` for one of these three,
you are looking at $\sigma$, not a mean. Every R-side function in the package
immediately rebinds it to `sigma`, so nothing else ever calls it `mu`. The
free-mean families (`skellam2`, `dlaplace2`, `dnorm2`) are the ones whose `mu`
genuinely is the mean.

**`sigmaexcess`, not `sigma_excess`.** `custom_family()` disallows dots and
underscores in dpar names, so `skellam2()`'s excess-spread parameter is spelled
`sigmaexcess`.

### Notation

| Code | Math | Meaning |
|---|---|---|
| `d` (the response) | $d = y_1 - y_2$ | the difference, on $\mathbb{Z}$ |
| `mu` (free-mean families) | $\mu$ | mean of $d$ |
| `mu` (fixed-mean families) | $\sigma$ | the spread — see above; **not** a mean |
| `sigma` | $\sigma$ | SD of $d$, log-linked |
| `sigmaexcess` | $\sigma_{\text{excess}}$ | `skellam2()`'s free spread, with $\sigma^2 = \lvert\mu\rvert + \sigma_{\text{excess}}^2$ |
| — | $\theta_1$, $\theta_2$ | the two Poisson rates a Skellam difference is built from |

Source 1 is whichever count you subtract from. Which of the two counts takes
that label is arbitrary, and nothing else in the package depends on it.

## Installation

```r
# install.packages("pak")
pak::pak("anhsmith/skellambrms")
```

Stan and a C++ toolchain are required. On Windows, install
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
Works with either rstan or cmdstanr as the brms backend.

**Package history.** These families were released as `skellambrms` (0.1.0–0.5.0).
At 0.6.0 the package was renamed `pairedcountbrms`, and a second, unrelated
suite was added to it: families that model the count *pair* jointly rather than
its difference. The two suites shared no code. After `pairedcountbrms` 0.8.0,
they were separated again: the difference families return here as `skellambrms`
0.6.0, under the name they were first released with, and the joint families
continue as [`bicountbrms`](https://github.com/anhsmith/bicountbrms), which was
0.9.0 at the split.

No family, dpar or Stan function name in *this* package changed in any of it, so
no stored fit needs refitting — brms resolves a fit's `log_lik_*` /
`posterior_predict_*` / `posterior_epred_*` methods off the attached search path
at call time. The only source change is `library(pairedcountbrms)` →
`library(skellambrms)`. (One joint family was renamed on the other side of the
split, `binegbin_joint` → `binegbin_cens`, and `bicountbrms` 0.10.0 then
removed the `_cens` names outright, with no deprecation layer; see
`bicountbrms`'s `NEWS.md` if you hold fits under either name.)

Because this repository was recreated at the old name, `anhsmith/skellambrms`
now serves this package directly rather than redirecting to `pairedcountbrms`.
That is the intended behaviour: a user who installed `skellambrms` wanted these
families, and would otherwise have been redirected to a package that no longer
contains them.

Documentation is at <https://anhsmith.github.io/skellambrms/>.

## Limitations

**`posterior_epred()` errors on truncated fits (a `brms` limitation).**
`brms::posterior_epred()` — and everything built on it, including `fitted()`
and `conditional_effects()` — errors on *any* truncated custom-family fit.
`brms`'s `posterior_epred.brmsprep()` checks truncation *before* family type
and routes truncated fits to `brms:::posterior_epred_trunc()`, which has no
fallback to a custom family's own `posterior_epred_<family>()` (it looks for a
non-existent `posterior_epred_trunc_custom()`). This is `brms`'s dispatch, not
this package's computation. **Workaround** — call the family method directly:

```r
prep  <- brms::prepare_predictions(fit)
epred <- posterior_epred_dnorm2(prep)   # or skellam1(), etc.
```

Each family's `posterior_epred_<family>()` accounts for `resp_trunc()` bounds
correctly when called this way. `posterior_predict()` is unaffected and works
for truncated fits of every family.

**The pair's level and correlation are not recoverable.** By construction — see
[Information discarded by reducing a pair to its difference](#information-discarded-by-reducing-a-pair-to-its-difference).
This is not a defect of the implementation but the definition of what a
difference model is.

## Testing

For every family the suite (`tests/testthat/`) checks:

- The R-side log-PMF and log-CCDF against a trusted external reference —
  `skellam::dskellam()`/`pskellam()` for Skellam, and a hand-derived,
  numerically stable log-space CDF-differencing reference for the discrete
  Laplace and discrete normal (`extraDistr::ddlaplace()` implements a
  *different* discrete Laplace and is unusable as a reference).
- PMF sums to $1$ across a parameter grid; numerical stability (no `NaN`/`Inf`)
  across a realistic-but-stressed range, deep into the tails.
- Stan log-PMF/log-CCDF agree with the R references via
  `rstan::expose_stan_functions()`.
- For Skellam, the exact and normal-approx log-CCDF branches agree at the
  threshold seam, and changing `normal_approx_threshold` moves the cutover.
- Structural (not rejection-based) validity, by inspecting `make_stancode()`:
  `skellam2()`'s $\theta_1,\theta_2\ge 0$ constraint, and the *absence* of
  coupling in `dlaplace2()`/`dnorm2()`.
- Free-mean families reduce exactly to their fixed-mean counterparts at
  $\mu=0$.
- End-to-end parameter recovery from simulated (and truncated) data, with
  divergence/Rhat checks.
- That `log_lik_<family>()` and the internal `_lpmf_r()`/`_lccdf_r()` helpers
  return one value **per posterior draw** for a single observation — the
  direction that silently broke `log_lik_dlaplace1()` (and `loo()`) before
  0.3.2, since R's `ifelse()` takes its length from its test argument.

### Running the tests

Everything that needs a Stan toolchain is behind `skip_on_cran()`, and every
file-level compilation behind `stan_tests_enabled()`, so a plain `R CMD check`
runs the analytic and R-side tests but **skips every model fit and every
compile**. To include them:

```bash
NOT_CRAN=true R CMD check --no-manual skellambrms_0.6.0.tar.gz
```

The recovery tests are **smoke gates**, not calibration statements. Each fits
once and checks convergence (zero divergences, $\hat R < 1.01$) and that the
true value lies within a wide (99%) posterior interval, which detects gross
mis-specification. Asserting instead that a truth falls within a *90%* interval
from a single fit fails 10% of the time by construction for a correct model, and
that is how a pair of assertions in this suite came to fail on every run for
several releases without anyone noticing. See
`tests/testthat/helper-recovery.R`, which also says where the matching
calibration instrument went.

**What runs automatically.** `.github/workflows/R-CMD-check.yaml` defines two
jobs on different cadences, since the checks differ by an order of magnitude in
cost:

| job | when | what |
|---|---|---|
| `check` | every push and pull request | builds, docs, examples, and the tests that need no Stan backend |
| `check-stan` | weekly, and on demand | sets `NOT_CRAN=true` and installs both backends, so every model fit actually runs |

Trigger `check-stan` by hand from the repository's Actions tab whenever you want
it. GitHub disables scheduled workflows after 60 days without a commit to the
repository, so the weekly run stops once 60 days pass with nothing pushed; the
manual trigger is the fallback.

## Function reference

Each family exports the family object, its `_stanvars()`, and a
`_lccdf_stanvars()` for truncation, plus `log_lik_`, `posterior_predict_` and
`posterior_epred_` interface functions.

| Function | Purpose |
|---|---|
| `skellam1()` / `skellam1_stanvars()` | Symmetric Skellam (mean $0$) |
| `skellam1_lccdf_stanvars()` | Truncation log-CCDF for `skellam1` |
| `skellam2()` / `skellam2_stanvars()` | Asymmetric Skellam (free mean) |
| `skellam2_lccdf_stanvars()` | Truncation log-CCDF for `skellam2` |
| `skellam2_dpars()` | Reports `mu`, `sigma`, `sigmasq`, `theta1`, `theta2` from a fitted `skellam2()` (via `get_dpar()`) |
| `dlaplace1()` / `dlaplace1_stanvars()` / `dlaplace1_lccdf_stanvars()` | Discrete Laplace (mean $0$) |
| `dlaplace2()` / `dlaplace2_stanvars()` / `dlaplace2_lccdf_stanvars()` | Discrete Laplace (free mean and scale) |
| `dnorm1()` / `dnorm1_stanvars()` / `dnorm1_lccdf_stanvars()` | Discrete normal (mean $0$) |
| `dnorm2()` / `dnorm2_stanvars()` / `dnorm2_lccdf_stanvars()` | Discrete normal (free mean and scale) |

The `log_lik_`, `posterior_predict_` and `posterior_epred_` functions are
located by `brms` via name convention and are not normally called directly
(except for the truncated-`posterior_epred` workaround above).

## References

Bland JM, Altman DG (1986) Statistical methods for assessing agreement between
two methods of clinical measurement. *The Lancet* 327:307–310.

Bürkner P-C (2017) brms: an R package for Bayesian multilevel models using Stan.
*Journal of Statistical Software* 80:1–28.

Skellam JG (1946) The Frequency Distribution of the Difference Between Two
Poisson Variates Belonging to Different Populations. *Journal of the Royal
Statistical Society* 109:296.

Karlis D, Ntzoufras I (2006) Bayesian Analysis of the Differences of Count
Data. *Statistics in Medicine* 25:1885–1905.

Holgate P (1964) Estimation for the Bivariate Poisson Distribution.
*Biometrika* 51:241–245. (The trivariate-reduction construction underlying the
joint families in `bicountbrms`, referred to above.)

Genest C, Nešlehová J (2007) A primer on copulas for count data.
*ASTIN Bulletin* 37:475–515.
