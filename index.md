# pairedcountbrms

Custom [brms](https://paul-buerkner.github.io/brms/) families for
modelling a **pair of counts** from two sources that are meant to
measure the same thing — two observers, two instruments, two reporting
channels — where the question is how much, and how systematically, they
disagree.

The package offers two complementary ways to pose that question, and it
is worth deciding which you want before reading further:

- **Difference families** model the single integer $`d = y_1 - y_2 \in
  \mathbb{Z}`$ directly (Skellam, discrete Laplace, discrete normal).
  Reach for these when *only the disagreement matters* and you want
  truncation support.
- **Joint bivariate-count families** model the pair $`(y_1, y_2)`$
  itself (bivariate Poisson and bivariate negative-binomial). Reach for
  these when the pair’s **overall level, correlation, or marginal
  overdispersion** matters, or when you need to **simulate one count
  given the other**.

The two suites answer different questions and are documented in their
own sections below. The [Overview](#overview) contrasts them; if you
already know which you need, skip to [Difference
families](#difference-families) or [Joint bivariate-count
families](#joint-bivariate-count-families).

## Overview

Many count-comparison problems produce a response that is an integer but
can be negative — the difference between two observers’ counts of the
same event, the goal difference between two teams, the gap between a
camera’s and a logbook’s tally of the same fishing set. Standard `brms`
count families (Poisson, negative binomial, …) model a non-negative
response and cannot take a $`\mathbb{Z}`$-valued one, nor a
jointly-modelled pair. This package fills both gaps.

|  | **Difference families** | **Joint families** |
|----|----|----|
| Response | $`d = y_1 - y_2`$ (one value per pair) | $`(y_1, y_2)`$ (the pair; $`y_2`$ via `vint()`) |
| Families | `skellam1/2`, `dlaplace1/2`, `dnorm1/2` | `bipois`, `binegbin`, `binegbin_joint` |
| Captures | location and spread of the disagreement | level, correlation, marginal overdispersion, **and** the disagreement |
| Discards | the pair’s level and correlation | nothing — but needs both counts observed (except `binegbin_joint`) |
| Truncation ([`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)) | **yes** | no |
| Conditional simulation $`y_1 \mid y_2`$ | n/a | **yes** |
| Overdispersed margins | n/a (models $`d`$, not the margins) | `binegbin` / `binegbin_joint` only |

**Choosing.** If your data are a set of differences and you care only
about whether and by how much the two sources disagree — with, say, a
lower truncation bound because $`d`$ could not physically have fallen
below some value — use a **difference** family. If you have the two
counts side by side and their correlation or absolute level is itself of
interest, or the margins are overdispersed, or you want to impute one
count from the other, use a **joint** family. The two views are
connected: for the bivariate Poisson, the difference $`y_1 - y_2`$ is
*exactly* Skellam-distributed (see below), so `bipois` contains
`skellam2` as its induced difference model.

## Difference families

A difference family models $`d = y_1 - y_2`$ with a distribution on
$`\mathbb{Z}`$. All three underlying distributions are parameterised on
a common **(mean, SD-scale)** convention so that fits are directly
comparable: a location $`\mu`$ (the mean of $`d`$) and a spread
$`\sigma`$ (its standard deviation, on the log scale). Each comes in two
flavours — mean fixed at $`0`$ (does the pair agree on average?) and
free mean (how large is the systematic bias?):

| Family | Mean | Spread | Mean–spread coupling |
|----|----|----|----|
| [`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md) | fixed at $`0`$ | $`\sigma`$ (log) | internally $`\theta_1=\theta_2=\sigma^2/2`$ |
| [`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md) | free $`\mu`$ (identity) | $`\sigma_{\text{excess}}`$ (log) | $`\sigma^2 = \lvert\mu\rvert + \sigma_{\text{excess}}^2`$ — the genuine Skellam constraint; $`\to`$[`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md) at $`\mu=0`$ |
| [`dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md) | fixed at $`0`$ | $`\sigma`$ (log) | none |
| [`dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md) | free $`\mu`$ (identity) | $`\sigma`$ (log) | none — deliberately decoupled |
| [`dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md) | fixed at $`0`$ | $`\sigma`$ (log) | none |
| [`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md) | free $`\mu`$ (identity) | $`\sigma`$ (log) | none — deliberately decoupled |

### Skellam

If $`Y_1 \sim \text{Poisson}(\theta_1)`$ and
$`Y_2 \sim \text{Poisson}(\theta_2)`$ independently, then
$`D = Y_1 - Y_2`$ has the **Skellam** distribution

``` math
P(D = k) = e^{-(\theta_1+\theta_2)}\left(\frac{\theta_1}{\theta_2}\right)^{k/2}
           I_{\lvert k\rvert}\!\left(2\sqrt{\theta_1\theta_2}\right),
\qquad k \in \mathbb{Z},
```

where $`I_\nu`$ is the modified Bessel function of the first kind, with
mean $`\theta_1-\theta_2`$ and variance $`\theta_1+\theta_2`$. The
package samples on $`(\mu,\sigma)`$ rather than $`(\theta_1,\theta_2)`$:

- **[`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)**
  fixes $`\mu=0`$, so $`\theta_1=\theta_2=\sigma^2/2`$ and
  $`\mathrm{Var}(D)=\sigma^2`$. One parameter, $`\sigma`$.
- **[`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)**
  frees the mean via $`\theta_1=(\sigma^2+\mu)/2`$,
  $`\theta_2=(\sigma^2-\mu)/2`$. Because $`\theta_1,\theta_2\ge 0`$
  requires $`\sigma^2\ge\lvert\mu\rvert`$ (variance
  $`\ge \lvert\text{mean}\rvert`$ — a sum of two non-negative rates can
  never be smaller than the size of their difference), the family sets
  $`\sigma^2 = \lvert\mu\rvert +
  \sigma_{\text{excess}}^2`$ with $`\sigma_{\text{excess}}\ge 0`$ free.
  This makes the constraint hold structurally for every $`\mu`$, and
  reduces exactly to
  [`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)
  at $`\mu=0`$. It is a *genuine* coupling: a Skellam difference with a
  large mean must also have large variance.

### Discrete Laplace and discrete normal

Both are obtained by discretising a continuous distribution $`F`$ onto
the integers by **CDF differencing**,

``` math
P(Z = z) = F\!\left(z + \tfrac12\right) - F\!\left(z - \tfrac12\right),
```

with $`F`$ the $`\mathrm{Laplace}(\mu, b)`$ or
$`\mathrm{Normal}(\mu, \sigma)`$ CDF. The scale is put on the same SD
footing as the Skellam families: for the Laplace, $`\mathrm{Var}=2b^2`$,
so $`b=\sigma/\sqrt2`$; for the normal, $`\sigma`$ is already the SD.
The `*1` versions fix $`\mu=0`$; the `*2` versions free it.

Unlike
[`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md),
[`dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
and
[`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md)
impose **no** coupling between $`\mu`$ and $`\sigma`$ — they are free,
independent parameters. That contrast is deliberate and is the reason to
have all three: fitting
[`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)
(bias and spread structurally coupled) against
[`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md)
/
[`dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
(uncoupled) tests whether your data’s disagreement obeys the Skellam
$`\sigma^2 \ge \lvert\text{mean}\rvert`$ relationship or not. The
discrete normal is the light-tailed reference; the discrete Laplace the
heavy-tailed one.

### Usage

Every family follows the same pattern: pass `family = <family>()` and
`stanvars = <family>_stanvars()` to
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html); add
`<family>_lccdf_stanvars()` (combined with `+`) to enable truncation.

``` r

library(brms)
library(pairedcountbrms)

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
call shape (their free-scale dpar is `sigma`,
e.g. `bf(d ~ 1 + x, sigma ~ 1)`). `neg_bound` is a column giving a
(possibly row-varying) lower truncation bound — e.g. `-y_2`, if $`d`$
could not have fallen more than $`y_2`$ below zero for that row. All
families accept arbitrary `brms` formula syntax: random effects, and
non-linear or covariate-dependent predictors on the spread dpar (and,
for the free-mean families, on `mu`).

### Truncation

Each difference family exports `<family>_lccdf_stanvars()`, which
defines a Stan function `<family>_lccdf` — the log complementary CDF
$`\log P(Z > y)`$. `brms`’s
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
finds it purely by name convention and uses it for the truncated
likelihood’s normalising constant, including a row-varying bound. No
wiring beyond adding the stanvar is needed.

For
[`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)/[`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)
the exact log-CCDF is an iterative tail-sum over the Bessel-function
PMF. Above a configurable `normal_approx_threshold` (default `100`, on
the underlying $`\mu_{\text{skellam}}`$ scale — $`\sigma^2/2`$ for
[`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md),
$`(\theta_1+\theta_2)/2`$ for
[`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md))
the exact sum is replaced by a normal approximation. This guards two
confirmed failure modes that occur when HMC warmup pushes the
(log-linked, hence unbounded) spread to an extreme: a crash from a
huge-order Bessel evaluation, and a slower blow-up in cost/memory when
many rows hit the exact loop inside one deep NUTS tree. See
[`?skellam1_lccdf_stanvars`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1_lccdf_stanvars.md)
for how to pick a threshold for your data’s scale. The discrete Laplace
and discrete normal families have closed-form log-CCDFs (via
`double_exponential_lcdf` and an `erfc`-based survival function), so
their `_lccdf_stanvars()` take no threshold — there is no large-argument
mode to guard.

## Joint bivariate-count families

A joint family models the pair $`(y_1, y_2)`$ directly, keeping the
information a difference family throws away — the pair’s level and its
correlation. All three are built by the same **trivariate reduction**:
three independent latent counts

``` math
N_{\text{shared}} \sim \mathcal{D}(\text{mu}), \quad
N_{10} \sim \mathcal{D}(\lambda_{\text{em}}), \quad
N_{01} \sim \mathcal{D}(\lambda_{\text{lb}}),
```

combined as

``` math
y_1 = y_{\text{em}} = N_{\text{shared}} + N_{10}, \qquad
y_2 = y_{\text{lb}} = N_{\text{shared}} + N_{01}.
```

The shared count $`N_{\text{shared}}`$ appears in both, inducing
positive correlation; the two private counts $`N_{10}, N_{01}`$ drive
the difference. (The `em`/`lb` labels — electronic-monitoring vs logbook
— are historical; read them as “count 1, the modelled response” and
“count 2, supplied via `vint()`”. `mu` is the shared *rate*, **not** the
mean of either response.) The likelihood marginalises the unobserved
$`N_{\text{shared}}`$ analytically:

``` math
P(y_{\text{em}}=x,\, y_{\text{lb}}=y)
  = \sum_{k=0}^{\min(x,y)}
    f_{\text{s}}(k)\, f_{10}(x-k)\, f_{01}(y-k),
```

where $`f_{\text{s}}, f_{10}, f_{01}`$ are the pmfs of the three latent
counts. Two consequences follow directly and are worth internalising:

- **The margins.** $`\mathrm{E}[y_{\text{em}}] = \text{mu} +
  \lambda_{\text{em}}`$ and likewise for $`y_{\text{lb}}`$;
  $`\mathrm{Cov}(y_{\text{em}}, y_{\text{lb}}) =
  \mathrm{Var}(N_{\text{shared}})`$, so the correlation is
  $`\mathrm{Var}(N_{\text{shared}}) / \sqrt{\mathrm{Var}(y_{\text{em}})\mathrm{Var}(y_{\text{lb}})}`$.
- **The difference.**
  $`d = y_{\text{em}} - y_{\text{lb}} = N_{10} - N_{01}`$ — the shared
  count *cancels*. So the difference depends only on the two private
  components, exactly the quantity a difference family models. For the
  Poisson case this difference is precisely
  $`\mathrm{Skellam}(\lambda_{\text{em}}, \lambda_{\text{lb}})`$, tying
  the two suites together.

### `bipois` vs `binegbin`: overdispersion

| Family | Latent law | Dispersion | Var of each latent | Use when |
|----|----|----|----|----|
| [`bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md) | Poisson | none | $`\mathrm{Var}=\text{mean}`$ | margins are **not** overdispersed |
| [`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md) | negative-binomial | scalar `shapes` (shared), `shapex` (private) | $`m + m^2/\phi`$ | margins **are** overdispersed (the usual case) |

With Poisson latents, each component has $`\mathrm{Var}=\text{mean}`$,
so `bipois` cannot represent overdispersed margins and underfits the
marginal (and difference) variance of real count data. `binegbin`
replaces each latent with a negative-binomial,
$`N \sim \text{NB2}(m, \phi)`$ — Stan’s `neg_binomial_2`, R’s
`dnbinom(size = φ, mu = m)` — with mean $`m`$ and variance
$`m + m^2/\phi`$. It carries the extra spread in two **scalar**
dispersion dpars: `shapes` $`=\phi_{\text{s}}`$ for the shared count,
`shapex` $`=\phi_{\text{x}}`$ shared across both private counts. The
moments become

``` math
\mathrm{Var}(y_{\text{em}}) = \Big(\text{mu}+\tfrac{\text{mu}^2}{\phi_{\text{s}}}\Big)
  + \Big(\lambda_{\text{em}}+\tfrac{\lambda_{\text{em}}^2}{\phi_{\text{x}}}\Big),
\qquad
\mathrm{Cov}(y_{\text{em}},y_{\text{lb}}) = \text{mu}+\tfrac{\text{mu}^2}{\phi_{\text{s}}},
```

and, with $`\lambda_{\text{em}}=\lambda_{\text{lb}}=\lambda`$,
$`\mathrm{Var}(d) = 2\big(\lambda + \lambda^2/\phi_{\text{x}}\big)`$. As
$`\phi_{\text{s}},\phi_{\text{x}}\to\infty`$ the negative-binomials
collapse to Poissons and `binegbin` $`\to`$`bipois`.

**Why scalar dispersion, not a random effect.** The obvious alternative
— a per-observation random effect (OLRE) on the private components — was
tried and rejected. With one pair per unit but three latent deviates per
unit, the excess deviates act as residual absorbers: their population SD
collapses toward the prior, and drawing fresh deviates fails to
regenerate the observed spread (recovered excess SD $`0.37`$ vs true
$`0.85`$; fresh-draw $`\mathrm{Var}(d)`$$`2.9`$ vs true $`19.2`$ in the
motivating case). A *conditional* posterior-predictive check hides this
entirely; only a **marginal** (fresh-draw) check exposes it. Scalar
`shapes`/`shapex` are identified from the aggregate mean–variance
mismatch across units instead, with no per-unit overfitting.

### `binegbin_joint`: partially-observed pairs (censoring)

`binegbin_joint` is `binegbin` extended to rows where $`y_{\text{em}}`$
was **not observed** — the second margin is censored, not matched. Each
row carries a second `vint()` integer, an `em_obs` $`\in\{0,1\}`$ flag:

- **`em_obs == 1` (matched):** the full `binegbin` joint lpmf on
  $`(y_{\text{em}}, y_{\text{lb}})`$ — identical, term for term.
- **`em_obs == 0` (LB-only):** the $`y_{\text{em}}`$-**integrated
  marginal** of that *same* joint,

``` math
P(y_{\text{lb}}=y) = \sum_{k=0}^{y} f_{\text{s}}(k)\, f_{01}(y-k),
```

i.e. the joint with $`N_{10}`$ summed out (the inner sum over
$`y_{\text{em}}`$ telescopes to $`1`$). This is a convolution of the
shared and LB-excess negative-binomials — **not** a separate
single-dispersion `neg_binomial_2` on $`y_{\text{lb}}`$, which would be
a different, incoherent model.

One [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) call
thus pools matched and censored rows under one likelihood.
$`\lambda_{\text{em}}`$ and the EM/LB bias are identified **only** by
the matched rows; the censored rows sharpen `mu`, `shapes`,
$`\lambda_{\text{lb}}`$, and any shared random-effect structure. See
[Limitations](#limitations) for the identifiability caveat this implies.

### Usage

``` r

library(brms)
library(pairedcountbrms)

# bipois(): joint bivariate Poisson (non-overdispersed margins)
fit_bp <- brm(
  bf(y_em | vint(y_lb) ~ 1,
     mu ~ 1 + (1 | vessel), lambdaem ~ 1, lambdalb ~ 1),
  data = dat, family = bipois(), stanvars = bipois_stanvars(), chains = 4
)

# binegbin(): joint bivariate negative-binomial (overdispersed margins)
fit_nb <- brm(
  bf(y_em | vint(y_lb) ~ 1,
     mu ~ 1 + (1 | vessel),
     nlf(lambdaem ~ lamx), nlf(lambdalb ~ lamx), lamx ~ 1,
     shapes ~ 1, shapex ~ 1, nl = TRUE),
  data = dat, family = binegbin(), stanvars = binegbin_stanvars(), chains = 4
)

# binegbin_joint(): same model, but y_em is unobserved where em_obs == 0
fit_cj <- brm(
  bf(y_em | vint(y_lb, em_obs) ~ 1,
     mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
     nlf(lambdaem ~ lamx + methd),
     nlf(lambdalb ~ lamx - methd),
     lamx ~ 1, methd ~ 1, shapes ~ 1, shapex ~ 1, nl = TRUE),
  data = dat, family = binegbin_joint(), stanvars = binegbin_joint_stanvars(),
  chains = 4
)
```

The `nlf(lambdaem ~ lamx)` / `nlf(lambdalb ~ lamx)` idiom ties the two
private rates to one value — a “no systematic bias” assumption,
$`\mathrm{E}[y_{\text{em}}]=\mathrm{E}[y_{\text{lb}}]`$. Splitting them
as `lamx + methd` / `lamx - methd` introduces a directional bias
parameter `methd` (half the log ratio of the two excess rates); giving
the two rates separate predictors (`nlf(lambdaem ~ lem)`,
`nlf(lambdalb ~ llb)`, `lem ~ 1`, `llb ~ 1`) is the fully unconstrained
version.

The second count (and, for `binegbin_joint`, the flag) travel via
`vint()` because
[`custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
declares a single response column — `vint(y_lb, em_obs)` binds
`vint1 = y_lb`, `vint2 = em_obs` in listed order.

### Prediction

[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
simulates $`y_{\text{em}}`$**conditional on the observed
$`y_{\text{lb}}`$** (which is fixed data, not itself re-simulated):

- **`bipois`:** the conditional split is closed-form,
  $`N_{\text{shared}}\mid y_{\text{lb}} \sim \text{Binomial}\big(y_{\text{lb}},\,
  \text{mu}/(\text{mu}+\lambda_{\text{lb}})\big)`$, then a fresh
  $`N_{10}`$.
- **`binegbin` / `binegbin_joint`:** a negative-binomial sum has no
  Binomial conditional, so the discrete law
  $`P(N_{\text{shared}}=k\mid y_{\text{lb}})
  \propto f_{\text{s}}(k)\,f_{01}(y_{\text{lb}}-k)`$ is sampled
  directly, then a fresh $`N_{10}`$ added.
- For `binegbin_joint`, `em_obs` is **ignored** at prediction time:
  every row (matched and censored alike) gets a conditional
  $`y_{\text{em}}`$ draw, so you can impute the unobserved margin
  fleet-wide.

## Parameterisation and naming notes

**The forced `"mu"` dpar.**
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
unconditionally requires one dpar to be named literally `"mu"`. For
[`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md),
[`dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md),
and
[`dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md)
— the fixed-mean families, whose mean is structurally $`0`$ and *not* a
parameter — that forced `"mu"` slot actually holds `sigma`. If you read
[`make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
output or call `get_dpar(prep, "mu")` for one of these three, you are
looking at $`\sigma`$, not a mean. Every R-side function in the package
immediately rebinds it to `sigma`, so nothing else ever calls it `mu`.
For the joint families the same slot holds the shared *rate*, again not
a response mean. The free-mean difference families
(`skellam2/dlaplace2/dnorm2`) are the only ones whose `mu` genuinely is
the mean.

**`sigmaexcess`, not `sigma_excess`.**
[`custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
disallows dots and underscores in dpar names, so
[`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)’s
excess-spread parameter is spelled `sigmaexcess`.

## Installation

``` r

# install.packages("pak")
pak::pak("anhsmith/pairedcountbrms")
```

Stan and a C++ toolchain are required. On Windows, install
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
Works with either rstan or cmdstanr as the brms backend.

Formerly released as `skellambrms` (versions 0.1.0–0.5.0). No family
names changed in the rename, so existing fits still resolve their
post-processing methods and nothing needs refitting — only
[`library(skellambrms)`](https://rdrr.io/r/base/library.html) becomes
[`library(pairedcountbrms)`](https://github.com/anhsmith/pairedcountbrms).

Documentation, including a worked getting-started vignette that
simulates, fits, and recovers `binegbin` parameters end to end, is at
<https://anhsmith.github.io/pairedcountbrms/>. Locally:

``` r

vignette("pairedcountbrms")
```

## Limitations

**[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
errors on truncated fits (a `brms` limitation).**
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
— and everything built on it, including
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
— errors on *any* truncated custom-family fit. `brms`’s
`posterior_epred.brmsprep()` checks truncation *before* family type and
routes truncated fits to `brms:::posterior_epred_trunc()`, which has no
fallback to a custom family’s own `posterior_epred_<family>()` (it looks
for a non-existent `posterior_epred_trunc_custom()`). This is `brms`’s
dispatch, not this package’s computation. **Workaround** — call the
family method directly:

``` r

prep  <- brms::prepare_predictions(fit)
epred <- posterior_epred_dnorm2(prep)   # or skellam1(), etc.
```

Each family’s `posterior_epred_<family>()` accounts for
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
bounds correctly when called this way.
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
is unaffected and works for truncated fits of every family.

**Joint families do not support truncation.**
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
does not apply to `bipois`/`binegbin`/`binegbin_joint`; no
`_lccdf_stanvars()` is provided for them.

**`posterior_epred` for the joint families.** `bipois` returns the exact
$`\mathrm{E}[y_{\text{em}}\mid y_{\text{lb}}]`$; `binegbin`’s is a point
*approximation* (no closed-form conditional mean for a negative-binomial
sum — use
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
for exact conditional simulation); `binegbin_joint` defines **no**
`posterior_epred` at all, so use
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
for it.

**`binegbin_joint` identifiability.** Because $`\lambda_{\text{em}}`$
and the EM/LB bias are informed only by the matched (`em_obs == 1`)
rows, a fit with few matched rows will learn the bias weakly even if the
total sample is large; the censored rows add power for the shared/level
parameters, not the bias. Fits that lean on the bias should be judged
against the matched subset, not the full $`n`$.

## Testing

For every **difference** family the suite (`tests/testthat/`) checks:

- The R-side log-PMF and log-CCDF against a trusted external reference —
  [`skellam::dskellam()`](https://rdrr.io/pkg/skellam/man/skellam.html)/`pskellam()`
  for Skellam, and a hand-derived, numerically stable log-space
  CDF-differencing reference for the discrete Laplace and discrete
  normal (`extraDistr::ddlaplace()` implements a *different* discrete
  Laplace and is unusable as a reference).
- PMF sums to $`1`$ across a parameter grid; numerical stability (no
  `NaN`/`Inf`) across a realistic-but-stressed range, deep into the
  tails.
- Stan log-PMF/log-CCDF agree with the R references via
  [`rstan::expose_stan_functions()`](https://mc-stan.org/rstan/reference/expose_stan_functions.html).
- For Skellam, the exact and normal-approx log-CCDF branches agree at
  the threshold seam, and changing `normal_approx_threshold` moves the
  cutover.
- Structural (not rejection-based) validity, by inspecting
  [`make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html):
  [`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)’s
  $`\theta_1,\theta_2\ge 0`$ constraint, and the *absence* of coupling
  in
  [`dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)/[`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md).
- Free-mean families reduce exactly to their fixed-mean counterparts at
  $`\mu=0`$.
- End-to-end parameter recovery from simulated (and truncated) data,
  with divergence/Rhat checks.
- That `log_lik_<family>()` and the internal `_lpmf_r()`/`_lccdf_r()`
  helpers return one value **per posterior draw** for a single
  observation — the direction that silently broke
  [`log_lik_dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md)
  (and [`loo()`](https://mc-stan.org/loo/reference/loo.html)) before
  0.3.2, since R’s [`ifelse()`](https://rdrr.io/r/base/ifelse.html)
  takes its length from its test argument.

For every **joint** family the suite validates the marginalised joint
log-PMF against an independent R brute-force reference across a
rate/shape grid and at edge cases; normalisation to $`1`$; the analytic
moment identities (mean, marginal variance, difference variance,
covariance); the Poisson-limit reduction `binegbin` $`\to`$`bipois`; and
end-to-end recovery with divergence/Rhat checks. For `binegbin_joint` it
additionally pins the **marginal identity** ($`\sum_{y_{\text{em}}}`$ of
the matched branch equals the LB-only branch), the **equivalence** to
`binegbin` on `em_obs == 1` rows (R and Stan), and the
**conditional-prediction identity** (`posterior_predict` draws match
joint/marginal), plus a censored end-to-end fit verifying
[`loo()`](https://mc-stan.org/loo/reference/loo.html)/[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
dispatch.

## Function reference

**Difference families** — each exports the family object, its
`_stanvars()`, and a `_lccdf_stanvars()` for truncation, plus
`log_lik_`, `posterior_predict_`, and `posterior_epred_` interface
functions.

| Function | Purpose |
|----|----|
| [`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md) / [`skellam1_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md) | Symmetric Skellam (mean $`0`$) |
| [`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1_lccdf_stanvars.md) | Truncation log-CCDF for `skellam1` |
| [`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md) / [`skellam2_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md) | Asymmetric Skellam (free mean) |
| [`skellam2_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2_lccdf_stanvars.md) | Truncation log-CCDF for `skellam2` |
| [`skellam2_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2_dpars.md) | Reports `mu`, `sigma`, `sigmasq`, `theta1`, `theta2` from a fitted [`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md) (via [`get_dpar()`](https://paulbuerkner.com/brms/reference/get_dpar.html)) |
| [`dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md) / [`dlaplace1_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md) / [`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1_lccdf_stanvars.md) | Discrete Laplace (mean $`0`$) |
| [`dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md) / [`dlaplace2_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md) / [`dlaplace2_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2_lccdf_stanvars.md) | Discrete Laplace (free mean and scale) |
| [`dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md) / [`dnorm1_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md) / [`dnorm1_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1_lccdf_stanvars.md) | Discrete normal (mean $`0`$) |
| [`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md) / [`dnorm2_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md) / [`dnorm2_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2_lccdf_stanvars.md) | Discrete normal (free mean and scale) |

**Joint families** — each exports the family object, its `_stanvars()`,
and `log_lik_` / `posterior_predict_` interface functions
(`posterior_epred_` for `bipois`/`binegbin` only).

| Function | Purpose |
|----|----|
| [`bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md) / [`bipois_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md) | Joint bivariate Poisson |
| [`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md) / [`binegbin_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md) | Joint bivariate negative-binomial (overdispersed margins) |
| [`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md) / [`binegbin_joint_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md) | Censoring-aware bivariate negative-binomial (partially-observed second margin) |

The `log_lik_`, `posterior_predict_`, and `posterior_epred_` functions
are located by `brms` via name convention and are not normally called
directly (except for the truncated-`posterior_epred` workaround above).

## References

Skellam JG (1946) The Frequency Distribution of the Difference Between
Two Poisson Variates Belonging to Different Populations. *Journal of the
Royal Statistical Society* 109:296.

Holgate P (1964) Estimation for the Bivariate Poisson Distribution.
*Biometrika* 51:241–245. (The trivariate-reduction construction
underlying the joint families.)

Karlis D, Ntzoufras I (2003) Analysis of Sports Data by Using Bivariate
Poisson Models. *Journal of the Royal Statistical Society: Series D (The
Statistician)* 52:381–393.

Karlis D, Ntzoufras I (2006) Bayesian Analysis of the Differences of
Count Data. *Statistics in Medicine* 25:1885–1905.

Karlis D, Michels R, Ötting M (2026) Modelling Handball Outcomes Using
Univariate and Bivariate Approaches. *Statistical Methods &
Applications* 35:263–284.
