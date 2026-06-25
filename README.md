# skellambrms

A [brms](https://paul-buerkner.github.io/brms/) custom family for the symmetric 
Skellam distribution — the distribution of the difference of two independent 
Poisson random variables with equal rates.

## Background

If X ~ Poisson(μ) and Y ~ Poisson(μ) independently, then X − Y ~ Skellam(μ, μ), 
with mean zero and variance 2μ. This arises naturally when modelling set-level 
differences between paired count sources (e.g. two independent observers of the 
same process).

The single parameter μ (modelled on the log scale) controls dispersion, not 
location. The expected value of the response is always zero regardless of 
covariates; covariates explain variation in *spread*.

## Installation

```r
# install.packages("pak")
pak::pak("anhsmith/skellambrms")
```

Stan and a C++ toolchain are required. On Windows, install 
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html). 
Works with either rstan or cmdstanr as the brms backend.

## Usage

```r
library(brms)
library(skellambrms)

fit <- brm(
  y ~ 1 + x + (1 | group),
  data     = dat,
  family   = skellam1(),
  stanvars = skellam1_stanvars(),
  chains   = 4
)
```

The family supports arbitrary brms formula syntax including random effects,
offsets, and non-linear formulas.

## Truncation

`skellam1_lccdf_stanvars()` adds the log-CCDF brms needs to support
`resp_trunc()`, including a row-varying bound:

```r
fit <- brm(
  bf(y | trunc(lb = neg_bound) ~ x),
  data     = dat,
  family   = skellam1(),
  stanvars = skellam1_stanvars() + skellam1_lccdf_stanvars(),
  chains   = 4
)
```

For `mu` above `normal_approx_threshold` (default `100`), the exact
Bessel-sum tail is replaced by a normal approximation, both for speed and
to guard against a confirmed crash when an unadapted HMC proposal pushes
`mu` to an extreme value during warmup. The default was calibrated to one
project's data (real `mu` topping out around 30); see
`?skellam1_lccdf_stanvars` for how to choose this for your own data before
relying on the default elsewhere.

## Functions

| Function | Purpose |
|---|---|
| `skellam1()` | Custom family object for use in `brm()` |
| `skellam1_stanvars()` | Stan code block for use in `brm()` |
| `skellam1_lccdf_stanvars()` | Stan log-CCDF block enabling `resp_trunc()` support |

## Scope

This package implements the symmetric case Skellam(μ, μ) only. The asymmetric 
case (μ₁ ≠ μ₂) is out of scope. `skellam1_lccdf_stanvars()` adds truncation
support to this same symmetric family — it is not a new family.

## Reference

Karlis, D. & Ntzoufras, I. (2006). Bayesian analysis of the differences of count 
data. *Statistics in Medicine*, 25, 1885–1905.
