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
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.htm). 
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

## Functions

| Function | Purpose |
|---|---|
| `skellam1()` | Custom family object for use in `brm()` |
| `skellam1_stanvars()` | Stan code block for use in `brm()` |

## Scope

This package implements the symmetric case Skellam(μ, μ) only. The asymmetric 
case (μ₁ ≠ μ₂) is out of scope.

## Reference

Karlis, D. & Ntzoufras, I. (2006). Bayesian analysis of the differences of count 
data. *Statistics in Medicine*, 25, 1885–1905.
