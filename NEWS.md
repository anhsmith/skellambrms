# skellambrms 0.3.2

* Fixed a silent `ifelse()` length-collapse bug in `log_lik_dlaplace1()`:
  `ifelse(test, yes, no)` takes its output length from `test`, not from
  the (vectorised) `yes`/`no` branches. Because `dlaplace1()` has no free
  mean, the CDF-differencing argument built from the observation `z` is a
  scalar, while `sigma` (and the derived `b`) varies across posterior
  draws — so the `ifelse()` test was evaluated at length 1 and the whole
  per-observation log-likelihood silently collapsed to length 1 instead of
  `ndraws`. This broke `brms::add_criterion(fit, "loo")` for every
  `dlaplace1()` fit (`is.matrix(unnormalized_log_weights) is not TRUE`),
  while sampling, `posterior_predict()`, and truncation were entirely
  unaffected — confirmed isolated to this one function's R-side length
  handling, not a data or convergence issue.
* Applied the same fix pre-emptively to five more internal R-side helpers
  sharing the identical `ifelse()` shape — `dlaplace1_lpmf_r()`,
  `dlaplace1_lccdf_r()`, `dlaplace2_lccdf_r()`, `dnorm1_lpmf_r()`, and
  `dnorm2_lpmf_r()` in `R/truncation.R`. None were triggering the bug at
  their current call sites (which happen to keep argument lengths
  matched), but all shared the same landmine.

# skellambrms 0.3.1

* Fixed `posterior_predict_<family>()` for all six families: previously
  ignored `trunc()`/`resp_trunc()` bounds entirely, drawing from the
  untruncated distribution and returning it verbatim (confirmed to produce
  out-of-bound draws, e.g. `posterior_predict_dnorm2()` returning values
  well below a `lb = -14` bound). Now performs correct inverse-CDF sampling
  within the truncation bounds, reusing each family's already-validated
  log-CCDF math (transcribed to R in the new internal `R/truncation.R`)
  rather than rejection sampling, which was confirmed empirically
  slow/low-acceptance for tight bounds — especially costly for
  `skellam2()`, whose per-evaluation cost (an iterative Bessel-function
  tail-sum) is comparatively high.
* Fixed `posterior_epred_<family>()` for all six families: previously
  returned the untruncated mean (`mu`, or `0` for `skellam1()`/
  `dlaplace1()`/`dnorm1()`) even when a truncation bound was tight enough
  to meaningfully shift the conditional expectation, with no warning. Now
  computes the correct truncated conditional expectation via deterministic
  numerical summation of the truncated PMF — not Monte Carlo, so the
  result is exact to a documented tolerance and fully reproducible.
* **Known limitation surfaced (not introduced) by this fix:**
  `brms::posterior_epred()` — and anything built on it, including
  `fitted()` and `conditional_effects()` — errors on any truncated fit of
  a custom family, for all six families here, under the currently
  installed `brms`. This is a `brms` limitation: its internal dispatcher
  checks whether a fit is truncated *before* checking whether the family
  is a custom one, and has no fallback to a custom family's own
  `posterior_epred_<family>()` on the truncated branch. Call
  `posterior_epred_<family>(brms::prepare_predictions(fit))` directly as a
  workaround. `brms::posterior_predict()` is unaffected and works
  correctly for truncated fits of every family. See the README's
  "Limitations" section for details.

# skellambrms 0.3.0

* **Breaking change:** `skellam1()` now samples on `sigma`, the SD of the
  difference (log-linked), rather than the underlying Skellam rate
  directly. `mu_skellam = sigma^2 / 2` is derived internally; the
  Bessel-sum likelihood itself is unchanged. A prior previously stated on
  `log(mu_skellam)` translates as
  `log(sigma) = 0.5*log(2) + 0.5*log(mu_skellam)` — e.g. an old
  `normal(1, 1.5)` becomes `normal(0.847, 0.75)`. This reparameterisation
  establishes a common (mean, SD-scale) convention shared by every family
  below.
* Added `skellam2()` / `skellam2_stanvars()` / `skellam2_lccdf_stanvars()`:
  the asymmetric Skellam (Koopman parameterisation), with a free mean
  (`mu`) and a free `sigmaexcess` (so that
  `sigma^2 = |mu| + sigmaexcess^2`, guaranteeing Skellam validity for
  every `mu` and `sigmaexcess >= 0` — a corrected constraint relative to
  the originally-specified `sigma = sqrt(mu^2 + sigmaexcess^2)`, which
  only guarantees the weaker `sigma >= |mu|` and admits invalid
  (negative-rate) parameter combinations for `|mu| < 1`). Reduces exactly
  to `skellam1()` at `mu = 0`.
* Added `dlaplace1()` / `dlaplace1_stanvars()` / `dlaplace1_lccdf_stanvars()`:
  a discrete Laplace distribution (location fixed at 0, free `sigma`),
  discretised from the continuous Laplace via CDF differencing.
* Added `dlaplace2()` / `dlaplace2_stanvars()` / `dlaplace2_lccdf_stanvars()`:
  the free-location/free-scale discrete Laplace, with no constraint
  coupling `mu` and `sigma` — a deliberate structural contrast with
  `skellam2()`, for comparing models where bias and spread are
  structurally coupled against ones where they vary independently.
* Added `dnorm1()` / `dnorm1_stanvars()` / `dnorm1_lccdf_stanvars()`: a
  discrete normal distribution (location fixed at 0, free `sigma`), via
  the same CDF-differencing pattern as `dlaplace1()`.
* Added `dnorm2()` / `dnorm2_stanvars()` / `dnorm2_lccdf_stanvars()`: the
  free-location/free-scale discrete normal, structurally analogous to
  `dlaplace2()`.
* Fixed a numerical-stability issue affecting `skellam1_lccdf_stanvars()`
  and `skellam2_lccdf_stanvars()`'s normal-approximation branch, and
  `dnorm1`/`dnorm2`'s `_lpmf`/`_lccdf`: Stan's built-in `normal_lccdf` is
  not safe to call directly in this context. This is a documented Stan
  limitation, not a guess — the Stan Functions Reference states
  `normal_lccdf` underflows to `-inf` for `(y-mu)/sigma > ~8.25`, and
  [stan-dev/math#1985](https://github.com/stan-dev/math/issues/1985)
  confirms `normal_lccdf` (unlike `normal_lcdf`) was never updated with
  the more accurate Mills-ratio approximation. Fixed via an exact
  `erfc()`-based closed form throughout, confirmed to match a trusted R
  reference to machine precision out to 30+ SDs.
* Fixed a Stan-compiler portability bug: `skellam2_lpmf`/`skellam2_lccdf`
  used `fabs()`, which compiles under `rstan`'s bundled Stan version but
  is not a valid identifier under `cmdstanr`'s (use `abs()`, which is
  type-generic and already used elsewhere in the same functions).
* Added `cmdstanr` to `Suggests` (previously only `rstan` was declared,
  so `R CMD check`'s isolated test environment could not see an
  already-installed `cmdstanr`).

# skellambrms 0.2.0

* Added `skellam1_lccdf_stanvars()`, providing the log-CCDF of the
  symmetric Skellam(mu, mu) distribution so that brms's `resp_trunc()`
  can be used with `skellam1()`, including row-varying truncation
  bounds. Still the symmetric Skellam(mu, mu) case only — this adds
  truncation support to the existing family, not a new family or the
  asymmetric case.
* The exact log-CCDF (an iterative Bessel-sum tail) switches to a normal
  approximation above a configurable `normal_approx_threshold` (default
  `100`), guarding against a confirmed `std::bad_alloc` crash and a
  confirmed multi-GB memory blowup when an unadapted HMC proposal pushes
  `mu` to an extreme value during warmup. See `?skellam1_lccdf_stanvars`
  for guidance on choosing this threshold for your own data.

# skellambrms 0.1.0

* Initial release: `skellam1()` and `skellam1_stanvars()`, a brms custom
  family for the symmetric Skellam(mu, mu) distribution.
