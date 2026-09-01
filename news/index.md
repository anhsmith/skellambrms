# Changelog

## skellambrms 0.6.0

- **The difference families are back under their original package
  name.** Between 0.5.0 and this release, this code was distributed as
  `pairedcountbrms` (0.6.0–0.8.0), which contained two unrelated suites:
  these six families, which model the difference `d = y1 - y2`, and a
  set that models the count *pair* jointly. They shared no code — not a
  helper, not a Stan block, not a dependency — and have been separated
  again. The joint families are now
  [`bicountbrms`](https://github.com/anhsmith/bicountbrms). Its 0.9.0
  was the matching half of this release, at the time of the split.
  Install its current version rather than 0.9.0: `bicountbrms` 0.10.0
  removed the `bipois_cens()`, `binegbin_cens()` and `binegbin_joint`
  names that 0.9.0 defined, and supplies `bipois()`, `binegbin()` and a
  `_partialobs()` constructor for each in their place.

  **No family, dpar or Stan function name changed**, so no fitted model
  needs refitting. brms resolves a fit’s `log_lik_*` /
  `posterior_predict_*` / `posterior_epred_*` methods off the attached
  search path at call time, so a stored fit works as soon as this
  package is attached. The only source change required is
  [`library(pairedcountbrms)`](https://rdrr.io/r/base/library.html) →
  [`library(skellambrms)`](https://github.com/anhsmith/skellambrms), and
  any `pairedcountbrms::` prefix.

  Note on the version number: `pairedcountbrms` also had a 0.6.0, a
  different release of a differently-named package. This release
  continues *this* package’s line from 0.5.0.

- **`anhsmith/skellambrms` no longer redirects.** From the 0.6.0 rename
  until now, GitHub served a permanent redirect from this path to
  `anhsmith/pairedcountbrms`. Recreating the repository here breaks that
  redirect deliberately: a user who installed `skellambrms` wanted the
  Skellam families, and following the redirect would now land them on a
  package that no longer has any. `pak::pak("anhsmith/skellambrms")`
  installs these families again, which is what it did originally.

- **Roxygen markdown is enabled, which changes every help page.**
  `DESCRIPTION` gained `Roxygen: list(markdown = TRUE)`, which this
  package had never set, so roxygen had been passing markdown through to
  Rd uninterpreted. Regenerating the 13 Rd files converted 502 literal
  backticks into 251 code spans. Bold and italic markers that had been
  printing as their own asterisks now render. The topics also gained
  cross-references, of which there had been none: each family links to
  its `_lccdf_stanvars()` companion and to the families it is meant to
  be fitted against.

  Markdown interpretation is retroactive, so every topic was rendered
  with [`tools::Rd2txt()`](https://rdrr.io/r/tools/Rd2HTML.html) before
  and after and diffed word by word. One expression had been altered
  silently. In
  [`?skellam1`](https://anhsmith.github.io/skellambrms/reference/skellam1.md),
  the prior-translation intercept `0.5*log(2) + 0.5*1` was read as
  emphasis between its two asterisks and rendered as `0.5log(2) + 0.51`,
  dropping both multiplication signs; it is now a code span. Under
  [`?skellam1_lccdf_stanvars`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md),
  the two guarded failure modes had been running together into one
  paragraph, and are now a list.

- **Every documented function states what it returns.** The six family
  topics document five functions each and had a single `\value{}`
  describing the constructor alone. Each now gives the return value and
  shape of `log_lik_*`, `posterior_predict_*` and `posterior_epred_*`
  separately, as well as those of the constructor and its `_stanvars()`
  companion. Two behaviours that were documented nowhere are stated
  there: `posterior_predict_*` draws subject to a row’s
  [`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
  bounds, and `posterior_epred_*` takes its mean over the truncated
  distribution on any row that is bounded. `Title:` and `Description:`
  also quote `'brms'` and `'Stan'`, which is the convention CRAN applies
  to software names.

- **Metadata for archiving.** `LICENSE` and `LICENSE.md` named
  “skellambrms authors”, the placeholder `usethis` writes, and now name
  the copyright holder. A `.zenodo.json` records the creator, ORCID,
  affiliation and licence for the Zenodo deposit, and is excluded from
  the built package. `inst/CITATION` contains a commented `doi` field to
  fill in once the release is archived, with a note that a manuscript
  should wait for the DOI rather than cite the repository URL.

- **Wording, across README, NEWS, roxygen and test comments.** A sweep
  against the prose conventions in `CLAUDE.md` replaced vague verbs of
  agency and possession with the relations they stood for, split
  balanced “X, and Y” constructions into separate sentences, and gave
  bare pronouns and quantifiers their nouns. Two statements were wrong
  rather than loose. The README said the normal-approximation threshold
  “guards two confirmed failure modes” when it guards against them, and
  it described `bipois()` as containing
  [`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  when the relation is that the difference implied by `bipois()` is
  [`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md).

- The calibration instrument `coverage_recovery()` and its
  `PAIREDCOUNTBRMS_COVERAGE` gate were built for the joint families’
  dispersion parameters and went with them to `bicountbrms`. The smoke
  gate this suite actually calls, `recovery_ok()`, stays, in
  `tests/testthat/helper-recovery.R`, which records where the
  calibration half went and why the gate here should be widened rather
  than tightened if a correct model starts failing it.

- `make_synthetic_prep()` (`tests/testthat/helper-prep.R`) loses its
  `vint1`/`vint2` arguments. These families declare a single response
  and take no supplementary integer data; only the joint families needed
  them.

- The package no longer ships a vignette. The one it had fitted a
  `binegbin()` model and went to `bicountbrms` with that family. The
  README documents every family here with runnable examples.

------------------------------------------------------------------------

Entries below record this code’s history under its former names.
Releases 0.4.0 and 0.5.0 added and extended the joint bivariate-count
families and say nothing about the six families in this package;
0.6.0–0.8.0 were made under the name `pairedcountbrms`, and concern the
joint families almost throughout. They are left as written rather than
retrospectively edited. The full record of those releases, including
everything omitted here, is in [`bicountbrms`’s
NEWS](https://github.com/anhsmith/bicountbrms/blob/master/NEWS.md).

Of the entries in that span, two concern these families:

- **0.7.0** re-parameterised
  [`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  onto `sigma` (see that release’s notes for the prior translation) and
  rewrote the recovery tests as smoke gates.
- **0.6.0** renamed the package to `pairedcountbrms`.

## skellambrms 0.3.2

- Fixed a silent [`ifelse()`](https://rdrr.io/r/base/ifelse.html)
  length-collapse bug in
  [`log_lik_dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md):
  `ifelse(test, yes, no)` takes its output length from `test`, not from
  the (vectorised) `yes`/`no` branches. Because
  [`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  has no free mean, the CDF-differencing argument built from the
  observation `z` is a scalar, while `sigma` (and the derived `b`)
  varies across posterior draws — so the
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) test was evaluated at
  length 1 and the whole per-observation log-likelihood silently
  collapsed to length 1 instead of `ndraws`. This broke
  `brms::add_criterion(fit, "loo")` for every
  [`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  fit (`is.matrix(unnormalized_log_weights) is not TRUE`), while
  sampling,
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
  and truncation were entirely unaffected — confirmed isolated to this
  one function’s R-side length handling, not a data or convergence
  issue.
- Applied the same fix pre-emptively to five more internal R-side
  helpers sharing the identical
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) shape —
  `dlaplace1_lpmf_r()`, `dlaplace1_lccdf_r()`, `dlaplace2_lccdf_r()`,
  `dnorm1_lpmf_r()`, and `dnorm2_lpmf_r()` in `R/truncation.R`. None
  were triggering the bug at their current call sites (which happen to
  keep argument lengths matched), but all shared the same landmine.

## skellambrms 0.3.1

- Fixed `posterior_predict_<family>()` for all six families: previously
  ignored
  [`trunc()`](https://rdrr.io/r/base/Round.html)/[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
  bounds entirely, drawing from the untruncated distribution and
  returning it verbatim (confirmed to produce out-of-bound draws,
  e.g. [`posterior_predict_dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  returning values well below a `lb = -14` bound). Now performs correct
  inverse-CDF sampling within the truncation bounds, reusing each
  family’s already-validated log-CCDF math (transcribed to R in the new
  internal `R/truncation.R`) rather than rejection sampling, which was
  confirmed empirically slow/low-acceptance for tight bounds —
  especially costly for
  [`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md),
  whose per-evaluation cost (an iterative Bessel-function tail-sum) is
  comparatively high.
- Fixed `posterior_epred_<family>()` for all six families: previously
  returned the untruncated mean (`mu`, or `0` for
  [`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)/
  [`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)/[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md))
  even when a truncation bound was tight enough to meaningfully shift
  the conditional expectation, with no warning. Now computes the correct
  truncated conditional expectation via deterministic numerical
  summation of the truncated PMF — not Monte Carlo, so the result is
  exact to a documented tolerance and fully reproducible.
- **Known limitation surfaced (not introduced) by this fix:**
  [`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  — and anything built on it, including
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  — errors on any truncated fit of a custom family, for all six families
  here, under the currently installed `brms`. This is a `brms`
  limitation: its internal dispatcher checks whether a fit is truncated
  *before* checking whether the family is a custom one, and has no
  fallback to a custom family’s own `posterior_epred_<family>()` on the
  truncated branch. Call
  `posterior_epred_<family>(brms::prepare_predictions(fit))` directly as
  a workaround.
  [`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  is unaffected and works correctly for truncated fits of every family.
  See the README’s “Limitations” section for details.

## skellambrms 0.3.0

- **Breaking change:**
  [`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  now samples on `sigma`, the SD of the difference (log-linked), rather
  than the underlying Skellam rate directly. `mu_skellam = sigma^2 / 2`
  is derived internally; the Bessel-sum likelihood itself is unchanged.
  A prior previously written on `log(mu_skellam)` translates as
  `log(sigma) = 0.5*log(2) + 0.5*log(mu_skellam)` — e.g. an old
  `normal(1, 1.5)` becomes `normal(0.847, 0.75)`. This
  reparameterisation establishes a common (mean, SD-scale) convention
  shared by every family below.
- Added
  [`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  /
  [`skellam2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
  /
  [`skellam2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_lccdf_stanvars.md):
  the asymmetric Skellam (Koopman parameterisation), with a free mean
  (`mu`) and a free `sigmaexcess` (so that
  `sigma^2 = |mu| + sigmaexcess^2`, guaranteeing Skellam validity for
  every `mu` and `sigmaexcess >= 0` — a corrected constraint relative to
  the originally-specified `sigma = sqrt(mu^2 + sigmaexcess^2)`, which
  only guarantees the weaker `sigma >= |mu|` and admits invalid
  (negative-rate) parameter combinations for `|mu| < 1`). Reduces
  exactly to
  [`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  at `mu = 0`.
- Added
  [`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  /
  [`dlaplace1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
  /
  [`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md):
  a discrete Laplace distribution (location fixed at 0, free `sigma`),
  discretised from the continuous Laplace via CDF differencing.
- Added
  [`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
  /
  [`dlaplace2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
  /
  [`dlaplace2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2_lccdf_stanvars.md):
  the free-location/free-scale discrete Laplace, with no constraint
  coupling `mu` and `sigma` — a deliberate structural contrast with
  [`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md),
  for comparing models where bias and spread are structurally coupled
  against ones where they vary independently.
- Added
  [`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
  /
  [`dnorm1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
  /
  [`dnorm1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1_lccdf_stanvars.md):
  a discrete normal distribution (location fixed at 0, free `sigma`),
  via the same CDF-differencing pattern as
  [`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md).
- Added
  [`dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  /
  [`dnorm2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
  /
  [`dnorm2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm2_lccdf_stanvars.md):
  the free-location/free-scale discrete normal, structurally analogous
  to
  [`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md).
- Fixed a numerical-stability issue affecting
  [`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
  and
  [`skellam2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_lccdf_stanvars.md)’s
  normal-approximation branch, and `dnorm1`/`dnorm2`’s `_lpmf`/`_lccdf`:
  Stan’s built-in `normal_lccdf` is not safe to call directly in this
  context. This is a documented Stan limitation, not a guess — the Stan
  Functions Reference states `normal_lccdf` underflows to `-inf` for
  `(y-mu)/sigma > ~8.25`, and
  [stan-dev/math#1985](https://github.com/stan-dev/math/issues/1985)
  confirms `normal_lccdf` (unlike `normal_lcdf`) was never updated with
  the more accurate Mills-ratio approximation. Fixed via an exact
  `erfc()`-based closed form throughout, confirmed to match a trusted R
  reference to machine precision out to 30+ SDs.
- Fixed a Stan-compiler portability bug:
  `skellam2_lpmf`/`skellam2_lccdf` used `fabs()`, which compiles under
  `rstan`’s bundled Stan version but is not a valid identifier under
  `cmdstanr`’s (use [`abs()`](https://rdrr.io/r/base/MathFun.html),
  which is type-generic and already used elsewhere in the same
  functions).
- Added `cmdstanr` to `Suggests` (previously only `rstan` was declared,
  so `R CMD check`’s isolated test environment could not see an
  already-installed `cmdstanr`).

## skellambrms 0.2.0

- Added
  [`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md),
  providing the log-CCDF of the symmetric Skellam(mu, mu) distribution
  so that brms’s
  [`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
  can be used with
  [`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md),
  including row-varying truncation bounds. Still the symmetric
  Skellam(mu, mu) case only — this adds truncation support to the
  existing family, not a new family or the asymmetric case.
- The exact log-CCDF (an iterative Bessel-sum tail) switches to a normal
  approximation above a configurable `normal_approx_threshold` (default
  `100`), guarding against a confirmed `std::bad_alloc` crash and a
  confirmed multi-GB memory blowup when an unadapted HMC proposal pushes
  `mu` to an extreme value during warmup. See
  [`?skellam1_lccdf_stanvars`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
  for guidance on choosing this threshold for your own data.

## skellambrms 0.1.0

- Initial release:
  [`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
  and
  [`skellam1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md),
  a brms custom family for the symmetric Skellam(mu, mu) distribution.
