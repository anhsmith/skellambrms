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
