# Log-CCDF of the asymmetric Skellam distribution, for truncated fits

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `skellam2_lccdf`, the log complementary CDF of the asymmetric
Skellam(theta1, theta2) distribution:
`skellam2_lccdf(y, mu, sigmaexcess)` = log P(delta \> y).
`skellam2_lccdf_stanvars()` is used in the same way as
[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md),
and the documentation of
[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
describes how
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
finds the function, the normal approximation and the exact sum. In
`skellam2_lccdf`, the threshold is compared with
`mu_skellam = (theta1 + theta2) / 2`, which equals the `mu_skellam` of
[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
when `theta1 = theta2`, and the exact sum changes from the upper to the
lower tail at `y = mu` rather than at `y = 0`.

## Usage

``` r
skellam2_lccdf_stanvars(normal_approx_threshold = 100)
```

## Arguments

- normal_approx_threshold:

  Numeric scalar, compared with `mu_skellam` as described above. Default
  `100`.

## Value

A
[`brms::stanvars`](https://paulbuerkner.com/brms/reference/stanvar.html)
object defining the `skellam2_lccdf` Stan function, for combining with
[`skellam2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
via `+`.

## See also

[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
for the family itself;
[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
for the normal approximation and the exact sum.
