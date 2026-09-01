# Report skellam2's derived quantities from a fitted model

Returns `mu`, `sigma`, `sigma^2`, `theta1`, and `theta2` (each a draws x
observations matrix) from a
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
`brmsfit`, computed in R via
[`brms::get_dpar()`](https://paulbuerkner.com/brms/reference/get_dpar.html)
rather than a Stan `generated quantities` block — see
"Generated-quantities note" in
[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
for why the latter isn't available for this family.

## Usage

``` r
skellam2_dpars(fit, newdata = NULL)
```

## Arguments

- fit:

  A `brmsfit` fitted with `family = skellam2()`.

- newdata:

  Optional new data, passed to
  [`brms::prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html).

## Value

A named list of draws x observations matrices: `mu`, `sigma`, `sigmasq`,
`theta1`, `theta2`.

## See also

[`skellam2()`](https://anhsmith.github.io/skellambrms/reference/skellam2.md)
for the family and the algebra these quantities are derived from.
