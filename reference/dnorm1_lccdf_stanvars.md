# Log-CCDF of the discrete normal distribution with location 0, for truncated fits

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `dnorm1_lccdf`, the log complementary CDF of the discrete
Normal(0, sigma) distribution: `dnorm1_lccdf(y, sigma)` = log P(Z \> y).
`dnorm1_lccdf_stanvars()` is used in the same way as
[`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md).
`dnorm1_lccdf` computes the upper tail with `erfc()` rather than with
the Stan function `normal_lccdf`, which returns negative infinity for
standardised arguments above about 8.25 (see
[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)).
`dnorm1_lccdf_stanvars()` takes no threshold argument.

## Usage

``` r
dnorm1_lccdf_stanvars()
```

## Value

A
[`brms::stanvars`](https://paulbuerkner.com/brms/reference/stanvar.html)
object defining the `dnorm1_lccdf` Stan function, for combining with
[`dnorm1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
via `+`.

## See also

[`dnorm1()`](https://anhsmith.github.io/skellambrms/reference/dnorm1.md)
for the family itself;
[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
for how
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
locates these functions by name.
