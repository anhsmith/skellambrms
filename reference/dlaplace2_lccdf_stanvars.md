# Log-CCDF of the discrete Laplace distribution with free location, for truncated fits

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `dlaplace2_lccdf`, the log complementary CDF of the discrete
Laplace(mu, sigma) distribution: `dlaplace2_lccdf(y, mu, sigma)` = log
P(Z \> y). `dlaplace2_lccdf_stanvars()` is used in the same way as
[`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md),
and likewise takes no threshold argument.

## Usage

``` r
dlaplace2_lccdf_stanvars()
```

## Value

A
[`brms::stanvars`](https://paulbuerkner.com/brms/reference/stanvar.html)
object defining the `dlaplace2_lccdf` Stan function, for combining with
[`dlaplace2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
via `+`.

## See also

[`dlaplace2()`](https://anhsmith.github.io/skellambrms/reference/dlaplace2.md)
for the family itself;
[`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md)
for the fixed-mean version.
