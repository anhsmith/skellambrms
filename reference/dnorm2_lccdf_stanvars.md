# Log-CCDF of the discrete normal distribution with free location, for truncated fits

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `dnorm2_lccdf`, the log complementary CDF of the discrete
Normal(mu, sigma) distribution: `dnorm2_lccdf(y, mu, sigma)` = log P(Z
\> y). `dnorm2_lccdf_stanvars()` is used in the same way as
[`dnorm1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1_lccdf_stanvars.md),
and likewise takes no threshold argument.

## Usage

``` r
dnorm2_lccdf_stanvars()
```

## Value

A
[`brms::stanvars`](https://paulbuerkner.com/brms/reference/stanvar.html)
object defining the `dnorm2_lccdf` Stan function, for combining with
[`dnorm2_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
via `+`.

## See also

[`dnorm2()`](https://anhsmith.github.io/skellambrms/reference/dnorm2.md)
for the family itself;
[`dnorm1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dnorm1_lccdf_stanvars.md)
for the fixed-mean version.
