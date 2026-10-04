# Log-CCDF of the discrete Laplace distribution with location 0, for truncated fits

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `dlaplace1_lccdf`, the log complementary CDF of the discrete
Laplace(0, sigma) distribution: `dlaplace1_lccdf(y, sigma)` = log P(Z \>
y). `dlaplace1_lccdf_stanvars()` is used in the same way as
[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md),
but takes no threshold argument, because the closed form
`log1m_exp(double_exponential_lcdf(...))` involves no Bessel function
and no iterative sum.

## Usage

``` r
dlaplace1_lccdf_stanvars()
```

## Value

A
[`brms::stanvars`](https://paulbuerkner.com/brms/reference/stanvar.html)
object defining the `dlaplace1_lccdf` Stan function, for combining with
[`dlaplace1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
via `+`.

## See also

[`dlaplace1()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1.md)
for the family itself;
[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md)
for how
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
locates these functions by name.
