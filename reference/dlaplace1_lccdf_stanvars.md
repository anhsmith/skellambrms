# Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc()

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `dlaplace1_lccdf`, the log complementary CDF of the discrete
Laplace(0, sigma) family – `dlaplace1_lccdf(y, sigma)` = log P(Z \> y).
Same role and calling convention as
[`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1_lccdf_stanvars.md).
Unlike the Skellam families' lccdf stanvars, this takes no threshold
argument: the closed-form `log1m_exp(double_exponential_lcdf(...))` has
no large-argument failure mode to guard against.

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
