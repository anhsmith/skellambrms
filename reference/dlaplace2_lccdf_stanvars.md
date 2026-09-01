# Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc() (free location and scale)

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `dlaplace2_lccdf`, the log complementary CDF of the discrete
Laplace(mu, sigma) family – `dlaplace2_lccdf(y, mu, sigma)` = log P(Z \>
y). Same role, calling convention, and no-threshold-argument rationale
as
[`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md).

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
