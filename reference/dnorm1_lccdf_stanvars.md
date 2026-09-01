# Truncated-discrete-normal log-CCDF for use with brms's resp_trunc()

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `dnorm1_lccdf`, the log complementary CDF of the discrete
Normal(0, sigma) family – `dnorm1_lccdf(y, sigma)` = log P(Z \> y). Same
role and calling convention as
[`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/dlaplace1_lccdf_stanvars.md),
but built directly on Stan's `normal_lccdf` (an upper-tail log-survival
function Stan exposes as a built-in for the normal), rather than a
`log1m_exp(lcdf(...))` composition – no threshold argument and no
large-argument failure mode to guard against.

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
