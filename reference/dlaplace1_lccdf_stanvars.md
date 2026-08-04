# Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc()

Returns a \`brms::stanvar()\` defining \`dlaplace1_lccdf\`, the log
complementary CDF of the discrete Laplace(0, sigma) family –
\`dlaplace1_lccdf(y, sigma)\` = log P(Z \> y). Same role and calling
convention as \`skellam1_lccdf_stanvars()\`. Unlike the Skellam
families' lccdf stanvars, this takes no threshold argument: the
closed-form \`log1m_exp(double_exponential_lcdf(...))\` has no
large-argument failure mode to guard against.

## Usage

``` r
dlaplace1_lccdf_stanvars()
```

## Value

A \`brms::stanvars\` object defining the \`dlaplace1_lccdf\` Stan
function, for combining with \`dlaplace1_stanvars()\` via \`+\`.
