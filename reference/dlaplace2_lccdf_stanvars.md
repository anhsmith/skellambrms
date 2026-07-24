# Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc() (free location and scale)

Returns a \`brms::stanvar()\` defining \`dlaplace2_lccdf\`, the log
complementary CDF of the discrete Laplace(mu, sigma) family –
\`dlaplace2_lccdf(y, mu, sigma)\` = log P(Z \> y). Same role, calling
convention, and no-threshold-argument rationale as
\`dlaplace1_lccdf_stanvars()\`.

## Usage

``` r
dlaplace2_lccdf_stanvars()
```

## Value

A \`brms::stanvars\` object defining the \`dlaplace2_lccdf\` Stan
function, for combining with \`dlaplace2_stanvars()\` via \`+\`.
