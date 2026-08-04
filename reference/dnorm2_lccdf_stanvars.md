# Truncated-discrete-normal log-CCDF for use with brms's resp_trunc() (free location and scale)

Returns a \`brms::stanvar()\` defining \`dnorm2_lccdf\`, the log
complementary CDF of the discrete Normal(mu, sigma) family –
\`dnorm2_lccdf(y, mu, sigma)\` = log P(Z \> y). Same role, calling
convention, and no-threshold-argument rationale as
\`dnorm1_lccdf_stanvars()\`.

## Usage

``` r
dnorm2_lccdf_stanvars()
```

## Value

A \`brms::stanvars\` object defining the \`dnorm2_lccdf\` Stan function,
for combining with \`dnorm2_stanvars()\` via \`+\`.
