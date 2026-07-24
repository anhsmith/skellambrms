# Discrete-normal custom family for brms (location 0, free scale)

Returns a brms custom family for the discrete normal distribution,
location fixed at 0, discretised from the continuous Normal(0, sigma)
via CDF differencing: \`P(Z=z) = F(z+0.5) - F(z-0.5)\`. One parameter,
sigma (link = "log"), the SD; the mean is always zero. Same
CDF-differencing pattern as \`dlaplace1()\`, using Stan's built-in
\`normal_lcdf\`/\`normal_lccdf\` directly – no Bessel function and no
iteration cap needed, but see the cancellation note below for a branch
this family's PMF does need.

Use in a brm() call as: brm(y ~ ..., family = dnorm1(), stanvars =
dnorm1_stanvars(), data = ...)

## Usage

``` r
dnorm1()

dnorm1_stanvars()

log_lik_dnorm1(i, prep)

posterior_predict_dnorm1(i, prep, ...)

posterior_epred_dnorm1(prep)
```

## Value

A brms custom_family object.

## Details

\*\*Naming note.\*\* Same forced naming as
\`skellam1()\`/\`dlaplace1()\`: \`brms::custom_family()\` requires a
dpar literally named \`"mu"\`; here it represents sigma (the SD), not a
mean. See \`?skellam1\` Details for the full rationale.

\*\*No scale conversion needed.\*\* Unlike \`dlaplace1()\`, where Stan's
\`double_exponential_lcdf\` expects the continuous Laplace's own scale
\`b\` (requiring \`b = sigma / sqrt(2)\` first), the continuous normal's
own SD parameter \*is\* sigma directly – \`sigma\` is passed straight to
\`normal_lcdf\`/\`normal_lccdf\` with no intermediate conversion.

\*\*Cancellation in the PMF, fixed by branching on z's sign.\*\* The
naive \`log_diff_exp(normal_lcdf(z+0.5), normal_lcdf(z-0.5))\` fails
once \`z\` is far enough into the positive tail that both
\`normal_lcdf\` calls round to the same double (both within machine
epsilon of \`log(1)=0\`) – confirmed to occur at only ~10 SDs out, well
inside this package's realistic-but-stressed test range for the other
families, and far sooner than the analogous direct-subtraction form in
\`dlaplace1()\` (the normal's thinner tail saturates near 1 much faster
per SD than the Laplace's). Fixed in \`dnorm1_lpmf\` (and the R-side
\`log_lik_dnorm1\`) by differencing two \*survival\* values
(\`normal_lccdf\`, both small and hence distinguishable) instead of two
\*CDF\* values when \`z \>= 0\` – the same exact-survival-form idea
\`dlaplace1_lccdf\`/\`dlaplace2_lccdf\` already use, applied here to the
PMF rather than the CCDF, since CDF differencing is itself the operation
that creates the cancellation risk in the first place.
