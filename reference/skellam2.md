# Asymmetric Skellam custom family for brms

Returns a brms custom family for the general (asymmetric) Skellam
distribution, Skellam(theta1, theta2) — the distribution of the
difference of two independent Poisson(theta1), Poisson(theta2) random
variables with possibly unequal rates. Two parameters: \`mu\` (link =
"identity"), the mean of the difference, and \`sigmaexcess\` (link =
"log", so \`\>= 0\`), from which \`sigma\`, the SD of the difference,
and the underlying rates \`theta1\`, \`theta2\` are derived as
transformed quantities (see Details).

Use in a brm() call as: brm(y ~ ..., family = skellam2(), stanvars =
skellam2_stanvars(), data = ...)

## Usage

``` r
skellam2()

skellam2_stanvars()

log_lik_skellam2(i, prep)

posterior_predict_skellam2(i, prep, ...)

posterior_epred_skellam2(prep)
```

## Value

A brms custom_family object.

## Details

\*\*Naming note.\*\* \`brms::custom_family()\` disallows underscores in
\`dpars\` (\`stop2("Dots or underscores are not allowed in
'dpars'.")\`), so the second parameter is spelled \`sigmaexcess\`, not
\`sigma_excess\` as in the package's design notes and Stan code comments
— the two names refer to the same quantity.

\*\*Constraint algebra — corrected from the original design.\*\* The
natural-seeming construction \`sigma = sqrt(mu^2 + sigmaexcess^2)\`
(Pythagorean in mu and sigmaexcess) only guarantees \`sigma \>=
\|mu\|\`. That is NOT the condition Skellam validity actually needs.
With \`theta1 = (sigma^2 + mu) / 2\` and \`theta2 = (sigma^2 - mu) / 2\`
(from \`theta1 + theta2 = sigma^2\` and \`theta1 - theta2 = mu\`),
\`theta1, theta2 \>= 0\` requires \`sigma^2 \>= \|mu\|\` — i.e. Var \>=
\|mean\|, the genuine Skellam constraint (sum of two nonnegative Poisson
rates is always \>= their difference's absolute value). \`sigma \>=
\|mu\|\` and \`sigma^2 \>= \|mu\|\` coincide only when \`\|mu\| \>= 1\`;
for \`\|mu\| \< 1\` they diverge, and \`sigma = sqrt(mu^2 +
sigmaexcess^2)\` can produce a \*negative\* theta1 or theta2 — confirmed
numerically, e.g. \`mu = 0.5, sigmaexcess = 0\` gives \`sigma = 0.5\`,
\`theta2 = -0.125\`. This package instead uses: sigma^2 = \|mu\| +
sigmaexcess^2 which guarantees \`sigma^2 \>= \|mu\|\` directly (the
right-hand side is \`\|mu\|\` plus a nonnegative term), for every \`mu\`
and every \`sigmaexcess \>= 0\`, with equality (the minimal-spread
boundary) at \`sigmaexcess = 0\`. \`theta1\` and \`theta2\` are then
both sums of nonnegative terms (verify: for \`mu \>= 0\`, \`theta1 =
mu + sigmaexcess^2/2 \>= 0\` and \`theta2 = sigmaexcess^2/2 \>= 0\`; for
\`mu \< 0\`, the roles swap) — strictly positive whenever \`sigmaexcess
\> 0\`, which the log link guarantees for any finite linear predictor.
At \`mu = 0\` this reduces exactly to skellam1's symmetric family
(\`sigma = sigmaexcess\`, \`theta1 = theta2 = sigmaexcess^2 / 2\`).

\*\*Generated-quantities note.\*\* This family does \*not\* expose
\`mu\`, \`sigma\`, \`sigma^2\`, \`theta1\`, \`theta2\` via a Stan
\`generated quantities\` block. Confirmed via \`make_stancode()\`: brms
declares a custom family's per-observation dpar vectors (\`mu\`,
\`sigmaexcess\` here) as local variables inside the generated model's
\`model\` block, not \`transformed parameters\` — out of Stan-scope for
\`generated quantities\`, regardless of \`loop = TRUE/FALSE\`.
Reconstructing them from brms's internal linear-predictor variable names
(\`Xc\`, \`b\`, \`Intercept\`, ...) would only work for simple
fixed-effects-only formulas, breaking silently for anything with random
effects or splines. \`skellam2_dpars()\` (below) reports the same five
quantities from R instead, via \`brms::get_dpar()\` — works for any
formula.
