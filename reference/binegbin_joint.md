# Censoring-aware joint EM/logbook bivariate-Negative-Binomial family for brms

Censoring-aware extension of \[binegbin()\]. Models the same trivariate-
reduction bivariate Negative-Binomial pair \`(y_em, y_lb)\` – \`y_em =
N_shared + N10\`, \`y_lb = N_shared + N01\`, with \`N_shared ~ NB2(mu,
shapes)\`, \`N10 ~ NB2(lambdaem, shapex)\`, \`N01 ~ NB2(lambdalb,
shapex)\` mutually independent given their rates – but allows the EM
margin (\`y_em\`) to be UNOBSERVED on some rows. Each row carries two
supplementary integers via \`vint()\`: \`y_lb\` (the always-observed
logbook count) and \`em_obs\` (a 0/1 flag marking whether \`y_em\` was
observed for that row).

On \`em_obs == 1\` (matched) rows the likelihood is the full joint
\[binegbin()\] lpmf on \`(y_em, y_lb)\`. On \`em_obs == 0\` (LB-only)
rows it is the \`y_em\`-integrated marginal of that same joint,
\`P(y_lb) = sum_k NB2(k \| mu, shapes) NB2(y_lb - k \| lambdalb,
shapex)\` – NOT a separate single-dispersion \`neg_binomial_2\` on
\`y_lb\`, which would be a different model inconsistent with the matched
decomposition. This lets one \`brm()\` call pool matched and LB-only
rows under one coherent likelihood: \`lambdaem\` and the EM/LB bias are
identified only by the matched rows, while the LB-only rows sharpen
\`mu\`, \`shapes\`, \`lambdalb\`, and the shared vessel/trip
random-effect structure.

Five dpars, identical to \[binegbin()\]: the three rates (\`mu\` =
shared rate, \`lambdaem\`/\`lambdalb\` = EM-/LB-only rates) plus two
dispersions – \`shapes\` for the shared component and \`shapex\` shared
across the two excess components. All five use \`link = "log"\` (see
\[binegbin()\]). To share the excess level across the two rates and
split them by a directional bias, supply them through non-linear
formulas \*without\* an explicit \`exp()\` – the log link applies it, so
\`nlf(lambdaem ~ lamx + methd)\` gives \`lambdaem = exp(lamx + methd)\`.

Use in a brm() call as: brm( bf(y_em \| vint(y_lb, em_obs) ~ 1, mu ~ 1 +
(1 \| vessel) + (1 \| vessel:trip_id), nlf(lambdaem ~ lamx + methd),
nlf(lambdalb ~ lamx - methd), lamx ~ 1, methd ~ 1, shapes ~ 1, shapex ~
1, nl = TRUE), family = binegbin_joint(), stanvars =
binegbin_joint_stanvars(), data = dat )

## Usage

``` r
binegbin_joint()

binegbin_joint_stanvars()

log_lik_binegbin_joint(i, prep)

posterior_predict_binegbin_joint(i, prep, ...)
```

## Value

A brms custom_family object.

## Details

\*\*Two \`vint()\` arguments, in declared order.\*\* brms appends
\`vint()\` integers to the generated lpmf call in the order they are
listed in the formula's \`vint()\` term, matching the \`vars\` declared
here (\`c("vint1\[n\]", "vint2\[n\]")\`): so \`vint(y_lb, em_obs)\`
binds \`vint1 = y_lb\` and \`vint2 = em_obs\`. brms generates \`target
+= binegbin_joint_lpmf(Y\[n\] \| mu\[n\], lambdaem\[n\], lambdalb\[n\],
shapes\[n\], shapex\[n\], vint1\[n\], vint2\[n\])\` – dpars in the order
declared here, then the two vint args. \`binegbin_joint_stan_funs\`
declares \`binegbin_joint_lpmf\` with exactly this signature; reordering
the dpars or the two \`vint()\` terms without matching the Stan
signature silently swaps which quantity governs which component or which
integer is the branch flag.

\*\*Forced \`mu\` naming, and the second count via \`vint()\`.\*\*
Identical conventions to \[binegbin()\]/\[bipois()\] – \`mu\` is brms's
mandatory dpar name, here bound to the shared component's rate, not a
mean of either response; \`y_lb\` (and \`em_obs\`) travel as
supplementary integer data through \`vint()\` because
\`custom_family()\` declares a single response column. See \[bipois()\]
for the full explanation.

\*\*Relationship to \[binegbin()\].\*\* On \`em_obs == 1\` rows this
family's lpmf equals the \[binegbin()\] lpmf exactly (same
marginalisation sum). The \`em_obs == 0\` branch is the
\`y_em\`-integrated marginal of that same bivariate model. The package
tests pin both identities (marginal identity; binegbin equivalence).
