# Joint EM/logbook bivariate-Negative-Binomial custom family for brms

Overdispersed sibling of \[bipois()\]. Returns a brms custom family for
the joint distribution of a matched count pair \`(y_em, y_lb)\` via
trivariate reduction with Negative-Binomial (rather than Poisson) latent
components: \`y_em = N_shared + N10\`, \`y_lb = N_shared + N01\`, with
\`N_shared ~ NB2(mu, shapes)\`, \`N10 ~ NB2(lambdaem, shapex)\`, \`N01 ~
NB2(lambdalb, shapex)\` mutually independent given their rates. \`NB2(m,
phi)\` has mean \`m\` and variance \`m + m^2/phi\` (Stan
\`neg_binomial_2\`; R \`dnbinom(size = phi, mu = m)\`).

Five dpars: the three rates (\`mu\` = shared rate,
\`lambdaem\`/\`lambdalb\` = EM-/LB-only rates) plus two dispersions –
\`shapes\` for the shared component and \`shapex\` shared across the two
excess components. All five use \`link = "log"\`. Supply the excess
rates through a non-linear formula without an explicit \`exp()\` (the
log link applies it): \`nlf(lambdaem ~ lamx)\` gives \`lambdaem =
exp(lamx)\`.

See the \`binegbin.R\` file header and the \`tnc001-belize-em\` project
docs (05-07 generative rationale; the OLRE-failure / NegBin-resolution
finding) for why NegBin components are used instead of an
observation-level random effect on \[bipois()\].

Use in a brm() call as: brm( bf(y_em \| vint(y_lb) ~ 1, mu ~ 1 + (1 \|
vessel), nlf(lambdaem ~ lamx), nlf(lambdalb ~ lamx), lamx ~ 1, shapes ~
1, shapex ~ 1, nl = TRUE), family = binegbin(), stanvars =
binegbin_stanvars(), data = dat )

## Usage

``` r
binegbin()

binegbin_stanvars()

log_lik_binegbin(i, prep)

posterior_predict_binegbin(i, prep, ...)

posterior_epred_binegbin(prep)
```

## Value

A brms custom_family object.

## Details

\*\*Forced \`mu\` naming, and \`y_lb\` via \`vint()\`.\*\* Identical
conventions to \[bipois()\] – \`mu\` is brms's mandatory dpar name, here
bound to the shared component's rate (\`lambda_shared\`), not a mean of
either response; \`y_lb\` travels as supplementary integer data through
\`vint()\` because \`custom_family()\` declares a single response
column. See \[bipois()\] for the full explanation.

\*\*Order of dpars matters for the generated Stan call.\*\* brms
generates \`target += binegbin_lpmf(Y\[n\] \| mu\[n\], lambdaem\[n\],
lambdalb\[n\], shapes\[n\], shapex\[n\], vint1\[n\])\` – dpars in the
order declared here, then vint args. \`binegbin_stan_funs\` declares
\`binegbin_lpmf\` with exactly this signature; reordering one without
the other silently swaps which rate or dispersion governs which
component.
