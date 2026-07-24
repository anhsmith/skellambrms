# The anatomy of a paired count

Every joint family in this package —
[`bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md),
[`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md),
[`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)
— is built by **trivariate reduction**. Three independent counts are
drawn, and the two observed counts share one of them:

``` math
y_{\text{em}} = N_{\text{shared}} + N_{10}
\qquad
y_{\text{lb}} = N_{\text{shared}} + N_{01}
```

That shared term is the whole point. It is what makes the pair
correlated, it is never observed, and it is marginalised out
analytically in the likelihood. Everything else is the part each source
saw alone.

The widget below lets you take that apart. The novel bit is that it
exposes **two coordinate systems at once, each driving the other**.

## Why two coordinate systems

The dpars a family actually takes are three rates — `mu` for the shared
component, `lambdaem` and `lambdalb` for the two excesses. That is what
the likelihood consumes, but it is awkward to reason in: raise the
overall level of counting and all three move together, so no single one
of them answers “how much was there”, “how much did the two sources
agree”, or “which source ran high”.

Those three questions have their own coordinates:

|  |  |  |
|----|----|----|
| $`M`$ | overall level | $`\mu + (\lambda_{\text{em}} + \lambda_{\text{lb}})/2`$ |
| $`f`$ | congruence — the share of $`M`$ both sources saw | $`\mu / M`$ |
| $`\beta`$ | method bias, bounded on $`[-1,1]`$ | $`(\lambda_{\text{em}} - \lambda_{\text{lb}})/(\lambda_{\text{em}} + \lambda_{\text{lb}})`$ |

The map between them is a bijection, so neither set is more “real” —
they are two descriptions of one object. Drag anything below and watch
the other five respond.

re-simulate

seed

reset

#### Interpretable coordinates

#### Native dpars (what `binegbin()` takes)

## Things to try

**Turn `f` up towards 1.** Both excess rates fall to zero and the two
bars converge — the sources agree completely. Now notice what happens to
$`\beta`$: it stops meaning anything. There is no excess left to be
biased, so the bias is *unidentified*, and the widget flags it and holds
the last value rather than snapping to zero. Zero would be a claim (the
methods are unbiased) that the state cannot support.

**Drag `lambdaem` alone.** $`M`$, $`f`$ and $`\beta`$ all move, because
changing one excess rate changes the overall level, the shared share,
*and* the imbalance simultaneously. This is exactly why the native
coordinates are awkward to reason in, and it is much easier to see than
to describe.

**Turn $`\beta`$ up and watch $`M`$.** The two bars separate, but the
$`M`$ rule does not move — it stays at their *average*, touching
neither. The average of the two excess rates is $`M(1-f)`$ for any
$`\beta`$, so $`M`$ is pinned to the midpoint whatever the bias. $`M`$
is a midpoint of what the sources *report*, not a property of the
underlying process.

**Compare $`\kappa_A`$ against $`\kappa_X`$.** They are near-orthogonal
channels. $`\kappa_A`$ moves the pair up and down *together* — it cannot
touch the difference, because the shared component cancels from it.
$`\kappa_X`$ pulls the pair *apart*, and drives the whole difference.

**Push $`\kappa_A`$ to 0.** The counts do not stop moving. That is the
Poisson floor, not determinism: there is no parameter setting that
stills the stacks.

## The same map, in R

The widget’s arithmetic is not a separate model — it is the package’s
own map, which you can call directly. This is the authoritative version,
and it is what the JavaScript above is checked against:

``` r

library(pairedcountbrms)

# The widget's defaults
binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = atanh(0), kappas = 0.6, kappax = 1.0)
#> $mu
#> [1] 8.04
#> 
#> $lambdaem
#> [1] 3.96
#> 
#> $lambdalb
#> [1] 3.96
#> 
#> $shapes
#> [1] 2.777778
#> 
#> $shapex
#> [1] 1
```

`delta` is the unbounded log-ratio bias; the widget’s $`\beta`$ is the
bounded $`\tanh\delta`$. Going back the other way:

``` r

d <- binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.3)
binegbin_dpars_to_mfd(d$mu, d$lambdaem, d$lambdalb)[c("M", "f", "delta", "beta")]
#> $M
#> [1] 12
#> 
#> $f
#> [1] 0.67
#> 
#> $delta
#> [1] 0.3
#> 
#> $beta
#> [1] 0.2913126
```

And the degenerate case the widget flags:

``` r

# Perfect congruence: no excess, so no bias to identify
binegbin_dpars_to_mfd(mu = 12, lambdaem = 0, lambdalb = 0)$delta
#> [1] NA
```

## Fitting in these coordinates

Nothing above fits anything — these are coordinate transforms. To *fit*
in $`(M, f, \delta)`$ you do not need a different family, because the
reparameterisation is reachable through a non-linear formula. All five
dpars are log-linked, so the link supplies the
[`exp()`](https://rdrr.io/r/base/Log.html):

``` r

bf(y_em | vint(y_lb) ~ 1, nl = TRUE) +
  nlf(mu       ~ eta + log_inv_logit(con)) +
  nlf(lambdaem ~ log(2) + eta + log_inv_logit(-con) + log_inv_logit( 2 * methd)) +
  nlf(lambdalb ~ log(2) + eta + log_inv_logit(-con) + log_inv_logit(-2 * methd)) +
  lf(eta ~ 1, con ~ 1, methd ~ 1)
```

with `eta` $`= \log M`$, `con` $`= \operatorname{logit} f`$, and `methd`
$`= \delta`$. Putting a random effect on `eta` asks whether groups
differ in overall level; putting one on `con` asks whether they differ
in congruence. Those are separable questions in these coordinates and
entangled ones in the native dpars, which is the practical reason to
bother.

## A caveat about the widget

The JavaScript is an independent implementation of the *generative*
model — it draws from the same trivariate reduction, but it is not the
likelihood, and no part of the package’s inference runs in your browser.
The rate map it uses is the one tested in `test-mfd.R`; the sampler is
illustrative. Treat the picture as intuition, and the R above as the
specification.
