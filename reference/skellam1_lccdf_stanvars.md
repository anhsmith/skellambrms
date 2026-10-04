# Log-CCDF of the symmetric Skellam distribution, for truncated fits

Returns a
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
defining `skellam1_lccdf`, the log complementary CDF of the symmetric
Skellam(mu_skellam, mu_skellam) distribution: `skellam1_lccdf(y, sigma)`
= log P(delta \> y), with `mu_skellam = sigma^2 / 2` derived inside the
function. For truncation with
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html),
brms finds the log-CCDF of a custom family by the name `<family>_lccdf`.
Adding this stanvar to
[`skellam1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
is therefore all that a truncated fit requires, including a fit whose
lower bound varies by row.

## Usage

``` r
skellam1_lccdf_stanvars(normal_approx_threshold = 100)
```

## Arguments

- normal_approx_threshold:

  Numeric scalar. Values of `mu_skellam` above this threshold use the
  normal approximation instead of the exact tail sum. Default `100`.

## Value

A
[`brms::stanvars`](https://paulbuerkner.com/brms/reference/stanvar.html)
object defining the `skellam1_lccdf` Stan function, for combining with
[`skellam1_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
via `+`.

## Details

For `mu_skellam` above `normal_approx_threshold`, `skellam1_lccdf` uses
a normal approximation with variance `2 * mu_skellam`. At or below the
threshold, it sums the PMF exactly, over whichever tail lies away from
the mean: for `y >= 0`, the upper tail upward from `y + 1`; for `y < 0`,
the lower tail downward from `y`, returning the log of one minus that
sum. Each term of the sum evaluates a modified Bessel function, and the
log link places no upper bound on `sigma`, so warmup can propose values
of `mu_skellam` at which the exact sum is slow to evaluate. The sum
stops after 500 terms, or earlier once a term is more than 40 log-units
below the running sum; neither limit can be changed.

A larger threshold applies the exact sum to more evaluations, and a
smaller one applies the normal approximation, which is less accurate in
the tails, to more. The threshold is on the `mu_skellam` scale, not the
`sigma` scale.

## See also

[`skellam1()`](https://anhsmith.github.io/skellambrms/reference/skellam1.md)
for the family itself;
[`skellam2_lccdf_stanvars()`](https://anhsmith.github.io/skellambrms/reference/skellam2_lccdf_stanvars.md)
for the free-mean counterpart.

## Examples

``` r
if (FALSE) { # \dontrun{
library(brms)

brm(
  bf(y | trunc(lb = neg_bound) ~ x),
  family   = skellam1(),
  stanvars = skellam1_stanvars() + skellam1_lccdf_stanvars(),
  data     = dat
)
} # }
```
