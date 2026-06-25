# Stan function block for the symmetric Skellam log-PMF. skellam1 is
# sampled on sigma (the SD of the difference, log-linked) rather than the
# underlying Skellam rate -- mu_skellam = sigma^2 / 2 is the single
# Skellam(mu_skellam, mu_skellam) rate parameter (theta1 = theta2 =
# mu_skellam); see skellam1() in family.R for the old-prior/new-prior
# correspondence on the sigma scale. The Bessel-sum likelihood itself is
# unchanged -- only the parameter it's expressed in terms of changed.
# Injected into the Stan model via stanvars.
skellam1_stan_funs <- "
  real skellam1_lpmf(int k, real sigma) {
    real mu_skellam = square(sigma) / 2;
    return -2 * mu_skellam + log_modified_bessel_first_kind(abs(k), 2 * mu_skellam);
  }
"

# Stan function block for the symmetric Skellam log-CCDF (truncation
# support via brms's resp_trunc()). Operates on mu_skellam = sigma^2 / 2,
# the same derived quantity as skellam1_lpmf above. The normal-
# approximation threshold is templated in at call time rather than fixed,
# since the right value depends on the data's plausible mu_skellam range
# -- see skellam1_lccdf_stanvars() for the rationale and how to choose it
# (the threshold is on the mu_skellam scale, not sigma, so existing
# calibration advice carries over unchanged under the reparameterisation).
# The iteration cap (500) and early-exit tolerance are fixed: they guard
# a confirmed std::bad_alloc crash in log_modified_bessel_first_kind and
# a confirmed multi-GB memory blowup at extreme mu_skellam, independent of
# where the normal-approximation threshold is set.
skellam1_lccdf_stan <- function(normal_approx_threshold = 100) {
  sprintf("
  real skellam1_lccdf(int y, real sigma) {
    // log P(delta > y). mu_skellam = sigma^2 / 2 is the underlying
    // Skellam(mu,mu) rate. Beyond the threshold, skip the exact
    // Bessel-sum tail and use a normal approximation (Skellam(mu,mu) has
    // variance 2*mu_skellam, so CLT applies).
    real mu_skellam = square(sigma) / 2;
    if (mu_skellam > %s) {
      real z = (y + 0.5) / sqrt(2 * mu_skellam);
      // NOT normal_lccdf(z | 0, 1): documented Stan limitation, not a
      // hunch -- the Stan Functions Reference states normal_lccdf
      // underflows to -inf for (y-mu)/sigma above ~8.25, and
      // stan-dev/math#1985 confirms normal_lccdf (unlike normal_lcdf)
      // was never updated with the more accurate Mills-ratio
      // approximation, so it still hits this floor. Confirmed directly
      // in this package's own testing: normal_lccdf(9.5|0,1) returns
      // -inf here vs the true -48.3. Exact closed form instead:
      // P(Z>z) = 0.5*erfc(z/sqrt(2)), confirmed to match R's pnorm() to
      // high precision out to z=30+ (erfc() itself does not share
      // normal_lccdf's accuracy gap).
      return log(0.5) + log(erfc(z / sqrt(2)));
    }
    real acc = negative_infinity();
    int k = y + 1;
    int hard_cap = y + 1 + 500;
    while (k < hard_cap) {
      real lp_k = -2 * mu_skellam + log_modified_bessel_first_kind(abs(k), 2 * mu_skellam);
      real new_acc = log_sum_exp(acc, lp_k);
      if (lp_k < new_acc - 40 && k > y + 5) {
        acc = new_acc;
        break;
      }
      acc = new_acc;
      k += 1;
    }
    return acc;
  }
", normal_approx_threshold)
}

# Stan function block for the asymmetric Skellam log-PMF (Koopman-style
# mean/SD parameterisation). Sampled dpars are mu (the mean, identity
# link) and sigmaexcess (lb=0, log link; brms forbids underscores in
# dpars names, hence the no-underscore spelling -- see skellam2() in
# family.R). sigma^2 = |mu| + sigmaexcess^2 is the key construction: it
# guarantees sigma^2 >= |mu|, the actual Skellam constraint (Var >=
# |mean|, equivalently theta1, theta2 both >= 0), by construction, for
# every mu and every sigmaexcess >= 0 -- see family.R for why the more
# obvious-looking sigma = sqrt(mu^2 + sigmaexcess^2) does NOT work (it
# only guarantees sigma >= |mu|, a weaker and, for |mu| < 1, insufficient
# condition). With sigmaexcess > 0 strictly (guaranteed by the log link
# for any finite linear predictor), theta1 and theta2 are both strictly
# positive for every mu, including mu = 0, where the family reduces
# exactly to skellam1's symmetric lpmf (sigma = sigmaexcess,
# theta1 = theta2 = sigmaexcess^2 / 2). Matches skellam::dskellam's own
# internal formula exactly (verified against its source).
skellam2_stan_funs <- "
  real skellam2_lpmf(int k, real mu, real sigmaexcess) {
    real sigmasq = abs(mu) + square(sigmaexcess);
    real theta1  = (sigmasq + mu) / 2;
    real theta2  = (sigmasq - mu) / 2;
    return -theta1 - theta2 + (k / 2.0) * log(theta1 / theta2)
           + log_modified_bessel_first_kind(abs(k), 2 * sqrt(theta1 * theta2));
  }
"

# Stan function block for the asymmetric Skellam log-CCDF. Same
# Bessel-tail-sum + normal-approximation-above-threshold pattern as
# skellam1_lccdf_stan, generalised to a nonzero mean: the threshold is
# checked against mu_skellam = (theta1 + theta2) / 2, exactly the
# quantity skellam1_lccdf_stan thresholds on (there, theta1 = theta2 =
# mu_skellam, so this reduces identically), and the normal approximation
# uses mean mu and SD sigma (rather than mean 0 and SD sqrt(2*mu_skellam)
# as in the symmetric case) -- Var(Skellam(theta1,theta2)) = theta1 +
# theta2 = sigma^2 regardless of asymmetry, so this is the direct
# generalisation, not a new derivation.
skellam2_lccdf_stan <- function(normal_approx_threshold = 100) {
  sprintf("
  real skellam2_lccdf(int y, real mu, real sigmaexcess) {
    // log P(delta > y).
    real sigmasq    = abs(mu) + square(sigmaexcess);
    real theta1     = (sigmasq + mu) / 2;
    real theta2     = (sigmasq - mu) / 2;
    real mu_skellam = (theta1 + theta2) / 2;
    if (mu_skellam > %s) {
      real sigma = sqrt(sigmasq);
      real z = (y + 0.5 - mu) / sigma;
      // Same documented normal_lccdf limitation and erfc-based exact fix
      // as skellam1_lccdf_stan above (stan-dev/math#1985) -- see there
      // for the citation and the confirmed failure point.
      return log(0.5) + log(erfc(z / sqrt(2)));
    }
    real acc = negative_infinity();
    int k = y + 1;
    int hard_cap = y + 1 + 500;
    while (k < hard_cap) {
      real lp_k = -theta1 - theta2 + (k / 2.0) * log(theta1 / theta2)
                  + log_modified_bessel_first_kind(abs(k), 2 * sqrt(theta1 * theta2));
      real new_acc = log_sum_exp(acc, lp_k);
      if (lp_k < new_acc - 40 && k > y + 5) {
        acc = new_acc;
        break;
      }
      acc = new_acc;
      k += 1;
    }
    return acc;
  }
", normal_approx_threshold)
}

# Stan function block for the discrete-Laplace log-PMF (location fixed at
# 0, free scale). Discretised from the continuous Laplace via CDF
# differencing -- P(Z=z) = F(z+0.5) - F(z-0.5), where F is the
# continuous Laplace(0, b) CDF -- using Stan's built-in CDF for this
# distribution, a closed form with no Bessel-style blowup risk and so,
# unlike skellam1/skellam2, no large-argument branch or iteration cap
# needed. Naming note: Stan calls the Laplace distribution
# "double_exponential", not "laplace" -- `laplace_lcdf` does not exist
# and fails to compile ("undeclared identifier"); the correct builtin is
# `double_exponential_lcdf`, confirmed to match the hand-derived R
# reference exactly.
#
# Sampled parameter is sigma (log-linked, same convention as the other
# three families); b is the scale double_exponential_lcdf expects.
# Algebra: Var(continuous Laplace(0,b)) = 2*b^2, so SD = b*sqrt(2).
# Treating sigma as exactly that continuous-Laplace SD (the
# discretisation perturbs the true discrete variance only slightly, and
# this keeps "sigma" on a common, directly-comparable scale across all
# four families, which is the point of the convention) gives
# sigma = b*sqrt(2), i.e.
#   b = sigma / sqrt(2)
# computed first thing in both functions below, before anything is
# passed to double_exponential_lcdf.
dlaplace1_stan_funs <- "
  real dlaplace1_lpmf(int z, real sigma) {
    real b = sigma / sqrt2();
    return log_diff_exp(double_exponential_lcdf(z + 0.5 | 0, b), double_exponential_lcdf(z - 0.5 | 0, b));
  }
"

# Stan function block for the discrete-Laplace log-CCDF (truncation
# support via brms's resp_trunc()). log P(Z > y) = log(1 - F(y+0.5)) =
# log1m_exp(double_exponential_lcdf(y+0.5 | 0, b)) -- single closed form,
# same b conversion as dlaplace1_lpmf above, no branching.
dlaplace1_lccdf_stan <- "
  real dlaplace1_lccdf(int y, real sigma) {
    real b = sigma / sqrt2();
    return log1m_exp(double_exponential_lcdf(y + 0.5 | 0, b));
  }
"

# Stan function block for the discrete-Laplace log-PMF with free location
# AND free scale -- no constraint coupling mu and sigma (a genuine
# structural difference from skellam2's sigma >= |mu| floor: the point
# of having both an asymmetric-Skellam and a free-location discrete-
# Laplace family is to compare a model where bias and spread are
# structurally coupled against one where they're independent, so no
# constraint is imposed here). Same CDF-differencing and b conversion as
# dlaplace1 (b = sigma / sqrt(2)), but mu is passed straight through to
# double_exponential_lcdf's own location argument rather than shifting z
# manually -- Stan's double_exponential_lcdf(y | mu, b), like
# normal_lcdf(y | mu, sigma), takes location and scale directly.
dlaplace2_stan_funs <- "
  real dlaplace2_lpmf(int z, real mu, real sigma) {
    real b = sigma / sqrt2();
    return log_diff_exp(double_exponential_lcdf(z + 0.5 | mu, b), double_exponential_lcdf(z - 0.5 | mu, b));
  }
"

# Stan function block for the discrete-Laplace log-CCDF, free location
# and scale. Same closed form as dlaplace1_lccdf_stan, generalised to a
# nonzero mu via double_exponential_lcdf's own location argument.
dlaplace2_lccdf_stan <- "
  real dlaplace2_lccdf(int y, real mu, real sigma) {
    real b = sigma / sqrt2();
    return log1m_exp(double_exponential_lcdf(y + 0.5 | mu, b));
  }
"

# Stan function block for the discrete-normal log-PMF (location fixed at
# 0, free scale). Discretised from the continuous Normal via CDF
# differencing -- P(Z=z) = F(z+0.5) - F(z-0.5). Unlike dlaplace1/
# dlaplace2, no scale conversion is needed first: the continuous
# normal's own SD *is* sigma (no b-style intermediate parameter), so
# sigma is passed straight through.
#
# NOT simply log_diff_exp(normal_lcdf(z+0.5|0,sigma), normal_lcdf(z-0.5|0,sigma))
# on both sides -- confirmed that naive form catastrophically cancels
# once z is far enough into the *positive* tail that normal_lcdf(z-0.5)
# and normal_lcdf(z+0.5) both round to the same double (both within
# machine epsilon of log(1)=0), at which point log_diff_exp(a,a) = -inf
# regardless of how the true (tiny but distinct) difference should come
# out. This isn't a remote edge case: it was confirmed to occur at only
# ~10 SDs out for sigma=1 -- inside this package's existing
# "realistic-but-stressed" test range for the other families (e.g.
# sigma up to 100, k within 10 SDs; see test-skellam2.R /
# test-dlaplace1.R). Fixed the same way dlaplace1/dlaplace2's lccdf uses
# the exact survival form for x>=0 rather than 1-F(x): for z on the far
# side of the mean (z>=0 here), difference two *survival* values instead
# of two *CDF* values (both near 1 and hence not distinguishable).
#
# The survival values themselves are NOT computed via normal_lccdf,
# despite that being Stan's built-in upper-tail log-survival function --
# this is a documented Stan limitation, not a guess: the Stan Functions
# Reference states normal_lccdf underflows to -inf for
# (y-mu)/sigma > ~8.25, and stan-dev/math#1985 confirms normal_lccdf
# (unlike normal_lcdf, which received it) was never updated with the
# more accurate Mills-ratio approximation, so it still hits this floor.
# Confirmed directly: normal_lccdf(9.5|0,1) returns -inf here vs the true
# -48.3. Used instead: P(Z>x) = 0.5*erfc(x/(sigma*sqrt(2))), confirmed to
# match R's pnorm() to high precision out to 30+ SDs (erfc() does not
# share normal_lccdf's accuracy gap).
dnorm1_stan_funs <- "
  real dnorm1_lpmf(int z, real sigma) {
    if (z >= 0) {
      real lo = (z - 0.5) / (sigma * sqrt2());
      real hi = (z + 0.5) / (sigma * sqrt2());
      return log(0.5) + log_diff_exp(log(erfc(lo)), log(erfc(hi)));
    }
    return log_diff_exp(normal_lcdf(z + 0.5 | 0, sigma), normal_lcdf(z - 0.5 | 0, sigma));
  }
"

# Stan function block for the discrete-normal log-CCDF (truncation
# support via brms's resp_trunc()). NOT normal_lccdf directly -- same
# documented Stan limitation as dnorm1_lpmf above (stan-dev/math#1985);
# see there for the citation. Exact closed form via erfc() instead.
dnorm1_lccdf_stan <- "
  real dnorm1_lccdf(int y, real sigma) {
    real z = (y + 0.5) / (sigma * sqrt2());
    return log(0.5) + log(erfc(z));
  }
"

# Stan function block for the discrete-normal log-PMF with free location
# AND free scale -- no constraint coupling mu and sigma, the same
# structural contrast with skellam2 already documented for dlaplace2 (see
# dlaplace2_stan_funs above): the point of having both an asymmetric-
# Skellam and free-location discrete families is to compare a model where
# bias and spread are structurally coupled against ones where they are
# not, so no constraint is imposed here either. mu enters by re-centring
# before the erfc() argument, since erfc() itself takes no location/scale.
#
# Same z>=mu vs z<mu branch as dnorm1_lpmf above (centred on mu rather
# than 0: the cancellation risk is about which side of the *mean* z
# falls on, not the sign of z itself), and the same erfc()-based exact
# survival form in place of the documented-broken normal_lccdf
# (stan-dev/math#1985, see dnorm1_lpmf above for the citation).
dnorm2_stan_funs <- "
  real dnorm2_lpmf(int z, real mu, real sigma) {
    if (z >= mu) {
      real lo = (z - mu - 0.5) / (sigma * sqrt2());
      real hi = (z - mu + 0.5) / (sigma * sqrt2());
      return log(0.5) + log_diff_exp(log(erfc(lo)), log(erfc(hi)));
    }
    return log_diff_exp(normal_lcdf(z + 0.5 | mu, sigma), normal_lcdf(z - 0.5 | mu, sigma));
  }
"

# Stan function block for the discrete-normal log-CCDF, free location and
# scale. Same erfc()-based exact form as dnorm1_lccdf_stan (not
# normal_lccdf -- stan-dev/math#1985, see dnorm1_lpmf above), generalised
# to a nonzero mu via re-centring before the erfc() argument.
dnorm2_lccdf_stan <- "
  real dnorm2_lccdf(int y, real mu, real sigma) {
    real z = (y + 0.5 - mu) / (sigma * sqrt2());
    return log(0.5) + log(erfc(z));
  }
"
