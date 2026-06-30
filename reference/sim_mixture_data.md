# Generate simulated dataset with two exposure mixtures

Generates a dataset for simulation studies involving two mixtures of
correlated exposures and their interaction. Exposures are generated from
a multivariate normal distribution with user-specified within-mixture
and between-mixture correlations, then discretized into quantiles. The
outcome is generated from family- and estimand-scale-specific models
where the parameters correspond to simultaneous one-quantile increases
in each mixture.

## Usage

``` r
sim_mixture_data(
  n,
  pA,
  pB,
  rho_within_A,
  rho_within_B,
  rho_between,
  psi1,
  psi2,
  psi12,
  sigma_eps = 1,
  q = 4,
  return_quantized = FALSE,
  seed = NULL,
  family = stats::gaussian(),
  estimand_scale = NULL,
  baseline_mean = 0,
  baseline_risk = 0.2,
  baseline_rate = 1,
  beta_C = 1
)
```

## Arguments

- n:

  Number of observations.

- pA:

  Number of components in mixture A.

- pB:

  Number of components in mixture B.

- rho_within_A:

  Correlation between components within mixture A.

- rho_within_B:

  Correlation between components within mixture B.

- rho_between:

  Correlation between components across mixtures.

- psi1:

  Effect of a one-quantile increase in mixture A.

- psi2:

  Effect of a one-quantile increase in mixture B.

- psi12:

  Interaction effect representing the change in the effect of mixture A
  when mixture B increases by one quantile.

- sigma_eps:

  Standard deviation of the outcome error term for Gaussian generation.

- q:

  Number of quantiles used to discretize exposures.

- return_quantized:

  Logical; if TRUE returns quantized exposures, otherwise returns the
  original continuous exposures.

- seed:

  Optional random seed for reproducibility.

- family:

  Outcome family. Supported families are
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html), and
  [`poisson()`](https://rdrr.io/r/stats/family.html).

- estimand_scale:

  Optional outcome-generating estimand scale. Supported values mirror
  the fitting API: `"mean_difference"`, `"risk_difference"`,
  `"odds_ratio"`, and `"rate_ratio"`. If omitted, defaults follow the
  `0.5.0` fitting conventions.

- baseline_mean:

  Baseline mean used for Gaussian generation.

- baseline_risk:

  Baseline risk used for binomial generation.

- baseline_rate:

  Baseline rate used for Poisson generation.

- beta_C:

  Effect of covariate `C` in the outcome-generating model.

## Value

A data frame containing the simulated outcome `Y`, mixture exposures
(`X1...XpA` and `W1...WpB`), and covariate `C`.

## Details

The exposure summaries are:

\$\$S_A = \textrm{rowMeans}(X_q), \qquad S_B =
\textrm{rowMeans}(W_q)\$\$

where `X_q` and `W_q` are the quantized mixture components. The outcome
is then generated from one of the following models:

- Gaussian mean-difference generation: \$\$Y = \mu_0 + \psi_1 S_A +
  \psi_2 S_B + \psi\_{12} S_A S_B + \beta_C C + \epsilon\$\$

- Binomial odds-ratio generation: \$\$\textrm{logit}\\P(Y = 1)\\ =
  \textrm{logit}(p_0) + \log(\psi_1) S_A + \log(\psi_2) S_B +
  \log(\psi\_{12}) S_A S_B + \beta_C C\$\$

- Binomial risk-difference generation: \$\$P(Y = 1) = p_0 + \psi_1 S_A +
  \psi_2 S_B + \psi\_{12} S_A S_B + \beta_C C\$\$

- Poisson rate-ratio generation: \$\$\log\\E(Y)\\ = \log(\lambda_0) +
  \log(\psi_1) S_A + \log(\psi_2) S_B + \log(\psi\_{12}) S_A S_B +
  \beta_C C\$\$

- Poisson mean-difference generation: \$\$E(Y) = \lambda_0 + \psi_1
  S_A + \psi_2 S_B + \psi\_{12} S_A S_B + \beta_C C\$\$

Ratio-valued inputs are accepted on their natural scale. For example,
`psi1 = 1.5` with `family = poisson()` and
`estimand_scale = "rate_ratio"` means that a one-quantile increase in
mixture A corresponds to a rate ratio of `1.5` before any interaction
modification by mixture B.

Continuous exposures are generated from a multivariate normal
distribution with a block correlation structure. Exposures are then
discretized into quantile categories ranging from `0` to `q - 1`. The
outcome is always generated from the quantized mixture summaries so that
the requested parameters `psi1`, `psi2`, and `psi12` correspond directly
to the intended one-quantile mixture effects on the chosen generating
scale.

Additive binary-risk and additive Poisson-mean generation are supported,
but the requested baseline and effect values must imply valid
probabilities or strictly positive means for all simulated observations.

## Examples

``` r
dat_gaussian <- sim_mixture_data(
  n = 1000,
  pA = 4,
  pB = 4,
  rho_within_A = 0.5,
  rho_within_B = 0.5,
  rho_between = 0.3,
  psi1 = 0.5,
  psi2 = 0.3,
  psi12 = 0.2,
  return_quantized = TRUE,
  seed = 123
)

dat_binomial <- sim_mixture_data(
  n = 1000,
  pA = 3,
  pB = 3,
  rho_within_A = 0.4,
  rho_within_B = 0.4,
  rho_between = 0.2,
  psi1 = 1.5,
  psi2 = 1.2,
  psi12 = 1.1,
  family = binomial(),
  baseline_risk = 0.2,
  beta_C = 0.25,
  seed = 123
)

dat_poisson <- sim_mixture_data(
  n = 1000,
  pA = 3,
  pB = 3,
  rho_within_A = 0.4,
  rho_within_B = 0.4,
  rho_between = 0.2,
  psi1 = 1.5,
  psi2 = 1.2,
  psi12 = 1.1,
  family = poisson(),
  baseline_rate = 1,
  beta_C = 0.15,
  seed = 123
)
```
