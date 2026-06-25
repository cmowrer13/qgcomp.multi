# Fit a quantile g-computation model for two exposure mixtures

Fits a two-mixture g-computation model with an optional interaction
between the mixture intervention variables. The function fits an outcome
regression, computes predicted potential outcomes under joint
interventions on the two mixtures, and then fits a marginal structural
model (MSM) to summarize the resulting intervention-response surface.
The current interface is designed for analyses with exactly two
mixtures, supplied through `mix1` and `mix2`. Uncertainty is estimated
using a nonparametric bootstrap.

## Usage

``` r
qgcomp.glm.multi(
  f,
  data,
  mix1,
  mix2,
  interaction = TRUE,
  family = gaussian(),
  q = 4,
  centering = "none",
  B = 200,
  id = NULL,
  MCsize = nrow(data),
  seed = NULL,
  progress = FALSE,
  parallel = FALSE,
  workers = NULL
)
```

## Arguments

- f:

  A model formula for the outcome regression. The formula should include
  the outcome and any baseline covariates. Mixture variables listed in
  `mix1` and `mix2` should also appear in the formula if they are to be
  included in the outcome model. The `mix1` and `mix2` arguments define
  how exposures are grouped into mixtures, but they do not add variables
  to the formula automatically.

- data:

  A data frame containing the outcome, exposure variables, and any
  covariates in the model.

- mix1:

  A character vector giving the names of the variables in the first
  exposure mixture.

- mix2:

  A character vector giving the names of the variables in the second
  exposure mixture.

- interaction:

  Logical; if `TRUE`, the package includes an interaction term in both
  the outcome regression and the fitted MSM. In the current
  implementation, the outcome model is augmented with a product between
  the sums of the components in `mix1` and `mix2`, and the MSM includes
  the `psi1 * psi2` interaction term. If `FALSE`, both models are fit
  without that interaction.

- family:

  A GLM family object (e.g.,
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`poisson()`](https://rdrr.io/r/stats/family.html)) specifying the
  outcome model.

- q:

  Integer greater than or equal to 2 giving the number of quantiles used
  to discretize the exposure variables, or `NULL` to skip quantization
  and fit the outcome model on the original exposure scale. When
  `q = NULL`, the fit-time intervention grid is defined by the pooled
  25th, 50th, and 75th percentile values within each mixture, and under
  each intervention every component in a mixture is set to the same
  pooled mixture-specific value.

- centering:

  Character string controlling how the marginal structural model
  intervention variables are coded when `q = NULL`. Must be one of
  `"none"` or `"median"`. With `"none"`, the MSM uses the raw
  intervention values. With `"median"`, the MSM uses intervention values
  centered at the pooled median within each mixture. Centering affects
  only the MSM fit, not the outcome regression. This argument is ignored
  when `q` is numeric.

- B:

  Integer greater than or equal to 2 giving the number of bootstrap
  replications used for standard error estimation.

- id:

  Optional character string giving the name of a cluster identifier
  variable. If supplied, bootstrap resampling is performed at the
  cluster level rather than the observation level.

- MCsize:

  Optional integer controlling the Monte Carlo sample size used in the
  g-computation step. If `MCsize` is greater than or equal to the
  current analysis sample size, the empirical covariate distribution is
  fully enumerated. Smaller values approximate the marginalization step
  using a random subsample, which can reduce computation time in large
  datasets. When `id` is supplied and `MCsize < nrow(data)`, the
  approximation is implemented by sampling `MCsize` clusters with
  replacement.

- seed:

  Optional integer master seed used to make the main fit, bootstrap
  resampling, and any Monte Carlo subsampling reproducible. When
  supplied, the function deterministically derives one full-fit seed and
  one bootstrap-worker seed per replication. If `NULL`, the current RNG
  state is used and not modified by `qgcomp.glm.multi()`.

- progress:

  Logical; if `TRUE`, display a compact single-line bootstrap progress
  bar in the console while the model is fitting. The display reports
  percent complete, attempted bootstrap replications, elapsed time, and
  an estimated time remaining based on the running average bootstrap
  iteration time. Off by default. Serial mode only; if
  `parallel = TRUE`, progress is disabled with an explicit warning.

- parallel:

  Logical; if `TRUE`, dispatch bootstrap replications through
  `future.apply`. Parallelism is limited to the bootstrap replication
  layer.

- workers:

  Optional integer worker count for the internal
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html)
  path. Leave `NULL` to use the active non-sequential `future` plan when
  one is already set, or otherwise let `qgcomp.glm.multi()` choose a
  temporary local worker count automatically. If a non-sequential
  `future` plan is already active, supplying an explicit `workers` value
  is treated as an unsupported combination and errors clearly.

## Value

An object of class `"qgcompmulti"` representing the fitted two-mixture
quantile g-computation model. Major components include:

- `data_info`:

  Outcome name, sample-size metadata, and indicators for quantization
  and clustered fitting.

- `mixtures`:

  The mixture definitions, quantization setting `q`, and original-scale
  centering choice when `q = NULL`.

- `analysis`:

  Model settings such as the GLM family, interaction status, bootstrap
  count, cluster identifier, Monte Carlo size, and any supplied random
  seed.

- `fits`:

  The fitted outcome regression and marginal structural model (MSM)
  objects.

- `prediction`:

  Stored fit-time prediction objects for later prediction, plotting, and
  diagnostic methods, including the intervention grid, the MSM-coded
  grid, the exact counterfactual surface, the corresponding MSM fitted
  surface, and a comparison object on the common fit-time grid.

- `bootstrap`:

  Retained bootstrap coefficient draws, bootstrap replication counts and
  lightweight failure data.

- `results`:

  The MSM coefficient vector, standard errors, covariance matrix, and
  coefficient table.

- `labels`:

  Internal coefficient names and human-readable labels used by the print
  and summary methods.

## Details

The fitted MSM has the form

\$\$ E\[Y^{x(q_1), w(q_2)}\] = \psi_0 + \psi_1 q_1 + \psi_2 q_2 +
\psi\_{12} q_1 q_2 \$\$

when `interaction = TRUE`. When `q` is an integer, `q1` and `q2` index
quantized joint intervention levels for mixtures 1 and 2. When
`q = NULL`, the same linear MSM is fit over original-scale intervention
values using a `3 x 3` grid formed from the pooled 25th, 50th, and 75th
percentile values within each mixture. In that setting, every component
in a mixture is set to the same pooled mixture-specific value under a
given intervention.

This function extends quantile g-computation to two exposure mixtures by
evaluating predicted outcomes over a two-dimensional intervention grid.
For each bootstrap replication, the observed data are resampled,
exposures are either quantized or left on their original scale depending
on `q`, the outcome model is fit, and predicted potential outcomes are
computed under joint interventions on the two mixtures. A marginal
structural model is then fit to those predicted counterfactual means to
obtain the reported mixture effect estimates.

When `q` is an integer, the intervention grid is `0, 1, ..., q - 1` for
each mixture. In that setting, `psi1` and `psi2` are interpreted as the
change in the marginal mean outcome associated with simultaneously
increasing every component in the corresponding mixture by one quantile,
holding the other mixture intervention level fixed at the lowest
quantile.

When `q = NULL`, the exposure variables are left on their original
analysis scales. The intervention grid is then defined by assigning
every component in a mixture to a common pooled percentile value from
that mixture. In this setting, the MSM coefficients are defined with
respect to that original-scale intervention coding, so their units
depend on the measurement scale of the underlying exposures. If
`centering = "median"`, the MSM is fit on centered intervention values
and the intercept corresponds to the pooled-median intervention for both
mixtures.

If `interaction = TRUE`, the current implementation adds an interaction
term to both the outcome regression and the MSM. This means the MSM is
summarizing a counterfactual surface implied by an interacting outcome
model rather than simply adding an interaction at the final summary
step.

The outcome model is fit using
[`glm()`](https://rdrr.io/r/stats/glm.html), so this function can be
used with Gaussian, binomial, Poisson, and other generalized linear
models supported by the supplied formula and family specification.
Predicted potential outcomes are computed on the response scale and the
MSM is fit using an identity link, regardless of the outcome-model link
function.

Interpretation therefore depends on the outcome type:

- For continuous outcomes, parameters represent mean differences.

- For binary outcomes, parameters represent differences in predicted
  probabilities, that is, risk differences.

- For count outcomes, parameters represent differences in expected
  counts.

If `progress = TRUE`, the bootstrap loop prints a compact single-line
status display in serial mode. The failed-replicate counter is shown
only after the first failed bootstrap iteration, so clean runs do not
carry extra visual noise. Parallel execution is intentionally limited to
one level in Version `0.4.0`, so requesting `progress = TRUE` together
with `parallel = TRUE` disables the progress display with an explicit
warning.

When `parallel = TRUE`, the bootstrap replications are dispatched with
`future.apply`. If no non-sequential `future` plan is already active and
`workers` is left `NULL`, the function uses a temporary local
[`future::multisession`](https://future.futureverse.org/reference/multisession.html)
path. Reproducibility is defined within a fixed backend and execution
mode. The function treats `seed` as a master seed and deterministically
expands it into one full-fit seed plus one bootstrap-worker seed per
replication.

For causal interpretation, the usual identifying conditions for
g-computation still apply: consistency, conditional exchangeability,
positivity for the interventions under study, and adequate specification
of the outcome model. The support diagnostic can help users inspect the
intervention grid, but it should not be read as a full positivity proof.

## Examples

``` r
dat <- sim_mixture_data(
  n = 500,
  pA = 4,
  pB = 4,
  rho_within_A = 0.3,
  rho_within_B = 0.3,
  rho_between = 0.2,
  psi1 = 0.5,
  psi2 = 0.3,
  psi12 = 0.2,
  return_quantized = FALSE,
  seed = 123
)

fit <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3", "X4"),
  mix2 = c("W1", "W2", "W3", "W4"),
  interaction = TRUE,
  q = 4,
  B = 100,
  MCsize = nrow(dat),
  seed = 13
)

fit
#> qgcompmulti fit
#> 
#> Call:
#> qgcomp.glm.multi(f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + 
#>     C, data = dat, mix1 = c("X1", "X2", "X3", "X4"), mix2 = c("W1", 
#>     "W2", "W3", "W4"), interaction = TRUE, q = 4, B = 100, MCsize = nrow(dat), 
#>     seed = 13)
#> 
#> Model:
#>   Outcome: Y
#>   Family: gaussian (identity)
#>   Observations used: 500
#>   Exposure mode: Quantized exposures (q = 4)
#>   MSM interaction: included
#>   Random seed: 13
#> 
#> Mixtures:
#>   Mixture 1: X1, X2, X3, X4
#>   Mixture 2: W1, W2, W3, W4
#> 
#> MSM coefficients:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> Intercept             -0.22367    0.18580 -1.2038  0.228680    
#> Mixture 1 main effect  0.59476    0.12982  4.5816 4.614e-06 ***
#> Mixture 2 main effect  0.38901    0.13640  2.8520  0.004344 ** 
#> Mixture interaction    0.15243    0.08085  1.8854  0.059379 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(fit)
#> Summary of qgcompmulti fit
#> 
#> Call:
#> qgcomp.glm.multi(f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + 
#>     C, data = dat, mix1 = c("X1", "X2", "X3", "X4"), mix2 = c("W1", 
#>     "W2", "W3", "W4"), interaction = TRUE, q = 4, B = 100, MCsize = nrow(dat), 
#>     seed = 13)
#> 
#> Model overview:
#>   Formula: Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C
#>   Outcome: Y
#>   Family: gaussian (identity)
#>   Observations used: 500
#>   Exposure mode: Quantized exposures (q = 4)
#>   MSM interaction: included
#>   Bootstrap replications: 100
#>   Monte Carlo size: 500
#>   Random seed: 13
#> 
#> Mixtures:
#>   Mixture 1: X1, X2, X3, X4
#>   Mixture 2: W1, W2, W3, W4
#> 
#> MSM coefficients:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> Intercept             -0.22367    0.18580 -1.2038  0.228680    
#> Mixture 1 main effect  0.59476    0.12982  4.5816 4.614e-06 ***
#> Mixture 2 main effect  0.38901    0.13640  2.8520  0.004344 ** 
#> Mixture interaction    0.15243    0.08085  1.8854  0.059379 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Outcome model context:
#>   Model class: glm
#>   Estimated parameters: 11
#>   AIC: 1425.626
#>   Null deviance: 1325.786
#>   Residual deviance: 482.984
coef(fit)
#> (Intercept)        psi1        psi2   psi1:psi2 
#>  -0.2236658   0.5947622   0.3890148   0.1524334 
confint(fit)
#>                    2.5 %    97.5 %
#> (Intercept) -0.587836841 0.1405052
#> psi1         0.340329402 0.8491951
#> psi2         0.121677162 0.6563524
#> psi1:psi2   -0.006030253 0.3108970

# Optional bootstrap-level parallelized computation
parallel_fit <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3", "X4"),
  mix2 = c("W1", "W2", "W3", "W4"),
  interaction = TRUE,
  q = 4,
  B = 100,
  MCsize = nrow(dat),
  seed = 13,
  parallel = TRUE,
  workers = 2
)

# Public prediction and plotting workflow
if (FALSE) { # \dontrun{
predict(fit)
predict(
  fit,
  type = "msm_contrast",
  from = c(psi1 = 0, psi2 = 0),
  to = c(psi1 = 3, psi2 = 3),
  interval = TRUE
)
plot(fit)
plot(fit, style = "contour")

# Diagnostics
support(fit)
diagnostics(fit, type = "bootstrap")
adequacy(fit)

# Sensitivity helpers
mcsize_sensitivity(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  MCsize_values = c(250, 500),
  q = 4,
  B = 100,
  seed = 13
)
} # }

dat_cont <- sim_mixture_data(
  n = 500,
  pA = 4,
  pB = 4,
  rho_within_A = 0.3,
  rho_within_B = 0.3,
  rho_between = 0.2,
  psi1 = 0.5,
  psi2 = 0.3,
  psi12 = 0.2,
  return_quantized = FALSE,
  seed = 321
)

fit_cont <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C,
  data = dat_cont,
  mix1 = c("X1", "X2", "X3", "X4"),
  mix2 = c("W1", "W2", "W3", "W4"),
  interaction = TRUE,
  q = NULL,
  centering = "median",
  B = 100,
  MCsize = nrow(dat_cont),
  seed = 13
)

fit_cont
#> qgcompmulti fit
#> 
#> Call:
#> qgcomp.glm.multi(f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + 
#>     C, data = dat_cont, mix1 = c("X1", "X2", "X3", "X4"), mix2 = c("W1", 
#>     "W2", "W3", "W4"), interaction = TRUE, q = NULL, centering = "median", 
#>     B = 100, MCsize = nrow(dat_cont), seed = 13)
#> 
#> Model:
#>   Outcome: Y
#>   Family: gaussian (identity)
#>   Observations used: 500
#>   Exposure mode: Original-scale exposures (centering = "median")
#>   MSM interaction: included
#>   Random seed: 13
#> 
#> Mixtures:
#>   Mixture 1: X1, X2, X3, X4
#>   Mixture 2: W1, W2, W3, W4
#> 
#> MSM coefficients:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> Intercept             1.701759   0.082863 20.5369 < 2.2e-16 ***
#> Mixture 1 main effect 0.827489   0.085637  9.6627 < 2.2e-16 ***
#> Mixture 2 main effect 0.615797   0.076362  8.0641 7.376e-16 ***
#> Mixture interaction   0.291838   0.094064  3.1026  0.001919 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(fit_cont)
#> Summary of qgcompmulti fit
#> 
#> Call:
#> qgcomp.glm.multi(f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + 
#>     C, data = dat_cont, mix1 = c("X1", "X2", "X3", "X4"), mix2 = c("W1", 
#>     "W2", "W3", "W4"), interaction = TRUE, q = NULL, centering = "median", 
#>     B = 100, MCsize = nrow(dat_cont), seed = 13)
#> 
#> Model overview:
#>   Formula: Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C
#>   Outcome: Y
#>   Family: gaussian (identity)
#>   Observations used: 500
#>   Exposure mode: Original-scale exposures (centering = "median")
#>   MSM interaction: included
#>   Bootstrap replications: 100
#>   Monte Carlo size: 500
#>   Random seed: 13
#> 
#> Mixtures:
#>   Mixture 1: X1, X2, X3, X4
#>   Mixture 2: W1, W2, W3, W4
#> 
#> MSM coefficients:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> Intercept             1.701759   0.082863 20.5369 < 2.2e-16 ***
#> Mixture 1 main effect 0.827489   0.085637  9.6627 < 2.2e-16 ***
#> Mixture 2 main effect 0.615797   0.076362  8.0641 7.376e-16 ***
#> Mixture interaction   0.291838   0.094064  3.1026  0.001919 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Outcome model context:
#>   Model class: glm
#>   Estimated parameters: 11
#>   AIC: 1496.806
#>   Null deviance: 1363.142
#>   Residual deviance: 556.876
coef(fit_cont)
#> (Intercept)        psi1        psi2   psi1:psi2 
#>   1.7017593   0.8274885   0.6157973   0.2918381 
confint(fit_cont)
#>                 2.5 %    97.5 %
#> (Intercept) 1.5393502 1.8641685
#> psi1        0.6596428 0.9953343
#> psi2        0.4661296 0.7654650
#> psi1:psi2   0.1074769 0.4761994
```
