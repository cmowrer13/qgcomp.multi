# Fit the core outcome-model and MSM components for two-mixture g-computation

Fits the outcome regression, computes predicted potential outcomes under
joint interventions on two exposure mixtures, and estimates the marginal
structural model (MSM) coefficients that summarize the resulting
two-dimensional intervention-response surface.

## Usage

``` r
qgcompmulti_msm_fit(
  f,
  data,
  mix1,
  mix2,
  interaction = TRUE,
  family = gaussian(),
  q = 4,
  centering = "none",
  id = NULL,
  MCsize = nrow(data),
  seed = NULL
)
```

## Arguments

- f:

  A model formula for the outcome regression. The formula should include
  the outcome and any baseline covariates. Mixture variables listed in
  `mix1` and `mix2` should also appear in the formula if they are to be
  included in the outcome model.

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
  implementation, the outcome model is augmented with a cross-product
  between the sums of the components in `mix1` and `mix2`, and the MSM
  includes the `psi1 * psi2` interaction term. If `FALSE`, both models
  are fit without that interaction.

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
  `"none"` or `"median"`. Centering affects only the MSM predictors and
  does not change the outcome regression fit. This argument is ignored
  when `q` is numeric.

- id:

  Optional character string giving the name of a cluster identifier
  variable. If supplied, Monte Carlo subsampling is performed at the
  cluster level rather than the observation level.

- MCsize:

  Optional integer controlling the Monte Carlo sample size used to
  approximate the marginalization step in g-computation. If `MCsize` is
  greater than or equal to the current analysis sample size, all
  observations are used. Smaller values compute predicted outcomes over
  a random subsample drawn from the empirical distribution. When `id` is
  supplied and `MCsize < nrow(data)`, the approximation is implemented
  by sampling `MCsize` clusters with replacement.

- seed:

  Optional integer random seed used to make Monte Carlo subsampling
  reproducible when `MCsize < nrow(data)`. If `NULL`, the current RNG
  state is used and not modified by `qgcompmulti_msm_fit()`.

## Value

A list with components:

- `outcome_fit`:

  The fitted outcome regression model object.

- `msm_fit`:

  The fitted marginal structural model object.

- `coefficients`:

  A named vector of MSM coefficients.

- `n_used`:

  The number of observations used in the g-computation prediction step.

- `intervention_grid`:

  The fit-time intervention grid on the intervention-value scale.

- `msm_grid`:

  The corresponding fit-time grid on the MSM coding scale.

- `counterfactual_surface`:

  The exact fit-time counterfactual mean surface implied by the fitted
  outcome model.

- `msm_surface`:

  The fitted MSM surface evaluated on the common fit-time grid.

- `surface_comparison`:

  A direct exact-versus-MSM comparison object on the common fit-time
  grid.

## Details

This is a lower-level fitting helper used internally by
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md).
It is exported because it can be useful for method development, testing,
and direct inspection of the fitted outcome model, MSM, and stored
fit-time surfaces. Most users will want
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
instead.

This function carries out the core g-computation step for the
two-mixture extension. It first fits a generalized linear model for the
outcome using [`glm()`](https://rdrr.io/r/stats/glm.html). It then
constructs a grid of intervention levels over the two mixtures, replaces
the observed mixture values with the intervention values, and computes
predicted outcomes for each observation (or for a Monte Carlo subsample
if `MCsize < nrow(data)`). These predicted outcomes are stacked into a
pseudo-dataset and used to fit a marginal structural model that
summarizes the dose-response surface. Note that predicted potential
outcomes are computed on the response scale and the marginal structural
model is fit using an identity link, regardless of the outcome model.

When `q` is an integer, the intervention grid is `0, 1, ..., q - 1` for
each mixture, corresponding to simultaneous quantile increases in all
components of that mixture.

When `q = NULL`, each mixture is instead set to common pooled percentile
values on the original exposure scale. The resulting MSM coefficients
are therefore scale-dependent and should be interpreted in the units of
the underlying exposures. If `centering = "median"`, the intercept
corresponds to the pooled-median intervention for both mixtures.

Because the number of intervention combinations grows as `q^2`, this
step can become computationally expensive for large datasets. The
`MCsize` argument reduces this burden by approximating the empirical
covariate distribution using a random subset while leaving the outcome
model fit unchanged.

## Examples

``` r
dat <- sim_mixture_data(
  n = 500,
  pA = 3,
  pB = 3,
  rho_within_A = 0.3,
  rho_within_B = 0.3,
  rho_between = 0.2,
  psi1 = 0.5,
  psi2 = 0.3,
  psi12 = 0.2,
  return_quantized = FALSE,
  seed = 123
)

qgcompmulti_msm_fit(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  MCsize = nrow(dat),
  seed = 13
)
#> $outcome_fit
#> 
#> Call:  glm(formula = outcome_formula, family = family, data = data)
#> 
#> Coefficients:
#>                        (Intercept)                                  X1  
#>                            1.69913                             0.17268  
#>                                 X2                                  X3  
#>                            0.31351                             0.28412  
#>                                 W1                                  W2  
#>                            0.22104                             0.26571  
#>                                 W3                                   C  
#>                            0.15212                             1.03490  
#> I((X1 + X2 + X3) * (W1 + W2 + W3))  
#>                            0.01312  
#> 
#> Degrees of Freedom: 499 Total (i.e. Null);  491 Residual
#> Null Deviance:       1370 
#> Residual Deviance: 522.2     AIC: 1461
#> 
#> $msm_fit
#> 
#> Call:  glm(formula = msm_formula, data = msmdat)
#> 
#> Coefficients:
#> (Intercept)         psi1         psi2    psi1:psi2  
#>      1.6757       0.7703       0.6389       0.1181  
#> 
#> Degrees of Freedom: 7999 Total (i.e. Null);  7996 Residual
#> Null Deviance:       15810 
#> Residual Deviance: 3.221e-24     AIC: -481900
#> 
#> $coefficients
#> (Intercept)        psi1        psi2   psi1:psi2 
#>   1.6757274   0.7703144   0.6388732   0.1180513 
#> 
#> $n_used
#> [1] 500
#> 
#> $intervention_grid
#>    grid_id psi1 psi2
#> 1        1    0    0
#> 2        2    1    0
#> 3        3    2    0
#> 4        4    3    0
#> 5        5    0    1
#> 6        6    1    1
#> 7        7    2    1
#> 8        8    3    1
#> 9        9    0    2
#> 10      10    1    2
#> 11      11    2    2
#> 12      12    3    2
#> 13      13    0    3
#> 14      14    1    3
#> 15      15    2    3
#> 16      16    3    3
#> 
#> $msm_grid
#>    grid_id psi1 psi2
#> 1        1    0    0
#> 2        2    1    0
#> 3        3    2    0
#> 4        4    3    0
#> 5        5    0    1
#> 6        6    1    1
#> 7        7    2    1
#> 8        8    3    1
#> 9        9    0    2
#> 10      10    1    2
#> 11      11    2    2
#> 12      12    3    2
#> 13      13    0    3
#> 14      14    1    3
#> 15      15    2    3
#> 16      16    3    3
#> 
#> $counterfactual_surface
#>    grid_id intervention_psi1 intervention_psi2 msm_psi1 msm_psi2 exact_mean
#> 1        1                 0                 0        0        0   1.675727
#> 2        2                 1                 0        1        0   2.446042
#> 3        3                 2                 0        2        0   3.216356
#> 4        4                 3                 0        3        0   3.986671
#> 5        5                 0                 1        0        1   2.314601
#> 6        6                 1                 1        1        1   3.202966
#> 7        7                 2                 1        2        1   4.091332
#> 8        8                 3                 1        3        1   4.979698
#> 9        9                 0                 2        0        2   2.953474
#> 10      10                 1                 2        1        2   3.959891
#> 11      11                 2                 2        2        2   4.966308
#> 12      12                 3                 2        3        2   5.972725
#> 13      13                 0                 3        0        3   3.592347
#> 14      14                 1                 3        1        3   4.716815
#> 15      15                 2                 3        2        3   5.841284
#> 16      16                 3                 3        3        3   6.965752
#> 
#> $msm_surface
#>    grid_id intervention_psi1 intervention_psi2 msm_psi1 msm_psi2 msm_mean
#> 1        1                 0                 0        0        0 1.675727
#> 2        2                 1                 0        1        0 2.446042
#> 3        3                 2                 0        2        0 3.216356
#> 4        4                 3                 0        3        0 3.986671
#> 5        5                 0                 1        0        1 2.314601
#> 6        6                 1                 1        1        1 3.202966
#> 7        7                 2                 1        2        1 4.091332
#> 8        8                 3                 1        3        1 4.979698
#> 9        9                 0                 2        0        2 2.953474
#> 10      10                 1                 2        1        2 3.959891
#> 11      11                 2                 2        2        2 4.966308
#> 12      12                 3                 2        3        2 5.972725
#> 13      13                 0                 3        0        3 3.592347
#> 14      14                 1                 3        1        3 4.716815
#> 15      15                 2                 3        2        3 5.841284
#> 16      16                 3                 3        3        3 6.965752
#> 
#> $surface_comparison
#>    grid_id intervention_psi1 intervention_psi2 msm_psi1 msm_psi2 exact_mean
#> 1        1                 0                 0        0        0   1.675727
#> 2        2                 1                 0        1        0   2.446042
#> 3        3                 2                 0        2        0   3.216356
#> 4        4                 3                 0        3        0   3.986671
#> 5        5                 0                 1        0        1   2.314601
#> 6        6                 1                 1        1        1   3.202966
#> 7        7                 2                 1        2        1   4.091332
#> 8        8                 3                 1        3        1   4.979698
#> 9        9                 0                 2        0        2   2.953474
#> 10      10                 1                 2        1        2   3.959891
#> 11      11                 2                 2        2        2   4.966308
#> 12      12                 3                 2        3        2   5.972725
#> 13      13                 0                 3        0        3   3.592347
#> 14      14                 1                 3        1        3   4.716815
#> 15      15                 2                 3        2        3   5.841284
#> 16      16                 3                 3        3        3   6.965752
#>    msm_mean     residual
#> 1  1.675727 2.109424e-14
#> 2  2.446042 2.087219e-14
#> 3  3.216356 1.998401e-14
#> 4  3.986671 2.042810e-14
#> 5  2.314601 2.087219e-14
#> 6  3.202966 2.087219e-14
#> 7  4.091332 1.953993e-14
#> 8  4.979698 1.953993e-14
#> 9  2.953474 2.042810e-14
#> 10 3.959891 2.042810e-14
#> 11 4.966308 1.953993e-14
#> 12 5.972725 2.042810e-14
#> 13 3.592347 1.998401e-14
#> 14 4.716815 1.953993e-14
#> 15 5.841284 1.865175e-14
#> 16 6.965752 1.865175e-14
#> 
```
