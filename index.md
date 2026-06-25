# qgcomp.multi

`qgcomp.multi` extends quantile g-computation to settings with **two
exposure mixtures** and their interaction. It is designed for analyses
where exposures come in distinct groups and the effect of shifting one
mixture may depend on the level of another.

The package fits an outcome regression model, computes predicted
counterfactual means under joint interventions on the two mixtures, and
then fits a marginal structural model (MSM) to summarize the resulting
intervention-response surface.

## What the package does

`qgcomp.multi` currently supports:

- two-mixture estimation through
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
- optional interaction between the two mixture intervention variables
- quantized analyses with `q >= 2`
- original-scale analyses with `q = NULL`
- MSM-based prediction, plotting, diagnostics, and adequacy checks
- sensitivity helpers for Monte Carlo size, bootstrap size, and
  quantization
- native multiple-imputation pooling through
  [`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
- optional bootstrap-level parallel execution
- broom-style `tidy()` and `glance()` methods for single-fit and pooled
  MI objects

## Installation

You can install the development version from GitHub:

``` r

devtools::install_github("cmowrer13/qgcomp.multi")
```

## Main quantized workflow

For a standard quantized analysis, fit a two-mixture model with
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md):

``` r

library(qgcomp.multi)

fit <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dataset,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 100,
  seed = 13
)

summary(fit)
coef(fit)
confint(fit)
```

In this setting, the MSM coefficients summarize the expected outcome
change associated with simultaneous one-quantile shifts in the component
exposures of each mixture.

## Original-scale fitting with `q = NULL`

If quantization is not scientifically appropriate for the analysis, set
`q = NULL` to fit the outcome model on the original analysis scale and
build the MSM over a pooled `3 x 3` intervention grid within each
mixture.

``` r

fit_cont <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dataset,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = NULL,
  centering = "median",
  B = 100,
  seed = 13
)

summary(fit_cont)
coef(fit_cont)
confint(fit_cont)
```

With `q = NULL`, the MSM coefficients are not one-quantile effects. They
are defined with respect to the original-scale intervention coding used
to fit the MSM, so interpretation should stay tied to the actual
intervention grid.

## Model interpretation

The main estimands are the coefficients of the fitted marginal
structural model:

- `psi1`: Mixture 1 main effect
- `psi2`: Mixture 2 main effect
- `psi1:psi2`: interaction between the two mixtures when
  `interaction = TRUE`

These are summaries of the intervention-response surface implied by the
fitted outcome model. For causal interpretation, the usual g-computation
conditions still matter: consistency, exchangeability, positivity, and
adequate outcome model specification.

## Beyond the basic fit

The package also includes tools for:

- MSM prediction and direct contrasts with
  [`predict()`](https://rdrr.io/r/stats/predict.html)
- surface visualization with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
- intervention support, bootstrap, and adequacy diagnostics
- sensitivity analysis with
  [`mcsize_sensitivity()`](https://cmowrer13.github.io/qgcomp.multi/reference/mcsize_sensitivity.md),
  [`b_sensitivity()`](https://cmowrer13.github.io/qgcomp.multi/reference/b_sensitivity.md),
  and
  [`q_sensitivity()`](https://cmowrer13.github.io/qgcomp.multi/reference/q_sensitivity.md)
- native multiply imputed analysis with
  [`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
- optional bootstrap-level parallelism for larger workflows

Those features are documented in more detail in the workflow article,
function reference, and the fuller README on GitHub.

## Where to go next

- Read the [workflow
  article](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-workflow.md)
  for a fuller end-to-end example
- Browse the [function
  reference](https://cmowrer13.github.io/qgcomp.multi/reference/index.md)
  for the complete API
- Review the [news
  page](https://cmowrer13.github.io/qgcomp.multi/news/index.md) for
  recent release history
- View the source repository on
  [GitHub](https://github.com/cmowrer13/qgcomp.multi)

## References

Keil AP, Buckley JP, O’Brien KM, Ferguson KK, Zhao S, White AJ. A
Quantile-Based g-Computation Approach to Addressing the Effects of
Exposure Mixtures. Environmental Health Perspectives.
2020;128(4):047004.
