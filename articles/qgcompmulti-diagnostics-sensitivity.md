# Diagnostics and Sensitivity Checks

## Introduction

This article is a short companion to the main workflow vignette. It
focuses on the post-fit tools that help you judge whether a
`qgcomp.multi` analysis is stable and interpretable:

- [`support()`](https://cmowrer13.github.io/qgcomp.multi/reference/support.md)
  for the intervention grid behind the stored surface
- [`diagnostics()`](https://cmowrer13.github.io/qgcomp.multi/reference/diagnostics.md)
  for bootstrap behavior
- [`adequacy()`](https://cmowrer13.github.io/qgcomp.multi/reference/adequacy.md)
  for how well the fitted MSM summarizes the exact fit-time surface
- sensitivity helpers for Monte Carlo size, bootstrap size, and
  quantization choice

These are not replacements for scientific judgment. They are structured
checks that help you see where the analysis is easy to defend and where
it needs more thought.

## Fit a small working example

``` r

library(qgcomp.multi)

dat <- sim_mixture_data(
  n = 300,
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

fit <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 10,
  MCsize = nrow(dat),
  seed = 13
)
```

## Intervention support

``` r

support(fit)
#> qgcompmulti intervention support diagnostic
#> 
#> Mode: quantized
#> Centering: none
#> Grid points: 16
#> Intervention psi1 range: [0, 3.00]
#> Intervention psi2 range: [0, 3.00]
```

The support diagnostic summarizes the intervention grid used to define
the stored surface. It is especially easy to read for quantized models,
where the grid is finite and discrete.

What it does **not** do is prove causal positivity. It tells you what
intervention values are represented in the stored fit-time surface, not
whether every underlying intervention is equally plausible in the source
population.

## Bootstrap behavior

``` r

diagnostics(fit, type = "bootstrap")
#> qgcompmulti bootstrap diagnostic
#> 
#> Requested replications: 10
#> Successful replications: 10
#> Failed replications: 0
#> Success rate: 100.000%
```

This reports how many bootstrap replications were requested, how many
were retained, and whether any failures were logged.

If failures appear in a real analysis, that is a sign to look more
closely at model specification, numerical stability, and covariate
sparsity.

## MSM adequacy

``` r

adequacy(fit)
#> qgcompmulti MSM adequacy diagnostic
#> 
#> Grid points: 16
#> Mean absolute error: 0.000
#> RMSE: 0.000
#> Maximum absolute error: 0.000
#> Mean signed error: -0.000
#> Correlation: 1.000
#> 
#> Adequacy compares the exact fit-time counterfactual surface to the fitted MSM surface on the response scale.
```

The adequacy diagnostic asks a narrower question than model truth:

> Does the fitted MSM do a good job approximating the
> intervention-response surface implied by the fitted outcome model?

Good adequacy means the MSM tracks the stored exact fit-time surface
well on the fit-time grid. It does not prove that the outcome model is
correct, and it does not prove that the exact surface is globally linear
between or beyond those grid points.

## Sensitivity helpers

The package includes dedicated helpers for practical repeated-fit
checks.

For Monte Carlo size:

``` r

mcsize_sensitivity(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  MCsize_values = c(150, 300),
  interaction = TRUE,
  q = 4,
  B = 10,
  seed = 13
)
```

For bootstrap size:

``` r

b_sensitivity(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  B_values = c(10, 20, 40),
  interaction = TRUE,
  q = 4,
  MCsize = nrow(dat),
  seed = 13
)
```

For quantization choice:

``` r

q_sensitivity(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  q_values = c(3, 4, 5),
  interaction = TRUE,
  B = 10,
  seed = 13
)
```

## How to interpret these checks

The post-fit checks answer different questions:

- [`support()`](https://cmowrer13.github.io/qgcomp.multi/reference/support.md)
  helps you inspect the intervention grid behind the stored surface
- `diagnostics(..., type = "bootstrap")` tells you whether the
  repeated-fit uncertainty calculation behaved cleanly
- [`adequacy()`](https://cmowrer13.github.io/qgcomp.multi/reference/adequacy.md)
  tells you whether the MSM is a reasonable summary of the exact
  fit-time surface on the stored grid
- the sensitivity helpers tell you whether conclusions are stable to
  practical analysis choices such as `MCsize`, `B`, and `q`

These checks are most useful together. A clean coefficient table is
easier to defend when the intervention grid is interpretable, bootstrap
behavior is stable, adequacy is acceptable, and the main conclusions do
not move much under reasonable sensitivity changes.

## Related Articles

- [Applied Workflow for Two-Mixture Quantile
  g-Computation](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-workflow.md)
- [Multiple Imputation and Parallel
  Workflows](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-mi-parallel.md)
- [Clustered Data and Repeated
  Measures](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-clustered-data.md)
