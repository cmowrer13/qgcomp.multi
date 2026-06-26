# Clustered Data and Repeated Measures

## Introduction

This article is a short companion to the main workflow vignette. It
focuses on one practical feature that matters whenever rows are not the
right resampling unit for inference:

- clustered sampling, where observations are grouped within higher-level
  units
- repeated measures, where multiple rows belong to the same individual

In these settings, `qgcomp.multi` can still be used through the usual
fitting interface. The key change is to supply the cluster or subject
identifier through the `id` argument.

When `id` is supplied, the package performs bootstrap resampling at the
cluster level rather than at the individual row level. This changes the
uncertainty calculation to account for within-cluster dependence.

## A small clustered example

The example below creates a small clustered dataset with 50 clusters and
3 observations per cluster. The same setup also works as a simple
repeated- measures template if each cluster ID corresponds to one
participant observed multiple times.

``` r

library(qgcomp.multi)

n_clusters <- 50
cluster_size <- 3
n_total <- n_clusters * cluster_size

dat <- sim_mixture_data(
  n = n_total,
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

dat$cluster_id <- rep(seq_len(n_clusters), each = cluster_size)

head(dat[c("cluster_id", "Y", "X1", "W1", "C")])
#>   cluster_id         Y          X1         W1          C
#> 1          1 1.2062091 -0.27889714 -0.6249504 -1.0141142
#> 2          1 0.7188395  0.04136179 -0.9336474 -0.7913139
#> 3          1 1.1034198 -1.53834328 -1.3698203  0.2995937
#> 4          2 2.0207057 -1.15308989 -0.3414860  1.6390519
#> 5          2 3.8536386 -0.32546695 -0.5866225  1.0846170
#> 6          2 0.2601378 -0.91541887  0.7952908 -0.6245675
```

## Fit the model with cluster-level resampling

Supply the cluster identifier through `id`.

``` r

fit_clustered <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 5,
  id = "cluster_id",
  MCsize = nrow(dat),
  seed = 13
)
```

This is still the ordinary quantized two-mixture workflow. The
difference is that bootstrap uncertainty is now based on resampling the
50 clusters instead of resampling the 150 rows one by one.

## Confirm that the fit recognized the clustering

The printed summary reports clustering metadata directly:

``` r

summary(fit_clustered)
#> Summary of qgcompmulti fit
#> 
#> Call:
#> qgcomp.glm.multi(f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C, data = dat, 
#>     mix1 = c("X1", "X2", "X3"), mix2 = c("W1", "W2", "W3"), interaction = TRUE, 
#>     q = 4, B = 5, id = "cluster_id", MCsize = nrow(dat), seed = 13)
#> 
#> Model overview:
#>   Formula: Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C
#>   Outcome: Y
#>   Family: gaussian (identity)
#>   Observations used: 150
#>   Exposure mode: Quantized exposures (q = 4)
#>   MSM interaction: included
#>   Bootstrap replications: 5
#>   Monte Carlo size: 150
#>   Random seed: 13
#>   Clusters: 50 (id = cluster_id)
#> 
#> Mixtures:
#>   Mixture 1: X1, X2, X3
#>   Mixture 2: W1, W2, W3
#> 
#> MSM coefficients:
#>                       Estimate Std. Error z value  Pr(>|z|)    
#> Intercept             -0.05731    0.31305 -0.1831 0.8547451    
#> Mixture 1 main effect  0.68283    0.18483  3.6943 0.0002205 ***
#> Mixture 2 main effect  0.29431    0.31141  0.9451 0.3446253    
#> Mixture interaction    0.12858    0.12852  1.0005 0.3170693    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Outcome model context:
#>   Model class: glm
#>   Estimated parameters: 9
#>   AIC: 434.522
#>   Null deviance: 431.637
#>   Residual deviance: 139.245
```

You can also inspect the stored fit metadata programmatically:

``` r

fit_clustered$data_info[c("has_clusters", "cluster_var", "n_clusters", "n_used")]
#> $has_clusters
#> [1] TRUE
#> 
#> $cluster_var
#> [1] "cluster_id"
#> 
#> $n_clusters
#> [1] 50
#> 
#> $n_used
#> [1] 150
```

For repeated-measures data, this is the same idea: set `id` to the
participant identifier so that all rows from the same participant move
together during the bootstrap resampling step.

## What the id argument changes

Supplying `id` changes the resampling unit used for bootstrap inference.

That is the right adjustment when:

- observations are clustered within a higher-level unit
- repeated measurements are recorded on the same individual
- rows within a cluster should not be treated as independently resampled
  units

In practical terms, this lets you keep the same
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
fitting interface while asking the bootstrap to respect the structure of
the data.

## What the id argument does not change

The `id` argument is about inference, not about changing the outcome
model into a mixed-effects or subject-specific model.

It does not:

- automatically add random effects
- change the formula interface
- solve every longitudinal modeling question on its own

Instead, it tells `qgcomp.multi` what the correct resampling unit is
when the rows are clustered or repeated within a subject.

## Brief notes on other workflows

The same clustered-resampling idea also works in other package
workflows.

- If you fit an original-scale model with `q = NULL`, you can still
  supply `id` in the same way.
- If you use
  [`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
  with clustered completed datasets, the wrapper also supports clustered
  bootstrap resampling.

For multiply imputed clustered workflows, the completed datasets must
share the same cluster IDs in the same order.

## Practical takeaways

The main question is whether the bootstrap should resample rows or
resample clusters. In `qgcomp.multi`, supplying `id` tells the package
to use the clustered resampling path while keeping the rest of the
two-mixture workflow the same.

## Related Articles

- [Applied Workflow for Two-Mixture Quantile
  g-Computation](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-workflow.md)
- [Multiple Imputation and Parallel
  Workflows](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-mi-parallel.md)
- [Diagnostics and Sensitivity
  Checks](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-diagnostics-sensitivity.md)
