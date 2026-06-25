# Multiple Imputation and Parallel Workflows

## Introduction

Version `0.4.0` adds two workflow-oriented capabilities that often
matter together in practice:

- native multiple-imputation fitting through
  [`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
- optional bootstrap-level parallel execution for ordinary and MI
  workflows

This article is a short companion to the main workflow vignette. It
focuses on when to use the native MI wrapper, what the pooled object
contains, and how the parallel option fits into larger repeated-fit
analyses.

## When to use the native MI wrapper

Use
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
when you want `qgcomp.multi` to:

- fit one
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
  model per imputation
- pool the MSM coefficients and covariance internally using Rubin’s
  rules

In Version `0.4.0`, the pooled object is focused on inference and
summary. It does not try to reproduce the full prediction and
diagnostics surface of a single-fit `qgcompmulti` object.

## A small completed-data example

``` r

library(qgcomp.multi)

dat <- sim_mixture_data(
  n = 120,
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

imp_list <- lapply(
  seq_len(2),
  function(i) {
    dat_i <- dat
    dat_i$X1 <- dat_i$X1 + rnorm(nrow(dat_i), sd = 0.02 * i)
    dat_i$W2 <- dat_i$W2 + rnorm(nrow(dat_i), sd = 0.02 * i)
    dat_i
  }
)
```

Fit the pooled model:

``` r

fit_mi <- qgcomp.glm.multi.mi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = imp_list,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 5,
  seed = 13
)
```

Inspect the pooled result:

``` r

summary(fit_mi)
#> Summary of qgcompmulti multiple-imputation fit
#> 
#> Call:
#> qgcomp.glm.multi.mi(f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C, 
#>     data = imp_list, mix1 = c("X1", "X2", "X3"), mix2 = c("W1", 
#>         "W2", "W3"), interaction = TRUE, q = 4, B = 5, seed = 13)
#> 
#> Model overview:
#>   Formula: Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C
#>   Outcome: Y
#>   Family: gaussian (identity)
#>   Observations used: 120
#>   Exposure mode: Quantized exposures (q = 4)
#>   MSM interaction: included
#>   Bootstrap replications per imputation: 5
#>   Monte Carlo size: 120
#> 
#> Multiple imputation overview:
#>   Imputed datasets: 2
#>   Input type: Completed data list
#>   Retained per-imputation fits: no
#>   Master seed: 13
#>   Stored fit-specific seeds: 2
#> 
#> Mixtures:
#>   Mixture 1: X1, X2, X3
#>   Mixture 2: W1, W2, W3
#> 
#> MSM coefficients:
#>                        Estimate Std. Error   t value         df  Pr(>|t|)    
#> Intercept             -0.341659   0.275525 -1.240028    2362459   0.21496    
#> Mixture 1 main effect  0.843965   0.181627  4.646708    3000888 3.373e-06 ***
#> Mixture 2 main effect  0.548387   0.263028  2.084901 1948154365   0.03708 *  
#> Mixture interaction    0.047446   0.150156  0.315982  298091684   0.75202    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Pooling diagnostics:
#>                               df        RIV   FMI
#> Intercept             2.3625e+06 6.5103e-04 7e-04
#> Mixture 1 main effect 3.0009e+06 5.7760e-04 6e-04
#> Mixture 2 main effect 1.9482e+09 2.2657e-05 0e+00
#> Mixture interaction   2.9809e+08 5.7923e-05 1e-04
coef(fit_mi)
#> (Intercept)        psi1        psi2   psi1:psi2 
#> -0.34165870  0.84396540  0.54838743  0.04744643
confint(fit_mi)
#>                   2.5 %    97.5 %
#> (Intercept) -0.88167781 0.1983604
#> psi1         0.48798382 1.1999470
#> psi2         0.03286187 1.0639130
#> psi1:psi2   -0.24685318 0.3417460
```

## Retaining per-imputation fits

If you want the full completed-data fits in addition to the pooled
result, set `keep_fits = TRUE`.

``` r

fit_mi_keep <- qgcomp.glm.multi.mi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = imp_list,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = NULL,
  centering = "median",
  B = 5,
  seed = 17,
  keep_fits = TRUE
)

fit_mi_keep$mi_info$keep_fits
fit_mi_keep$fits$imputation_fits[[1]]
```

This is useful when you want pooled inference as the main result but
still need to inspect one or more single-fit `qgcompmulti` objects
directly.

## `mids` objects are also supported

The wrapper also accepts
[`mice::mids`](https://amices.org/mice/reference/mids.html) objects when
the optional `mice` package is installed.

``` r

library(mice)

dat_miss <- dat
dat_miss$X1[sample.int(nrow(dat_miss), 8)] <- NA_real_
dat_miss$W2[sample.int(nrow(dat_miss), 8)] <- NA_real_

mids_obj <- mice(dat_miss, m = 3, maxit = 1, printFlag = FALSE)

fit_mi_mids <- qgcomp.glm.multi.mi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = mids_obj,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 5,
  seed = 19
)
```

## Optional bootstrap-level parallelism

Version `0.4.0` supports one level of optional parallelism.

For ordinary fits, bootstrap replications can be dispatched in parallel
inside
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md):

``` r

fit_parallel <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 20,
  seed = 29,
  parallel = TRUE,
  workers = 2
)
```

For native MI fits, the imputation loop stays serial, but the bootstrap
work inside each completed-data fit can still be parallelized:

``` r

fit_mi_parallel <- qgcomp.glm.multi.mi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = imp_list,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 10,
  seed = 31,
  parallel = TRUE,
  workers = 2
)
```

## Practical notes

Several practical limits are worth stating plainly:

- the MI loop itself is still serial in Version `0.4.0`
- reproducibility is defined within a fixed backend and execution mode,
  not as exact equality between serial and parallel runs
- `progress = TRUE` is supported only in serial mode and is disabled
  with a clear warning when `parallel = TRUE`
- pooled prediction and pooled diagnostics are not yet the main target
  of the pooled MI object

In practice, the native MI wrapper is best used when pooled MSM
inference is the main goal and you want the package to manage the
repeated fitting and Rubin pooling directly.

## Related Articles

- [Applied Workflow for Two-Mixture Quantile
  g-Computation](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-workflow.md)
- [Diagnostics and Sensitivity
  Checks](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-diagnostics-sensitivity.md)
