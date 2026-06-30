# Fit a pooled multiple-imputation quantile g-computation model

Fits one ordinary
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
model per completed imputation and pools the marginal structural model
(MSM) coefficient inference internally using Rubin's rules. This is the
primary native multiple-imputation workflow for `qgcomp.multi`.

## Usage

``` r
qgcomp.glm.multi.mi(
  f,
  data,
  mix1,
  mix2,
  interaction = TRUE,
  family = gaussian(),
  estimand_scale = NULL,
  q = 4,
  centering = "none",
  B = 200,
  id = NULL,
  MCsize = NULL,
  seed = NULL,
  keep_fits = FALSE,
  progress = FALSE,
  parallel = FALSE,
  workers = NULL
)
```

## Arguments

- f:

  A model formula for the outcome regression.

- data:

  Either a [`mice::mids`](https://amices.org/mice/reference/mids.html)
  object or a non-empty list of completed data frames.

- mix1:

  A character vector giving the names of the variables in the first
  exposure mixture.

- mix2:

  A character vector giving the names of the variables in the second
  exposure mixture.

- interaction:

  Logical; if `TRUE`, include the mixture interaction in each
  per-imputation fit and in the pooled MSM.

- family:

  A GLM family object for the per-imputation outcome model.

- estimand_scale:

  Optional character string naming the fitted MSM estimand scale.
  Supported values depend on `family` and its link and are forwarded
  into each underlying
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
  fit. If `NULL`, the ordinary-fit family default is used within each
  completed-data fit.

- q:

  Integer greater than or equal to `2` giving the number of quantiles
  used to discretize the exposure variables, or `NULL` to fit on the
  original exposure scale.

- centering:

  Character string controlling how the MSM intervention variables are
  coded when `q = NULL`. Must be one of `"none"` or `"median"`.

- B:

  Integer greater than or equal to `2` giving the number of bootstrap
  replications used within each completed-data fit.

- id:

  Optional character string giving the name of a cluster identifier
  variable for clustered bootstrap resampling.

- MCsize:

  Optional integer controlling the Monte Carlo sample size used within
  each completed-data fit. If `NULL`, the wrapper uses the completed
  dataset row count.

- seed:

  Optional integer master seed. When supplied, the wrapper
  deterministically derives one distinct fit-specific seed per
  imputation and stores both the master seed and the fit-specific seeds
  in the pooled result.

- keep_fits:

  Logical; if `TRUE`, retain the full per-imputation `qgcompmulti`
  fitted objects inside the pooled result. Defaults to `FALSE` to keep
  the pooled object lean.

- progress:

  Logical; if `TRUE`, request the existing bootstrap progress display
  within each per-imputation fit. Progress remains serial-only, so
  `progress = TRUE` with `parallel = TRUE` is normalized once at the
  wrapper level and disabled with a single warning.

- parallel:

  Logical; if `TRUE`, allow each per-imputation fit to use the
  established bootstrap-level parallel engine from
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md).
  The multiple-imputation loop itself remains serial.

- workers:

  Optional integer worker count passed through to the per-imputation
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
  calls when `parallel = TRUE`. Leave `NULL` to use the active
  non-sequential `future` plan when one is already set, or otherwise let
  the ordinary-fit path choose a temporary local worker count
  automatically. If a non-sequential `future` plan is already active,
  supplying an explicit `workers` value is treated as an unsupported
  combination and errors cleanly

## Value

An object of class `"qgcompmulti_mi"` containing pooled MSM
coefficients, pooled covariance and standard errors, Rubin pooling
components, MI metadata, and optionally the retained per-imputation
fitted objects.

## Details

This wrapper is an orchestrator over repeated ordinary fits rather than
a second modeling engine. After normalizing the imputations into a
completed data-frame list, the function validates that the completed
datasets are structurally compatible, fits one
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
model per imputation, and pools the resulting MSM coefficient vectors
and covariance matrices using the internal Rubin-pooling helpers defined
for the `"qgcompmulti_mi"` class. This is the recommended route when you
want `qgcomp.multi` to handle the completed-data normalization, repeated
fitting, and Rubin pooling step itself. Users who need to stay inside a
broader `mice` pipeline can still fit
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
through
[`mice::with()`](https://amices.org/mice/reference/with.mids.html) and
then rely on the package's `tidy()`, `glance()`, and
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html)
compatibility work when that better matches the analysis pipeline.

Optional parallel execution stays at one level only. The wrapper itself
still iterates over imputations serially, but each ordinary
per-imputation fit may dispatch its bootstrap replications through the
same `future.apply` backend already used for single-dataset fits.
Reproducibility for these workflows is defined within a fixed backend
and execution mode: the wrapper treats `seed` as a master seed,
deterministically expands it into independent per-imputation fit seeds,
and then each ordinary fit expands its own fit seed into full-fit and
worker-level bootstrap seeds. Clustered workflows are also supported:
when `id` is supplied, each completed dataset must share identical
cluster IDs in the same order, and the ordinary per-imputation fits
continue to resample clusters rather than rows.

Pooled results stay focused on MSM inference and a compact summary
layer. The pooled object does not yet try to replicate the full
prediction and diagnostics surface of ordinary single-fit `qgcompmulti`
objects. Pooled confidence intervals use Wald intervals on the fitting
scale; pooled bootstrap interval methods and pooled confidence regions
are out of scope for this release. Use the pooled extractors,
[`summary()`](https://rdrr.io/r/base/summary.html), `tidy()`, and
`glance()` for reporting, and use `keep_fits = TRUE` only when you
explicitly need the stored per-imputation single-fit objects.

When `data` is a `mids` object, the optional `mice` package must be
installed. Inputs supplied as a plain list must already be fully
completed: required analysis variables cannot still contain missing
values.

## Examples

``` r
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
  seq_len(3),
  function(i) {
    dat_i <- dat
    dat_i$X1 <- dat_i$X1 + rnorm(nrow(dat_i), sd = 0.02 * i)
    dat_i$W2 <- dat_i$W2 + rnorm(nrow(dat_i), sd = 0.02 * i)
    dat_i
  }
)

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

fit_mi$results$coef_table
#>                Estimate Std. Error    t value       df     Pr(>|t|)
#> (Intercept) -0.29854692  0.3195680 -0.9342202 375.1619 0.3507912013
#> psi1         0.81187841  0.2313854  3.5087708 336.6317 0.0005113559
#> psi2         0.51048473  0.2967411  1.7203033 469.4994 0.0860362498
#> psi1:psi2    0.07315005  0.1770865  0.4130752 281.5038 0.6798659465
coef(fit_mi)
#> (Intercept)        psi1        psi2   psi1:psi2 
#> -0.29854692  0.81187841  0.51048473  0.07315005 
confint(fit_mi)
#>                   2.5 %    97.5 %
#> (Intercept) -0.92691591 0.3298221
#> psi1         0.35673495 1.2670219
#> psi2        -0.07262036 1.0935898
#> psi1:psi2   -0.27543176 0.4217319

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
#> [1] TRUE
fit_mi_keep$mixtures$centering
#> [1] "median"

if (FALSE) { # \dontrun{
# Bootstrap-level parallelism inside each per-imputation fit
fit_mi_parallel <- qgcomp.glm.multi.mi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = imp_list,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  interaction = TRUE,
  q = 4,
  B = 5,
  seed = 19,
  parallel = TRUE,
  workers = 2
)

# `mids` inputs remain supported when the optional `mice` package is
# installed
if (requireNamespace("mice", quietly = TRUE)) {
  dat_miss <- dat
  dat_miss$X1[sample.int(nrow(dat_miss), 8)] <- NA_real_
  mids_obj <- mice::mice(dat_miss, m = 3, maxit = 1, printFlag = FALSE)

  fit_mi_mids <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = mids_obj,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 5,
    seed = 23
  )
}
} # }
```
