# Sensitivity to Monte Carlo size in qgcompmulti fits

Re-fits
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
across multiple `MCsize` values while preserving the rest of the
analysis specification.

## Usage

``` r
mcsize_sensitivity(
  f,
  data,
  mix1,
  mix2,
  MCsize_values,
  interaction = TRUE,
  family = gaussian(),
  q = 4,
  centering = "none",
  B = 200,
  id = NULL,
  seed = NULL,
  keep_fits = TRUE
)
```

## Arguments

- f, data, mix1, mix2, interaction, family, q, centering, B, id, seed:

  Arguments passed through to
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md).

- MCsize_values:

  Integer vector of Monte Carlo sizes to compare.

- keep_fits:

  Logical; if `TRUE`, retain the full fitted objects.

## Value

An object of class `"qgcompmulti_mcsize_sensitivity"`.

## Details

This helper is intended to answer a practical computational question:
are the fitted results reasonably stable as the Monte Carlo
approximation size changes?

`MCsize` sensitivity is implemented as a repeated-fit workflow rather
than a diagnostic of one existing fit. The helper keeps the model
formula, outcome family, mixture definitions, interaction setting,
quantization choice, bootstrap count, clustering, and seed fixed while
varying only `MCsize`.

The resulting object is designed for stability assessment rather than
automatic tuning. Users should look for whether the fitted coefficients,
adequacy summaries, and bootstrap behavior are reasonably consistent
across the requested `MCsize` values.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim_mixture_data(
  n = 400,
  pA = 3,
  pB = 3,
  rho_within_A = 0.3,
  rho_within_B = 0.3,
  rho_between = 0.2,
  psi1 = 0.5,
  psi2 = 0.3,
  psi12 = 0.2,
  seed = 123
)

mcsize_sensitivity(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  MCsize_values = c(100, 200, 400),
  q = 4,
  B = 100,
  seed = 13
)
} # }
```
