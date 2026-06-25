# Sensitivity to quantization choice in qgcompmulti fits

Re-fits
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
across multiple integer `q` values while preserving the rest of the
analysis specification.

## Usage

``` r
q_sensitivity(
  f,
  data,
  mix1,
  mix2,
  q_values,
  interaction = TRUE,
  family = gaussian(),
  centering = "none",
  B = 200,
  id = NULL,
  MCsize = nrow(data),
  seed = NULL,
  keep_fits = TRUE
)
```

## Arguments

- f, data, mix1, mix2, interaction, family, centering, B, id, MCsize,
  seed:

  Arguments passed through to
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md).

- q_values:

  Integer vector of quantization choices to compare.

- keep_fits:

  Logical; if `TRUE`, retain the full fitted objects.

## Value

An object of class `"qgcompmulti_q_sensitivity"`.

## Details

This helper is intended for robustness assessment, not for coefficient
ranking across different quantization choices.

`q` sensitivity is implemented as a repeated-fit workflow in which all
settings other than `q` are held fixed.

Users should be cautious when comparing raw coefficient magnitudes
across different values of `q`. A larger `q` implies a smaller
one-quantile intervention step, so smaller coefficients may be expected
mechanically even when the broader qualitative pattern of the fitted
surface is stable. For that reason, the printed sensitivity object
includes an explicit comparability note.

The helper therefore supports sensitivity assessment, not a claim that
one choice of `q` produces “stronger” or “weaker” effects based only on
raw one-step coefficient magnitudes.

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

q_sensitivity(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  q_values = c(3, 4, 5),
  B = 100,
  seed = 13
)
} # }
```
