# Sensitivity to bootstrap iteration count in qgcompmulti fits

Re-fits
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
across multiple `B` values while preserving the rest of the analysis
specification.

## Usage

``` r
b_sensitivity(
  f,
  data,
  mix1,
  mix2,
  B_values,
  interaction = TRUE,
  family = gaussian(),
  q = 4,
  centering = "none",
  id = NULL,
  MCsize = nrow(data),
  seed = NULL,
  keep_fits = TRUE
)
```

## Arguments

- f, data, mix1, mix2, interaction, family, q, centering, id, MCsize,
  seed:

  Arguments passed through to
  [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md).

- B_values:

  Integer vector of bootstrap iteration counts to compare.

- keep_fits:

  Logical; if `TRUE`, retain the full fitted objects.

## Value

An object of class `"qgcompmulti_b_sensitivity"`.

## Details

This helper is intended to answer a practical inferential question: are
the reported bootstrap-based standard errors and related summaries
reasonably stable as the number of bootstrap replications changes?

`B` sensitivity is implemented as a repeated-fit workflow in which all
settings other than `B` are held fixed.

When `seed` is supplied, this helper treats it as a master seed and
deterministically derives one distinct fit-specific seed for each
requested value of `B`. This preserves reproducibility of the overall
sensitivity workflow while avoiding reuse of the same bootstrap
resamples across the different refits.

The resulting object is designed for stability assessment rather than
automatic tuning. Users should focus especially on whether
bootstrap-based standard errors, bootstrap retention counts, and broad
qualitative conclusions are reasonably consistent across the requested
`B` values.

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

b_sensitivity(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  B_values = c(50, 100, 200),
  q = 4,
  MCsize = nrow(dat),
  seed = 13
)
} # }
```
