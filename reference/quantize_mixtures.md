# Quantize exposure variables for two-mixture quantile g-computation

Converts continuous exposure variables belonging to two mixtures into
discrete quantile categories. Each exposure is partitioned into `q`
approximately equal-sized groups using empirical quantiles, and replaced
with integer values ranging from 0 to `q - 1`. This transformation
places mixture components on a common scale, enabling interpretation of
model parameters as effects of simultaneous one-quantile increases in
each component.

## Usage

``` r
quantize_mixtures(data, mix1, mix2, q = 4)
```

## Arguments

- data:

  A data frame containing the exposure variables and any additional
  variables (e.g., outcome or covariates).

- mix1:

  A character vector giving the names of the variables in the first
  exposure mixture.

- mix2:

  A character vector giving the names of the variables in the second
  exposure mixture.

- q:

  Integer greater than or equal to 2 giving the number of quantiles used
  to discretize each exposure variable, or `NULL` to skip quantization
  and return the mixture variables unchanged after validation.

## Value

A data frame of the same dimensions as `data`. When `q` is numeric, the
variables listed in `mix1` and `mix2` have been replaced by their
quantized versions. Quantized variables take integer values in
`{0, 1, ..., q - 1}`. When `q = NULL`, the returned data are unchanged.

## Details

Quantization is a key step in quantile g-computation, as it ensures that
exposure variables measured on different scales are made comparable. By
transforming each exposure into quantile categories, a one-unit change
in a quantized exposure corresponds to a shift of one quantile category.
In the fitted two-mixture MSM, this supports the interpretation of the
mixture coefficients as summaries of joint one-quantile shifts in all
components of a mixture.

Quantiles assign observations to groups of approximately equal size. The
resulting categories are then shifted to start at 0 (rather than 1) to
align with the intervention levels used in the marginal structural
model.

When `q = NULL`, this function performs no discretization and simply
returns `data` unchanged after checking that the mixture definitions are
valid.

## Examples

``` r
dat <- sim_mixture_data(
  n = 100,
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

dat_q <- quantize_mixtures(
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  q = 4
)

head(dat_q)
#>            Y X1 X2 X3 W1 W2 W3           C
#> 1  2.3346820  3  1  0  3  2  2  1.07401226
#> 2  0.3655440  3  1  2  2  2  0 -0.02734697
#> 3 -0.3153138  0  1  0  1  1  0 -0.03333034
#> 4  0.2138929  1  2  1  3  2  0 -1.51606762
#> 5  0.9256759  1  0  2  2  1  3  0.79038534
#> 6  0.7903596  0  0  1  0  2  0 -0.21073418

dat_cont <- quantize_mixtures(
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  q = NULL
  )
```
