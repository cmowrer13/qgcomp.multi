# Confidence regions for qgcompmulti coefficients

Builds a bootstrap-covariance chi-squared confidence region for selected
marginal structural model (MSM) coefficients from a fitted
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
object.

## Usage

``` r
confregion(object, ...)

# S3 method for class 'qgcompmulti'
confregion(
  object,
  parm = c("psi1", "psi2"),
  level = 0.95,
  method = "bootstrap_chisq",
  npoints = 200L,
  ...
)

# S3 method for class 'qgcompmulti_mi'
confregion(object, ...)
```

## Arguments

- object:

  A fitted object.

- ...:

  Additional arguments passed to methods.

- parm:

  Character or integer vector identifying the MSM coefficients to
  include in the region. The default is `c("psi1", "psi2")`. Use `NULL`
  to include all MSM coefficients.

- level:

  Confidence level for the chi-squared cutoff. Must be a single number
  strictly between 0 and 1.

- method:

  Region construction method. Version `0.5.0` supports only
  `"bootstrap_chisq"`, the chi-squared ellipsoid based on the bootstrap
  covariance matrix of the selected MSM coefficients.

- npoints:

  Number of boundary points stored for two-parameter plotting. Ignored
  for regions with more than two selected coefficients.

## Value

A structured `qgcompmulti_region` object containing the selected
coefficient names, fitting-scale center, restricted bootstrap covariance
matrix, confidence level, chi-squared cutoff, degrees of freedom, scale
metadata, and two-dimensional plot data when applicable.

## Details

Confidence regions are computed on the MSM fitting coefficient scale.
For mean-difference and risk-difference estimands this is the identity
scale. For odds-ratio and rate-ratio estimands this is the log
coefficient scale. Version `0.5.0` does not geometrically transform
confidence ellipsoids onto nonlinear ratio display scales.

The region is defined as

\$\$(\theta - \hat\theta)' \hat\Sigma^{-1} (\theta - \hat\theta) \leq
\chi^2\_{p, level},\$\$

where `p` is the number of selected coefficients and the covariance
matrix is the bootstrap covariance matrix restricted to those
coefficients.

Pooled multiple-imputation confidence regions are intentionally out of
scope for Version `0.5.0`.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim_mixture_data(
  n = 400, pA = 3, pB = 3, rho_within_A = 0.3, rho_within_B = 0.3,
  rho_between = 0.2, psi1 = 0.5, psi2 = 0.3, psi12 = 0.2, seed = 123
)
fit <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat, mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"), B = 100, seed = 13
)
region <- confregion(fit, parm = c("psi1", "psi2"))
region
plot(region)
} # }
```
