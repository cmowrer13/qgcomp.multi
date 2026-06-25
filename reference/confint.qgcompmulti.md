# Wald confidence intervals for qgcompmulti coefficients

Returns Wald confidence intervals for the marginal structural model
(MSM) coefficients reported by a fitted `qgcompmulti` object.

## Usage

``` r
# S3 method for class 'qgcompmulti'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A fitted `qgcompmulti` object.

- parm:

  Optional specification of which coefficients to include. May be `NULL`
  for all coefficients, an integer vector of coefficient positions, or a
  character vector of coefficient names.

- level:

  Confidence level for the returned intervals. Must be strictly between
  0 and 1.

- ...:

  Not used

## Value

A numeric matrix with one row per selected coefficient and two columns
giving the lower and upper Wald confidence limits.

## Details

Confidence intervals are Wald intervals based on the MSM coefficient
estimates and the stored covariance matrix returned by
[`vcov.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md).
The rows of the returned matrix align with
[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md).

## See also

[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`vcov.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
