# Confidence intervals for qgcompmulti coefficients

Returns confidence intervals for the marginal structural model (MSM)
coefficients reported by a fitted `qgcompmulti` object. Wald,
percentile, and basic bootstrap intervals are supported for single-fit
objects. For ratio estimands, intervals are computed on the MSM fitting
scale and then transformed to the user-facing ratio scale.

## Usage

``` r
# S3 method for class 'qgcompmulti'
confint(object, parm = NULL, level = 0.95, method = NULL, ...)
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

- method:

  Optional interval method. `NULL` uses the fitted object's stored
  default interval method. Supported values are `"wald"`,
  `"percentile"`, and `"basic"`. The basic method is the
  reverse-percentile bootstrap interval.

- ...:

  Not used.

## Value

A numeric matrix with one row per selected coefficient and two columns
giving the lower and upper confidence limits on the display scale
recorded in `object$analysis$estimand_scale`.

## Details

Wald confidence intervals are based on the MSM coefficient estimates and
the stored bootstrap covariance matrix returned by
[`vcov.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md).
Percentile intervals use empirical quantiles of the stored bootstrap
coefficient draws. Basic intervals use the reverse-percentile
construction, reflecting the bootstrap quantiles around the full-sample
estimate.

All interval calculations are carried out on the MSM fitting scale. When
the active estimand is `"odds_ratio"` or `"rate_ratio"`, the returned
limits are exponentiated for display. The rows of the returned matrix
align with
[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
but the numeric scale may differ for ratio estimands. Pooled
multiple-imputation objects intentionally remain Wald-only.

## See also

[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`vcov.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
