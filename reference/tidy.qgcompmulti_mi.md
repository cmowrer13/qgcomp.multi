# Tidy coefficient summaries for pooled qgcompmulti multiple-imputation fits

Returns a broom-style coefficient table for the pooled marginal
structural model (MSM) coefficients stored in a fitted `qgcompmulti_mi`
object. The core `estimate` and `std.error` columns stay on the
Rubin-pooled fitting scale, while display columns expose the active
estimand scale for reporting.

## Usage

``` r
# S3 method for class 'qgcompmulti_mi'
tidy(x, conf.int = FALSE, conf.level = 0.95, method = NULL, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti_mi` object.

- conf.int:

  Logical; if `TRUE`, add confidence interval columns.

- conf.level:

  Confidence level for interval columns when `conf.int = TRUE`.

- method:

  Optional interval method. Pooled multiple-imputation coefficient
  reporting supports only `"wald"` in Version 0.5.0.

- ...:

  Not used.

## Value

A data frame with one row per pooled MSM coefficient. Columns
`estimate`, `std.error`, `statistic`, `df`, and `p.value` are on the
fitting scale. Columns `display.estimate`, `display.conf.low`, and
`display.conf.high` are on the active estimand scale when present.

## Details

Pooled multiple-imputation inference is carried out on the fitting
scale. For odds-ratio and rate-ratio estimands, the display columns
transform only the final pooled estimates and interval limits. Rubin
pooling, standard errors, test statistics, and degrees of freedom remain
on the fitting scale.

## See also

[`broom::glance()`](https://broom.tidymodels.org/reference/reexports.html),
[`coef.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`vcov.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`confint.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`summary.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti_mi.md),
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
