# Tidy coefficient summaries for qgcompmulti fits

Returns a broom-style coefficient table for the marginal structural
model (MSM) coefficients stored in a fitted `qgcompmulti` object. The
core `estimate` and `std.error` columns stay on the MSM fitting scale so
downstream tools such as
[`mice::pool()`](https://amices.org/mice/reference/pool.html) can use
them coherently. Display columns expose the active estimand scale for
ordinary reporting.

## Usage

``` r
# S3 method for class 'qgcompmulti'
tidy(x, conf.int = FALSE, conf.level = 0.95, method = NULL, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti` object.

- conf.int:

  Logical; if `TRUE`, add confidence interval columns.

- conf.level:

  Confidence level for interval columns when `conf.int = TRUE`.

- method:

  Optional interval method for confidence intervals. `NULL` uses the
  fitted object's stored default interval method. Supported values are
  `"wald"`, `"percentile"`, and `"basic"`.

- ...:

  Not used.

## Value

A data frame with one row per MSM coefficient. Columns `estimate`,
`std.error`, `statistic`, and `p.value` are on the fitting scale.
Columns `display.estimate`, `display.conf.low`, and `display.conf.high`
are on the active estimand scale when present.

## Details

Keeping `estimate` and `std.error` on the fitting scale preserves
machine pooling behavior for multiple-imputation workflows that call
[`mice::pool()`](https://amices.org/mice/reference/pool.html). For
odds-ratio and rate-ratio estimands, use `display.estimate` and the
display confidence interval columns for user-facing ratio summaries.

## See also

[`broom::glance()`](https://broom.tidymodels.org/reference/reexports.html),
[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`vcov.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`confint.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/confint.qgcompmulti.md),
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
