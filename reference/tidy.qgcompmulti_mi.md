# Tidy coefficient summaries for pooled qgcompmulti multiple-imputation fits

Returns a broom-style coefficient table for the pooled marginal
structural model (MSM) coefficients stored in a fitted `qgcompmulti_mi`
object.

## Usage

``` r
# S3 method for class 'qgcompmulti_mi'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti_mi` object.

- conf.int:

  Logical; if `TRUE`, add confidence interval columns.

- conf.level:

  Confidence level for interval columns when `conf.int = TRUE`.

- ...:

  Not used.

## Value

A data frame with one row per pooled MSM coefficient and columns `term`,
`estimate`, `std.error`, `statistic`, `df`, and `p.value`. When
`conf.int = TRUE`, the returned data frame also includes `conf.low` and
`conf.high`.

## Details

This method stays coefficient-centric. It reports the pooled MSM
coefficients and Rubin-style inferential quantities that belong
naturally in a reporting table, while leaving fit-level metadata to
[`broom::glance()`](https://broom.tidymodels.org/reference/reexports.html)
and leaving per-imputation fitted-model inspection to the stored
`"qgcompmulti"` objects when `keep_fits = TRUE`.

## See also

[`broom::glance()`](https://broom.tidymodels.org/reference/reexports.html),
[`coef.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`vcov.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`confint.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`summary.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti_mi.md),
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
