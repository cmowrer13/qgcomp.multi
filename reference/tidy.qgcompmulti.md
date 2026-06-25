# Tidy coefficient summaries for qgcompmulti fits

Returns a broom-style coefficient table for the marginal structural
model (MSM) coefficients stored in a fitted `qgcompmulti` object. The
output is coefficient-centric and uses the internal machine-readable MSM
term names.

## Usage

``` r
# S3 method for class 'qgcompmulti'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti` object.

- conf.int:

  Logical; if `TRUE`, add Wald confidence interval columns.

- conf.level:

  Confidence level for interval columns when `conf.int = TRUE`.

- ...:

  Not used.

## Value

A data frame with one row per MSM coefficient and columns `term`,
`estimate`, `std.error`, `statistic`, and `p.value`. When
`conf.int = TRUE`, the returned data frame also includes `conf.low` and
`conf.high`.

## Details

This method is intentionally coefficient-centric. It does not add
fit-level metadata columns or presentation-oriented labels because those
belong in
[`broom::glance()`](https://broom.tidymodels.org/reference/reexports.html)
or the existing print and summary methods.

## See also

[`broom::glance()`](https://broom.tidymodels.org/reference/reexports.html),
[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`vcov.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`confint.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/confint.qgcompmulti.md),
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
