# Glance summaries for pooled qgcompmulti multiple-imputation fits

Returns a one-row metadata summary for a fitted `qgcompmulti_mi` object.
The output is intentionally compact and inference-focused.

## Usage

``` r
# S3 method for class 'qgcompmulti_mi'
glance(x, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti_mi` object.

- ...:

  Not used.

## Value

A one-row data frame containing compact pooled-fit metadata for the
`qgcompmulti_mi` object, including sample-size information, the number
of imputations, input-type metadata, mixture coding metadata, Monte
Carlo size, interaction status, retained-fit status, and clustering
metadata.

## Details

`glance()` stays compact for pooled MI objects. It summarizes how the
pooled fit was constructed, but it does not try to stand in for pooled
prediction output, pooled diagnostics, or a full scientific narrative
about model adequacy.

## See also

[`broom::tidy()`](https://broom.tidymodels.org/reference/reexports.html),
[`summary.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti_mi.md),
[`print.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti_mi.md),
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
