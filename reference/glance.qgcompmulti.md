# Glance summaries for qgcompmulti fits

Returns a one-row metadata summary for a fitted `qgcompmulti` object.
The output stays compact and fit-oriented rather than repeating
coefficient-level results or diagnostic summaries.

## Usage

``` r
# S3 method for class 'qgcompmulti'
glance(x, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti` object.

- ...:

  Not used.

## Value

A one-row data frame containing compact fit metadata for the
`qgcompmulti` object, including sample-size information, outcome-family
metadata, estimand-scale metadata, mixture coding metadata, bootstrap
replication counts, Monte Carlo size, interaction status, and clustering
metadata.

## Details

For original-scale fits (`q = NULL`), `q` is returned as `NA_integer_`
and the `centering` column records how the MSM intervention variables
were coded. For quantized fits, `centering` still reports the stored
setting from the fit object, but it is not inferentially active.

## See also

[`broom::tidy()`](https://broom.tidymodels.org/reference/reexports.html),
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md),
[`print.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
