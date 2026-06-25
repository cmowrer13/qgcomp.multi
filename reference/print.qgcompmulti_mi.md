# Print a pooled qgcompmulti multiple-imputation fit

Prints a moderately detailed overview of a pooled `qgcompmulti_mi`
object for ordinary interactive use. The display mirrors the single-fit
`qgcompmulti` print method while adding the key multiple-imputation
metadata needed to interpret the pooled inference.

## Usage

``` r
# S3 method for class 'qgcompmulti_mi'
print(x, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti_mi` object.

- ...:

  Not used.

## Value

The input object, invisibly.

## See also

[`summary.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti_mi.md),
[`coef.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`confint.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
