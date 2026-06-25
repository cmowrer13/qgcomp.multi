# Print a qgcompmulti fit

Prints a moderately detailed overview of a fitted `qgcompmulti` object
for ordinary interactive use. The display highlights the marginal
structural model (MSM) results while keeping the full outcome-model
context for
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md).

## Usage

``` r
# S3 method for class 'qgcompmulti'
print(x, ...)
```

## Arguments

- x:

  A fitted `qgcompmulti` object.

- ...:

  Not used.

## Value

The input object, invisibly.

## Details

`print.qgcompmulti()` is designed for ordinary interactive use. It gives
a moderately detailed overview of the fitted model, including the call,
analysis mode, mixture definitions, and labeled MSM coefficient table.
For a fuller report, use
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md).

## See also

[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md),
[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`confint.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/confint.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
