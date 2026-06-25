# Summarize a pooled qgcompmulti multiple-imputation fit

Builds a structured summary object for a fitted `qgcompmulti_mi` model.
The summary follows the same general presentation as
[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md)
while adding the pooled multiple-imputation metadata and Rubin
diagnostics that matter for downstream reporting.

Prints a structured summary of a pooled `qgcompmulti_mi` fit,
emphasizing the pooled MSM results while reporting the shared analysis
settings and compact Rubin diagnostics as supporting context.

## Usage

``` r
# S3 method for class 'qgcompmulti_mi'
summary(object, ...)

# S3 method for class 'summary.qgcompmulti_mi'
print(x, ...)
```

## Arguments

- object:

  A fitted `qgcompmulti_mi` object.

- ...:

  Not used.

- x:

  A `"summary.qgcompmulti_mi"` object produced by
  `summary.qgcompmulti_mi()`.

## Value

An object of class `"summary.qgcompmulti_mi"` with components for the
matched call, original formula, model overview, MI overview, mixture
definitions, labeled pooled MSM coefficient table, pooled Rubin
diagnostics, and stored labels.

The input object, invisibly.

## See also

[`print.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti_mi.md),
[`coef.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`confint.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md),
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
