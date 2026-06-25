# Summarize a qgcompmulti fit

Builds a structured summary object for a fitted `qgcompmulti` model. The
summary foregrounds the MSM results because those are the main reported
estimands, but those coefficients should still be interpreted in light
of the intervention coding and the adequacy of the MSM as a summary of
the fitted outcome-model surface.

Prints a structured summary of a fitted `qgcompmulti` model, emphasizing
the marginal structural model results while reporting the outcome model
as supporting context.

## Usage

``` r
# S3 method for class 'qgcompmulti'
summary(object, ...)

# S3 method for class 'summary.qgcompmulti'
print(x, ...)
```

## Arguments

- object:

  A fitted `qgcompmulti` object.

- ...:

  Not used.

- x:

  A `"summary.qgcompmulti"` object produced by `summary.qgcompmulti()`.

## Value

An object of class `"summary.qgcompmulti"` with components for the
matched call, original formula, model overview, mixture definitions,
labeled MSM coefficient table, outcome-model context, and stored labels.

The input object, invisibly.

## Details

The printed summary highlights the marginal structural model (MSM)
coefficient results while reporting the outcome-model family, formula,
exposure-mode setting, sample-size metadata, and compact outcome-model
fit statistics as supporting context.

## See also

[`print.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti.md),
[`coef.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md),
[`confint.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/confint.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
