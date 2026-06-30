# Core extractor methods for pooled qgcompmulti multiple-imputation objects

Extract the pooled marginal structural model (MSM) coefficient results
from a fitted `qgcompmulti_mi` object using standard R generics.

## Usage

``` r
# S3 method for class 'qgcompmulti_mi'
coef(object, ...)

# S3 method for class 'qgcompmulti_mi'
vcov(object, ...)

# S3 method for class 'qgcompmulti_mi'
confint(object, parm = NULL, level = 0.95, method = NULL, ...)
```

## Arguments

- object:

  A fitted `qgcompmulti_mi` object.

- ...:

  Not used.

- parm:

  Optional specification of which coefficients to include. May be `NULL`
  for all coefficients, an integer vector of coefficient positions, or a
  character vector of coefficient names.

- level:

  Confidence level for the returned intervals. Must be strictly between
  0 and 1.

- method:

  Optional interval method. Pooled multiple-imputation coefficient
  reporting supports only `"wald"` in Version 0.5.0. `NULL` uses the
  fitted object's stored default, which is `"wald"` for pooled MI fits.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector of pooled MSM coefficients on the fitting scale.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the pooled
covariance matrix aligned with `coef(object)`.

[`confint()`](https://rdrr.io/r/stats/confint.html) returns pooled
Wald-style confidence intervals. The Wald calculation uses Rubin-pooled
coefficients, standard errors, and term-specific degrees of freedom on
the fitting scale. For odds-ratio and rate-ratio estimands, returned
interval limits are exponentiated for display.

## Details

These extractors operate on the pooled multiple-imputation result, not
on the individual completed-data fits. They therefore return
Rubin-pooled MSM coefficient summaries for the inferential target
defined by
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md),
rather than per-imputation coefficient tables, pooled prediction
objects, or pooled diagnostics. If you need to inspect the stored
completed-data fits directly, fit with `keep_fits = TRUE` and extract
the individual `"qgcompmulti"` objects from
`object$fits$imputation_fits`.

## See also

[`summary.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti_mi.md),
[`print.qgcompmulti_mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti_mi.md),
[`broom::tidy()`](https://broom.tidymodels.org/reference/reexports.html),
[`broom::glance()`](https://broom.tidymodels.org/reference/reexports.html),
[`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
