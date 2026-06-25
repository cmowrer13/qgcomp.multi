# Core extractor methods for qgcompmulti objects

Extract core fitted-model components from a `qgcompmulti` object using
standard R generic functions. These methods give the primary marginal
structural model (MSM) results and the metadata needed for model
inspection.

## Usage

``` r
# S3 method for class 'qgcompmulti'
coef(object, ...)

# S3 method for class 'qgcompmulti'
formula(x, ...)

# S3 method for class 'qgcompmulti'
nobs(object, ...)

# S3 method for class 'qgcompmulti'
df.residual(object, ...)

# S3 method for class 'qgcompmulti'
residuals(object, ...)

# S3 method for class 'qgcompmulti'
vcov(object, ...)
```

## Arguments

- object:

  A fitted `qgcompmulti` object for
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html), and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html).

- ...:

  Not used.

- x:

  A fitted `qgcompmulti` object for
  [`formula()`](https://rdrr.io/r/stats/formula.html).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector of MSM coefficients.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the full
covariance matrix aligned with `coef(object)`.

[`formula()`](https://rdrr.io/r/stats/formula.html) returns the original
formula supplied to
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md).

[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
observations actually used in the model-fitting workflow.

[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) returns a
coefficient-level residual degrees-of-freedom approximation based on the
stored MSM parameter count and the number of observations used. This is
primarily intended to support
[`mice::pool()`](https://amices.org/mice/reference/pool.html)
complete-data degrees-of-freedom handling for single-fit `qgcompmulti`
analyses.

[`residuals()`](https://rdrr.io/r/stats/residuals.html) returns the
residual vector from the stored outcome-model fit. This preserves access
to ordinary GLM residual behavior for downstream tools that expect a
conventional fitted-model residual interface.

## See also

[`summary.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md),
[`print.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti.md),
[`confint.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/confint.qgcompmulti.md),
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
