# Diagnostic summaries for a qgcompmulti fit

Returns structured diagnostics for intervention support, bootstrap
behavior, and MSM adequacy from a fitted
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
object.

## Usage

``` r
diagnostics(object, type = c("all", "support", "bootstrap", "adequacy"), ...)
```

## Arguments

- object:

  A fitted `"qgcompmulti"` object.

- type:

  Character string indicating which diagnostic to return. Supported
  values are `"all"`, `"support"`, `"bootstrap"`, and `"adequacy"`.

- ...:

  Unused.

## Value

A structured diagnostic object, or a named list of all diagnostics when
`type = "all"`.

## Details

These diagnostics are deliberately kept outside the main fitted-model
summary so that users can inspect model fit, support, and computational
behavior explicitly rather than having a very dense
[`summary()`](https://rdrr.io/r/base/summary.html) method.

The three main diagnostic families are:

- **Support diagnostics**, which summarize the fit-time intervention
  grid. For `q = NULL`, these diagnostics are especially important
  because the intervention grid is built from pooled percentile values
  within each mixture, and under each intervention every component in a
  mixture is set to the same pooled mixture-specific value.

- **Bootstrap diagnostics**, which summarize how many bootstrap
  replications were requested, retained, or failed, together with any
  lightweight failure metadata stored in the fitted object.

- **MSM adequacy diagnostics**, which compare the exact fit-time
  counterfactual surface to the fitted MSM surface. For
  transformed-scale fits, the primary adequacy comparison is made on the
  MSM fitting scale.

Support diagnostics should not be read as a full positivity proof. They
are designed to help users see which intervention values define the
stored surface and to highlight when original-scale pooled interventions
may deserve extra scrutiny.

MSM adequacy is one of the most method-specific diagnostics in the
package. It asks whether the fitted MSM is a reasonable low-dimensional
summary of the exact fit-time surface implied by the fitted outcome
model. It is not a test of whether the outcome model is true, and it is
not a general check of causal identification.

Adequacy is evaluated on the stored fit-time grid. A good adequacy
result therefore means that the MSM tracks the exact fit-time surface
well on that grid. It does not imply that the surface is globally linear
between or beyond those intervention points.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- sim_mixture_data(
  n = 400,
  pA = 3,
  pB = 3,
  rho_within_A = 0.3,
  rho_within_B = 0.3,
  rho_between = 0.2,
  psi1 = 0.5,
  psi2 = 0.3,
  psi12 = 0.2,
  seed = 123
)

fit <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  q = 4,
  B = 100,
  seed = 13
)

# Full diagnostic bundle
diagnostics(fit)

# Focused diagnostics
support(fit)
diagnostics(fit, type = "bootstrap")
adequacy(fit)
} # }
```
