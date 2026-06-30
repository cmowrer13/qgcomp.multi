# Changelog

## qgcomp.multi 0.5.0

### Added

- Added estimand-scale choices for ordinary and multiple-imputation
  fits, including odds-ratio and rate-ratio summaries for supported
  family/link combinations.
- Added percentile and basic bootstrap confidence intervals for
  single-fit coefficient summaries and MSM-based prediction intervals.
- Added
  [`confregion()`](https://cmowrer13.github.io/qgcomp.multi/reference/confregion.md)
  for single-fit bootstrap-covariance chi-squared confidence regions for
  selected MSM coefficients.
- Added an effect-scales vignette explaining fitting scale, display
  scale, response-scale prediction, interval methods, and current MI
  limits.

### Changed

- Updated binomial-logit and Poisson-log defaults to use ratio estimands
  unless the user requests additive alternatives.
- Clarified reporting so ratio estimands show both the fitting-scale
  coefficient and the exponentiated ratio where appropriate.
- Kept prediction and surface plotting on the response scale by default
  while allowing direct MSM contrasts on the estimand scale.

### Documentation

- Expanded function documentation, README, workflow material, and
  pkgdown configuration for the v0.5.0 inference features.

## qgcomp.multi 0.4.0

### Added

- Added native multiple-imputation fitting and pooling through
  [`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md).
- Added support for multiply imputed inputs from completed-data lists
  and [`mice::mids`](https://amices.org/mice/reference/mids.html)
  objects.
- Added broom-style
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html) methods
  for single-fit and pooled MI objects.
- Added optional bootstrap-level parallel execution for ordinary fits
  and native multiple-imputation workflows.

### Changed

- Expanded the public workflow beyond single-fit estimation to include
  pooled MI inference and runtime-oriented parallel controls.
- Formalized the pooled MI result structure around a dedicated
  fitted-object class for summary and inference tasks.

### Documentation

- Expanded the README and package documentation to cover native
  multiple-imputation workflows, parallel execution, and broom-based
  reporting.

## qgcomp.multi 0.3.0

### Added

- Added MSM-based prediction support, including direct contrasts and
  exact fit-time surface extraction.
- Added surface plotting for fitted intervention-response surfaces.
- Added intervention support, bootstrap, and adequacy diagnostics as
  first-class public features.
- Added sensitivity helpers for Monte Carlo size, bootstrap size, and
  quantization choice.
- Added the main applied workflow vignette for end-to-end package use.

### Changed

- Extended fitted `qgcompmulti` objects to retain fit-time prediction
  information needed for prediction, plotting, and diagnostic methods.
- Reorganized the package around a clearer post-fit workflow built on
  stored intervention grids and surfaces.

### Documentation

- Expanded examples and method documentation for prediction, plotting,
  diagnostics, and sensitivity analysis.

## qgcomp.multi 0.2.0

### Added

- Introduced the structured `qgcompmulti` fitted-object class.
- Added core fitted-model methods including
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`formula()`](https://rdrr.io/r/stats/formula.html), and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html).
- Added a cleaner object contract to support future extensions in
  prediction, diagnostics, and multiple-imputation workflows.

### Changed

- Moved the package away from an earlier proof-of-concept output
  structure toward a more conventional model-object interface.
- Standardized the organization of fitted-model results and metadata.

### Documentation

- Updated the package documentation to reflect the new fitted-object
  structure and model-method workflow.
