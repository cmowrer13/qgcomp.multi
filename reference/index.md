# Package index

## Model fitting

- [`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
  : Fit a quantile g-computation model for two exposure mixtures
- [`qgcomp.glm.multi.mi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.mi.md)
  : Fit a pooled multiple-imputation quantile g-computation model
- [`qgcompmulti_msm_fit()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti_msm_fit.md)
  : Fit the core outcome-model and MSM components for two-mixture
  g-computation

## Prediction and intervention surfaces

- [`predict(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/predict.qgcompmulti.md)
  : Predict from a qgcompmulti fit
- [`fitted_surface()`](https://cmowrer13.github.io/qgcomp.multi/reference/fitted_surface.md)
  : Extract stored fit-time surfaces from a qgcompmulti fit
- [`intervention_grid()`](https://cmowrer13.github.io/qgcomp.multi/reference/intervention_grid.md)
  : Extract stored intervention grids from a qgcompmulti fit
- [`plot(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/plot.qgcompmulti.md)
  : Plot a qgcompmulti fit

## Diagnostics

- [`support()`](https://cmowrer13.github.io/qgcomp.multi/reference/support.md)
  : Extract intervention support diagnostics from a qgcompmulti fit
- [`diagnostics()`](https://cmowrer13.github.io/qgcomp.multi/reference/diagnostics.md)
  : Diagnostic summaries for a qgcompmulti fit
- [`adequacy()`](https://cmowrer13.github.io/qgcomp.multi/reference/adequacy.md)
  : Extract MSM adequacy diagnostics from a qgcompmulti fit

## Sensitivity analysis

- [`mcsize_sensitivity()`](https://cmowrer13.github.io/qgcomp.multi/reference/mcsize_sensitivity.md)
  : Sensitivity to Monte Carlo size in qgcompmulti fits
- [`b_sensitivity()`](https://cmowrer13.github.io/qgcomp.multi/reference/b_sensitivity.md)
  : Sensitivity to bootstrap iteration count in qgcompmulti fits
- [`q_sensitivity()`](https://cmowrer13.github.io/qgcomp.multi/reference/q_sensitivity.md)
  : Sensitivity to quantization choice in qgcompmulti fits

## Reporting and extractors

- [`summary(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md)
  [`print(`*`<summary.qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti.md)
  : Summarize a qgcompmulti fit
- [`summary(`*`<qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti_mi.md)
  [`print(`*`<summary.qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/summary.qgcompmulti_mi.md)
  : Summarize a pooled qgcompmulti multiple-imputation fit
- [`print(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti.md)
  : Print a qgcompmulti fit
- [`print(`*`<qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/print.qgcompmulti_mi.md)
  : Print a pooled qgcompmulti multiple-imputation fit
- [`confint(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/confint.qgcompmulti.md)
  : Wald confidence intervals for qgcompmulti coefficients
- [`tidy(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/tidy.qgcompmulti.md)
  : Tidy coefficient summaries for qgcompmulti fits
- [`tidy(`*`<qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/tidy.qgcompmulti_mi.md)
  : Tidy coefficient summaries for pooled qgcompmulti
  multiple-imputation fits
- [`glance(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/glance.qgcompmulti.md)
  : Glance summaries for qgcompmulti fits
- [`glance(`*`<qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/glance.qgcompmulti_mi.md)
  : Glance summaries for pooled qgcompmulti multiple-imputation fits
- [`coef(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md)
  [`formula(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md)
  [`nobs(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md)
  [`df.residual(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md)
  [`residuals(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md)
  [`vcov(`*`<qgcompmulti>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-extractors.md)
  : Core extractor methods for qgcompmulti objects
- [`coef(`*`<qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md)
  [`vcov(`*`<qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md)
  [`confint(`*`<qgcompmulti_mi>`*`)`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcompmulti-mi-extractors.md)
  : Core extractor methods for pooled qgcompmulti multiple-imputation
  objects

## Data preparation and simulation

- [`quantize_mixtures()`](https://cmowrer13.github.io/qgcomp.multi/reference/quantize_mixtures.md)
  : Quantize exposure variables for two-mixture quantile g-computation
- [`sim_mixture_data()`](https://cmowrer13.github.io/qgcomp.multi/reference/sim_mixture_data.md)
  : Generate simulated dataset with two exposure mixtures
- [`sim_biomarker_data()`](https://cmowrer13.github.io/qgcomp.multi/reference/sim_biomarker_data.md)
  : Simulate a biomarker dataset for two-mixture inflammation analyses
