# Effect Scales, Intervals, and Prediction Scales

## Introduction

Version `0.5.0` adds explicit estimand scales for the marginal
structural model (MSM) coefficients. This matters most for binary and
count outcomes. A binary outcome can be summarized on either an
odds-ratio or risk-difference scale. A count outcome can be summarized
on either a rate-ratio or mean-difference scale.

The package now separates the following three ideas:

- the fitting scale used for coefficients, standard errors, covariance
  matrices, pooling, and confidence regions
- the display scale used in printed summaries, confidence intervals, and
  tidy output
- the response scale used for predicted risks, expected counts, and
  fitted surface plots

Those scales are related, but they are not interchangeable.

## Default estimand scales

The defaults follow the outcome family and link:

| Outcome family | Default estimand | Additive alternative |
|----|----|----|
| [`gaussian()`](https://rdrr.io/r/stats/family.html) | `"mean_difference"` | none |
| `binomial(link = "logit")` | `"odds_ratio"` | `"risk_difference"` |
| `poisson(link = "log")` | `"rate_ratio"` | `"mean_difference"` |

For binomial and Poisson fits with other links, ratio-scale MSM fitting
is not available. The package falls back to additive response-scale
summaries.

You can always make the choice explicit:

``` r

fit_or <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat_binary,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  family = binomial(),
  estimand_scale = "odds_ratio",
  B = 100,
  seed = 13
)

fit_rd <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat_binary,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  family = binomial(),
  estimand_scale = "risk_difference",
  B = 100,
  seed = 13
)
```

## Fitting scale versus display scale

For additive estimands, the fitting and display scales are the same. For
odds ratios and rate ratios, they differ.

A ratio-scale MSM is fit on a log coefficient scale. That is the scale
where standard errors, covariance matrices, Rubin pooling, and
confidence regions make statistical sense. The package therefore keeps
core machine-readable quantities on the fitting scale:

``` r

coef(fit_or)
vcov(fit_or)
tidy(fit_or)
```

For reporting, the same fitted object also gives display-scale values:

``` r

fit_or
summary(fit_or)
confint(fit_or)
tidy(fit_or, conf.int = TRUE)
```

In a ratio-scale fit, the printed coefficient table includes both the
exponentiated ratio and the fitting-scale coefficient. That is
deliberate. The ratio is easier to read, but the fitting-scale
coefficient is the one that has a standard error with the usual
large-sample interpretation.

## Binary outcomes

For binary outcomes fit with `family = binomial()`, the default v0.5.0
estimand is an odds ratio. The data generator accepts ratio parameters
on their natural scale:

``` r

dat_binary <- sim_mixture_data(
  n = 600,
  pA = 3,
  pB = 3,
  rho_within_A = 0.3,
  rho_within_B = 0.3,
  rho_between = 0.2,
  psi1 = 1.4,
  psi2 = 1.2,
  psi12 = 1.1,
  family = binomial(),
  baseline_risk = 0.2,
  seed = 123
)
```

Then fit the odds-ratio MSM:

``` r

fit_or <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat_binary,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  family = binomial(),
  B = 100,
  seed = 13
)
summary(fit_or)
```

If the scientific question is additive risk, request that scale
directly:

``` r

fit_rd <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat_binary,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  family = binomial(),
  estimand_scale = "risk_difference",
  B = 100,
  seed = 13
)
```

Risk-difference fits stay on the identity scale. Their coefficients and
display values are therefore the same scale.

## Count outcomes

For Poisson-log fits, the default v0.5.0 estimand is a rate ratio:

``` r

dat_count <- sim_mixture_data(
  n = 600,
  pA = 3,
  pB = 3,
  rho_within_A = 0.3,
  rho_within_B = 0.3,
  rho_between = 0.2,
  psi1 = 1.3,
  psi2 = 1.15,
  psi12 = 1.05,
  family = poisson(),
  baseline_rate = 1,
  seed = 123
)

fit_rr <- qgcomp.glm.multi(
  f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
  data = dat_count,
  mix1 = c("X1", "X2", "X3"),
  mix2 = c("W1", "W2", "W3"),
  family = poisson(),
  B = 100,
  seed = 13
)

summary(fit_rr)
```

If expected-count differences are the target, use
`estimand_scale = "mean_difference"`.

## Confidence intervals

For single-fit objects,
[`confint()`](https://rdrr.io/r/stats/confint.html) supports three
methods:

- `"wald"`, based on the coefficient estimate and stored bootstrap
  covariance
- `"percentile"`, based on empirical bootstrap quantiles
- `"basic"`, the reverse-percentile bootstrap interval

All three methods are computed on the MSM fitting scale. For odds-ratio
and rate-ratio estimands, the returned limits are exponentiated for
display:

``` r

confint(fit_or, method = "wald")
confint(fit_or, method = "percentile")
confint(fit_or, method = "basic")
```

Pooled multiple-imputation objects remain Wald-only. The pooled
coefficients and covariance are Rubin-pooled on the fitting scale, and
then display columns can be exponentiated for ratio summaries.

## Predictions and plots stay on the response scale

[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) answer a
different question from [`coef()`](https://rdrr.io/r/stats/coef.html).
They show predicted risks, expected counts, or means under intervention
regimes. Those are response-scale quantities:

``` r

predict(fit_or)
plot(fit_or)
```

This is true even when `fit_or` reports odds-ratio coefficients. A
surface of predicted risks is usually what users need to understand the
fitted intervention response.

For a direct comparison of two regimes on the fitted estimand scale,
request an MSM contrast:

``` r

predict(
  fit_or,
  type = "msm_contrast",
  from = c(psi1 = 0, psi2 = 0),
  to = c(psi1 = 3, psi2 = 3),
  contrast_scale = "estimand",
  interval = TRUE,
  method = "basic"
)
```

For an odds-ratio fit, that contrast is reported as an odds ratio. For a
rate-ratio fit, it is reported as a rate ratio.

## Confidence regions

[`confregion()`](https://cmowrer13.github.io/qgcomp.multi/reference/confregion.md)
builds confidence regions for selected MSM coefficients in single-fit
objects:

``` r

region <- confregion(fit_or, parm = c("psi1", "psi2"))
region
plot(region)
```

The region is based on the bootstrap covariance matrix and a chi-squared
cutoff. Its geometry is always on the fitting coefficient scale. For
odds-ratio and rate-ratio fits, that means the log coefficient scale.
The package does not draw nonlinear exponentiated ratio-scale ellipses
in v0.5.0.

Pooled MI confidence regions are not supported in this release.

## Practical notes

A few rules prevent most interpretation mistakes:

- Use [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) when downstream software
  needs fitting-scale coefficients and covariance matrices.
- Use print, [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and display
  columns from `tidy()` when writing ratio summaries for applied
  audiences.
- Use [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) when the
  question is about predicted risks, expected counts, or means.
- Use `predict(..., type = "msm_contrast", contrast_scale = "estimand")`
  when the question is a specific regime comparison on the odds-ratio or
  rate-ratio scale.
- Use
  [`confregion()`](https://cmowrer13.github.io/qgcomp.multi/reference/confregion.md)
  for single-fit coefficient regions, and remember that those regions
  live on the fitting scale.

## Related articles

- [Applied Workflow for Two-Mixture Quantile
  g-Computation](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-workflow.md)
- [Multiple Imputation and Parallel
  Workflows](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-mi-parallel.md)
- [Clustered Data and Repeated
  Measures](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-clustered-data.md)
- [Diagnostics and Sensitivity
  Checks](https://cmowrer13.github.io/qgcomp.multi/articles/qgcompmulti-diagnostics-sensitivity.md)
