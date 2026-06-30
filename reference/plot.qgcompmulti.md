# Plot a qgcompmulti fit

Produces a base-graphics display of the fitted marginal structural model
(MSM) surface from a
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
fit. The default display is a heatmap over the stored fit-time MSM grid.
Contour plotting is also supported, along with a slice-based interval
display for MSM predictions.

## Usage

``` r
# S3 method for class 'qgcompmulti'
plot(
  x,
  style = c("heatmap", "contour"),
  scale = "response",
  grid = NULL,
  interval = FALSE,
  slice = NULL,
  level = 0.95,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A fitted `"qgcompmulti"` object.

- style:

  Character string specifying the surface display style. Supported
  values are `"heatmap"` and `"contour"`.

- scale:

  Character string specifying the plotted outcome scale. Version `0.5.0`
  supports only `"response"` for plotting. Transformed-scale surface
  plots are intentionally out of scope; use
  [`predict.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/predict.qgcompmulti.md)
  with `type = "msm_contrast"` and `contrast_scale = "estimand"` for
  direct ratio-scale regime comparisons.

- grid:

  Optional user-specified MSM grid with columns `psi1` and `psi2`. If
  omitted, the stored fit-time MSM grid is used. User-supplied values
  are interpreted on the MSM coding scale.

- interval:

  Logical; if `TRUE`, produces a slice-based line display with bootstrap
  intervals instead of a 2D surface plot.

- slice:

  Optional list with elements `var` and `value` describing the fixed MSM
  coordinate for interval plotting. For example,
  `list(var = "psi2", value = 1)` fixes `psi2` and varies `psi1` over
  the stored MSM-grid support.

- level:

  Confidence level used when `interval = TRUE`.

- xlab:

  Optional x-axis label.

- ylab:

  Optional y-axis label.

- main:

  Optional plot title.

- ...:

  Additional graphical parameters passed to the underlying base graphics
  call.

## Value

Invisibly returns the structured prediction result used to draw the
plot.

## Details

Plotting in `qgcomp.multi` is intentionally prediction-driven: the
plotting method delegates the scientific computation to
[`predict.qgcompmulti()`](https://cmowrer13.github.io/qgcomp.multi/reference/predict.qgcompmulti.md)
and only handles rendering. This helps keep the plotted quantities
aligned with the documented public prediction interface.

The default `plot(fit)` call visualizes the fitted MSM surface over the
stored fit-time MSM grid. For quantized fits, the heatmap axes
correspond to the quantile-index intervention levels.

For stored `q = NULL` heatmaps, the plotted axis labels are mapped back
to the pooled percentile intervention values used at fit time, even when
the MSM itself was fit on centered coordinates. This means there are two
scales in play:

- the **intervention-value scale**, which is what users see on the
  stored heatmap axes for `q = NULL`; and

- the **MSM coding scale**, which is what `grid` and `slice$value` use
  for user-supplied prediction or plotting requests.

`interval = TRUE` switches to a slice-based line display rather than
trying to overlay uncertainty on the full two-dimensional surface. This
keeps the uncertainty display readable and aligned with the currently
supported MSM interval calculations.

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

# Default heatmap of the fitted MSM surface
plot(fit)

# Contour rendering of the same stored MSM surface
plot(fit, style = "contour")

# Slice-based interval display
plot(
  fit,
  interval = TRUE,
  slice = list(var = "psi2", value = 1)
)
} # }
```
