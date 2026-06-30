# Plot a qgcompmulti confidence region

Draws the two-parameter confidence-region ellipsoid produced by
[`confregion()`](https://cmowrer13.github.io/qgcomp.multi/reference/confregion.md).

## Usage

``` r
# S3 method for class 'qgcompmulti_region'
plot(
  x,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  col = "#2f6f9f",
  lwd = 2,
  pch = 19,
  point_col = "#1f2933",
  xlim = NULL,
  ylim = NULL,
  asp = 1,
  ...
)
```

## Arguments

- x:

  A `qgcompmulti_region` object returned by
  [`confregion()`](https://cmowrer13.github.io/qgcomp.multi/reference/confregion.md).

- xlab:

  Optional x-axis label.

- ylab:

  Optional y-axis label.

- main:

  Optional plot title.

- col:

  Line color for the ellipsoid boundary.

- lwd:

  Line width for the ellipsoid boundary.

- pch:

  Plotting character for the point estimate.

- point_col:

  Color for the point estimate.

- xlim:

  Optional x-axis limits.

- ylim:

  Optional y-axis limits.

- asp:

  Numeric aspect ratio passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).
  The default `1` preserves the ellipsoid geometry on the fitting
  coefficient scale.

- ...:

  Additional graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.

## Details

Direct plotting is supported only for two-parameter confidence regions.
For higher-dimensional regions, inspect the returned object fields
directly. Ratio-estimand regions are plotted on the log coefficient
scale used for fitting and covariance estimation.
