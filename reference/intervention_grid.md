# Extract stored intervention grids from a qgcompmulti fit

Returns the fit-time intervention grid and the corresponding MSM-coded
grid retained in a fitted
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
object.

## Usage

``` r
intervention_grid(object, type = c("all", "intervention", "msm"), ...)
```

## Arguments

- object:

  A fitted `"qgcompmulti"` object.

- type:

  Character string indicating which stored grid to return. Supported
  values are `"all"`, `"intervention"`, and `"msm"`.

- ...:

  Unused.

## Value

Either a list containing both stored grids or a single data frame,
depending on `type`.
