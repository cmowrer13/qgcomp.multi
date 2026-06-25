# Extract stored fit-time surfaces from a qgcompmulti fit

Returns the fit-time exact counterfactual surface, the fit-time MSM
fitted surface, or the stored exact-versus-MSM comparison object
retained in a fitted
[`qgcomp.glm.multi()`](https://cmowrer13.github.io/qgcomp.multi/reference/qgcomp.glm.multi.md)
object.

## Usage

``` r
fitted_surface(object, type = c("all", "exact", "msm", "comparison"), ...)
```

## Arguments

- object:

  A fitted `"qgcompmulti"` object.

- type:

  Character string indicating which stored surface object to return.
  Supported values are `"all"`, `"exact"`, `"msm"`, and `"comparison"`.

- ...:

  Unused.

## Value

Either a list of stored surfaces or one stored surface data frame,
depending on `type`.
