# Extract intervention support diagnostics from a qgcompmulti fit

Extract intervention support diagnostics from a qgcompmulti fit

## Usage

``` r
support(object, ...)
```

## Arguments

- object:

  A fitted `"qgcompmulti"` object.

- ...:

  Unused.

## Value

An object of class `"qgcompmulti_support_diagnostic"`.

## Details

The support diagnostic is particularly informative for `q = NULL`, where
the stored intervention grid is defined by pooled percentile values
rather than quantile indices. In that setting, every component in a
mixture is set to the same pooled mixture-specific intervention value,
so users should inspect the resulting grid for scientific plausibility
as well as numerical range.
