# Extract MSM adequacy diagnostics from a qgcompmulti fit

Extract MSM adequacy diagnostics from a qgcompmulti fit

## Usage

``` r
adequacy(object, ...)
```

## Arguments

- object:

  A fitted `"qgcompmulti"` object.

- ...:

  Unused.

## Value

An object of class `"qgcompmulti_adequacy_diagnostic"`.

## Details

The adequacy diagnostic evaluates how closely the fitted MSM reproduces
the exact fit-time counterfactual surface stored in the fitted object.
This is a diagnostic of MSM approximation on the stored fit-time grid,
not a test of whether the outcome model is correctly specified relative
to the true data-generating process.
