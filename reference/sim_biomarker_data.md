# Simulate a biomarker dataset for two-mixture inflammation analyses

Generates a narrative environmental-health teaching dataset with two
correlated exposure mixtures, a continuous inflammation outcome, and a
small amount of intentional model imperfection. The simulated biomarkers
are positive, right-skewed concentrations that naturally motivate log
transformation before model fitting.

## Usage

``` r
sim_biomarker_data(n = 800, include_log = FALSE, seed = NULL)
```

## Arguments

- n:

  Number of observations to generate.

- include_log:

  Logical; if `TRUE`, also include log-transformed biomarker variables
  named `ln_mep`, `ln_mibp`, `ln_mehpp`, `ln_bpa`, `ln_bp3`, and
  `ln_ppb`.

- seed:

  Optional integer random seed. When supplied, the generated dataset is
  reproducible and the caller's RNG state is restored on exit.

## Value

A data frame containing `inflammation_score`, six raw biomarker
concentrations, and the covariates `age_years`, `bmi_kg_m2`, and
`smoker`. If `include_log = TRUE`, log-transformed biomarker variables
are included as well.

## Details

Mixture 1 represents urinary phthalate metabolites:

- `mep_ng_ml`

- `mibp_ng_ml`

- `mehpp_ng_ml`

Mixture 2 represents urinary phenol biomarkers:

- `bpa_ng_ml`

- `bp3_ng_ml`

- `ppb_ng_ml`

The outcome `inflammation_score` is generated from weighted mixture
summaries on the log-concentration scale, a moderate positive
interaction between the two mixtures, common epidemiologic covariates,
and a small amount of deliberate nonlinearity. In the current
implementation, that nonlinearity comes from a centered cubic term in
the phthalate mixture score and a mild threshold term in the phenol
mixture score. This design makes the fitted MSM informative without
making the exact fit-time surface unrealistically perfect.

This helper is intended for examples, documentation, and teaching. The
exposure concentrations are simulated as correlated log-normal variables
with moderate within-mixture correlation, mild between-mixture
correlation, and deliberately different raw scales across chemicals.

Those scale differences are useful for demonstrating:

- standard quantized fitting,

- original-scale fitting with `q = NULL`,

- why log transformation can be scientifically natural for biomarker
  concentrations, and

- why median centering can make the `q = NULL` intercept more
  interpretable.

## Examples

``` r
dat <- sim_biomarker_data(n = 200, seed = 123)
head(dat)
#>   inflammation_score mep_ng_ml mibp_ng_ml mehpp_ng_ml bpa_ng_ml bp3_ng_ml
#> 1          2.8556787 162.80984  47.955682    84.08285  1.024531 19.225320
#> 2          1.9842073 106.39914  35.957787    33.15359  2.247171 24.590296
#> 3         -0.1851046  22.32440   5.774262    17.97875  1.038477 11.498077
#> 4          3.8436727  57.07891  12.301756    47.24134  4.087454  8.034531
#> 5          3.9717934  21.43953  24.179712    23.57163  1.380142 15.827428
#> 6          1.5624084  28.84729   3.866986    11.30165  1.280904 18.101201
#>   ppb_ng_ml age_years bmi_kg_m2 smoker
#> 1 2.6612586        54  22.02635      0
#> 2 1.3953884        36  24.34083      0
#> 3 0.9726333        57  31.93167      0
#> 4 1.8595875        36  25.62586      0
#> 5 5.8318242        54  29.25188      0
#> 6 0.5612947        60  35.12220      0

dat_log <- sim_biomarker_data(n = 200, include_log = TRUE, seed = 123)
head(dat_log[c("ln_mep", "ln_bpa", "inflammation_score")])
#>     ln_mep     ln_bpa inflammation_score
#> 1 5.092583 0.02423504          2.8556787
#> 2 4.667197 0.80967187          1.9842073
#> 3 3.105680 0.03775492         -0.1851046
#> 4 4.044435 1.40792233          3.8436727
#> 5 3.065236 0.32218666          3.9717934
#> 6 3.362016 0.24756595          1.5624084
```
