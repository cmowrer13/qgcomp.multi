
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qgcomp.multi

<!-- badges: start -->
<!-- badges: end -->

qgcomp.multi is an extension of Quantile g-Computation, developed by
Keil et al. to estimate the joint effect of multiple exposure mixtures
and their interaction using marginal structural models.

## Installation

You can install the development version of qgcomp.multi like so:

``` r
devtools::install_github("cmowrer13/qgcomp.multi")
```

## Example

This is a basic example to fit a Quantile g-Computation model with two
exposure mixtures:

``` r
library(qgcomp.multi)

fit <- qgcomp.glm.multi.boot(f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C, 
        data = dataset, 
        mix1 = c("X1", "X2", "X3"),
        mix2 = c("W1", "W2", "W3"))
```
