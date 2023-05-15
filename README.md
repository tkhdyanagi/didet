
<!-- README.md is generated from README.Rmd. Please edit that file -->

# didet: An Effective Treatment Approach to Difference-in-Differences with General Treatment Patterns

<!-- badges: start -->
<!-- badges: end -->

The **didet** package provides tools for differences-in-differences
(DiD) estimation when the treatment variable of interest may be
non-binary and its value may change in each time period. The package
enables to implement the doubly robust estimation of the average
treatment effect for movers using the concept of effective treatment.
For more details, see [Yanagi (2023) “An effective treatment approach to
difference-in-differences with general treatment
patterns”](https://doi.org/10.48550/arXiv.2212.13226).

## Installation

Get the package from GitHub:

``` r
# install.packages("devtools") # if necessary
devtools::install_github("tkhdyanagi/didet", build_vignettes = TRUE)
```

## Vignettes

For more details, see the package vignettes with:

``` r
library("didet")

# Getting Started with the didet Package
vignette("didet")

# Review of An Effective Treatment Approach to Difference-in-Differences with General Treatment Patterns
vignette("review", package = "didet")
```

## References

- Yanagi, T., 2023. An effective treatment approach to
  difference-in-differences with general treatment patterns. arXiv
  preprint arXiv:2212.13226.
  [Link](https://doi.org/10.48550/arXiv.2212.13226)
