
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ILSAstats

<!-- badges: start -->

<!-- badges: end -->

Calculates point estimates and standard errors using replicate weights
and plausible values for International Large-Scale Assessments (ILSA),
including: means, proportions, quantiles, correlations, singlelevel
regressions, and multilevel regressions.

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/ILSAstats)](https://CRAN.R-project.org/package=ILSAstats) -->

<!-- [![](https://img.shields.io/github/r-package/v/dopatendo/ILSAstats)](https://github.com/dopatendo/ILSAstats) -->

<!-- [![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) -->

<!-- ![Static Badge](https://img.shields.io/badge/dependencies-haven-brightgreen) -->

<!-- [![](https://img.shields.io/badge/doi-10.32614/CRAN.package.ILSAstats-green.svg)](https://doi.org/10.32614/CRAN.package.ILSAstats) -->

<!-- ![![](http://cranlogs.r-pkg.org/badges/grand-total/ILSAstats?color=blue)](https://cran.r-project.org/package=ILSAstats)-->

<!-- badges: end -->

## Installation

You can install the stable version of `ILSAstats` directly from CRAN:

``` r
install.packages("ILSAstats")
```

Or, if you wish to install the development version of `ILSAstats`:

``` r
remotes::install_github("dopatendo/ILSAstats")
```

## Create replicate weights

Before working with replicate weights statistics we need to be sure to
create the replicate weights (if they are not provided within the data).
For an example on how to create them see `vignette("repweights")`.

## Setup

When doing multiple estimations we may need to use the same options more
than one time. For an example on how to create a setup object for
`ILSAstats` see `vignette("setup")`.

## Means

For an example on estimating means using replicate weights and plausible
values see `vignette("means")`.

## Proportions

For an example on estimating proportions using replicate weights and
plausible values see `vignette("proportions")`.

## Quantiles

For an example on estimating quantiles using replicate weights and
plausible values see `vignette("quantiles")`.

## Correlations

For an example on estimating correlations using replicate weights and
plausible values see `vignette("correlations")`.

## League tables

For an example on estimating league tables automatically see
`vignette("leaguetables")`.
