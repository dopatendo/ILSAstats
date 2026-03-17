# ILSAstats

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
For an example on how to create them see
[`vignette("repweights")`](https://dopatendo.github.io/ILSAstats/articles/repweights.md).

## Setup

When doing multiple estimations we may need to use the same options more
than one time. For an example on how to create a setup object for
`ILSAstats` see
[`vignette("setup")`](https://dopatendo.github.io/ILSAstats/articles/setup.md).

## Means

For an example on estimating means using replicate weights and plausible
values see
[`vignette("means")`](https://dopatendo.github.io/ILSAstats/articles/means.md).

## Proportions

For an example on estimating proportions using replicate weights and
plausible values see
[`vignette("proportions")`](https://dopatendo.github.io/ILSAstats/articles/proportions.md).

## Quantiles

For an example on estimating quantiles using replicate weights and
plausible values see
[`vignette("quantiles")`](https://dopatendo.github.io/ILSAstats/articles/quantiles.md).

## Correlations

For an example on estimating correlations using replicate weights and
plausible values see
[`vignette("correlations")`](https://dopatendo.github.io/ILSAstats/articles/correlations.md).

## League tables

For an example on estimating league tables automatically see
[`vignette("leaguetables")`](https://dopatendo.github.io/ILSAstats/articles/leaguetables.md).
