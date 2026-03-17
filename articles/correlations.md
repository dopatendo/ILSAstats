# Correlations

We can estimate the correlations for any pair of variables using the
function
[`reprho()`](https://dopatendo.github.io/ILSAstats/reference/reprho.md),
as any other “rep” function of `ILSAstats`, we need to specify the data
(`df`), the total weights (`wt`), the replicate weights (`repwt`), and
the method (`method`).

Besides these basic options, other arguments can be used:

- `x`: a string with the name of the variable (or variables) to be used
  in the analysis.
- `pv`: a string containing the name of plausible value variables
  related to a construct.
- `pv2`: a string containing the name of plausible value variables
  related to **another** construct (diferente from the one in `pv`).
- `relatedpvs`: a logical value indicating if when using two plausible
  value constructs, there should be related or not. If `TRUE`
  correlations between plausible values will be estimated in pairs (1
  with 1, 2 with 2, etc.). If `FALSE` correlations will be estimated
  using all posible combination between both plausible value variables.
- `rho`: a string indicating the correlation coefficient to be computed,
  options are `"pearson"` (the default), `"spearman"`, and
  `"polychoric"`.
- `group`: a string containing the name of the variable that contains
  the groups of countries. If used all statistics will be estimated
  separately for each group, and groups will be treated as
  **independent** from each other, e.g., countries.
- `exclude`: a string containing which groups should be excluded from
  aggregate estimations.
- `aggregates`: a string containing the aggregate statistics that should
  be estimated. Options include: `"pooled"` for also estimating all
  groups (without exclusions) as a single group; and `"composite"` for
  averaging all the estimations for each single group (without
  exclusions).

## Weights and setup

For
[`reprho()`](https://dopatendo.github.io/ILSAstats/reference/reprho.md),
first we need to create the replicate weights. Using the included
`repdata` data, and using the `"LANA"` method:

``` r
RW <- repcreate(df = repdata,
                wt = "wt",
                jkzone = "jkzones",
                jkrep = "jkrep",
                method = "LANA")
```

To make it easier to specify some arguments, it is advised that we
create also a `"repsetup"` object. We will create three setups for this
example: one without groups, one with groups and without exclusions, and
one with groups and exclusions (excluding group 2):

``` r
# No groups
STNG <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "LANA")

# With groups
STGR <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "LANA",
                 group = "GROUP")

# With groups and exclusions
STGE <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "LANA",
                 group = "GROUP", exclude = "GR2")
```

## Two non-PV variables

For example, if we want to estimate the correlation between `"item01"`
and `"SES"`, we can use either of the setups to get the overall or group
results (notice that if we do not specify the type of correlation it
will default to `"pearson"`):

``` r
# No groups
reprho(x = c("SES","item01"), setup = STNG)
```

    ## More than one correlation method provided. Only first method will be used (pearson).

    ##   variable1 variable2      rho      se    n   tvalue  pvalue
    ## 1       SES    item01 -0.03101 0.01495 4483 -2.07395 0.03814

``` r
# With groups
reprho(x = c("item01","SES","item01"), setup = STGR, rho = "pearson")
```

    ##    variable1 variable2     group      rho      se    n        tvalue  pvalue
    ## 1     item01       SES    Pooled -0.03101 0.01495 4483 -2.073950e+00 0.03814
    ## 2     item01       SES Composite -0.02555 0.01444   NA -1.769220e+00      NA
    ## 3     item01       SES       GR1  0.01100 0.02541 1507  4.326500e-01 0.66533
    ## 4     item01       SES       GR2 -0.01740 0.02759 1491 -6.305000e-01 0.52847
    ## 5     item01       SES       GR3 -0.07026 0.02169 1485 -3.239700e+00 0.00122
    ## 6     item01    item01    Pooled  1.00000 0.00000 4483  1.313835e+15 0.00000
    ## 7     item01    item01 Composite  1.00000 0.00000   NA  9.421999e+14      NA
    ## 8     item01    item01       GR1  1.00000 0.00000 1507  5.796038e+14 0.00000
    ## 9     item01    item01       GR2  1.00000 0.00000 1491  4.870536e+14 0.00000
    ## 10    item01    item01       GR3  1.00000 0.00000 1485  5.826273e+14 0.00000
    ## 11       SES    item01    Pooled -0.03101 0.01495 4483 -2.073950e+00 0.03814
    ## 12       SES    item01 Composite -0.02555 0.01444   NA -1.769220e+00      NA
    ## 13       SES    item01       GR1  0.01100 0.02541 1507  4.326500e-01 0.66533
    ## 14       SES    item01       GR2 -0.01740 0.02759 1491 -6.305000e-01 0.52847
    ## 15       SES    item01       GR3 -0.07026 0.02169 1485 -3.239700e+00 0.00122

``` r
# With groups and exclusions
reprho(x = c("SES","item01"), setup = STGE, rho = "pearson")
```

    ##   variable1 variable2     group      rho      se    n   tvalue  pvalue
    ## 1       SES    item01    Pooled -0.03800 0.01641 2992 -2.31519 0.02067
    ## 2       SES    item01 Composite -0.02963 0.01670   NA -1.77390      NA
    ## 3       SES    item01       GR1  0.01100 0.02541 1507  0.43265 0.66533
    ## 4       SES    item01       GR2 -0.01740 0.02759 1491 -0.63050 0.52847
    ## 5       SES    item01       GR3 -0.07026 0.02169 1485 -3.23970 0.00122

We can notice that using no groups we would get the same results for the
pooled estimates if we use groups and no exclusions. But, when we
exclude group 2, the pooled and the composite estimate changes.

## Correlation methods

Besides the Pearson correlation, also the Spearman correlation and
polychoric correlation can be obtained, for this we use the `rho`
argument:

``` r
# Pearson
reprho(x = c("item02","item01"), setup = STGR, rho = "pearson")
```

    ##   variable1 variable2     group     rho      se    n   tvalue pvalue
    ## 1    item02    item01    Pooled 0.55841 0.01523 4029 36.65308      0
    ## 2    item02    item01 Composite 0.55884 0.01519   NA 36.79579     NA
    ## 3    item02    item01       GR1 0.55507 0.02625 1344 21.14770      0
    ## 4    item02    item01       GR2 0.54681 0.02690 1336 20.32532      0
    ## 5    item02    item01       GR3 0.57464 0.02575 1349 22.31239      0

``` r
# Spearman
reprho(x = c("item02","item01"), setup = STGR, rho = "spearman")
```

    ##   variable1 variable2     group     rho      se    n   tvalue pvalue
    ## 1    item02    item01    Pooled 0.60639 0.01121 4029 54.10745      0
    ## 2    item02    item01 Composite 0.60613 0.01127   NA 53.79030     NA
    ## 3    item02    item01       GR1 0.59916 0.01765 1344 33.93823      0
    ## 4    item02    item01       GR2 0.60349 0.02162 1336 27.91022      0
    ## 5    item02    item01       GR3 0.61574 0.01907 1349 32.29232      0

``` r
# Polychoric
reprho(x = c("item02","item01"), setup = STGR, rho = "polychoric")
```

    ##   variable1 variable2     group     rho      se    n   tvalue pvalue
    ## 1    item02    item01    Pooled 0.70979 0.01237 4029 57.40152      0
    ## 2    item02    item01 Composite 0.71015 0.01261   NA 56.31031     NA
    ## 3    item02    item01       GR1 0.70916 0.02066 1344 34.31764      0
    ## 4    item02    item01       GR2 0.70233 0.02233 1336 31.45694      0
    ## 5    item02    item01       GR3 0.71897 0.02249 1349 31.96448      0

We can notice that using no groups we would get the same results for the
pooled estimates if we use groups and no exclusions. But, when we
exclude group 2, the pooled and the composite estimate changes.

## Multiple non-PV variables

We can also estimate correlations between more than two variables:

``` r
reprho(x = c("SES","item01","item02"), setup = STGR, rho = "pearson")
```

    ##    variable1 variable2     group      rho      se    n   tvalue  pvalue
    ## 1        SES    item01    Pooled -0.03101 0.01495 4483 -2.07395 0.03814
    ## 2        SES    item01 Composite -0.02555 0.01444   NA -1.76922      NA
    ## 3        SES    item01       GR1  0.01100 0.02541 1507  0.43265 0.66533
    ## 4        SES    item01       GR2 -0.01740 0.02759 1491 -0.63050 0.52847
    ## 5        SES    item01       GR3 -0.07026 0.02169 1485 -3.23970 0.00122
    ## 6        SES    item02    Pooled -0.02674 0.01481 4501 -1.80551 0.07106
    ## 7        SES    item02 Composite -0.01879 0.01408   NA -1.33494      NA
    ## 8        SES    item02       GR1  0.01774 0.02317 1496  0.76558 0.44405
    ## 9        SES    item02       GR2 -0.00231 0.02877 1494 -0.08025 0.93605
    ## 10       SES    item02       GR3 -0.07181 0.02046 1511 -3.50960 0.00046
    ## 11    item01    item02    Pooled  0.55841 0.01523 4029 36.65308 0.00000
    ## 12    item01    item02 Composite  0.55884 0.01519   NA 36.79579      NA
    ## 13    item01    item02       GR1  0.55507 0.02625 1344 21.14770 0.00000
    ## 14    item01    item02       GR2  0.54681 0.02690 1336 20.32532 0.00000
    ## 15    item01    item02       GR3  0.57464 0.02575 1349 22.31239 0.00000

## One PV variable

To estimate the correlations between a non plausible value variables and
a plausible value variable we need to state the normal variables in `x`
and plausible values variables in argument `pv`:

``` r
# One variable
reprho(x = c("SES"),
       pv = c(paste0("Math",1:5)),
       setup = STGR, rho = "pearson")
```

    ##   variable1 variable2     group     rho      se    n  tvalue  pvalue
    ## 1       SES       PVs    Pooled 0.37372 0.05132 5000 7.28203 0.00000
    ## 2       SES       PVs Composite 0.24538 0.03991   NA 6.14846      NA
    ## 3       SES       PVs       GR1 0.27224 0.07236 1667 3.76237 0.00017
    ## 4       SES       PVs       GR2 0.21660 0.06682 1666 3.24165 0.00121
    ## 5       SES       PVs       GR3 0.24730 0.06808 1667 3.63272 0.00029

``` r
# More than one variable
reprho(x = c("SES","item01"),
       pv = c(paste0("Math",1:5)),
       setup = STGR, rho = "pearson")
```

    ##    variable1 variable2     group      rho      se    n   tvalue  pvalue
    ## 1        SES       PVs    Pooled  0.37372 0.05132 5000  7.28203 0.00000
    ## 2        SES       PVs Composite  0.24538 0.03991   NA  6.14846      NA
    ## 3        SES       PVs       GR1  0.27224 0.07236 1667  3.76237 0.00017
    ## 4        SES       PVs       GR2  0.21660 0.06682 1666  3.24165 0.00121
    ## 5        SES       PVs       GR3  0.24730 0.06808 1667  3.63272 0.00029
    ## 6     item01       PVs    Pooled -0.01759 0.01840 4483 -0.95618 0.33903
    ## 7     item01       PVs Composite -0.00954 0.01742   NA -0.54757      NA
    ## 8     item01       PVs       GR1 -0.01989 0.03014 1507 -0.66008 0.50931
    ## 9     item01       PVs       GR2 -0.00214 0.03088 1491 -0.06918 0.94485
    ## 10    item01       PVs       GR3 -0.00659 0.02948 1485 -0.22340 0.82325

## Multiple PV variables

It is also possible to correlate two plausible value variables using
argument `pv` and `pv2`, when doing so `x` should be `NULL`. For
example, to correlate math and reading achievement in `repdata`:

``` r
reprho(pv = c(paste0("Math",1:5)),
       pv2 = c(paste0("Reading",1:5)),
       setup = STGR, rho = "pearson")
```

    ##       group      rho      se    n   tvalue  pvalue
    ## 1    Pooled  0.14209 0.03848 5000  3.69259 0.00022
    ## 2 Composite -0.14197 0.03827   NA -3.70955      NA
    ## 3       GR1 -0.07250 0.06419 1667 -1.12932 0.25893
    ## 4       GR2 -0.22991 0.07829 1666 -2.93673 0.00336
    ## 5       GR3 -0.12351 0.05416 1667 -2.28063 0.02270

Please notice that by default
[`reprho()`](https://dopatendo.github.io/ILSAstats/reference/reprho.md)
assumes that both plausible value variables are related and correlates
the first plausible value of each variable, then the seconde one of each
and so on. For our example it will estimate 5 correlations
(Math1-Reading1, Math2-Reading2, Math3-Reading3, Math4-Reading4, and
Math5-Reading5) and average them.

Nevertheless, it is also possible to calculate de correlation between no
related plausible values, therefore instead of 5 estimations,
[`reprho()`](https://dopatendo.github.io/ILSAstats/reference/reprho.md)
will make 25 estimations, with all the possible combinations between
math and reading. For doing that, we can use the argument
`relatedpvs=FALSE`:

``` r
reprho(pv = c(paste0("Math",1:5)),
       pv2 = c(paste0("Reading",1:5)),
       relatedpvs = FALSE,
       setup = STGR, rho = "pearson")
```

    ##       group      rho      se    n   tvalue  pvalue
    ## 1    Pooled  0.15895 0.03956 5000  4.01755 0.00006
    ## 2 Composite -0.11932 0.03486   NA -3.42305      NA
    ## 3       GR1 -0.04893 0.06081 1667 -0.80465 0.42114
    ## 4       GR2 -0.20752 0.06326 1666 -3.28045 0.00106
    ## 5       GR3 -0.10153 0.05689 1667 -1.78448 0.07453

## Aggregates

When using groups we can always omit the pooled and composite
calculations if we need to, by default both estimates will be
calculated.

``` r
# Default
reprho(pv = paste0("Math",1:5), 
       pv2 = c(paste0("Reading",1:5)), 
       setup = STGR, rho = "pearson")
```

    ##       group      rho      se    n   tvalue  pvalue
    ## 1    Pooled  0.14209 0.03848 5000  3.69259 0.00022
    ## 2 Composite -0.14197 0.03827   NA -3.70955      NA
    ## 3       GR1 -0.07250 0.06419 1667 -1.12932 0.25893
    ## 4       GR2 -0.22991 0.07829 1666 -2.93673 0.00336
    ## 5       GR3 -0.12351 0.05416 1667 -2.28063 0.02270

``` r
# Only pooled
reprho(pv = paste0("Math",1:5), 
       pv2 = c(paste0("Reading",1:5)), 
       setup = STGR, rho = "pearson",
       aggregates = "pooled")
```

    ##    group      rho      se    n   tvalue  pvalue
    ## 1 Pooled  0.14209 0.03848 5000  3.69259 0.00022
    ## 2    GR1 -0.07250 0.06419 1667 -1.12932 0.25893
    ## 3    GR2 -0.22991 0.07829 1666 -2.93673 0.00336
    ## 4    GR3 -0.12351 0.05416 1667 -2.28063 0.02270

``` r
# Only composite
reprho(pv = paste0("Math",1:5), 
       pv2 = c(paste0("Reading",1:5)), 
       setup = STGR, rho = "pearson",
       aggregates = "composite")
```

    ##       group      rho      se    n   tvalue  pvalue
    ## 1 Composite -0.14197 0.03827   NA -3.70955      NA
    ## 2       GR1 -0.07250 0.06419 1667 -1.12932 0.25893
    ## 3       GR2 -0.22991 0.07829 1666 -2.93673 0.00336
    ## 4       GR3 -0.12351 0.05416 1667 -2.28063 0.02270

``` r
# No aggregates
reprho(pv = paste0("Math",1:5), 
       pv2 = c(paste0("Reading",1:5)), 
       setup = STGR, rho = "pearson",
       aggregates = NULL)
```

    ##   group      rho      se    n   tvalue  pvalue
    ## 1   GR1 -0.07250 0.06419 1667 -1.12932 0.25893
    ## 2   GR2 -0.22991 0.07829 1666 -2.93673 0.00336
    ## 3   GR3 -0.12351 0.05416 1667 -2.28063 0.02270
