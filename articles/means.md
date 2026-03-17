# Means

We can estimate the arithmetic mean of any variable using the function
[`repmean()`](https://dopatendo.github.io/ILSAstats/reference/repmean.md),
as any other “rep” function of `ILSAstats`, we need to specify the data
(`df`), the total weights (`wt`), the replicate weights (`repwt`), and
the method (`method`).

Besides these basic options, other arguments can be used:

- `x`: a string with the name of the variable (or variables) to be used
  in the analysis.
- `PV`: a logical value indicating if `x` are plausible values or not.
  If `FALSE` the statistics of each variable will be computed
  independently. If `TRUE` all variables in `x` will be treated as
  plausible values of the same construct, and all their statistics will
  be combined.
- `var`: a string of length 1 indicating the type of variance to be
  estimated. Options are `"unbiased"` or `"ML"`.
- `group`: a string containing the name of the variable that contains
  the groups of countries. If used all statistics will be estimated
  separately for each group, and groups will be treated as
  **independent** from each other, e.g., countries.
- `by`: a string containing the name of a second grouping variable. If
  used, all statistics will be estimated separately for each category,
  and categories will be treated as **non-independent** from each other,
  e.g., boys and girls.
- `exclude`: a string containing which groups should be excluded from
  aggregate estimations.
- `aggregates`: a string containing the aggregate statistics that should
  be estimated. Options include: `"pooled"` for also estimating all
  groups (without exclusions) as a single group; and `"composite"` for
  averaging all the estimations for each single group (without
  exclusions).

## Weights and setup

For
[`repmean()`](https://dopatendo.github.io/ILSAstats/reference/repmean.md),
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

## Single variable

For example, if we want to estimate the mean of variable `"SES"`, we can
use either of the setups to get the overall or group results:

``` r
# No groups
repmean(x = "SES", setup = STNG)
```

    ##      N     mean      se      sd  sdse     var   varse
    ## 1 5000 49.96674 0.01442 1.00784 0.008 1.01574 0.01613

``` r
# With groups
repmean(x = "SES", setup = STGR)
```

    ##       group    N     mean      se      sd    sdse     var   varse
    ## 1    Pooled 5000 49.96674 0.01442 1.00784 0.00800 1.01574 0.01613
    ## 2 Composite   NA 49.96646 0.01221 0.93833 0.00926 0.88046 0.01739
    ## 3       GR1 1667 49.52108 0.02037 0.93659 0.01468 0.87720 0.02750
    ## 4       GR2 1666 49.95496 0.01965 0.93962 0.01700 0.88289 0.03195
    ## 5       GR3 1667 50.42334 0.02325 0.93878 0.01636 0.88130 0.03072

``` r
# With groups and exclusions
repmean(x = "SES", setup = STGE)
```

    ##       group    N     mean      se      sd    sdse     var   varse
    ## 1    Pooled 3334 49.97266 0.01717 1.04050 0.01099 1.08263 0.02287
    ## 2 Composite   NA 49.97221 0.01546 0.93768 0.01099 0.87925 0.02062
    ## 3       GR1 1667 49.52108 0.02037 0.93659 0.01468 0.87720 0.02750
    ## 4       GR2 1666 49.95496 0.01965 0.93962 0.01700 0.88289 0.03195
    ## 5       GR3 1667 50.42334 0.02325 0.93878 0.01636 0.88130 0.03072

We can notice that using no groups we would get the same results for the
pooled estimates if we use groups and no exclusions. But, when we
exclude group 2, the pooled and the composite estimate changes.

## Multiple variables

We can also estimate multiple variables at once, for example `"SES"` and
`"Math1"`:

``` r
# No groups
repmean(x = c("SES","Math1"), setup = STNG)
```

    ##   variable    N     mean      se      sd    sdse     var   varse
    ## 1      SES 5000 49.96674 0.01442 1.00784 0.00800 1.01574 0.01613
    ## 2    Math1 5000  0.00191 0.01718 1.02336 0.00928 1.04726 0.01899

``` r
# With groups
repmean(x = c("SES","Math1"), setup = STGR)
```

    ##    variable     group    N     mean      se      sd    sdse     var   varse
    ## 1       SES    Pooled 5000 49.96674 0.01442 1.00784 0.00800 1.01574 0.01613
    ## 2       SES Composite   NA 49.96646 0.01221 0.93833 0.00926 0.88046 0.01739
    ## 3       SES       GR1 1667 49.52108 0.02037 0.93659 0.01468 0.87720 0.02750
    ## 4       SES       GR2 1666 49.95496 0.01965 0.93962 0.01700 0.88289 0.03195
    ## 5       SES       GR3 1667 50.42334 0.02325 0.93878 0.01636 0.88130 0.03072
    ## 6     Math1    Pooled 5000  0.00191 0.01718 1.02336 0.00928 1.04726 0.01899
    ## 7     Math1 Composite   NA  0.00149 0.01283 0.89950 0.00852 0.80920 0.01532
    ## 8     Math1       GR1 1667 -0.60164 0.02192 0.90847 0.01330 0.82532 0.02417
    ## 9     Math1       GR2 1666  0.01102 0.02382 0.88556 0.01502 0.78422 0.02660
    ## 10    Math1       GR3 1667  0.59510 0.02080 0.90446 0.01583 0.81805 0.02864

## Plausible values

When treating with plausible values, we need to specify the names of all
plausible values of a construct, and use the argument `"PV"` so all
estimates will be combined (if not all variables will be estimated
separately). For example, for estimating the mean achievement in math
for this sample we would use:

``` r
# No groups
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STNG)
```

    ##      N    mean      se      sd    sdse     var   varse
    ## 1 5000 0.00524 0.01634 1.02282 0.01099 1.04619 0.02247

``` r
# With groups
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR)
```

    ##       group    N     mean      se      sd    sdse     var   varse
    ## 1    Pooled 5000  0.00524 0.01634 1.02282 0.01099 1.04619 0.02247
    ## 2 Composite   NA  0.00482 0.01970 0.90173 0.01192 0.81329 0.02150
    ## 3       GR1 1667 -0.59327 0.03943 0.90148 0.02093 0.81283 0.03766
    ## 4       GR2 1666  0.01857 0.02395 0.89383 0.02075 0.79905 0.03708
    ## 5       GR3 1667  0.58917 0.03693 0.90989 0.02028 0.82799 0.03699

## Aggregates

When using groups we can always omit the pooled and composite
calculations if we need to, by default both estimates will be
calculated.

``` r
# Default
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR)
```

    ##       group    N     mean      se      sd    sdse     var   varse
    ## 1    Pooled 5000  0.00524 0.01634 1.02282 0.01099 1.04619 0.02247
    ## 2 Composite   NA  0.00482 0.01970 0.90173 0.01192 0.81329 0.02150
    ## 3       GR1 1667 -0.59327 0.03943 0.90148 0.02093 0.81283 0.03766
    ## 4       GR2 1666  0.01857 0.02395 0.89383 0.02075 0.79905 0.03708
    ## 5       GR3 1667  0.58917 0.03693 0.90989 0.02028 0.82799 0.03699

``` r
# Only pooled
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR, aggregates = "pooled")
```

    ##    group    N     mean      se      sd    sdse     var   varse
    ## 1 Pooled 5000  0.00524 0.01634 1.02282 0.01099 1.04619 0.02247
    ## 2    GR1 1667 -0.59327 0.03943 0.90148 0.02093 0.81283 0.03766
    ## 3    GR2 1666  0.01857 0.02395 0.89383 0.02075 0.79905 0.03708
    ## 4    GR3 1667  0.58917 0.03693 0.90989 0.02028 0.82799 0.03699

``` r
# Only composite
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR, aggregates = "composite")
```

    ##       group    N     mean      se      sd    sdse     var   varse
    ## 1 Composite   NA  0.00482 0.01970 0.90173 0.01192 0.81329 0.02150
    ## 2       GR1 1667 -0.59327 0.03943 0.90148 0.02093 0.81283 0.03766
    ## 3       GR2 1666  0.01857 0.02395 0.89383 0.02075 0.79905 0.03708
    ## 4       GR3 1667  0.58917 0.03693 0.90989 0.02028 0.82799 0.03699

``` r
# No aggregates
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR, aggregates = NULL)
```

    ##   group    N     mean      se      sd    sdse     var   varse
    ## 1   GR1 1667 -0.59327 0.03943 0.90148 0.02093 0.81283 0.03766
    ## 2   GR2 1666  0.01857 0.02395 0.89383 0.02075 0.79905 0.03708
    ## 3   GR3 1667  0.58917 0.03693 0.90989 0.02028 0.82799 0.03699

## Difference between non-independent groups

For estimating the mean of not-independent groups, we can use the
argument `by`. For example, for estimating the mean achievement in math
between `GENDER==0` and `GENDER==1`, we can:

``` r
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STNG,by = "GENDER")
```

    ## $ALL
    ##      N    mean      se      sd    sdse     var   varse
    ## 1 5000 0.00524 0.01634 1.02282 0.01099 1.04619 0.02247
    ## 
    ## $`GENDER==0`
    ##      N    mean      se      sd   sdse     var   varse   mean_1 meandiff_1
    ## 1 2483 0.45427 0.03749 0.93268 0.0166 0.86999 0.03097 -0.43552    0.88979
    ##   meandiffse_1 tvalue_1 df_1 pvalue_1
    ## 1        0.079 11.26338 4998        0
    ## 
    ## $`GENDER==1`
    ##      N     mean      se      sd    sdse     var   varse  mean_0 meandiff_0
    ## 1 2517 -0.43552 0.04587 0.90835 0.02226 0.82533 0.04067 0.45427   -0.88979
    ##   meandiffse_0  tvalue_0 df_0 pvalue_0
    ## 1        0.079 -11.26338 4998        0

This will provide us with the overall statistics, the statistics for
each category within `GENDER`, and statistics of the difference between
means between all categories.

Of course, we can also estimate this using groups:

``` r
repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR, aggregates = NULL,
        by = "GENDER")
```

    ## $ALL
    ##   group    N     mean      se      sd    sdse     var   varse
    ## 1   GR1 1667 -0.59327 0.03943 0.90148 0.02093 0.81283 0.03766
    ## 2   GR2 1666  0.01857 0.02395 0.89383 0.02075 0.79905 0.03708
    ## 3   GR3 1667  0.58917 0.03693 0.90989 0.02028 0.82799 0.03699
    ## 
    ## $`GENDER==0`
    ##   group   N     mean      se      sd    sdse     var   varse   mean_1
    ## 1   GR1 819 -0.17761 0.04462 0.79348 0.03070 0.62996 0.04868 -0.99226
    ## 2   GR2 846  0.48279 0.03513 0.78112 0.03694 0.61083 0.05735 -0.45673
    ## 3   GR3 818  1.05458 0.06162 0.78595 0.03401 0.61826 0.05327  0.14141
    ##   meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
    ## 1    0.81465      0.09313  8.74764 1665        0
    ## 2    0.93952      0.07148 13.14385 1664        0
    ## 3    0.91317      0.08944 10.20971 1665        0
    ## 
    ## $`GENDER==1`
    ##   group   N     mean      se      sd    sdse     var   varse   mean_0
    ## 1   GR1 848 -0.99226 0.07126 0.81277 0.03932 0.66135 0.06357 -0.17761
    ## 2   GR2 820 -0.45673 0.04786 0.73724 0.03835 0.54426 0.05668  0.48279
    ## 3   GR3 849  0.14141 0.05126 0.78658 0.03239 0.61908 0.05123  1.05458
    ##   meandiff_0 meandiffse_0  tvalue_0 df_0 pvalue_0
    ## 1   -0.81465      0.09313  -8.74764 1665        0
    ## 2   -0.93952      0.07148 -13.14385 1664        0
    ## 3   -0.91317      0.08944 -10.20971 1665        0

## Difference between independent groups

For estimating the mean of independent groups, we can use the argument
[`repmeandif()`](https://dopatendo.github.io/ILSAstats/reference/repmeandif.md).
The only argument of this function will be an object produced by
`repmean` (using or not using `by`):

``` r
m1 <- repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR, aggregates = NULL,
        by = NULL)

m2 <- repmean(x = paste0("Math",1:5), PV = TRUE, setup = STGR, aggregates = NULL,
        by = "GENDER")
```

``` r
repmeandif(m1)
```

    ##   group1 group2      dif      se    tvalue   df pvalue
    ## 1    GR1    GR1  0.00000 0.05576   0.00000 3332      1
    ## 2    GR1    GR2 -0.61185 0.04614 -13.26073 3331      0
    ## 3    GR1    GR3 -1.18245 0.05402 -21.88912 3332      0
    ## 4    GR2    GR1  0.61185 0.04614  13.26073 3331      0
    ## 5    GR2    GR2  0.00000 0.03388   0.00000 3330      1
    ## 6    GR2    GR3 -0.57060 0.04401 -12.96524 3331      0
    ## 7    GR3    GR1  1.18245 0.05402  21.88912 3332      0
    ## 8    GR3    GR2  0.57060 0.04401  12.96524 3331      0
    ## 9    GR3    GR3  0.00000 0.05222   0.00000 3332      1

``` r
repmeandif(m2)
```

    ## $ALL
    ##   group1 group2      dif      se    tvalue   df pvalue
    ## 1    GR1    GR1  0.00000 0.05576   0.00000 3332      1
    ## 2    GR1    GR2 -0.61185 0.04614 -13.26073 3331      0
    ## 3    GR1    GR3 -1.18245 0.05402 -21.88912 3332      0
    ## 4    GR2    GR1  0.61185 0.04614  13.26073 3331      0
    ## 5    GR2    GR2  0.00000 0.03388   0.00000 3330      1
    ## 6    GR2    GR3 -0.57060 0.04401 -12.96524 3331      0
    ## 7    GR3    GR1  1.18245 0.05402  21.88912 3332      0
    ## 8    GR3    GR2  0.57060 0.04401  12.96524 3331      0
    ## 9    GR3    GR3  0.00000 0.05222   0.00000 3332      1
    ## 
    ## $`GENDER==0`
    ##   group1 group2      dif      se    tvalue   df pvalue
    ## 1    GR1    GR1  0.00000 0.06310   0.00000 1636      1
    ## 2    GR1    GR2 -0.66040 0.05679 -11.62881 1663      0
    ## 3    GR1    GR3 -1.23219 0.07608 -16.19598 1635      0
    ## 4    GR2    GR1  0.66040 0.05679  11.62881 1663      0
    ## 5    GR2    GR2  0.00000 0.04968   0.00000 1690      1
    ## 6    GR2    GR3 -0.57180 0.07093  -8.06147 1662      0
    ## 7    GR3    GR1  1.23219 0.07608  16.19598 1635      0
    ## 8    GR3    GR2  0.57180 0.07093   8.06147 1662      0
    ## 9    GR3    GR3  0.00000 0.08715   0.00000 1634      1
    ## 
    ## $`GENDER==1`
    ##   group1 group2      dif      se    tvalue   df pvalue
    ## 1    GR1    GR1  0.00000 0.10078   0.00000 1694      1
    ## 2    GR1    GR2 -0.53552 0.08584  -6.23858 1666      0
    ## 3    GR1    GR3 -1.13367 0.08778 -12.91490 1695      0
    ## 4    GR2    GR1  0.53552 0.08584   6.23858 1666      0
    ## 5    GR2    GR2  0.00000 0.06768   0.00000 1638      1
    ## 6    GR2    GR3 -0.59814 0.07013  -8.52902 1667      0
    ## 7    GR3    GR1  1.13367 0.08778  12.91490 1695      0
    ## 8    GR3    GR2  0.59814 0.07013   8.52902 1667      0
    ## 9    GR3    GR3  0.00000 0.07249   0.00000 1696      1

## Confidence intervals

Some ILSAs use confidence intervals instead of point estimates for
means. For adding confidence intervals to an object produced by
[`repmean()`](https://dopatendo.github.io/ILSAstats/reference/repmean.md),
we can use
[`repmeanCI()`](https://dopatendo.github.io/ILSAstats/reference/repmeanCI.md),
selecting a confidence level (by default it is 0.05):

``` r
repmeanCI(m1, alpha = 0.05)
```

    ##   group    N     mean      se   CIdown     CIup      sd    sdse     var   varse
    ## 1   GR1 1667 -0.59327 0.03943 -0.67056 -0.51599 0.90148 0.02093 0.81283 0.03766
    ## 2   GR2 1666  0.01857 0.02395 -0.02838  0.06552 0.89383 0.02075 0.79905 0.03708
    ## 3   GR3 1667  0.58917 0.03693  0.51680  0.66155 0.90989 0.02028 0.82799 0.03699

``` r
repmeanCI(m2, alpha = 0.05)
```

    ## $ALL
    ##   group    N     mean      se   CIdown     CIup      sd    sdse     var   varse
    ## 1   GR1 1667 -0.59327 0.03943 -0.67056 -0.51599 0.90148 0.02093 0.81283 0.03766
    ## 2   GR2 1666  0.01857 0.02395 -0.02838  0.06552 0.89383 0.02075 0.79905 0.03708
    ## 3   GR3 1667  0.58917 0.03693  0.51680  0.66155 0.90989 0.02028 0.82799 0.03699
    ## 
    ## $`GENDER==0`
    ##   group   N     mean      se   CIdown     CIup      sd    sdse     var   varse
    ## 1   GR1 819 -0.17761 0.04462 -0.26505 -0.09016 0.79348 0.03070 0.62996 0.04868
    ## 2   GR2 846  0.48279 0.03513  0.41394  0.55164 0.78112 0.03694 0.61083 0.05735
    ## 3   GR3 818  1.05458 0.06162  0.93380  1.17537 0.78595 0.03401 0.61826 0.05327
    ##     mean_1 meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
    ## 1 -0.99226    0.81465      0.09313  8.74764 1665        0
    ## 2 -0.45673    0.93952      0.07148 13.14385 1664        0
    ## 3  0.14141    0.91317      0.08944 10.20971 1665        0
    ## 
    ## $`GENDER==1`
    ##   group   N     mean      se   CIdown     CIup      sd    sdse     var   varse
    ## 1   GR1 848 -0.99226 0.07126 -1.13193 -0.85259 0.81277 0.03932 0.66135 0.06357
    ## 2   GR2 820 -0.45673 0.04786 -0.55053 -0.36294 0.73724 0.03835 0.54426 0.05668
    ## 3   GR3 849  0.14141 0.05126  0.04094  0.24187 0.78658 0.03239 0.61908 0.05123
    ##     mean_0 meandiff_0 meandiffse_0  tvalue_0 df_0 pvalue_0
    ## 1 -0.17761   -0.81465      0.09313  -8.74764 1665        0
    ## 2  0.48279   -0.93952      0.07148 -13.14385 1664        0
    ## 3  1.05458   -0.91317      0.08944 -10.20971 1665        0

We can also not add the confidence intervals, just obtain them
separately:

``` r
repmeanCI(m1, alpha = 0.05, add = FALSE)
```

    ##   group   CIdown     CIup
    ## 1   GR1 -0.67056 -0.51599
    ## 2   GR2 -0.02838  0.06552
    ## 3   GR3  0.51680  0.66155
