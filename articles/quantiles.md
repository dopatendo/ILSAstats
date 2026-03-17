# Quantiles

We can estimate the quantiles of any variable using the function
[`repquant()`](https://dopatendo.github.io/ILSAstats/reference/repquant.md),
as any other “rep” function of `ILSAstats`, we need to specify the data
(`df`), the total weights (`wt`), the replicate weights (`repwt`), and
the method (`method`).

Besides these basic options, other arguments can be used:

- `x`: a string with the name of the variable (or variables) to be used
  in the analysis. If multiple variables are specified, they will be
  treated as plausible values.
- `qtl`: a numeric vector indicating the desired quantiles (between 0
  and 1).
- `group`: a string containing the name of the variable that contains
  the groups of countries. If used all statistics will be estimated
  separately for each group.
- `by`: a string containing the name of a second grouping variable. If
  used, all statistics will be estimated separately for each category,
  and categories will be treated as **non-independent** from each other,
  e.g., boys and girls.
- `exclude`: a string containing which groups should be excluded from
  aggregate estimations.

## Weights and setup

For
[`repquant()`](https://dopatendo.github.io/ILSAstats/reference/repquant.md),
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

For example, if we want to estimate the quantiles of variable `"SES"`,
we can use either of the setups to get the overall or group results:

``` r
# No groups
repquant(x = "SES", setup = STNG, qtl = c(.25, .50, .75))
```

    ##   variable      P25   P25se     P50   P50se      P75   P75se
    ## 1      SES 49.28631 0.01865 49.9748 0.01744 50.65044 0.02438

``` r
# With groups
repquant(x = "SES", setup = STGR, qtl = c(.25, .50, .75))
```

    ##   variable     group      P25   P25se      P50   P50se      P75   P75se
    ## 1      SES    Pooled 49.28631 0.01865 49.97480 0.01744 50.65044 0.02438
    ## 2      SES Composite 49.33738 0.01965 49.96956 0.01021 50.61569 0.02126
    ## 3      SES       GR1 48.89882 0.04066 49.55360 0.01464 50.16223 0.03980
    ## 4      SES       GR2 49.30108 0.02581 49.94809 0.02105 50.62155 0.03729
    ## 5      SES       GR3 49.81225 0.03397 50.40699 0.01673 51.06328 0.03305

``` r
# With groups and exclusions
repquant(x = "SES", setup = STGE, qtl = c(.25, .50, .75))
```

    ##   variable     group      P25   P25se      P50   P50se      P75   P75se
    ## 1      SES    Pooled 49.27329 0.02817 50.00339 0.03331 50.66175 0.01781
    ## 2      SES Composite 49.35554 0.02649 49.98030 0.01111 50.61276 0.02587
    ## 3      SES       GR1 48.89882 0.04066 49.55360 0.01464 50.16223 0.03980
    ## 4      SES       GR2 49.30108 0.02581 49.94809 0.02105 50.62155 0.03729
    ## 5      SES       GR3 49.81225 0.03397 50.40699 0.01673 51.06328 0.03305

We can notice that using no groups we would get the same results for the
pooled estimates if we use groups and no exclusions. But, when we
exclude group 2, the pooled and the composite estimate changes.

## Plausible values

When treating with plausible values, we need to specify the names of all
plausible values of a construct, and use the argument `"PV"` so all
estimates will be combined (if not all variables will be estimated
separately). For example, for estimating the mean achievement in math
for this sample we would use:

``` r
# No groups
repquant(x = paste0("Math",1:5), setup = STNG, qtl = c(.25, .50, .75))
```

    ## More than one variable provided. 'x' treated as PVs.

    ##   variable      P25   P25se     P50   P50se     P75  P75se
    ## 1      PVs -0.68384 0.02715 0.00226 0.02279 0.69452 0.0256

``` r
# With groups
repquant(x = paste0("Math",1:5), setup = STGR, qtl = c(.25, .50, .75))
```

    ## More than one variable provided. 'x' treated as PVs.

    ##   variable     group      P25   P25se      P50   P50se     P75   P75se
    ## 1      PVs    Pooled -0.68384 0.02715  0.00226 0.02279 0.69452 0.02560
    ## 2      PVs Composite -0.58998 0.03166 -0.00377 0.02864 0.60439 0.02730
    ## 3      PVs       GR1 -1.18328 0.05265 -0.58728 0.05141 0.01332 0.05311
    ## 4      PVs       GR2 -0.57888 0.04231 -0.00614 0.04634 0.61342 0.03974
    ## 5      PVs       GR3 -0.00778 0.06679  0.58212 0.05092 1.18644 0.04805
