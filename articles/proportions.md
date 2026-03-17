# Proportions

We can estimate the proportions of the categories of any variable using
the function
[`repprop()`](https://dopatendo.github.io/ILSAstats/reference/repprop.md),
as any other “rep” function of `ILSAstats`, we need to specify the data
(`df`), the total weights (`wt`), the replicate weights (`repwt`), and
the method (`method`).

Besides these basic options, other arguments can be used:

- `x`: a string with the name of the variable (or variables) to be used
  in the analysis. If multiple variables are specified, they will be
  treated as plausible values.
- `group`: a string containing the name of the variable that contains
  the groups of countries. If used all statistics will be estimated
  separately for each group.
- `exclude`: a string containing which groups should be excluded from
  aggregate estimations.
- `aggregates`: a string containing the aggregate statistics that should
  be estimated. Options include: `"pooled"` for also estimating all
  groups (without exclusions) as a single group; and `"composite"` for
  averaging all the estimations for each single group (without
  exclusions).

## Weights and setup

For
[`repprop()`](https://dopatendo.github.io/ILSAstats/reference/repprop.md),
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

For example, if we want to estimate the proportions of each category of
variable `"GENDER"`, we can use either of the setups to get the overall
or group results:

``` r
# No groups
repprop(x = "GENDER", setup = STNG)
```

    ## $`GENDER==0`
    ##      N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
    ## 1 2483 0.49535 0.00773 0.50465   -0.00929      0.01546  -0.6012
    ## 
    ## $`GENDER==1`
    ##      N    prop      se  prop_0 propdiff_0 propdiffse_0 tvalue_0
    ## 1 2517 0.50465 0.00773 0.49535    0.00929      0.01546   0.6012

``` r
# With groups
repprop(x = "GENDER", setup = STGR)
```

    ## $`GENDER==0`
    ##       group    N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
    ## 1    Pooled 2483 0.49535 0.00773 0.50465   -0.00929      0.01546 -0.60120
    ## 2 Composite   NA 0.49533 0.00774 0.50467   -0.00933      0.01548 -0.60308
    ## 3       GR1  819 0.48976 0.01298 0.51024   -0.02048      0.02596 -0.78907
    ## 4       GR2  846 0.50590 0.01248 0.49410    0.01180      0.02496  0.47286
    ## 5       GR3  818 0.49034 0.01465 0.50966   -0.01932      0.02931 -0.65928
    ## 
    ## $`GENDER==1`
    ##       group    N    prop      se  prop_0 propdiff_0 propdiffse_0 tvalue_0
    ## 1    Pooled 2517 0.50465 0.00773 0.49535    0.00929      0.01546  0.60120
    ## 2 Composite   NA 0.50467 0.00774 0.49533    0.00933      0.01548  0.60308
    ## 3       GR1  848 0.51024 0.01298 0.48976    0.02048      0.02596  0.78907
    ## 4       GR2  820 0.49410 0.01248 0.50590   -0.01180      0.02496 -0.47286
    ## 5       GR3  849 0.50966 0.01465 0.49034    0.01932      0.02931  0.65928

``` r
# With groups and exclusions
repprop(x = "GENDER", setup = STGE)
```

    ## $`GENDER==0`
    ##       group    N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
    ## 1    Pooled 1637 0.49005 0.01071 0.50995   -0.01990      0.02143 -0.92882
    ## 2 Composite   NA 0.49005 0.00979 0.50995   -0.01990      0.01957 -1.01671
    ## 3       GR1  819 0.48976 0.01298 0.51024   -0.02048      0.02596 -0.78907
    ## 4       GR2  846 0.50590 0.01248 0.49410    0.01180      0.02496  0.47286
    ## 5       GR3  818 0.49034 0.01465 0.50966   -0.01932      0.02931 -0.65928
    ## 
    ## $`GENDER==1`
    ##       group    N    prop      se  prop_0 propdiff_0 propdiffse_0 tvalue_0
    ## 1    Pooled 1697 0.50995 0.01071 0.49005    0.01990      0.02143  0.92882
    ## 2 Composite   NA 0.50995 0.00979 0.49005    0.01990      0.01957  1.01671
    ## 3       GR1  848 0.51024 0.01298 0.48976    0.02048      0.02596  0.78907
    ## 4       GR2  820 0.49410 0.01248 0.50590   -0.01180      0.02496 -0.47286
    ## 5       GR3  849 0.50966 0.01465 0.49034    0.01932      0.02931  0.65928

We can notice that using no groups we would get the same results for the
pooled estimates if we use groups and no exclusions. But, when we
exclude group 2, the pooled and the composite estimate changes.

## Proportion table

By default,
[`repprop()`](https://dopatendo.github.io/ILSAstats/reference/repprop.md)
will provide a list where each element corresponds to the statistics of
each category. Nevertheless, we can summarize these results into a
single data frame using
[`repprop.table()`](https://dopatendo.github.io/ILSAstats/reference/repprop.table.md):

``` r
p1 <- repprop(x = "GENDER", setup = STGR)

repprop.table(x = p1)
```

    ##        group category    prop      se
    ## 1     Pooled        0 0.49535 0.00773
    ## 2  Composite        0 0.49533 0.00774
    ## 3        GR1        0 0.48976 0.01298
    ## 4        GR2        0 0.50590 0.01248
    ## 5        GR3        0 0.49034 0.01465
    ## 6     Pooled        1 0.50465 0.00773
    ## 7  Composite        1 0.50467 0.00774
    ## 8        GR1        1 0.51024 0.01298
    ## 9        GR2        1 0.49410 0.01248
    ## 10       GR3        1 0.50966 0.01465

By default,
[`repprop.table()`](https://dopatendo.github.io/ILSAstats/reference/repprop.table.md)
will provide a long table, but we can also get wide tables were groups
are separated by columns or by rows, and combined or not the proportion
estimates with the standard errors:

``` r
# Groups by rows, separate SE
repprop.table(x = p1, type = "wide1")
```

    ## $prop
    ##       group  prop.0  prop.1
    ## 1    Pooled 0.49535 0.50465
    ## 2 Composite 0.49533 0.50467
    ## 3       GR1 0.48976 0.51024
    ## 4       GR2 0.50590 0.49410
    ## 5       GR3 0.49034 0.50966
    ## 
    ## $se
    ##       group    se.0    se.1
    ## 1    Pooled 0.00773 0.00773
    ## 2 Composite 0.00774 0.00774
    ## 3       GR1 0.01298 0.01298
    ## 4       GR2 0.01248 0.01248
    ## 5       GR3 0.01465 0.01465

``` r
# Groups by rows, non-separate SE
repprop.table(x = p1, type = "wide1", separateSE = FALSE)
```

    ##       group  prop.0    se.0  prop.1    se.1
    ## 1    Pooled 0.49535 0.00773 0.50465 0.00773
    ## 2 Composite 0.49533 0.00774 0.50467 0.00774
    ## 3       GR1 0.48976 0.01298 0.51024 0.01298
    ## 4       GR2 0.50590 0.01248 0.49410 0.01248
    ## 5       GR3 0.49034 0.01465 0.50966 0.01465

``` r
# Groups by columns, separate SE
repprop.table(x = p1, type = "wide2")
```

    ## $prop
    ##   category    Pooled Composite       GR1       GR2       GR3
    ## 1        0 0.4953532 0.4953333 0.4897593 0.5059013 0.4903392
    ## 3        1 0.5046468 0.5046667 0.5102407 0.4940987 0.5096608
    ## 
    ## $se
    ##   category      Pooled   Composite        GR1        GR2        GR3
    ## 2        0 0.007729275 0.007738154 0.01297812 0.01248006 0.01465359
    ## 4        1 0.007729275 0.007738154 0.01297812 0.01248006 0.01465359

``` r
# Groups by columns, non-separate SE
repprop.table(x = p1, type = "wide2", separateSE = FALSE)
```

    ##   category statistic      Pooled   Composite        GR1        GR2        GR3
    ## 1        0      prop 0.495353153 0.495333268 0.48975930 0.50590131 0.49033920
    ## 2        0        se 0.007729275 0.007738154 0.01297812 0.01248006 0.01465359
    ## 3        1      prop 0.504646847 0.504666732 0.51024070 0.49409869 0.50966080
    ## 4        1        se 0.007729275 0.007738154 0.01297812 0.01248006 0.01465359

## Plausible values

When treating with plausible values, we need to specify the names of all
plausible values of a construct, so all estimates will be combined. For
example, for estimating the proportions of each proficiency level in
math we would use:

mean achievement in math for this sample we would use:

``` r
# No groups
repprop(x = paste0("CatMath",1:5), setup = STNG)|>
repprop.table(type = "wide2")
```

    ## More than one variable provided. 'x' treated as PVs.

    ## $prop
    ##   category    Pooled
    ## 1        1 0.1604335
    ## 3        2 0.3388836
    ## 5        3 0.3343512
    ## 7        4 0.1663317
    ## 
    ## $se
    ##   category      Pooled
    ## 2        1 0.008401597
    ## 4        2 0.009643329
    ## 6        3 0.010323985
    ## 8        4 0.006739155

``` r
# With groups
repprop(x = paste0("CatMath",1:5), setup = STGR)|>
repprop.table(type = "wide2")
```

    ## More than one variable provided. 'x' treated as PVs.

    ## $prop
    ##   category    Pooled Composite        GR1       GR2        GR3
    ## 1        1 0.1604335 0.1605950 0.31764327 0.1230013 0.04114025
    ## 3        2 0.3388836 0.3388807 0.42791017 0.3791503 0.20958161
    ## 5        3 0.3343512 0.3342311 0.21615471 0.3616510 0.42488755
    ## 7        4 0.1663317 0.1662933 0.03829185 0.1361973 0.32439059
    ## 
    ## $se
    ##   category      Pooled   Composite         GR1        GR2        GR3
    ## 2        1 0.008401597 0.008045957 0.020893302 0.01068229 0.00565645
    ## 4        2 0.009643329 0.009302842 0.016750950 0.01374063 0.01759223
    ## 6        3 0.010323985 0.010593748 0.015851628 0.02033506 0.01858114
    ## 8        4 0.006739155 0.008037791 0.006531979 0.01193100 0.01991078

## Aggregates

When using groups we can always omit the pooled and composite
calculations if we need to, by default both estimates will be
calculated.

``` r
# Default
repprop(x = paste0("CatMath",1:5), setup = STGR)|>
repprop.table(type = "wide2",separateSE = FALSE)
```

    ## More than one variable provided. 'x' treated as PVs.

    ##   category statistic      Pooled   Composite         GR1        GR2        GR3
    ## 1        1      prop 0.160433532 0.160594953 0.317643271 0.12300134 0.04114025
    ## 2        1        se 0.008401597 0.008045957 0.020893302 0.01068229 0.00565645
    ## 3        2      prop 0.338883595 0.338880700 0.427910174 0.37915032 0.20958161
    ## 4        2        se 0.009643329 0.009302842 0.016750950 0.01374063 0.01759223
    ## 5        3      prop 0.334351157 0.334231088 0.216154706 0.36165101 0.42488755
    ## 6        3        se 0.010323985 0.010593748 0.015851628 0.02033506 0.01858114
    ## 7        4      prop 0.166331717 0.166293259 0.038291849 0.13619733 0.32439059
    ## 8        4        se 0.006739155 0.008037791 0.006531979 0.01193100 0.01991078

``` r
# Only pooled
repprop(x = paste0("CatMath",1:5), setup = STGR, aggregates = "pooled")|>
repprop.table(type = "wide2",separateSE = FALSE)
```

    ## More than one variable provided. 'x' treated as PVs.

    ##   category statistic      Pooled         GR1        GR2        GR3
    ## 1        1      prop 0.160433532 0.317643271 0.12300134 0.04114025
    ## 2        1        se 0.008401597 0.020893302 0.01068229 0.00565645
    ## 3        2      prop 0.338883595 0.427910174 0.37915032 0.20958161
    ## 4        2        se 0.009643329 0.016750950 0.01374063 0.01759223
    ## 5        3      prop 0.334351157 0.216154706 0.36165101 0.42488755
    ## 6        3        se 0.010323985 0.015851628 0.02033506 0.01858114
    ## 7        4      prop 0.166331717 0.038291849 0.13619733 0.32439059
    ## 8        4        se 0.006739155 0.006531979 0.01193100 0.01991078

``` r
# Only composite
repprop(x = paste0("CatMath",1:5), setup = STGR, aggregates = "composite")|>
repprop.table(type = "wide2",separateSE = FALSE)
```

    ## More than one variable provided. 'x' treated as PVs.

    ##   category statistic   Composite         GR1        GR2        GR3
    ## 1        1      prop 0.160594953 0.317643271 0.12300134 0.04114025
    ## 2        1        se 0.008045957 0.020893302 0.01068229 0.00565645
    ## 3        2      prop 0.338880700 0.427910174 0.37915032 0.20958161
    ## 4        2        se 0.009302842 0.016750950 0.01374063 0.01759223
    ## 5        3      prop 0.334231088 0.216154706 0.36165101 0.42488755
    ## 6        3        se 0.010593748 0.015851628 0.02033506 0.01858114
    ## 7        4      prop 0.166293259 0.038291849 0.13619733 0.32439059
    ## 8        4        se 0.008037791 0.006531979 0.01193100 0.01991078

``` r
# No aggregates
repprop(x = paste0("CatMath",1:5), setup = STGR, aggregates = NULL)|>
repprop.table(type = "wide2",separateSE = FALSE)
```

    ## More than one variable provided. 'x' treated as PVs.

    ##   category statistic         GR1        GR2        GR3
    ## 1        1      prop 0.317643271 0.12300134 0.04114025
    ## 2        1        se 0.020893302 0.01068229 0.00565645
    ## 3        2      prop 0.427910174 0.37915032 0.20958161
    ## 4        2        se 0.016750950 0.01374063 0.01759223
    ## 5        3      prop 0.216154706 0.36165101 0.42488755
    ## 6        3        se 0.015851628 0.02033506 0.01858114
    ## 7        4      prop 0.038291849 0.13619733 0.32439059
    ## 8        4        se 0.006531979 0.01193100 0.01991078
