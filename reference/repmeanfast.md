# Mean with Replicate Weights

Estimates only the mean with replicate weights for a variable or a group
of variables and for one or more populations. For a detailed explanation
on how the standard errors are estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
repmeanfast(
  x,
  PV = FALSE,
  setup = NULL,
  repindex,
  wt,
  df,
  groups = NULL,
  by = NULL,
  exclude = NULL,
  aggregates = c("pooled", "composite"),
  simplify = TRUE
)
```

## Arguments

- x:

  a string vector specifying variable names (within `df`) for analysis.

- PV:

  a logical value indicating if the variables in `x` are plausible
  values.

- setup:

  an optional list produced by
  [`repsetup`](https://dopatendo.github.io/ILSAstats/reference/repsetup.md).

- repindex:

  a `repweights.index` object generate with
  [`repcreate`](https://dopatendo.github.io/ILSAstats/reference/repcreate.md)`(..., index = TRUE)`.

- wt:

  a string specifying the name of the column (within `df`) with the
  total weights.

- df:

  a data frame.

- by:

  a string specifying a second variable (within `df`) for grouping.
  Categories used in `by` are not considered independent, e.g., gender
  within a country. If used, the output will be a list with the same
  length as the unique values of `by`. This can only be used for
  analyses with one variable or a group of PVs.

- exclude:

  a vector indicating which groups (in the same format as `group`)
  should be excluded from the pooled and composite estimates.

- aggregates:

  a string vector indicating which aggregates should be included,
  options are `"pooled"` and `"composite"`, both options can be used at
  the same time. If `NULL` no aggregate will be estimated.

- simplify:

  a logical value indicating if only the summary statistics should be
  printed. If `FALSE` estimations for all replicated will be provided
  and no aggregates will be estimated. Default is `TRUE`.

## Value

a data frame or a list.

## Examples

``` r
# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----

# One variable - weights within df
repmean(x = c("item01"),
        PV = FALSE,
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
        method = "ICILS",var = "ML",zones = "jkzones")
#> 50 weights found.
#>      N nzones    mean      se      sd    sdse     var   varse
#> 1 4483     50 3.61672 0.00951 0.66577 0.01106 0.44325 0.01474

# One variable - weights as a separate data frame
repmean(x = c("item01"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones")
#>      N nzones    mean      se      sd    sdse     var   varse
#> 1 4483     50 3.61672 0.00951 0.66577 0.01106 0.44325 0.01474

# Multiple variables
repmean(x = c("item01","item02","item03"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones")
#>   variable    N nzones    mean      se      sd    sdse     var   varse
#> 1   item01 4483     50 3.61672 0.00951 0.66577 0.01106 0.44325 0.01474
#> 2   item02 4501     50 3.22749 0.01325 0.81223 0.00730 0.65972 0.01186
#> 3   item03 4519     50 3.82527 0.00797 0.50518 0.01429 0.25521 0.01446

# One PV variable
repmean(x = paste0("Math",1:5),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones")
#>      N nzones    mean      se      sd    sdse     var   varse
#> 1 5000     50 0.00524 0.01634 1.02276 0.01099 1.04605 0.02247

### Groups ----

# One variable
repmean(x = c("item01"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#>       group    N nzones    mean      se      sd    sdse     var   varse
#> 1    Pooled 2992     50 3.62342 0.01165 0.65015 0.01420 0.42270 0.01845
#> 2 Composite   NA     NA 3.62334 0.01168 0.64994 0.01419 0.42250 0.01840
#> 3       GR1 1507     50 3.63936 0.01542 0.64103 0.02036 0.41092 0.02603
#> 4       GR2 1491     50 3.60337 0.01979 0.69566 0.02093 0.48394 0.02912
#> 5       GR3 1485     50 3.60732 0.01755 0.65886 0.01977 0.43409 0.02602

# Multiple variables
repmean(x = c("item01","item02","item03"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#>    variable     group    N nzones    mean      se      sd    sdse     var
#> 1    item01    Pooled 2992     50 3.62342 0.01165 0.65015 0.01420 0.42270
#> 2    item01 Composite   NA     NA 3.62334 0.01168 0.64994 0.01419 0.42250
#> 3    item01       GR1 1507     50 3.63936 0.01542 0.64103 0.02036 0.41092
#> 4    item01       GR2 1491     50 3.60337 0.01979 0.69566 0.02093 0.48394
#> 5    item01       GR3 1485     50 3.60732 0.01755 0.65886 0.01977 0.43409
#> 6    item02    Pooled 3007     50 3.23185 0.01599 0.80732 0.00993 0.65177
#> 7    item02 Composite   NA     NA 3.23198 0.01630 0.80691 0.00980 0.65112
#> 8    item02       GR1 1496     50 3.25676 0.02233 0.80355 0.01363 0.64569
#> 9    item02       GR2 1494     50 3.21879 0.02134 0.82188 0.01488 0.67549
#> 10   item02       GR3 1511     50 3.20721 0.02374 0.81028 0.01409 0.65656
#> 11   item03    Pooled 3022     50 3.82110 0.01033 0.51400 0.01838 0.26420
#> 12   item03 Composite   NA     NA 3.82110 0.00932 0.51393 0.01780 0.26420
#> 13   item03       GR1 1513     50 3.82045 0.01437 0.52302 0.02578 0.27355
#> 14   item03       GR2 1497     50 3.83362 0.01383 0.48693 0.02640 0.23710
#> 15   item03       GR3 1509     50 3.82176 0.01187 0.50483 0.02454 0.25486
#>      varse
#> 1  0.01845
#> 2  0.01840
#> 3  0.02603
#> 4  0.02912
#> 5  0.02602
#> 6  0.01604
#> 7  0.01581
#> 8  0.02190
#> 9  0.02450
#> 10 0.02281
#> 11 0.01894
#> 12 0.01835
#> 13 0.02705
#> 14 0.02569
#> 15 0.02480

# One PV variable
repmean(x = paste0("Math",1:5),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#>       group    N nzones     mean      se      sd    sdse     var   varse
#> 1    Pooled 3334     50 -0.00146 0.02130 1.08175 0.01625 1.17026 0.03519
#> 2 Composite   NA     NA -0.00205 0.02701 0.90551 0.01456 0.82008 0.02638
#> 3       GR1 1667     50 -0.59327 0.03943 0.90130 0.02092 0.81251 0.03764
#> 4       GR2 1666     50  0.01857 0.02397 0.89365 0.02078 0.79873 0.03712
#> 5       GR3 1667     50  0.58917 0.03692 0.90971 0.02026 0.82766 0.03697

### Groups and By ----

# One variable
repmean(x = c("item01"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> $ALL
#>       group    N nzones    mean      se      sd    sdse     var   varse
#> 1    Pooled 2992     50 3.62342 0.01165 0.65015 0.01420 0.42270 0.01845
#> 2 Composite   NA     NA 3.62334 0.01168 0.64994 0.01419 0.42250 0.01840
#> 3       GR1 1507     50 3.63936 0.01542 0.64103 0.02036 0.41092 0.02603
#> 4       GR2 1491     50 3.60337 0.01979 0.69566 0.02093 0.48394 0.02912
#> 5       GR3 1485     50 3.60732 0.01755 0.65886 0.01977 0.43409 0.02602
#> 
#> $`GENDER==0`
#>       group    N nzones    mean      se      sd    sdse     var   varse  mean_1
#> 1    Pooled 1458     50 3.63697 0.01813 0.63362 0.02024 0.40148 0.02560 3.61060
#> 2 Composite   NA     NA 3.63698 0.01755 0.63356 0.01984 0.40142 0.02507 3.61028
#> 3       GR1  734     50 3.63042 0.02467 0.63818 0.02781 0.40727 0.03545 3.64776
#> 4       GR2  760     50 3.59046 0.02919 0.72030 0.02785 0.51883 0.04025 3.61665
#> 5       GR3  724     50 3.64354 0.02497 0.62894 0.02832 0.39557 0.03545 3.57281
#>   meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1    0.02637      0.02540  1.03818 2990  0.29927
#> 2    0.02670      0.02492  1.07127   NA       NA
#> 3   -0.01734      0.03674 -0.47194 1505  0.63704
#> 4   -0.02618      0.04066 -0.64399 1489  0.51968
#> 5    0.07074      0.03369  2.09948 1483  0.03594
#> 
#> $`GENDER==1`
#>       group    N nzones    mean      se      sd    sdse     var   varse  mean_0
#> 1    Pooled 1534     50 3.61060 0.01636 0.66517 0.01895 0.44245 0.02521 3.63697
#> 2 Composite   NA     NA 3.61028 0.01663 0.66397 0.01904 0.44127 0.02507 3.63698
#> 3       GR1  773     50 3.64776 0.02332 0.64359 0.02991 0.41420 0.03840 3.63042
#> 4       GR2  731     50 3.61665 0.02746 0.66909 0.02575 0.44769 0.03444 3.59046
#> 5       GR3  761     50 3.57281 0.02372 0.68436 0.02356 0.46834 0.03224 3.64354
#>   meandiff_0 meandiffse_0 tvalue_0 df_0 pvalue_0
#> 1   -0.02637      0.02540 -1.03818 2990  0.29927
#> 2   -0.02670      0.02492 -1.07127   NA       NA
#> 3    0.01734      0.03674  0.47194 1505  0.63704
#> 4    0.02618      0.04066  0.64399 1489  0.51968
#> 5   -0.07074      0.03369 -2.09948 1483  0.03594
#> 

# One PV variable
repmean(x = paste0("Math",1:5),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",zones = "jkzones",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> $ALL
#>       group    N nzones     mean      se      sd    sdse     var   varse
#> 1    Pooled 3334     50 -0.00146 0.02130 1.08175 0.01625 1.17026 0.03519
#> 2 Composite   NA     NA -0.00205 0.02701 0.90551 0.01456 0.82008 0.02638
#> 3       GR1 1667     50 -0.59327 0.03943 0.90130 0.02092 0.81251 0.03764
#> 4       GR2 1666     50  0.01857 0.02397 0.89365 0.02078 0.79873 0.03712
#> 5       GR3 1667     50  0.58917 0.03692 0.90971 0.02026 0.82766 0.03697
#> 
#> $`GENDER==0`
#>       group    N nzones     mean      se      sd    sdse     var   varse
#> 1    Pooled 1637     50  0.43947 0.04440 1.00184 0.01593 1.00371 0.03191
#> 2 Composite   NA     NA  0.43849 0.03804 0.78939 0.02289 0.62360 0.03604
#> 3       GR1  819     50 -0.17761 0.04462 0.79316 0.03068 0.62944 0.04864
#> 4       GR2  846     50  0.48279 0.03513 0.78081 0.03692 0.61035 0.05731
#> 5       GR3  818     50  1.05458 0.06162 0.78563 0.03399 0.61776 0.05321
#>     mean_1 meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1 -0.42518    0.86465      0.08929  9.68332 3332        0
#> 2 -0.42542    0.86391      0.06456 13.38179   NA       NA
#> 3 -0.99226    0.81465      0.09314  8.74621 1665        0
#> 4 -0.45673    0.93952      0.07151 13.13918 1664        0
#> 5  0.14141    0.91317      0.08942 10.21235 1665        0
#> 
#> $`GENDER==1`
#>       group    N nzones     mean      se      sd    sdse     var   varse
#> 1    Pooled 1697     50 -0.42518 0.05160 0.98059 0.02236 0.96175 0.04405
#> 2 Composite   NA     NA -0.42542 0.04389 0.79936 0.02544 0.63971 0.04078
#> 3       GR1  848     50 -0.99226 0.07127 0.81245 0.03932 0.66083 0.06353
#> 4       GR2  820     50 -0.45673 0.04790 0.73694 0.03836 0.54382 0.05667
#> 5       GR3  849     50  0.14141 0.05124 0.78627 0.03229 0.61859 0.05114
#>     mean_0 meandiff_0 meandiffse_0  tvalue_0 df_0 pvalue_0
#> 1  0.43947   -0.86465      0.08929  -9.68332 3332        0
#> 2  0.43849   -0.86391      0.06456 -13.38179   NA       NA
#> 3 -0.17761   -0.81465      0.09314  -8.74621 1665        0
#> 4  0.48279   -0.93952      0.07151 -13.13918 1664        0
#> 5  1.05458   -0.91317      0.08942 -10.21235 1665        0
#> 
```
