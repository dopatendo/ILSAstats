# Mean, Variance and Standard Deviation with Replicate Weights

Estimates the mean, variance and standard deviation with replicate
weights for a variable or a group of variables and for one or more
populations. For a detailed explanation on how the standard errors are
estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
repmean(
  x,
  PV = FALSE,
  setup = NULL,
  repwt = NULL,
  repindex = NULL,
  wt,
  df,
  method,
  var = c("unbiased", "ML", "none"),
  group = NULL,
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

- repwt:

  a string indicating the common names for the replicate weights columns
  (within `df`), or a data frame with the replicate weights.

- repindex:

  a `repweights.index` object generated with
  [`repcreate`](https://dopatendo.github.io/ILSAstats/reference/repcreate.md)`(..., index = TRUE)`.
  Using this argument instead of `repwt` will speed up the estimations
  considerably.

- wt:

  a string specifying the name of the column (within `df`) with the
  total weights.

- df:

  a data frame.

- method:

  a string indicating the name of the replication method. Available
  options are: `"JK2-full"`, `"JK2-half"`, `"FAY-0.5"`, and
  `"JK2-half-1PV"`.  
    
  Additionally, ILSA names can be used, defaulting into:

  - `"TIMSS"`, `"PIRLS"`, or `"LANA"` for `"JK2-full"`;

  - `"ICILS"`, `"ICCS"`, or `"CIVED"` for `"JK2-half"`;

  - `"PISA"` or `"TALIS"` for `"FAY-0.5"`;

  - and `"oldTIMSS"`, `"oldPIRLS"`, or `"RLII"` for `"JK2-half-1PV"`.

  Note that `"oldTIMSS"` and `"oldPIRLS"` refer to the method used for
  TIMSS and PIRLS before 2015, where within imputation variance is
  estimated using only 1 plausible value.

- var:

  a string indicating the method to use for the variance: `"unbiased"`
  calculates the unbiased estimate (n-1); `"ML"` calculates the maximum
  likelihood estimate.

- group:

  a string specifying the variable name (within `df`) to be used for
  grouping. Categories in `group` are treated as independent, e.g.,
  countries.

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
  and no aggregates will be estimated. Default is `TRUE`. This argument
  will be ignored if `by` is used.

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
        method = "ICILS",var = "ML")
#> 50 weights found.
#>      N    mean      se      sd    sdse     var   varse
#> 1 4483 3.61672 0.00951 0.66582 0.01105 0.44331 0.01472

# One variable - weights as a separate data frame
repmean(x = c("item01"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML")
#>      N    mean      se      sd    sdse     var   varse
#> 1 4483 3.61672 0.00951 0.66582 0.01105 0.44331 0.01472

# Multiple variables
repmean(x = c("item01","item02","item03"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML")
#>   variable    N    mean      se      sd    sdse     var   varse
#> 1   item01 4483 3.61672 0.00951 0.66582 0.01105 0.44331 0.01472
#> 2   item02 4501 3.22749 0.01325 0.81229 0.00729 0.65981 0.01184
#> 3   item03 4519 3.82527 0.00797 0.50522 0.01426 0.25525 0.01443

# One PV variable
repmean(x = paste0("Math",1:5),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML")
#>      N    mean      se      sd    sdse     var   varse
#> 1 5000 0.00524 0.01634 1.02282 0.01099 1.04619 0.02247

### Groups ----

# One variable
repmean(x = c("item01"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#>       group    N    mean      se      sd    sdse     var   varse
#> 1    Pooled 2992 3.62342 0.01165 0.65023 0.01427 0.42279 0.01855
#> 2 Composite   NA 3.62334 0.01168 0.65009 0.01427 0.42269 0.01851
#> 3       GR1 1507 3.63936 0.01542 0.64117 0.02043 0.41110 0.02611
#> 4       GR2 1491 3.60337 0.01979 0.69582 0.02070 0.48416 0.02881
#> 5       GR3 1485 3.60732 0.01755 0.65900 0.01994 0.43429 0.02625

# Multiple variables
repmean(x = c("item01","item02","item03"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#>    variable     group    N    mean      se      sd    sdse     var   varse
#> 1    item01    Pooled 2992 3.62342 0.01165 0.65023 0.01427 0.42279 0.01855
#> 2    item01 Composite   NA 3.62334 0.01168 0.65009 0.01427 0.42269 0.01851
#> 3    item01       GR1 1507 3.63936 0.01542 0.64117 0.02043 0.41110 0.02611
#> 4    item01       GR2 1491 3.60337 0.01979 0.69582 0.02070 0.48416 0.02881
#> 5    item01       GR3 1485 3.60732 0.01755 0.65900 0.01994 0.43429 0.02625
#> 6    item02    Pooled 3007 3.23185 0.01599 0.80741 0.00995 0.65191 0.01606
#> 7    item02 Composite   NA 3.23198 0.01630 0.80709 0.00984 0.65141 0.01587
#> 8    item02       GR1 1496 3.25676 0.02233 0.80373 0.01355 0.64598 0.02177
#> 9    item02       GR2 1494 3.21879 0.02134 0.82206 0.01485 0.67579 0.02444
#> 10   item02       GR3 1511 3.20721 0.02374 0.81046 0.01427 0.65685 0.02310
#> 11   item03    Pooled 3022 3.82110 0.01033 0.51406 0.01834 0.26426 0.01890
#> 12   item03 Composite   NA 3.82110 0.00932 0.51404 0.01776 0.26432 0.01831
#> 13   item03       GR1 1513 3.82045 0.01437 0.52313 0.02572 0.27367 0.02698
#> 14   item03       GR2 1497 3.83362 0.01383 0.48703 0.02639 0.23720 0.02567
#> 15   item03       GR3 1509 3.82176 0.01187 0.50494 0.02450 0.25497 0.02476

# One PV variable
repmean(x = paste0("Math",1:5),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#>       group    N     mean      se      sd    sdse     var   varse
#> 1    Pooled 3334 -0.00146 0.02130 1.08185 0.01624 1.17050 0.03516
#> 2 Composite   NA -0.00205 0.02701 0.90569 0.01462 0.82041 0.02648
#> 3       GR1 1667 -0.59327 0.03943 0.90148 0.02101 0.81283 0.03779
#> 4       GR2 1666  0.01857 0.02397 0.89383 0.02093 0.79905 0.03740
#> 5       GR3 1667  0.58917 0.03692 0.90989 0.02034 0.82799 0.03711

### Groups and By ----

# One variable
repmean(x = c("item01"),
        PV = FALSE,
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS",var = "ML",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> $ALL
#>       group    N    mean      se      sd    sdse     var   varse
#> 1    Pooled 2992 3.62342 0.01165 0.65023 0.01427 0.42279 0.01855
#> 2 Composite   NA 3.62334 0.01168 0.65009 0.01427 0.42269 0.01851
#> 3       GR1 1507 3.63936 0.01542 0.64117 0.02043 0.41110 0.02611
#> 4       GR2 1491 3.60337 0.01979 0.69582 0.02070 0.48416 0.02881
#> 5       GR3 1485 3.60732 0.01755 0.65900 0.01994 0.43429 0.02625
#> 
#> $`GENDER==0`
#>       group    N    mean      se      sd    sdse     var   varse  mean_1
#> 1    Pooled 1458 3.63697 0.01813 0.63377 0.02041 0.40166 0.02583 3.61060
#> 2 Composite   NA 3.63698 0.01755 0.63385 0.02006 0.40179 0.02533 3.61028
#> 3       GR1  734 3.63042 0.02467 0.63847 0.02812 0.40764 0.03585 3.64776
#> 4       GR2  760 3.59046 0.02919 0.72062 0.02742 0.51929 0.03963 3.61665
#> 5       GR3  724 3.64354 0.02497 0.62923 0.02860 0.39594 0.03581 3.57281
#>   meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1    0.02637      0.02540  1.03818 2990  0.29927
#> 2    0.02670      0.02492  1.07127   NA       NA
#> 3   -0.01734      0.03674 -0.47194 1505  0.63704
#> 4   -0.02618      0.04066 -0.64399 1489  0.51968
#> 5    0.07074      0.03369  2.09948 1483  0.03594
#> 
#> $`GENDER==1`
#>       group    N    mean      se      sd    sdse     var   varse  mean_0
#> 1    Pooled 1534 3.61060 0.01636 0.66531 0.01903 0.44264 0.02532 3.63697
#> 2 Composite   NA 3.61028 0.01663 0.66426 0.01914 0.44166 0.02521 3.63698
#> 3       GR1  773 3.64776 0.02332 0.64386 0.02989 0.41456 0.03837 3.63042
#> 4       GR2  731 3.61665 0.02746 0.66940 0.02560 0.44809 0.03425 3.59046
#> 5       GR3  761 3.57281 0.02372 0.68466 0.02391 0.46875 0.03272 3.64354
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
        method = "ICILS",var = "ML",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> $ALL
#>       group    N     mean      se      sd    sdse     var   varse
#> 1    Pooled 3334 -0.00146 0.02130 1.08185 0.01624 1.17050 0.03516
#> 2 Composite   NA -0.00205 0.02701 0.90569 0.01462 0.82041 0.02648
#> 3       GR1 1667 -0.59327 0.03943 0.90148 0.02101 0.81283 0.03779
#> 4       GR2 1666  0.01857 0.02397 0.89383 0.02093 0.79905 0.03740
#> 5       GR3 1667  0.58917 0.03692 0.90989 0.02034 0.82799 0.03711
#> 
#> $`GENDER==0`
#>       group    N     mean      se      sd    sdse     var   varse   mean_1
#> 1    Pooled 1637  0.43947 0.04440 1.00205 0.01584 1.00412 0.03173 -0.42518
#> 2 Composite   NA  0.43849 0.03804 0.78972 0.02293 0.62411 0.03611 -0.42542
#> 3       GR1  819 -0.17761 0.04462 0.79348 0.03057 0.62996 0.04847 -0.99226
#> 4       GR2  846  0.48279 0.03513 0.78112 0.03688 0.61083 0.05727 -0.45673
#> 5       GR3  818  1.05458 0.06162 0.78595 0.03419 0.61826 0.05353  0.14141
#>   meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1    0.86465      0.08929  9.68332 3332        0
#> 2    0.86391      0.06456 13.38179   NA       NA
#> 3    0.81465      0.09314  8.74621 1665        0
#> 4    0.93952      0.07151 13.13918 1664        0
#> 5    0.91317      0.08942 10.21235 1665        0
#> 
#> $`GENDER==1`
#>       group    N     mean      se      sd    sdse     var   varse   mean_0
#> 1    Pooled 1697 -0.42518 0.05160 0.98078 0.02238 0.96213 0.04410  0.43947
#> 2 Composite   NA -0.42542 0.04389 0.79967 0.02551 0.64022 0.04090  0.43849
#> 3       GR1  848 -0.99226 0.07127 0.81277 0.03952 0.66135 0.06388 -0.17761
#> 4       GR2  820 -0.45673 0.04790 0.73724 0.03839 0.54426 0.05673  0.48279
#> 5       GR3  849  0.14141 0.05124 0.78658 0.03226 0.61908 0.05110  1.05458
#>   meandiff_0 meandiffse_0  tvalue_0 df_0 pvalue_0
#> 1   -0.86465      0.08929  -9.68332 3332        0
#> 2   -0.86391      0.06456 -13.38179   NA       NA
#> 3   -0.81465      0.09314  -8.74621 1665        0
#> 4   -0.93952      0.07151 -13.13918 1664        0
#> 5   -0.91317      0.08942 -10.21235 1665        0
#> 
```
