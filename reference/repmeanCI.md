# Confidence Intervals for Replicated Means

Calculates the confidence intervals for a
[`repmean`](https://dopatendo.github.io/ILSAstats/reference/repmean.md)
object.

## Usage

``` r
repmeanCI(x, alpha = 0.05, add = TRUE)
```

## Arguments

- x:

  an object produced by
  [`repmean`](https://dopatendo.github.io/ILSAstats/reference/repmean.md).

- alpha:

  a numeric value indicating confidence level.

- add:

  a logical value indicating if the confidence intervals should be added
  to the object or not. Defaults is `TRUE`.

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


### Groups ----

# One variable
reme <- repmean(x = c("item01"),
                PV = FALSE,
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#>       group    N    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 2992 3.62342 0.01165 3.60058 3.64626 0.65023 0.01427 0.42279
#> 2 Composite   NA 3.62334 0.01168 3.60044 3.64623 0.65009 0.01427 0.42269
#> 3       GR1 1507 3.63936 0.01542 3.60914 3.66957 0.64117 0.02043 0.41110
#> 4       GR2 1491 3.60337 0.01979 3.56457 3.64216 0.69582 0.02070 0.48416
#> 5       GR3 1485 3.60732 0.01755 3.57292 3.64171 0.65900 0.01994 0.43429
#>     varse
#> 1 0.01855
#> 2 0.01851
#> 3 0.02611
#> 4 0.02881
#> 5 0.02625


# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#>       group    N     mean      se   CIdown     CIup      sd    sdse     var
#> 1    Pooled 3334 -0.00146 0.02130 -0.04321  0.04028 1.08185 0.01624 1.17050
#> 2 Composite   NA -0.00205 0.02701 -0.05498  0.05088 0.90569 0.01462 0.82041
#> 3       GR1 1667 -0.59327 0.03943 -0.67056 -0.51599 0.90148 0.02101 0.81283
#> 4       GR2 1666  0.01857 0.02397 -0.02841  0.06555 0.89383 0.02093 0.79905
#> 5       GR3 1667  0.58917 0.03692  0.51682  0.66153 0.90989 0.02034 0.82799
#>     varse
#> 1 0.03516
#> 2 0.02648
#> 3 0.03779
#> 4 0.03740
#> 5 0.03711

### Groups and By ----

# One variable
reme <- repmean(x = c("item01"),
                PV = FALSE,
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#> $ALL
#>       group    N    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 2992 3.62342 0.01165 3.60058 3.64626 0.65023 0.01427 0.42279
#> 2 Composite   NA 3.62334 0.01168 3.60044 3.64623 0.65009 0.01427 0.42269
#> 3       GR1 1507 3.63936 0.01542 3.60914 3.66957 0.64117 0.02043 0.41110
#> 4       GR2 1491 3.60337 0.01979 3.56457 3.64216 0.69582 0.02070 0.48416
#> 5       GR3 1485 3.60732 0.01755 3.57292 3.64171 0.65900 0.01994 0.43429
#>     varse
#> 1 0.01855
#> 2 0.01851
#> 3 0.02611
#> 4 0.02881
#> 5 0.02625
#> 
#> $`GENDER==0`
#>       group    N    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 1458 3.63697 0.01813 3.60143 3.67251 0.63377 0.02041 0.40166
#> 2 Composite   NA 3.63698 0.01755 3.60258 3.67138 0.63385 0.02006 0.40179
#> 3       GR1  734 3.63042 0.02467 3.58206 3.67878 0.63847 0.02812 0.40764
#> 4       GR2  760 3.59046 0.02919 3.53324 3.64768 0.72062 0.02742 0.51929
#> 5       GR3  724 3.64354 0.02497 3.59460 3.69249 0.62923 0.02860 0.39594
#>     varse  mean_1 meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1 0.02583 3.61060    0.02637      0.02540  1.03818 2990  0.29927
#> 2 0.02533 3.61028    0.02670      0.02492  1.07127   NA       NA
#> 3 0.03585 3.64776   -0.01734      0.03674 -0.47194 1505  0.63704
#> 4 0.03963 3.61665   -0.02618      0.04066 -0.64399 1489  0.51968
#> 5 0.03581 3.57281    0.07074      0.03369  2.09948 1483  0.03594
#> 
#> $`GENDER==1`
#>       group    N    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 1534 3.61060 0.01636 3.57853 3.64267 0.66531 0.01903 0.44264
#> 2 Composite   NA 3.61028 0.01663 3.57768 3.64288 0.66426 0.01914 0.44166
#> 3       GR1  773 3.64776 0.02332 3.60205 3.69347 0.64386 0.02989 0.41456
#> 4       GR2  731 3.61665 0.02746 3.56283 3.67046 0.66940 0.02560 0.44809
#> 5       GR3  761 3.57281 0.02372 3.52632 3.61930 0.68466 0.02391 0.46875
#>     varse  mean_0 meandiff_0 meandiffse_0 tvalue_0 df_0 pvalue_0
#> 1 0.02532 3.63697   -0.02637      0.02540 -1.03818 2990  0.29927
#> 2 0.02521 3.63698   -0.02670      0.02492 -1.07127   NA       NA
#> 3 0.03837 3.63042    0.01734      0.03674  0.47194 1505  0.63704
#> 4 0.03425 3.59046    0.02618      0.04066  0.64399 1489  0.51968
#> 5 0.03272 3.64354   -0.07074      0.03369 -2.09948 1483  0.03594
#> 

# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#> $ALL
#>       group    N     mean      se   CIdown     CIup      sd    sdse     var
#> 1    Pooled 3334 -0.00146 0.02130 -0.04321  0.04028 1.08185 0.01624 1.17050
#> 2 Composite   NA -0.00205 0.02701 -0.05498  0.05088 0.90569 0.01462 0.82041
#> 3       GR1 1667 -0.59327 0.03943 -0.67056 -0.51599 0.90148 0.02101 0.81283
#> 4       GR2 1666  0.01857 0.02397 -0.02841  0.06555 0.89383 0.02093 0.79905
#> 5       GR3 1667  0.58917 0.03692  0.51682  0.66153 0.90989 0.02034 0.82799
#>     varse
#> 1 0.03516
#> 2 0.02648
#> 3 0.03779
#> 4 0.03740
#> 5 0.03711
#> 
#> $`GENDER==0`
#>       group    N     mean      se   CIdown     CIup      sd    sdse     var
#> 1    Pooled 1637  0.43947 0.04440  0.35244  0.52649 1.00205 0.01584 1.00412
#> 2 Composite   NA  0.43849 0.03804  0.36393  0.51305 0.78972 0.02293 0.62411
#> 3       GR1  819 -0.17761 0.04462 -0.26507 -0.09015 0.79348 0.03057 0.62996
#> 4       GR2  846  0.48279 0.03513  0.41394  0.55164 0.78112 0.03688 0.61083
#> 5       GR3  818  1.05458 0.06162  0.93380  1.17537 0.78595 0.03419 0.61826
#>     varse   mean_1 meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1 0.03173 -0.42518    0.86465      0.08929  9.68332 3332        0
#> 2 0.03611 -0.42542    0.86391      0.06456 13.38179   NA       NA
#> 3 0.04847 -0.99226    0.81465      0.09314  8.74621 1665        0
#> 4 0.05727 -0.45673    0.93952      0.07151 13.13918 1664        0
#> 5 0.05353  0.14141    0.91317      0.08942 10.21235 1665        0
#> 
#> $`GENDER==1`
#>       group    N     mean      se   CIdown     CIup      sd    sdse     var
#> 1    Pooled 1697 -0.42518 0.05160 -0.52632 -0.32405 0.98078 0.02238 0.96213
#> 2 Composite   NA -0.42542 0.04389 -0.51144 -0.33940 0.79967 0.02551 0.64022
#> 3       GR1  848 -0.99226 0.07127 -1.13195 -0.85256 0.81277 0.03952 0.66135
#> 4       GR2  820 -0.45673 0.04790 -0.55061 -0.36286 0.73724 0.03839 0.54426
#> 5       GR3  849  0.14141 0.05124  0.04099  0.24183 0.78658 0.03226 0.61908
#>     varse   mean_0 meandiff_0 meandiffse_0  tvalue_0 df_0 pvalue_0
#> 1 0.04410  0.43947   -0.86465      0.08929  -9.68332 3332        0
#> 2 0.04090  0.43849   -0.86391      0.06456 -13.38179   NA       NA
#> 3 0.06388 -0.17761   -0.81465      0.09314  -8.74621 1665        0
#> 4 0.05673  0.48279   -0.93952      0.07151 -13.13918 1664        0
#> 5 0.05110  1.05458   -0.91317      0.08942 -10.21235 1665        0
#> 
```
