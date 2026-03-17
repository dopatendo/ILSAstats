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
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#>       group    N nzones    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 2992     50 3.62342 0.01165 3.60058 3.64626 0.65015 0.01420 0.42270
#> 2 Composite   NA     NA 3.62334 0.01168 3.60044 3.64623 0.64994 0.01419 0.42250
#> 3       GR1 1507     50 3.63936 0.01542 3.60914 3.66957 0.64103 0.02036 0.41092
#> 4       GR2 1491     50 3.60337 0.01979 3.56457 3.64216 0.69566 0.02093 0.48394
#> 5       GR3 1485     50 3.60732 0.01755 3.57292 3.64171 0.65886 0.01977 0.43409
#>     varse
#> 1 0.01845
#> 2 0.01840
#> 3 0.02603
#> 4 0.02912
#> 5 0.02602


# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#>       group    N nzones     mean      se   CIdown     CIup      sd    sdse
#> 1    Pooled 3334     50 -0.00146 0.02130 -0.04321  0.04028 1.08175 0.01625
#> 2 Composite   NA     NA -0.00205 0.02701 -0.05498  0.05088 0.90551 0.01456
#> 3       GR1 1667     50 -0.59327 0.03943 -0.67056 -0.51599 0.90130 0.02092
#> 4       GR2 1666     50  0.01857 0.02397 -0.02841  0.06555 0.89365 0.02078
#> 5       GR3 1667     50  0.58917 0.03692  0.51682  0.66153 0.90971 0.02026
#>       var   varse
#> 1 1.17026 0.03519
#> 2 0.82008 0.02638
#> 3 0.81251 0.03764
#> 4 0.79873 0.03712
#> 5 0.82766 0.03697

### Groups and By ----

# One variable
reme <- repmean(x = c("item01"),
                PV = FALSE,
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#> $ALL
#>       group    N nzones    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 2992     50 3.62342 0.01165 3.60058 3.64626 0.65015 0.01420 0.42270
#> 2 Composite   NA     NA 3.62334 0.01168 3.60044 3.64623 0.64994 0.01419 0.42250
#> 3       GR1 1507     50 3.63936 0.01542 3.60914 3.66957 0.64103 0.02036 0.41092
#> 4       GR2 1491     50 3.60337 0.01979 3.56457 3.64216 0.69566 0.02093 0.48394
#> 5       GR3 1485     50 3.60732 0.01755 3.57292 3.64171 0.65886 0.01977 0.43409
#>     varse
#> 1 0.01845
#> 2 0.01840
#> 3 0.02603
#> 4 0.02912
#> 5 0.02602
#> 
#> $`GENDER==0`
#>       group    N nzones    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 1458     50 3.63697 0.01813 3.60143 3.67251 0.63362 0.02024 0.40148
#> 2 Composite   NA     NA 3.63698 0.01755 3.60258 3.67138 0.63356 0.01984 0.40142
#> 3       GR1  734     50 3.63042 0.02467 3.58206 3.67878 0.63818 0.02781 0.40727
#> 4       GR2  760     50 3.59046 0.02919 3.53324 3.64768 0.72030 0.02785 0.51883
#> 5       GR3  724     50 3.64354 0.02497 3.59460 3.69249 0.62894 0.02832 0.39557
#>     varse  mean_1 meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1 0.02560 3.61060    0.02637      0.02540  1.03818 2990  0.29927
#> 2 0.02507 3.61028    0.02670      0.02492  1.07127   NA       NA
#> 3 0.03545 3.64776   -0.01734      0.03674 -0.47194 1505  0.63704
#> 4 0.04025 3.61665   -0.02618      0.04066 -0.64399 1489  0.51968
#> 5 0.03545 3.57281    0.07074      0.03369  2.09948 1483  0.03594
#> 
#> $`GENDER==1`
#>       group    N nzones    mean      se  CIdown    CIup      sd    sdse     var
#> 1    Pooled 1534     50 3.61060 0.01636 3.57853 3.64267 0.66517 0.01895 0.44245
#> 2 Composite   NA     NA 3.61028 0.01663 3.57768 3.64288 0.66397 0.01904 0.44127
#> 3       GR1  773     50 3.64776 0.02332 3.60205 3.69347 0.64359 0.02991 0.41420
#> 4       GR2  731     50 3.61665 0.02746 3.56283 3.67046 0.66909 0.02575 0.44769
#> 5       GR3  761     50 3.57281 0.02372 3.52632 3.61930 0.68436 0.02356 0.46834
#>     varse  mean_0 meandiff_0 meandiffse_0 tvalue_0 df_0 pvalue_0
#> 1 0.02521 3.63697   -0.02637      0.02540 -1.03818 2990  0.29927
#> 2 0.02507 3.63698   -0.02670      0.02492 -1.07127   NA       NA
#> 3 0.03840 3.63042    0.01734      0.03674  0.47194 1505  0.63704
#> 4 0.03444 3.59046    0.02618      0.04066  0.64399 1489  0.51968
#> 5 0.03224 3.64354   -0.07074      0.03369 -2.09948 1483  0.03594
#> 

# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeanCI(reme)
#> $ALL
#>       group    N nzones     mean      se   CIdown     CIup      sd    sdse
#> 1    Pooled 3334     50 -0.00146 0.02130 -0.04321  0.04028 1.08175 0.01625
#> 2 Composite   NA     NA -0.00205 0.02701 -0.05498  0.05088 0.90551 0.01456
#> 3       GR1 1667     50 -0.59327 0.03943 -0.67056 -0.51599 0.90130 0.02092
#> 4       GR2 1666     50  0.01857 0.02397 -0.02841  0.06555 0.89365 0.02078
#> 5       GR3 1667     50  0.58917 0.03692  0.51682  0.66153 0.90971 0.02026
#>       var   varse
#> 1 1.17026 0.03519
#> 2 0.82008 0.02638
#> 3 0.81251 0.03764
#> 4 0.79873 0.03712
#> 5 0.82766 0.03697
#> 
#> $`GENDER==0`
#>       group    N nzones     mean      se   CIdown     CIup      sd    sdse
#> 1    Pooled 1637     50  0.43947 0.04440  0.35244  0.52649 1.00184 0.01593
#> 2 Composite   NA     NA  0.43849 0.03804  0.36393  0.51305 0.78939 0.02289
#> 3       GR1  819     50 -0.17761 0.04462 -0.26507 -0.09015 0.79316 0.03068
#> 4       GR2  846     50  0.48279 0.03513  0.41394  0.55164 0.78081 0.03692
#> 5       GR3  818     50  1.05458 0.06162  0.93380  1.17537 0.78563 0.03399
#>       var   varse   mean_1 meandiff_1 meandiffse_1 tvalue_1 df_1 pvalue_1
#> 1 1.00371 0.03191 -0.42518    0.86465      0.08929  9.68332 3332        0
#> 2 0.62360 0.03604 -0.42542    0.86391      0.06456 13.38179   NA       NA
#> 3 0.62944 0.04864 -0.99226    0.81465      0.09314  8.74621 1665        0
#> 4 0.61035 0.05731 -0.45673    0.93952      0.07151 13.13918 1664        0
#> 5 0.61776 0.05321  0.14141    0.91317      0.08942 10.21235 1665        0
#> 
#> $`GENDER==1`
#>       group    N nzones     mean      se   CIdown     CIup      sd    sdse
#> 1    Pooled 1697     50 -0.42518 0.05160 -0.52632 -0.32405 0.98059 0.02236
#> 2 Composite   NA     NA -0.42542 0.04389 -0.51144 -0.33940 0.79936 0.02544
#> 3       GR1  848     50 -0.99226 0.07127 -1.13195 -0.85256 0.81245 0.03932
#> 4       GR2  820     50 -0.45673 0.04790 -0.55061 -0.36286 0.73694 0.03836
#> 5       GR3  849     50  0.14141 0.05124  0.04099  0.24183 0.78627 0.03229
#>       var   varse   mean_0 meandiff_0 meandiffse_0  tvalue_0 df_0 pvalue_0
#> 1 0.96175 0.04405  0.43947   -0.86465      0.08929  -9.68332 3332        0
#> 2 0.63971 0.04078  0.43849   -0.86391      0.06456 -13.38179   NA       NA
#> 3 0.66083 0.06353 -0.17761   -0.81465      0.09314  -8.74621 1665        0
#> 4 0.54382 0.05667  0.48279   -0.93952      0.07151 -13.13918 1664        0
#> 5 0.61859 0.05114  1.05458   -0.91317      0.08942 -10.21235 1665        0
#> 
```
