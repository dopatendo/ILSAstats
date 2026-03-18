# Mean Difference of Independent Samples with Replicate Weights

Estimates the mean difference for a single variable with replicate
weights. For a detailed explanation on how the standard errors are
estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
repmeandif(x)
```

## Arguments

- x:

  a data frame produced by
  [`repmean`](https://dopatendo.github.io/ILSAstats/reference/repmean.md)
  for a single variable or an object produced by
  [`leaguetable`](https://dopatendo.github.io/ILSAstats/reference/leaguetable.md).

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

repmeandif(reme)
#>       group1    group2      dif      se   tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.01648  0.00000 5982 1.00000
#> 2     Pooled Composite  0.00009 0.01650  0.00545   NA      NA
#> 3     Pooled       GR1 -0.01593 0.01933 -0.82411 4497 0.40992
#> 4     Pooled       GR2  0.02005 0.02297  0.87288 4481 0.38278
#> 5     Pooled       GR3  0.01611 0.02107  0.76459 4475 0.44455
#> 6  Composite    Pooled -0.00009 0.01650 -0.00545   NA      NA
#> 7  Composite Composite  0.00000 0.01652  0.00000   NA      NA
#> 8  Composite       GR1 -0.01602 0.01934 -0.82834   NA      NA
#> 9  Composite       GR2  0.01997 0.02298  0.86902   NA      NA
#> 10 Composite       GR3  0.01602 0.02108  0.75996   NA      NA
#> 11       GR1    Pooled  0.01593 0.01933  0.82411 4497 0.40992
#> 12       GR1 Composite  0.01602 0.01934  0.82834   NA      NA
#> 13       GR1       GR1  0.00000 0.02180  0.00000 3012 1.00000
#> 14       GR1       GR2  0.03599 0.02509  1.43444 2996 0.15155
#> 15       GR1       GR3  0.03204 0.02336  1.37158 2990 0.17030
#> 16       GR2    Pooled -0.02005 0.02297 -0.87288 4481 0.38278
#> 17       GR2 Composite -0.01997 0.02298 -0.86902   NA      NA
#> 18       GR2       GR1 -0.03599 0.02509 -1.43444 2996 0.15155
#> 19       GR2       GR2  0.00000 0.02799  0.00000 2980 1.00000
#> 20       GR2       GR3 -0.00395 0.02645 -0.14934 2974 0.88130
#> 21       GR3    Pooled -0.01611 0.02107 -0.76459 4475 0.44455
#> 22       GR3 Composite -0.01602 0.02108 -0.75996   NA      NA
#> 23       GR3       GR1 -0.03204 0.02336 -1.37158 2990 0.17030
#> 24       GR3       GR2  0.00395 0.02645  0.14934 2974 0.88130
#> 25       GR3       GR3  0.00000 0.02482  0.00000 2968 1.00000


# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
#>       group1    group2      dif      se    tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.03012   0.00000 6666 1.00000
#> 2     Pooled Composite  0.00059 0.03440   0.01715   NA      NA
#> 3     Pooled       GR1  0.59181 0.04482  13.20415 4999 0.00000
#> 4     Pooled       GR2 -0.02003 0.03207  -0.62457 4998 0.53228
#> 5     Pooled       GR3 -0.59064 0.04262 -13.85828 4999 0.00000
#> 6  Composite    Pooled -0.00059 0.03440  -0.01715   NA      NA
#> 7  Composite Composite  0.00000 0.03819   0.00000   NA      NA
#> 8  Composite       GR1  0.59122 0.04779  12.37121   NA      NA
#> 9  Composite       GR2 -0.02062 0.03611  -0.57103   NA      NA
#> 10 Composite       GR3 -0.59122 0.04574 -12.92567   NA      NA
#> 11       GR1    Pooled -0.59181 0.04482 -13.20415 4999 0.00000
#> 12       GR1 Composite -0.59122 0.04779 -12.37121   NA      NA
#> 13       GR1       GR1  0.00000 0.05576   0.00000 3332 1.00000
#> 14       GR1       GR2 -0.61185 0.04614 -13.26073 3331 0.00000
#> 15       GR1       GR3 -1.18245 0.05402 -21.88912 3332 0.00000
#> 16       GR2    Pooled  0.02003 0.03207   0.62457 4998 0.53228
#> 17       GR2 Composite  0.02062 0.03611   0.57103   NA      NA
#> 18       GR2       GR1  0.61185 0.04614  13.26073 3331 0.00000
#> 19       GR2       GR2  0.00000 0.03390   0.00000 3330 1.00000
#> 20       GR2       GR3 -0.57060 0.04402 -12.96229 3331 0.00000
#> 21       GR3    Pooled  0.59064 0.04262  13.85828 4999 0.00000
#> 22       GR3 Composite  0.59122 0.04574  12.92567   NA      NA
#> 23       GR3       GR1  1.18245 0.05402  21.88912 3332 0.00000
#> 24       GR3       GR2  0.57060 0.04402  12.96229 3331 0.00000
#> 25       GR3       GR3  0.00000 0.05221   0.00000 3332 1.00000

### Groups and By ----

# One variable
reme <- repmean(x = c("item01"),
                PV = FALSE,
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
#> $ALL
#>       group1    group2      dif      se   tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.01648  0.00000 5982 1.00000
#> 2     Pooled Composite  0.00009 0.01650  0.00545   NA      NA
#> 3     Pooled       GR1 -0.01593 0.01933 -0.82411 4497 0.40992
#> 4     Pooled       GR2  0.02005 0.02297  0.87288 4481 0.38278
#> 5     Pooled       GR3  0.01611 0.02107  0.76459 4475 0.44455
#> 6  Composite    Pooled -0.00009 0.01650 -0.00545   NA      NA
#> 7  Composite Composite  0.00000 0.01652  0.00000   NA      NA
#> 8  Composite       GR1 -0.01602 0.01934 -0.82834   NA      NA
#> 9  Composite       GR2  0.01997 0.02298  0.86902   NA      NA
#> 10 Composite       GR3  0.01602 0.02108  0.75996   NA      NA
#> 11       GR1    Pooled  0.01593 0.01933  0.82411 4497 0.40992
#> 12       GR1 Composite  0.01602 0.01934  0.82834   NA      NA
#> 13       GR1       GR1  0.00000 0.02180  0.00000 3012 1.00000
#> 14       GR1       GR2  0.03599 0.02509  1.43444 2996 0.15155
#> 15       GR1       GR3  0.03204 0.02336  1.37158 2990 0.17030
#> 16       GR2    Pooled -0.02005 0.02297 -0.87288 4481 0.38278
#> 17       GR2 Composite -0.01997 0.02298 -0.86902   NA      NA
#> 18       GR2       GR1 -0.03599 0.02509 -1.43444 2996 0.15155
#> 19       GR2       GR2  0.00000 0.02799  0.00000 2980 1.00000
#> 20       GR2       GR3 -0.00395 0.02645 -0.14934 2974 0.88130
#> 21       GR3    Pooled -0.01611 0.02107 -0.76459 4475 0.44455
#> 22       GR3 Composite -0.01602 0.02108 -0.75996   NA      NA
#> 23       GR3       GR1 -0.03204 0.02336 -1.37158 2990 0.17030
#> 24       GR3       GR2  0.00395 0.02645  0.14934 2974 0.88130
#> 25       GR3       GR3  0.00000 0.02482  0.00000 2968 1.00000
#> 
#> $`GENDER==0`
#>       group1    group2      dif      se   tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.02565  0.00000 2914 1.00000
#> 2     Pooled Composite -0.00001 0.02524 -0.00040   NA      NA
#> 3     Pooled       GR1  0.00655 0.03062  0.21391 2190 0.83064
#> 4     Pooled       GR2  0.04650 0.03437  1.35292 2216 0.17622
#> 5     Pooled       GR3 -0.00658 0.03086 -0.21322 2180 0.83117
#> 6  Composite    Pooled  0.00001 0.02524  0.00040   NA      NA
#> 7  Composite Composite  0.00000 0.02482  0.00000   NA      NA
#> 8  Composite       GR1  0.00656 0.03028  0.21664   NA      NA
#> 9  Composite       GR2  0.04652 0.03406  1.36583   NA      NA
#> 10 Composite       GR3 -0.00656 0.03052 -0.21494   NA      NA
#> 11       GR1    Pooled -0.00655 0.03062 -0.21391 2190 0.83064
#> 12       GR1 Composite -0.00656 0.03028 -0.21664   NA      NA
#> 13       GR1       GR1  0.00000 0.03489  0.00000 1466 1.00000
#> 14       GR1       GR2  0.03996 0.03822  1.04553 1492 0.29595
#> 15       GR1       GR3 -0.01312 0.03511 -0.37368 1456 0.70869
#> 16       GR2    Pooled -0.04650 0.03437 -1.35292 2216 0.17622
#> 17       GR2 Composite -0.04652 0.03406 -1.36583   NA      NA
#> 18       GR2       GR1 -0.03996 0.03822 -1.04553 1492 0.29595
#> 19       GR2       GR2  0.00000 0.04129  0.00000 1518 1.00000
#> 20       GR2       GR3 -0.05308 0.03842 -1.38157 1482 0.16731
#> 21       GR3    Pooled  0.00658 0.03086  0.21322 2180 0.83117
#> 22       GR3 Composite  0.00656 0.03052  0.21494   NA      NA
#> 23       GR3       GR1  0.01312 0.03511  0.37368 1456 0.70869
#> 24       GR3       GR2  0.05308 0.03842  1.38157 1482 0.16731
#> 25       GR3       GR3  0.00000 0.03532  0.00000 1446 1.00000
#> 
#> $`GENDER==1`
#>       group1    group2      dif      se   tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.02314  0.00000 3066 1.00000
#> 2     Pooled Composite  0.00032 0.02333  0.01372   NA      NA
#> 3     Pooled       GR1 -0.03716 0.02849 -1.30432 2305 0.19226
#> 4     Pooled       GR2 -0.00605 0.03196 -0.18930 2263 0.84988
#> 5     Pooled       GR3  0.03779 0.02882  1.31124 2293 0.18991
#> 6  Composite    Pooled -0.00032 0.02333 -0.01372   NA      NA
#> 7  Composite Composite  0.00000 0.02352  0.00000   NA      NA
#> 8  Composite       GR1 -0.03747 0.02864 -1.30831   NA      NA
#> 9  Composite       GR2 -0.00637 0.03210 -0.19844   NA      NA
#> 10 Composite       GR3  0.03747 0.02897  1.29341   NA      NA
#> 11       GR1    Pooled  0.03716 0.02849  1.30432 2305 0.19226
#> 12       GR1 Composite  0.03747 0.02864  1.30831   NA      NA
#> 13       GR1       GR1  0.00000 0.03298  0.00000 1544 1.00000
#> 14       GR1       GR2  0.03111 0.03602  0.86369 1502 0.38790
#> 15       GR1       GR3  0.07495 0.03326  2.25346 1532 0.02437
#> 16       GR2    Pooled  0.00605 0.03196  0.18930 2263 0.84988
#> 17       GR2 Composite  0.00637 0.03210  0.19844   NA      NA
#> 18       GR2       GR1 -0.03111 0.03602 -0.86369 1502 0.38790
#> 19       GR2       GR2  0.00000 0.03883  0.00000 1460 1.00000
#> 20       GR2       GR3  0.04384 0.03628  1.20838 1490 0.22709
#> 21       GR3    Pooled -0.03779 0.02882 -1.31124 2293 0.18991
#> 22       GR3 Composite -0.03747 0.02897 -1.29341   NA      NA
#> 23       GR3       GR1 -0.07495 0.03326 -2.25346 1532 0.02437
#> 24       GR3       GR2 -0.04384 0.03628 -1.20838 1490 0.22709
#> 25       GR3       GR3  0.00000 0.03355  0.00000 1520 1.00000
#> 

# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
#> $ALL
#>       group1    group2      dif      se    tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.03012   0.00000 6666 1.00000
#> 2     Pooled Composite  0.00059 0.03440   0.01715   NA      NA
#> 3     Pooled       GR1  0.59181 0.04482  13.20415 4999 0.00000
#> 4     Pooled       GR2 -0.02003 0.03207  -0.62457 4998 0.53228
#> 5     Pooled       GR3 -0.59064 0.04262 -13.85828 4999 0.00000
#> 6  Composite    Pooled -0.00059 0.03440  -0.01715   NA      NA
#> 7  Composite Composite  0.00000 0.03819   0.00000   NA      NA
#> 8  Composite       GR1  0.59122 0.04779  12.37121   NA      NA
#> 9  Composite       GR2 -0.02062 0.03611  -0.57103   NA      NA
#> 10 Composite       GR3 -0.59122 0.04574 -12.92567   NA      NA
#> 11       GR1    Pooled -0.59181 0.04482 -13.20415 4999 0.00000
#> 12       GR1 Composite -0.59122 0.04779 -12.37121   NA      NA
#> 13       GR1       GR1  0.00000 0.05576   0.00000 3332 1.00000
#> 14       GR1       GR2 -0.61185 0.04614 -13.26073 3331 0.00000
#> 15       GR1       GR3 -1.18245 0.05402 -21.88912 3332 0.00000
#> 16       GR2    Pooled  0.02003 0.03207   0.62457 4998 0.53228
#> 17       GR2 Composite  0.02062 0.03611   0.57103   NA      NA
#> 18       GR2       GR1  0.61185 0.04614  13.26073 3331 0.00000
#> 19       GR2       GR2  0.00000 0.03390   0.00000 3330 1.00000
#> 20       GR2       GR3 -0.57060 0.04402 -12.96229 3331 0.00000
#> 21       GR3    Pooled  0.59064 0.04262  13.85828 4999 0.00000
#> 22       GR3 Composite  0.59122 0.04574  12.92567   NA      NA
#> 23       GR3       GR1  1.18245 0.05402  21.88912 3332 0.00000
#> 24       GR3       GR2  0.57060 0.04402  12.96229 3331 0.00000
#> 25       GR3       GR3  0.00000 0.05221   0.00000 3332 1.00000
#> 
#> $`GENDER==0`
#>       group1    group2      dif      se    tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.06279   0.00000 3272 1.00000
#> 2     Pooled Composite  0.00098 0.05847   0.01676   NA      NA
#> 3     Pooled       GR1  0.61707 0.06295   9.80254 2454 0.00000
#> 4     Pooled       GR2 -0.04332 0.05662  -0.76510 2481 0.44428
#> 5     Pooled       GR3 -0.61512 0.07595  -8.09901 2453 0.00000
#> 6  Composite    Pooled -0.00098 0.05847  -0.01676   NA      NA
#> 7  Composite Composite  0.00000 0.05380   0.00000   NA      NA
#> 8  Composite       GR1  0.61610 0.05864  10.50648   NA      NA
#> 9  Composite       GR2 -0.04430 0.05178  -0.85554   NA      NA
#> 10 Composite       GR3 -0.61610 0.07242  -8.50732   NA      NA
#> 11       GR1    Pooled -0.61707 0.06295  -9.80254 2454 0.00000
#> 12       GR1 Composite -0.61610 0.05864 -10.50648   NA      NA
#> 13       GR1       GR1  0.00000 0.06311   0.00000 1636 1.00000
#> 14       GR1       GR2 -0.66040 0.05679 -11.62881 1663 0.00000
#> 15       GR1       GR3 -1.23219 0.07608 -16.19598 1635 0.00000
#> 16       GR2    Pooled  0.04332 0.05662   0.76510 2481 0.44428
#> 17       GR2 Composite  0.04430 0.05178   0.85554   NA      NA
#> 18       GR2       GR1  0.66040 0.05679  11.62881 1663 0.00000
#> 19       GR2       GR2  0.00000 0.04968   0.00000 1690 1.00000
#> 20       GR2       GR3 -0.57180 0.07093  -8.06147 1662 0.00000
#> 21       GR3    Pooled  0.61512 0.07595   8.09901 2453 0.00000
#> 22       GR3 Composite  0.61610 0.07242   8.50732   NA      NA
#> 23       GR3       GR1  1.23219 0.07608  16.19598 1635 0.00000
#> 24       GR3       GR2  0.57180 0.07093   8.06147 1662 0.00000
#> 25       GR3       GR3  0.00000 0.08715   0.00000 1634 1.00000
#> 
#> $`GENDER==1`
#>       group1    group2      dif      se    tvalue   df  pvalue
#> 1     Pooled    Pooled  0.00000 0.07298   0.00000 3392 1.00000
#> 2     Pooled Composite  0.00024 0.06774   0.00354   NA      NA
#> 3     Pooled       GR1  0.56707 0.08799   6.44471 2543 0.00000
#> 4     Pooled       GR2  0.03155 0.07040   0.44815 2515 0.65408
#> 5     Pooled       GR3 -0.56659 0.07272  -7.79139 2544 0.00000
#> 6  Composite    Pooled -0.00024 0.06774  -0.00354   NA      NA
#> 7  Composite Composite  0.00000 0.06207   0.00000   NA      NA
#> 8  Composite       GR1  0.56683 0.08370   6.77216   NA      NA
#> 9  Composite       GR2  0.03131 0.06496   0.48199   NA      NA
#> 10 Composite       GR3 -0.56683 0.06746  -8.40246   NA      NA
#> 11       GR1    Pooled -0.56707 0.08799  -6.44471 2543 0.00000
#> 12       GR1 Composite -0.56683 0.08370  -6.77216   NA      NA
#> 13       GR1       GR1  0.00000 0.10080   0.00000 1694 1.00000
#> 14       GR1       GR2 -0.53552 0.08587  -6.23640 1666 0.00000
#> 15       GR1       GR3 -1.13367 0.08778 -12.91490 1695 0.00000
#> 16       GR2    Pooled -0.03155 0.07040  -0.44815 2515 0.65408
#> 17       GR2 Composite -0.03131 0.06496  -0.48199   NA      NA
#> 18       GR2       GR1  0.53552 0.08587   6.23640 1666 0.00000
#> 19       GR2       GR2  0.00000 0.06773   0.00000 1638 1.00000
#> 20       GR2       GR3 -0.59814 0.07014  -8.52780 1667 0.00000
#> 21       GR3    Pooled  0.56659 0.07272   7.79139 2544 0.00000
#> 22       GR3 Composite  0.56683 0.06746   8.40246   NA      NA
#> 23       GR3       GR1  1.13367 0.08778  12.91490 1695 0.00000
#> 24       GR3       GR2  0.59814 0.07014   8.52780 1667 0.00000
#> 25       GR3       GR3  0.00000 0.07246   0.00000 1696 1.00000
#> 
```
