# Proportions with Replicate Weights

Estimates proportions using replicate weights for a variable or a group
of plausible values variables and for one or more populations. For a
detailed explanation on how the standard errors are estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
repprop(
  x,
  categories = NULL,
  setup = NULL,
  repwt,
  wt,
  df,
  method,
  group = NULL,
  exclude = NULL,
  aggregates = c("pooled", "composite")
)
```

## Arguments

- x:

  a string vector specifying variable names (within `df`) for analysis.

- categories:

  a vector indicating all possible response categories. If `NULL`,
  categories will be derived from the data.

- setup:

  an optional list produced by
  [`repsetup`](https://dopatendo.github.io/ILSAstats/reference/repsetup.md).

- repwt:

  a string indicating the common names for the replicate weights columns
  (within `df`), or a data frame with the replicate weights.

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

- group:

  a string specifying the variable name (within `df`) to be used for
  grouping. Categories in `group` are treated as independent, e.g.,
  countries.

- exclude:

  a vector indicating which groups (in the same format as `group`)
  should be excluded from the pooled and composite estimates.

- aggregates:

  a string vector indicating which aggregates should be included,
  options are `"pooled"` and `"composite"`, both options can be used at
  the same time. If `NULL` no aggregate will be estimated.

## Value

a list.

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
repprop(x = c("item01"),
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
        method = "ICILS")
#> $`item01==1`
#>    N    prop      se  prop_2 propdiff_2 propdiffse_2  tvalue_2  prop_3
#> 1 69 0.01532 0.00184 0.05747   -0.04215      0.00361 -11.69008 0.22237
#>   propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1   -0.20705      0.00691 -29.94822 0.70484   -0.68952      0.00723 -95.3045
#> 
#> $`item01==2`
#>     N    prop     se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_3
#> 1 258 0.05747 0.0032 0.01532    0.04215      0.00361 11.69008 0.22237
#>   propdiff_3 propdiffse_3 tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1    -0.1649      0.00786 -20.9785 0.70484   -0.64736      0.00812 -79.7547
#> 
#> $`item01==3`
#>     N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_2
#> 1 996 0.22237 0.00648 0.01532    0.20705      0.00691 29.94822 0.05747
#>   propdiff_2 propdiffse_2 tvalue_2  prop_4 propdiff_4 propdiffse_4  tvalue_4
#> 1     0.1649      0.00786  20.9785 0.70484   -0.48247      0.01256 -38.41063
#> 
#> $`item01==4`
#>      N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_2
#> 1 3160 0.70484 0.00664 0.01532    0.68952      0.00723  95.3045 0.05747
#>   propdiff_2 propdiffse_2 tvalue_2  prop_3 propdiff_3 propdiffse_3 tvalue_3
#> 1    0.64736      0.00812  79.7547 0.22237    0.48247      0.01256 38.41063
#> 

# One variable - weights weights as a separate data frame
repprop(x = c("item01"),
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")
#> $`item01==1`
#>    N    prop      se  prop_2 propdiff_2 propdiffse_2  tvalue_2  prop_3
#> 1 69 0.01532 0.00184 0.05747   -0.04215      0.00361 -11.69008 0.22237
#>   propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1   -0.20705      0.00691 -29.94822 0.70484   -0.68952      0.00723 -95.3045
#> 
#> $`item01==2`
#>     N    prop     se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_3
#> 1 258 0.05747 0.0032 0.01532    0.04215      0.00361 11.69008 0.22237
#>   propdiff_3 propdiffse_3 tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1    -0.1649      0.00786 -20.9785 0.70484   -0.64736      0.00812 -79.7547
#> 
#> $`item01==3`
#>     N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_2
#> 1 996 0.22237 0.00648 0.01532    0.20705      0.00691 29.94822 0.05747
#>   propdiff_2 propdiffse_2 tvalue_2  prop_4 propdiff_4 propdiffse_4  tvalue_4
#> 1     0.1649      0.00786  20.9785 0.70484   -0.48247      0.01256 -38.41063
#> 
#> $`item01==4`
#>      N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_2
#> 1 3160 0.70484 0.00664 0.01532    0.68952      0.00723  95.3045 0.05747
#>   propdiff_2 propdiffse_2 tvalue_2  prop_3 propdiff_3 propdiffse_3 tvalue_3
#> 1    0.64736      0.00812  79.7547 0.22237    0.48247      0.01256 38.41063
#> 

# Multiple variables - PVs are assumed
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")
#> More than one variable provided. 'x' treated as PVs.
#> $`PVs==1`
#>      prop      se prop_2 propdiff_2 propdiffse_2  tvalue_2  prop_3 propdiff_3
#> 1 0.16165 0.00928  0.336   -0.17435      0.01632 -10.68332 0.33836   -0.17671
#>   propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1      0.01631 -10.83619 0.16399   -0.00234      0.00931 -0.25091
#> 
#> $`PVs==2`
#>    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_3 propdiff_3
#> 1 0.336 0.00897 0.16165    0.17435      0.01632 10.68332 0.33836   -0.00237
#>   propdiffse_3 tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1       0.0127 -0.18627 0.16399    0.17201      0.01299 13.24406
#> 
#> $`PVs==3`
#>      prop     se  prop_1 propdiff_1 propdiffse_1 tvalue_1 prop_2 propdiff_2
#> 1 0.33836 0.0087 0.16165    0.17671      0.01631 10.83619  0.336    0.00237
#>   propdiffse_2 tvalue_2  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1       0.0127  0.18627 0.16399    0.17438      0.01225 14.23312
#> 
#> $`PVs==4`
#>      prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1 prop_2 propdiff_2
#> 1 0.16399 0.00572 0.16165    0.00234      0.00931  0.25091  0.336   -0.17201
#>   propdiffse_2  tvalue_2  prop_3 propdiff_3 propdiffse_3  tvalue_3
#> 1      0.01299 -13.24406 0.33836   -0.17438      0.01225 -14.23312
#> 

### Groups ----

# One variable - weights within df
repprop(x = c("item01"),
        group = "GROUP",
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
        method = "ICILS")
#> $`item01==1`
#>       group  N    prop      se  prop_2 propdiff_2 propdiffse_2  tvalue_2
#> 1    Pooled 69 0.01532 0.00184 0.05747   -0.04215      0.00361 -11.69008
#> 2 Composite NA 0.01531 0.00188 0.05749   -0.04218      0.00363 -11.62939
#> 3       GR1 20 0.01336 0.00360 0.05008   -0.03671      0.00647  -5.67067
#> 4       GR2 31 0.02077 0.00341 0.06000   -0.03923      0.00705  -5.56380
#> 5       GR3 18 0.01180 0.00270 0.06240   -0.05060      0.00517  -9.77994
#>    prop_3 propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4
#> 1 0.22237   -0.20705      0.00691 -29.94822 0.70484   -0.68952      0.00723
#> 2 0.22240   -0.20709      0.00716 -28.93243 0.70480   -0.68948      0.00764
#> 3 0.22040   -0.20704      0.01261 -16.41898 0.71616   -0.70280      0.01235
#> 4 0.21432   -0.19355      0.01171 -16.53349 0.70491   -0.68414      0.01423
#> 5 0.23248   -0.22068      0.01285 -17.17754 0.69332   -0.68152      0.01305
#>    tvalue_4
#> 1 -95.30450
#> 2 -90.23505
#> 3 -56.91751
#> 4 -48.06334
#> 5 -52.21059
#> 
#> $`item01==2`
#>       group   N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
#> 1    Pooled 258 0.05747 0.00320 0.01532    0.04215      0.00361 11.69008
#> 2 Composite  NA 0.05749 0.00335 0.01531    0.04218      0.00363 11.62939
#> 3       GR1  79 0.05008 0.00503 0.01336    0.03671      0.00647  5.67067
#> 4       GR2  88 0.06000 0.00662 0.02077    0.03923      0.00705  5.56380
#> 5       GR3  91 0.06240 0.00563 0.01180    0.05060      0.00517  9.77994
#>    prop_3 propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4
#> 1 0.22237   -0.16490      0.00786 -20.97850 0.70484   -0.64736      0.00812
#> 2 0.22240   -0.16491      0.00804 -20.49973 0.70480   -0.64730      0.00869
#> 3 0.22040   -0.17032      0.01378 -12.36156 0.71616   -0.66608      0.01307
#> 4 0.21432   -0.15431      0.01381 -11.17465 0.70491   -0.64491      0.01667
#> 5 0.23248   -0.17008      0.01421 -11.97115 0.69332   -0.63092      0.01520
#>    tvalue_4
#> 1 -79.75470
#> 2 -74.49372
#> 3 -50.96709
#> 4 -38.69789
#> 5 -41.50947
#> 
#> $`item01==3`
#>       group   N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
#> 1    Pooled 996 0.22237 0.00648 0.01532    0.20705      0.00691 29.94822
#> 2 Composite  NA 0.22240 0.00668 0.01531    0.20709      0.00716 28.93243
#> 3       GR1 330 0.22040 0.01161 0.01336    0.20704      0.01261 16.41898
#> 4       GR2 322 0.21432 0.01117 0.02077    0.19355      0.01171 16.53349
#> 5       GR3 344 0.23248 0.01195 0.01180    0.22068      0.01285 17.17754
#>    prop_2 propdiff_2 propdiffse_2 tvalue_2  prop_4 propdiff_4 propdiffse_4
#> 1 0.05747    0.16490      0.00786 20.97850 0.70484   -0.48247      0.01256
#> 2 0.05749    0.16491      0.00804 20.49973 0.70480   -0.48240      0.01309
#> 3 0.05008    0.17032      0.01378 12.36156 0.71616   -0.49576      0.02217
#> 4 0.06000    0.15431      0.01381 11.17465 0.70491   -0.49059      0.02267
#> 5 0.06240    0.17008      0.01421 11.97115 0.69332   -0.46084      0.02317
#>    tvalue_4
#> 1 -38.41063
#> 2 -36.84529
#> 3 -22.35972
#> 4 -21.63820
#> 5 -19.88514
#> 
#> $`item01==4`
#>       group    N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
#> 1    Pooled 3160 0.70484 0.00664 0.01532    0.68952      0.00723 95.30450
#> 2 Composite   NA 0.70480 0.00701 0.01531    0.68948      0.00764 90.23505
#> 3       GR1 1078 0.71616 0.01133 0.01336    0.70280      0.01235 56.91751
#> 4       GR2 1050 0.70491 0.01276 0.02077    0.68414      0.01423 48.06334
#> 5       GR3 1032 0.69332 0.01230 0.01180    0.68152      0.01305 52.21059
#>    prop_2 propdiff_2 propdiffse_2 tvalue_2  prop_3 propdiff_3 propdiffse_3
#> 1 0.05747    0.64736      0.00812 79.75470 0.22237    0.48247      0.01256
#> 2 0.05749    0.64730      0.00869 74.49372 0.22240    0.48240      0.01309
#> 3 0.05008    0.66608      0.01307 50.96709 0.22040    0.49576      0.02217
#> 4 0.06000    0.64491      0.01667 38.69789 0.21432    0.49059      0.02267
#> 5 0.06240    0.63092      0.01520 41.50947 0.23248    0.46084      0.02317
#>   tvalue_3
#> 1 38.41063
#> 2 36.84529
#> 3 22.35972
#> 4 21.63820
#> 5 19.88514
#> 

# One variable - weights weights as a separate data frame
repprop(x = c("item01"),
        group = "GROUP",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")
#> $`item01==1`
#>       group  N    prop      se  prop_2 propdiff_2 propdiffse_2  tvalue_2
#> 1    Pooled 69 0.01532 0.00184 0.05747   -0.04215      0.00361 -11.69008
#> 2 Composite NA 0.01531 0.00188 0.05749   -0.04218      0.00363 -11.62939
#> 3       GR1 20 0.01336 0.00360 0.05008   -0.03671      0.00647  -5.67067
#> 4       GR2 31 0.02077 0.00341 0.06000   -0.03923      0.00705  -5.56380
#> 5       GR3 18 0.01180 0.00270 0.06240   -0.05060      0.00517  -9.77994
#>    prop_3 propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4
#> 1 0.22237   -0.20705      0.00691 -29.94822 0.70484   -0.68952      0.00723
#> 2 0.22240   -0.20709      0.00716 -28.93243 0.70480   -0.68948      0.00764
#> 3 0.22040   -0.20704      0.01261 -16.41898 0.71616   -0.70280      0.01235
#> 4 0.21432   -0.19355      0.01171 -16.53349 0.70491   -0.68414      0.01423
#> 5 0.23248   -0.22068      0.01285 -17.17754 0.69332   -0.68152      0.01305
#>    tvalue_4
#> 1 -95.30450
#> 2 -90.23505
#> 3 -56.91751
#> 4 -48.06334
#> 5 -52.21059
#> 
#> $`item01==2`
#>       group   N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
#> 1    Pooled 258 0.05747 0.00320 0.01532    0.04215      0.00361 11.69008
#> 2 Composite  NA 0.05749 0.00335 0.01531    0.04218      0.00363 11.62939
#> 3       GR1  79 0.05008 0.00503 0.01336    0.03671      0.00647  5.67067
#> 4       GR2  88 0.06000 0.00662 0.02077    0.03923      0.00705  5.56380
#> 5       GR3  91 0.06240 0.00563 0.01180    0.05060      0.00517  9.77994
#>    prop_3 propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4
#> 1 0.22237   -0.16490      0.00786 -20.97850 0.70484   -0.64736      0.00812
#> 2 0.22240   -0.16491      0.00804 -20.49973 0.70480   -0.64730      0.00869
#> 3 0.22040   -0.17032      0.01378 -12.36156 0.71616   -0.66608      0.01307
#> 4 0.21432   -0.15431      0.01381 -11.17465 0.70491   -0.64491      0.01667
#> 5 0.23248   -0.17008      0.01421 -11.97115 0.69332   -0.63092      0.01520
#>    tvalue_4
#> 1 -79.75470
#> 2 -74.49372
#> 3 -50.96709
#> 4 -38.69789
#> 5 -41.50947
#> 
#> $`item01==3`
#>       group   N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
#> 1    Pooled 996 0.22237 0.00648 0.01532    0.20705      0.00691 29.94822
#> 2 Composite  NA 0.22240 0.00668 0.01531    0.20709      0.00716 28.93243
#> 3       GR1 330 0.22040 0.01161 0.01336    0.20704      0.01261 16.41898
#> 4       GR2 322 0.21432 0.01117 0.02077    0.19355      0.01171 16.53349
#> 5       GR3 344 0.23248 0.01195 0.01180    0.22068      0.01285 17.17754
#>    prop_2 propdiff_2 propdiffse_2 tvalue_2  prop_4 propdiff_4 propdiffse_4
#> 1 0.05747    0.16490      0.00786 20.97850 0.70484   -0.48247      0.01256
#> 2 0.05749    0.16491      0.00804 20.49973 0.70480   -0.48240      0.01309
#> 3 0.05008    0.17032      0.01378 12.36156 0.71616   -0.49576      0.02217
#> 4 0.06000    0.15431      0.01381 11.17465 0.70491   -0.49059      0.02267
#> 5 0.06240    0.17008      0.01421 11.97115 0.69332   -0.46084      0.02317
#>    tvalue_4
#> 1 -38.41063
#> 2 -36.84529
#> 3 -22.35972
#> 4 -21.63820
#> 5 -19.88514
#> 
#> $`item01==4`
#>       group    N    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1
#> 1    Pooled 3160 0.70484 0.00664 0.01532    0.68952      0.00723 95.30450
#> 2 Composite   NA 0.70480 0.00701 0.01531    0.68948      0.00764 90.23505
#> 3       GR1 1078 0.71616 0.01133 0.01336    0.70280      0.01235 56.91751
#> 4       GR2 1050 0.70491 0.01276 0.02077    0.68414      0.01423 48.06334
#> 5       GR3 1032 0.69332 0.01230 0.01180    0.68152      0.01305 52.21059
#>    prop_2 propdiff_2 propdiffse_2 tvalue_2  prop_3 propdiff_3 propdiffse_3
#> 1 0.05747    0.64736      0.00812 79.75470 0.22237    0.48247      0.01256
#> 2 0.05749    0.64730      0.00869 74.49372 0.22240    0.48240      0.01309
#> 3 0.05008    0.66608      0.01307 50.96709 0.22040    0.49576      0.02217
#> 4 0.06000    0.64491      0.01667 38.69789 0.21432    0.49059      0.02267
#> 5 0.06240    0.63092      0.01520 41.50947 0.23248    0.46084      0.02317
#>   tvalue_3
#> 1 38.41063
#> 2 36.84529
#> 3 22.35972
#> 4 21.63820
#> 5 19.88514
#> 

# Multiple variables - PVs are assumed
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        group = "GROUP",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")
#> More than one variable provided. 'x' treated as PVs.
#> $`PVs==1`
#>       group    prop      se  prop_2 propdiff_2 propdiffse_2  tvalue_2  prop_3
#> 1    Pooled 0.16165 0.00928 0.33600   -0.17435      0.01632 -10.68332 0.33836
#> 2 Composite 0.16182 0.00844 0.33599   -0.17417      0.01343 -12.96714 0.33823
#> 3       GR1 0.32388 0.02222 0.42585   -0.10197      0.02856  -3.57039 0.21327
#> 4       GR2 0.12176 0.01087 0.37793   -0.25617      0.01790 -14.30941 0.36914
#> 5       GR3 0.03982 0.00533 0.20421   -0.16439      0.02208  -7.44425 0.43229
#>   propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4  tvalue_4
#> 1   -0.17671      0.01631 -10.83619 0.16399   -0.00234      0.00931  -0.25091
#> 2   -0.17642      0.01731 -10.18944 0.16395   -0.00213      0.01225  -0.17424
#> 3    0.11061      0.03939   2.80766 0.03701    0.28687      0.02675  10.72467
#> 4   -0.24738      0.02823  -8.76186 0.13117   -0.00941      0.01327  -0.70874
#> 5   -0.39247      0.01868 -21.01544 0.32368   -0.28386      0.02142 -13.25115
#> 
#> $`PVs==2`
#>       group    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_3
#> 1    Pooled 0.33600 0.00897 0.16165    0.17435      0.01632 10.68332 0.33836
#> 2 Composite 0.33599 0.00886 0.16182    0.17417      0.01343 12.96714 0.33823
#> 3       GR1 0.42585 0.01370 0.32388    0.10197      0.02856  3.57039 0.21327
#> 4       GR2 0.37793 0.01291 0.12176    0.25617      0.01790 14.30941 0.36914
#> 5       GR3 0.20421 0.01876 0.03982    0.16439      0.02208  7.44425 0.43229
#>   propdiff_3 propdiffse_3 tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1   -0.00237      0.01270 -0.18627 0.16399    0.17201      0.01299 13.24406
#> 2   -0.00224      0.01577 -0.14209 0.16395    0.17204      0.01518 11.33583
#> 3    0.21258      0.02731  7.78343 0.03701    0.38884      0.01723 22.56245
#> 4    0.00879      0.02957  0.29717 0.13117    0.24676      0.01710 14.43422
#> 5   -0.22809      0.02487 -9.17118 0.32368   -0.11948      0.03852 -3.10172
#> 
#> $`PVs==3`
#>       group    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_2
#> 1    Pooled 0.33836 0.00870 0.16165    0.17671      0.01631 10.83619 0.33600
#> 2 Composite 0.33823 0.01091 0.16182    0.17642      0.01731 10.18944 0.33599
#> 3       GR1 0.21327 0.01933 0.32388   -0.11061      0.03939 -2.80766 0.42585
#> 4       GR2 0.36914 0.02045 0.12176    0.24738      0.02823  8.76186 0.37793
#> 5       GR3 0.43229 0.01669 0.03982    0.39247      0.01868 21.01544 0.20421
#>   propdiff_2 propdiffse_2 tvalue_2  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1    0.00237      0.01270  0.18627 0.16399    0.17438      0.01225 14.23312
#> 2    0.00224      0.01577  0.14209 0.16395    0.17428      0.01652 10.54814
#> 3   -0.21258      0.02731 -7.78343 0.03701    0.17626      0.01746 10.09498
#> 4   -0.00879      0.02957 -0.29717 0.13117    0.23797      0.02869  8.29586
#> 5    0.22809      0.02487  9.17118 0.32368    0.10861      0.03646  2.97904
#> 
#> $`PVs==4`
#>       group    prop      se  prop_1 propdiff_1 propdiffse_1  tvalue_1  prop_2
#> 1    Pooled 0.16399 0.00572 0.16165    0.00234      0.00931   0.25091 0.33600
#> 2 Composite 0.16395 0.00877 0.16182    0.00213      0.01225   0.17424 0.33599
#> 3       GR1 0.03701 0.00724 0.32388   -0.28687      0.02675 -10.72467 0.42585
#> 4       GR2 0.13117 0.01085 0.12176    0.00941      0.01327   0.70874 0.37793
#> 5       GR3 0.32368 0.02286 0.03982    0.28386      0.02142  13.25115 0.20421
#>   propdiff_2 propdiffse_2  tvalue_2  prop_3 propdiff_3 propdiffse_3  tvalue_3
#> 1   -0.17201      0.01299 -13.24406 0.33836   -0.17438      0.01225 -14.23312
#> 2   -0.17204      0.01518 -11.33583 0.33823   -0.17428      0.01652 -10.54814
#> 3   -0.38884      0.01723 -22.56245 0.21327   -0.17626      0.01746 -10.09498
#> 4   -0.24676      0.01710 -14.43422 0.36914   -0.23797      0.02869  -8.29586
#> 5    0.11948      0.03852   3.10172 0.43229   -0.10861      0.03646  -2.97904
#> 

# Multiple variables - excluding one group
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        group = "GROUP",
        exclude = "GR2",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")
#> More than one variable provided. 'x' treated as PVs.
#> $`PVs==1`
#>       group    prop      se  prop_2 propdiff_2 propdiffse_2  tvalue_2  prop_3
#> 1    Pooled 0.18171 0.01250 0.31492   -0.13321      0.02131  -6.24989 0.32289
#> 2 Composite 0.18185 0.01143 0.31503   -0.13318      0.01805  -7.37807 0.32278
#> 3       GR1 0.32388 0.02222 0.42585   -0.10197      0.02856  -3.57039 0.21327
#> 4       GR2 0.12176 0.01087 0.37793   -0.25617      0.01790 -14.30941 0.36914
#> 5       GR3 0.03982 0.00533 0.20421   -0.16439      0.02208  -7.44425 0.43229
#>   propdiff_3 propdiffse_3  tvalue_3  prop_4 propdiff_4 propdiffse_4  tvalue_4
#> 1   -0.14118      0.02421  -5.83148 0.18049    0.00122      0.01221   0.09984
#> 2   -0.14093      0.02180  -6.46528 0.18035    0.00150      0.01713   0.08768
#> 3    0.11061      0.03939   2.80766 0.03701    0.28687      0.02675  10.72467
#> 4   -0.24738      0.02823  -8.76186 0.13117   -0.00941      0.01327  -0.70874
#> 5   -0.39247      0.01868 -21.01544 0.32368   -0.28386      0.02142 -13.25115
#> 
#> $`PVs==2`
#>       group    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_3
#> 1    Pooled 0.31492 0.01105 0.18171    0.13321      0.02131  6.24989 0.32289
#> 2 Composite 0.31503 0.01162 0.18185    0.13318      0.01805  7.37807 0.32278
#> 3       GR1 0.42585 0.01370 0.32388    0.10197      0.02856  3.57039 0.21327
#> 4       GR2 0.37793 0.01291 0.12176    0.25617      0.01790 14.30941 0.36914
#> 5       GR3 0.20421 0.01876 0.03982    0.16439      0.02208  7.44425 0.43229
#>   propdiff_3 propdiffse_3 tvalue_3  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1   -0.00797      0.01645 -0.48471 0.18049    0.13443      0.01811  7.42144
#> 2   -0.00776      0.01847 -0.41990 0.18035    0.13468      0.02110  6.38315
#> 3    0.21258      0.02731  7.78343 0.03701    0.38884      0.01723 22.56245
#> 4    0.00879      0.02957  0.29717 0.13117    0.24676      0.01710 14.43422
#> 5   -0.22809      0.02487 -9.17118 0.32368   -0.11948      0.03852 -3.10172
#> 
#> $`PVs==3`
#>       group    prop      se  prop_1 propdiff_1 propdiffse_1 tvalue_1  prop_2
#> 1    Pooled 0.32289 0.01381 0.18171    0.14118      0.02421  5.83148 0.31492
#> 2 Composite 0.32278 0.01277 0.18185    0.14093      0.02180  6.46528 0.31503
#> 3       GR1 0.21327 0.01933 0.32388   -0.11061      0.03939 -2.80766 0.42585
#> 4       GR2 0.36914 0.02045 0.12176    0.24738      0.02823  8.76186 0.37793
#> 5       GR3 0.43229 0.01669 0.03982    0.39247      0.01868 21.01544 0.20421
#>   propdiff_2 propdiffse_2 tvalue_2  prop_4 propdiff_4 propdiffse_4 tvalue_4
#> 1    0.00797      0.01645  0.48471 0.18049    0.14240      0.02171  6.55819
#> 2    0.00776      0.01847  0.41990 0.18035    0.14244      0.02021  7.04717
#> 3   -0.21258      0.02731 -7.78343 0.03701    0.17626      0.01746 10.09498
#> 4   -0.00879      0.02957 -0.29717 0.13117    0.23797      0.02869  8.29586
#> 5    0.22809      0.02487  9.17118 0.32368    0.10861      0.03646  2.97904
#> 
#> $`PVs==4`
#>       group    prop      se  prop_1 propdiff_1 propdiffse_1  tvalue_1  prop_2
#> 1    Pooled 0.18049 0.00979 0.18171   -0.00122      0.01221  -0.09984 0.31492
#> 2 Composite 0.18035 0.01199 0.18185   -0.00150      0.01713  -0.08768 0.31503
#> 3       GR1 0.03701 0.00724 0.32388   -0.28687      0.02675 -10.72467 0.42585
#> 4       GR2 0.13117 0.01085 0.12176    0.00941      0.01327   0.70874 0.37793
#> 5       GR3 0.32368 0.02286 0.03982    0.28386      0.02142  13.25115 0.20421
#>   propdiff_2 propdiffse_2  tvalue_2  prop_3 propdiff_3 propdiffse_3  tvalue_3
#> 1   -0.13443      0.01811  -7.42144 0.32289   -0.14240      0.02171  -6.55819
#> 2   -0.13468      0.02110  -6.38315 0.32278   -0.14244      0.02021  -7.04717
#> 3   -0.38884      0.01723 -22.56245 0.21327   -0.17626      0.01746 -10.09498
#> 4   -0.24676      0.01710 -14.43422 0.36914   -0.23797      0.02869  -8.29586
#> 5    0.11948      0.03852   3.10172 0.43229   -0.10861      0.03646  -2.97904
#> 
```
