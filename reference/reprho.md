# Correlations with Replicate Weights

Estimates correlation coefficients using replicate weights. For a
detailed explanation on how the standard errors are estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
reprho(
  x = NULL,
  pv = NULL,
  pv2 = NULL,
  relatedpvs = TRUE,
  setup = NULL,
  repwt,
  wt,
  df,
  rho = c("pearson", "spearman", "polychoric"),
  method,
  group = NULL,
  exclude = NULL,
  aggregates = c("pooled", "composite")
)
```

## Arguments

- x:

  a string vector specifying variable names (within `df`) for analysis.
  If `pv` is `NULL`, this function estimates correlations between all
  variables in the vector. If `pv2` is NOT `NULL`, then `x` should be
  set to `NULL`.

- pv:

  a string vector indicating the variable names for all plausible values
  of a construct. If not `NULL`, this function estimates correlations
  only between `x` and the plausible values construct.

- pv2:

  a string vector indicating the variable names for all plausible values
  of a second construct (distinct from `pv`).

- relatedpvs:

  a logical value indicating if `pv` and `pv2` are drawn from the same
  model, and have the same number of plausible values. If `TRUE`
  (default), a total of \\n\\ estimations will be done, where \\n\\ is
  the number of plausible values of each. If `FALSE`, a total of \\n_1
  \times n_2\\ estimations will be done, where \\n_1\\ is the number of
  plausible values in `pv` and \\n_2\\ is the number of plausible values
  in `pv2`.

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

- rho:

  a string indicating the correlation coefficient to be computed:
  `"pearson"`, `"polychoric"`, or `"spearman"` (lower or uppercase).

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
  should be excluded from the estimation of pooled and composite
  estimates.

- aggregates:

  a string vector indicating which aggregates should be included,
  options are `"pooled"` and `"composite"`, both options can be used at
  the same time. If `NULL` no aggregate will be estimated.

## Value

a data frame.

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

# Non PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = NULL,
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")
#>   variable1 variable2      rho      se    n    tvalue pvalue
#> 1    GENDER     Math1 -0.46153 0.01068 5000 -43.22308      0
#> 2    GENDER     Math2 -0.43280 0.01101 5000 -39.30958      0
#> 3    GENDER     Math3 -0.46048 0.01106 5000 -41.65126      0
#> 4     Math1     Math2  0.84659 0.00387 5000 218.53970      0
#> 5     Math1     Math3  0.84121 0.00405 5000 207.59920      0
#> 6     Math2     Math3  0.81903 0.00429 5000 190.84450      0

# X var and PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = paste0("Reading",1:5),
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")
#>   variable1 variable2     rho      se    n  tvalue pvalue
#> 1    GENDER       PVs 0.38427 0.06145 5000 6.25331  0e+00
#> 2     Math1       PVs 0.15704 0.03603 5000 4.35889  1e-05
#> 3     Math2       PVs 0.18971 0.02750 5000 6.89934  0e+00
#> 4     Math3       PVs 0.15134 0.03067 5000 4.93366  0e+00

# PVs and PVs (related)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")
#>       rho      se    n  tvalue  pvalue
#> 1 0.14209 0.03848 5000 3.69231 0.00022

# PVs and PVs (UNrelated)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       relatedpvs = FALSE,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")
#>       rho      se    n  tvalue pvalue
#> 1 0.15895 0.03957 5000 4.01729  6e-05


### Groups ----

# Non PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = NULL,
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")
#>    variable1 variable2     group      rho      se    n    tvalue pvalue
#> 1     GENDER     Math1    Pooled -0.46153 0.01068 5000 -43.22308      0
#> 2     GENDER     Math1 Composite -0.52508 0.00907   NA -57.87475     NA
#> 3     GENDER     Math1       GR1 -0.48826 0.01665 1667 -29.31893      0
#> 4     GENDER     Math1       GR2 -0.55951 0.01413 1666 -39.61081      0
#> 5     GENDER     Math1       GR3 -0.52747 0.01625 1667 -32.46554      0
#> 6     GENDER     Math2    Pooled -0.43280 0.01101 5000 -39.30958      0
#> 7     GENDER     Math2 Composite -0.48714 0.00938   NA -51.90755     NA
#> 8     GENDER     Math2       GR1 -0.43410 0.01692 1667 -25.65883      0
#> 9     GENDER     Math2       GR2 -0.52950 0.01602 1666 -33.05124      0
#> 10    GENDER     Math2       GR3 -0.49781 0.01580 1667 -31.49878      0
#> 11    GENDER     Math3    Pooled -0.46048 0.01106 5000 -41.65126      0
#> 12    GENDER     Math3 Composite -0.53220 0.00925   NA -57.52529     NA
#> 13    GENDER     Math3       GR1 -0.49936 0.01801 1667 -27.72432      0
#> 14    GENDER     Math3       GR2 -0.55671 0.01455 1666 -38.26024      0
#> 15    GENDER     Math3       GR3 -0.54053 0.01530 1667 -35.32149      0
#> 16     Math1     Math2    Pooled  0.84659 0.00387 5000 218.53970      0
#> 17     Math1     Math2 Composite  0.80378 0.00537   NA 149.72312     NA
#> 18     Math1     Math2       GR1  0.79602 0.01112 1667  71.57672      0
#> 19     Math1     Math2       GR2  0.80766 0.00868 1666  93.00955      0
#> 20     Math1     Math2       GR3  0.80766 0.00776 1667 104.01320      0
#> 21     Math1     Math3    Pooled  0.84121 0.00405 5000 207.59920      0
#> 22     Math1     Math3 Composite  0.79142 0.00529   NA 149.67382     NA
#> 23     Math1     Math3       GR1  0.78923 0.00946 1667  83.40460      0
#> 24     Math1     Math3       GR2  0.77839 0.00992 1666  78.45221      0
#> 25     Math1     Math3       GR3  0.80665 0.00798 1667 101.10920      0
#> 26     Math2     Math3    Pooled  0.81903 0.00429 5000 190.84450      0
#> 27     Math2     Math3 Composite  0.76579 0.00634   NA 120.73227     NA
#> 28     Math2     Math3       GR1  0.76268 0.01141 1667  66.86056      0
#> 29     Math2     Math3       GR2  0.75316 0.01059 1666  71.08708      0
#> 30     Math2     Math3       GR3  0.78151 0.01094 1667  71.42797      0

# X var and PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = paste0("Reading",1:5),
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")
#>    variable1 variable2     group      rho      se    n   tvalue  pvalue
#> 1     GENDER       PVs    Pooled  0.38427 0.06145 5000  6.25331 0.00000
#> 2     GENDER       PVs Composite  0.45172 0.04264   NA 10.59468      NA
#> 3     GENDER       PVs       GR1  0.45084 0.07945 1667  5.67435 0.00000
#> 4     GENDER       PVs       GR2  0.45976 0.07022 1666  6.54767 0.00000
#> 5     GENDER       PVs       GR3  0.44455 0.07154 1667  6.21433 0.00000
#> 6      Math1       PVs    Pooled  0.15704 0.03603 5000  4.35889 0.00001
#> 7      Math1       PVs Composite -0.12633 0.03093   NA -4.08410      NA
#> 8      Math1       PVs       GR1 -0.05015 0.05796 1667 -0.86526 0.38702
#> 9      Math1       PVs       GR2 -0.21636 0.04946 1666 -4.37447 0.00001
#> 10     Math1       PVs       GR3 -0.11249 0.05297 1667 -2.12346 0.03386
#> 11     Math2       PVs    Pooled  0.18971 0.02750 5000  6.89934 0.00000
#> 12     Math2       PVs Composite -0.06939 0.02406   NA -2.88452      NA
#> 13     Math2       PVs       GR1  0.00558 0.03784 1667  0.14741 0.88283
#> 14     Math2       PVs       GR2 -0.15962 0.04469 1666 -3.57132 0.00037
#> 15     Math2       PVs       GR3 -0.05414 0.04218 1667 -1.28339 0.19953
#> 16     Math3       PVs    Pooled  0.15134 0.03067 5000  4.93366 0.00000
#> 17     Math3       PVs Composite -0.15367 0.02588   NA -5.93729      NA
#> 18     Math3       PVs       GR1 -0.09128 0.04686 1667 -1.94806 0.05158
#> 19     Math3       PVs       GR2 -0.23715 0.04104 1666 -5.77796 0.00000
#> 20     Math3       PVs       GR3 -0.13259 0.04636 1667 -2.86008 0.00429

# PVs and PVs (related)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")
#>       group      rho      se    n   tvalue  pvalue
#> 1    Pooled  0.14209 0.03848 5000  3.69231 0.00022
#> 2 Composite -0.14197 0.03827   NA -3.70947      NA
#> 3       GR1 -0.07250 0.06419 1667 -1.12933 0.25892
#> 4       GR2 -0.22991 0.07829 1666 -2.93648 0.00337
#> 5       GR3 -0.12351 0.05415 1667 -2.28080 0.02269

# PVs and PVs (UNrelated)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       relatedpvs = FALSE,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")
#>       group      rho      se    n   tvalue  pvalue
#> 1    Pooled  0.15895 0.03957 5000  4.01729 0.00006
#> 2 Composite -0.11932 0.03486   NA -3.42286      NA
#> 3       GR1 -0.04893 0.06081 1667 -0.80465 0.42113
#> 4       GR2 -0.20752 0.06327 1666 -3.27987 0.00106
#> 5       GR3 -0.10153 0.05689 1667 -1.78453 0.07452
```
