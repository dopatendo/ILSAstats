# Setup for Analysis with Replicate Weights

Creates a list with common arguments used for analysis with replicate
weights.

## Usage

``` r
repsetup(
  repwt = NULL,
  repindex = NULL,
  wt,
  df,
  method,
  group = NULL,
  exclude = NULL
)

repsetupILSA(
  study,
  year,
  repwt = NULL,
  repindex = NULL,
  df,
  group = NULL,
  exclude = NULL
)
```

## Arguments

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

- group:

  a string specifying the variable name (within `df`) to be used for
  grouping. Categories in `group` are treated as independent, e.g.,
  countries.

- exclude:

  a vector indicating which groups (in the same format as `group`)
  should be excluded from the pooled and composite estimates.

- study:

  a string indicating the study name. For checking available studies use
  `ILSAinfo$weights`.

- year:

  a numeric value indicating the study year. For checking available
  years use `ILSAinfo$weights`.

## Value

a list to be used in other functions.

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
stp1 <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "ICILS")
stp1
#> Error in get(x$repwt): object 'RW' not found

### Groups ----
stp2 <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "ICILS",
                 group = "GROUP", exclude = "GR2")
stp2
#> Error in get(x$repwt): object 'RW' not found


### repmean ----

repmean(x = "Math1",setup = stp1)
#>      N    mean      se      sd    sdse     var   varse
#> 1 5000 0.00191 0.01718 1.02336 0.00928 1.04726 0.01899

repmean(x = "Math1",setup = stp2)
#>       group    N     mean      se      sd    sdse     var   varse
#> 1    Pooled 3334 -0.00267 0.02192 1.08611 0.01133 1.17964 0.02462
#> 2 Composite   NA -0.00327 0.01511 0.90647 0.01033 0.82169 0.01874
#> 3       GR1 1667 -0.60164 0.02193 0.90847 0.01330 0.82532 0.02416
#> 4       GR2 1666  0.01102 0.02383 0.88556 0.01507 0.78422 0.02668
#> 5       GR3 1667  0.59510 0.02079 0.90446 0.01582 0.81805 0.02865
```
