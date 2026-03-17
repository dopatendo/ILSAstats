# Setup for Analysis with Replicate Weights

Creates a list with common arguments used for analysis with replicate
weights.

## Usage

``` r
repsetup(repwt, wt, df, method, group = NULL, exclude = NULL)

repsetupILSA(study, year, repwt, df, group = NULL, exclude = NULL)
```

## Arguments

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
#> Error in get(setup$repwt): object 'RW' not found

repmean(x = "Math1",setup = stp2)
#> Error in get(setup$repwt): object 'RW' not found
```
