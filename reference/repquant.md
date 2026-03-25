# Quantiles with Replicate Weights

Estimates quantiles with replicate weights for a variable or a group of
variables and for one or more populations. For a detailed explanation on
how the standard errors are estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
repquant(
  x,
  qtl = c(0.05, 0.25, 0.75, 0.95),
  setup = NULL,
  repwt,
  wt,
  df,
  method,
  group = NULL,
  by = NULL,
  exclude = NULL
)
```

## Arguments

- x:

  a string vector specifying variable names (within `df`) for analysis.

- qtl:

  a numeric vector indicating the desired quantiles (between 0 and 1).

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

- by:

  a string specifying a second variable (within `df`) for grouping.
  Categories used in `by` are not considered independent, e.g., gender
  within a country. If used, the output will be a list with the same
  length as the unique values of `by`. This can only be used for
  analyses with one variable or a group of PVs.

- exclude:

  a vector indicating which groups (in the same format as `group`)
  should be excluded from the pooled and composite estimates.

## Value

a data frame or a list.

## Examples

``` r
RWT <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----

# One variable - weights within df
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RWT),
        method = "ICILS")
#>   variable P05 P05se P25 P25se P75 P75se P95 P95se
#> 1   item01   2     0   3     0   4     0   4     0

# One variable - weights as a separate data frame
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS")
#>   variable P05 P05se P25 P25se P75 P75se P95 P95se
#> 1   item01   2     0   3     0   4     0   4     0

# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS")
#> More than one variable provided. 'x' treated as PVs.
#>   variable      P05   P05se      P25   P25se     P75   P75se     P95  P95se
#> 1      PVs -1.68118 0.04039 -0.68384 0.02698 0.69452 0.02582 1.68116 0.0399

### Groups ----

# One variable
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#>   variable     group P05 P05se P25 P25se P75 P75se P95 P95se
#> 1   item01    Pooled   2     0   3     0   4     0   4     0
#> 2   item01 Composite   2     0   3     0   4     0   4     0
#> 3   item01       GR1   2     0   3     0   4     0   4     0
#> 4   item01       GR2   2     0   3     0   4     0   4     0
#> 5   item01       GR3   2     0   3     0   4     0   4     0


# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> More than one variable provided. 'x' treated as PVs.
#>   variable     group      P05   P05se      P25   P25se     P75   P75se     P95
#> 1      PVs    Pooled -1.78116 0.05793 -0.74310 0.04099 0.73814 0.03017 1.76034
#> 2      PVs Composite -1.50158 0.04746 -0.59553 0.04241 0.59988 0.03536 1.48907
#> 3      PVs       GR1 -2.10212 0.07188 -1.18328 0.05214 0.01332 0.05149 0.88008
#> 4      PVs       GR2 -1.44166 0.06880 -0.57888 0.04137 0.61342 0.03943 1.50918
#> 5      PVs       GR3 -0.90104 0.06199 -0.00778 0.06690 1.18644 0.04848 2.09806
#>     P95se
#> 1 0.04780
#> 2 0.06260
#> 3 0.09722
#> 4 0.06852
#> 5 0.07890

### Groups and By ----

# One variable
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite

# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> More than one variable provided. 'x' treated as PVs.
```
