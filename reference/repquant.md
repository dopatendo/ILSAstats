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
        PV = FALSE,
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RWT),
        method = "ICILS")
#> Error in repquant(x = c("item01"), qtl = c(0.05, 0.25, 0.75, 0.95), PV = FALSE,     repwt = "REPWT", wt = "wt", df = cbind(repdata, RWT), method = "ICILS"): unused argument (PV = FALSE)

# One variable - weights as a separate data frame
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = FALSE,
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS")
#> Error in repquant(x = c("item01"), qtl = c(0.05, 0.25, 0.75, 0.95), PV = FALSE,     repwt = RWT, wt = "wt", df = repdata, method = "ICILS"): unused argument (PV = FALSE)

# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS")
#> Error in repquant(x = paste0("Math", 1:5), qtl = c(0.05, 0.25, 0.75, 0.95),     PV = TRUE, repwt = RWT, wt = "wt", df = repdata, method = "ICILS"): unused argument (PV = TRUE)

### Groups ----

# One variable
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = FALSE,
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> Error in repquant(x = c("item01"), qtl = c(0.05, 0.25, 0.75, 0.95), PV = FALSE,     repwt = RWT, wt = "wt", df = repdata, method = "ICILS", group = "GROUP",     exclude = "GR2"): unused argument (PV = FALSE)


# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> Error in repquant(x = paste0("Math", 1:5), qtl = c(0.05, 0.25, 0.75, 0.95),     PV = TRUE, repwt = RWT, wt = "wt", df = repdata, method = "ICILS",     group = "GROUP", exclude = "GR2"): unused argument (PV = TRUE)

### Groups and By ----

# One variable
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = FALSE,
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> Error in repquant(x = c("item01"), qtl = c(0.05, 0.25, 0.75, 0.95), PV = FALSE,     repwt = RWT, wt = "wt", df = repdata, method = "ICILS", group = "GROUP",     by = "GENDER", exclude = "GR2"): unused argument (PV = FALSE)

# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
#> Error in repquant(x = paste0("Math", 1:5), qtl = c(0.05, 0.25, 0.75, 0.95),     PV = TRUE, repwt = RWT, wt = "wt", df = repdata, method = "ICILS",     group = "GROUP", by = "GENDER", exclude = "GR2"): unused argument (PV = TRUE)
```
