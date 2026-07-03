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
                method = "ICILS",var = "ML",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
#> Error in `[<-.data.frame`(`*tmp*`, G1, "se", value = c(0.0116803773523922, 0.0116803773523922)): replacement has 2 rows, data has 25


# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
#> Error in `[<-.data.frame`(`*tmp*`, G1, "se", value = c(0.0270078280190179, 0.0270078280190179)): replacement has 2 rows, data has 25

### Groups and By ----

# One variable
reme <- repmean(x = c("item01"),
                PV = FALSE,
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
#> Error in `[<-.data.frame`(`*tmp*`, G1, "se", value = c(0.0116803773523922, 0.0116803773523922)): replacement has 2 rows, data has 25

# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
#> Error in `[<-.data.frame`(`*tmp*`, G1, "se", value = c(0.0270078280190179, 0.0270078280190179)): replacement has 2 rows, data has 25
```
