# Linear Models with Replicate Weights

Fits a linear model using [lm](https://rdrr.io/r/stats/lm.html) for
replicate weights. For a detailed explanation on how the standard errors
are estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
replm(
  formula,
  pvs = NULL,
  relatedpvs = TRUE,
  quiet = FALSE,
  summarize = TRUE,
  setup = NULL,
  df,
  wt,
  repwt,
  group = NULL,
  exclude = NULL,
  na.action = getOption("na.action"),
  method,
  aggregates = c("pooled", "composite")
)
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class): a symbolic description of the model to be fitted.

- pvs:

  if plausible values are not used, this should be `NULL`. Otherwise it
  is a list indicating which variables from `formula` should be replaced
  by which plausible values variables. For more details check the
  examples.

- relatedpvs:

  a logical value indicating if `pvs` are drawn from the same model. If
  `TRUE` (default), a total of \\n\\ estimations will be done, where
  \\n\\ is the number of plausible values for each plausible value
  variable. If `FALSE`, a total of \\n_1 \times n_2 \times n\_...\\
  estimations will be done, where \\n_i\\ is the number of plausible
  values in each plausible value variable.

- quiet:

  a logical value indicating if progress status should be shown while
  estimating models by group. Default is `FALSE`.

- summarize:

  a logical value indicating if `lm` objects should be converted to
  `summary.lm` or `summary.glm` objects and stripped from certain
  elements to reduce the size of the output object. Default is `TRUE`.

- setup:

  an optional list produced by
  [`repsetup`](https://dopatendo.github.io/ILSAstats/reference/repsetup.md).

- df:

  a data frame.

- wt:

  a string specifying the name of the column (within `df`) with the
  total weights.

- repwt:

  a string indicating the common names for the replicate weights columns
  (within `df`), or a data frame with the replicate weights.

- group:

  a string specifying the variable name (within `df`) to be used for
  grouping. Categories in `group` are treated as independent, e.g.,
  countries.

- exclude:

  a vector indicating which groups (in the same format as `group`)
  should be excluded from the pooled and composite estimates.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. The default is set by the `na.action` setting of
  [`options`](https://rdrr.io/r/base/options.html), and is
  [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.
  The ‘factory-fresh’ default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html). Another possible
  value is `NULL`, no action. Value
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) can be useful.

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

- aggregates:

  a string vector indicating which aggregates should be included,
  options are `"pooled"` and `"composite"`, both options can be used at
  the same time. If `NULL` no aggregate will be estimated.

## Value

a list.

## Examples

``` r
# Less data for shorter example
repdata2 <- repdata[1:1000,]

RW <- repcreate(df = repdata2, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----

# Simple regression - weights within df
replm(formula = Math1 ~ 1 + GENDER,
      wt = "wt", # Name of total weight column within df
      repwt = "REPWT", # Common names of replicate weights within df
      df = cbind(repdata2,RW), # Data frame
      method = "ICILS") # the name of the method aka the study name
#> 50 weights found.
#> 51 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> replm(formula = Math1 ~ 1 + GENDER, df = cbind(repdata2, RW), 
#>     wt = "wt", repwt = "REPWT", method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.5270     0.0364  14.470        0
#> GENDER       -1.0007     0.0525 -19.055        0
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model:
#> 
#> Call:
#> replm(formula = Math1 ~ 1 + GENDER, df = cbind(repdata2, RW), 
#>     wt = "wt", repwt = "REPWT", method = "ICILS")
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#>     NA     NA     NA     NA     NA 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.52703    0.03995   13.19   <2e-16 ***
#> GENDER      -1.00073    0.05629  -17.78   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.09 on 998 degrees of freedom
#> Multiple R-squared:  0.2405, Adjusted R-squared:  0.2398 
#> F-statistic: 316.1 on 1 and 998 DF,  p-value: < 2.2e-16
#> 

# Simple regression - weights as a separate data frame
replm(formula = Math1 ~ 1 + GENDER,
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       method = "ICILS") # the name of the method aka the study name
#> 51 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> replm(formula = Math1 ~ 1 + GENDER, df = repdata2, wt = "wt", 
#>     repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.5270     0.0364  14.470        0
#> GENDER       -1.0007     0.0525 -19.055        0
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model:
#> 
#> Call:
#> replm(formula = Math1 ~ 1 + GENDER, df = repdata2, wt = "wt", 
#>     repwt = RW, method = "ICILS")
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#>     NA     NA     NA     NA     NA 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.52703    0.03995   13.19   <2e-16 ***
#> GENDER      -1.00073    0.05629  -17.78   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.09 on 998 degrees of freedom
#> Multiple R-squared:  0.2405, Adjusted R-squared:  0.2398 
#> F-statistic: 316.1 on 1 and 998 DF,  p-value: < 2.2e-16
#> 

# Multiple regression
replm(formula = Math1 ~ 1 + GENDER + Reading1,
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       method = "ICILS") # the name of the method aka the study name
#> 51 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> replm(formula = Math1 ~ 1 + GENDER + Reading1, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.6057     0.0381  15.902        0
#> GENDER       -1.1401     0.0531 -21.473        0
#> Reading1      0.2532     0.0277   9.150        0
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model:
#> 
#> Call:
#> replm(formula = Math1 ~ 1 + GENDER + Reading1, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#>     NA     NA     NA     NA     NA 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.60574    0.03945  15.355   <2e-16 ***
#> GENDER      -1.14005    0.05638 -20.221   <2e-16 ***
#> Reading1     0.25321    0.02831   8.944   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.05 on 997 degrees of freedom
#> Multiple R-squared:  0.2969, Adjusted R-squared:  0.2955 
#> F-statistic: 210.5 on 2 and 997 DF,  p-value: < 2.2e-16
#> 

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 

replm(formula = Math ~ 1 + GENDER + Reading1, # Math1 now is "Math"
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       method = "ICILS") # the name of the method aka the study name
#> 153 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> replm(formula = Math ~ 1 + GENDER + Reading1, pvs = pvs, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.6043     0.0381  15.863        0
#> GENDER       -1.1397     0.0571 -19.964        0
#> Reading1      0.3023     0.0575   5.259        0
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model for first plausible value combination:
#> 
#> Call:
#> replm(formula = Math ~ 1 + GENDER + Reading1, pvs = pvs, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#>     NA     NA     NA     NA     NA 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.60574    0.03945  15.355   <2e-16 ***
#> GENDER      -1.14005    0.05638 -20.221   <2e-16 ***
#> Reading1     0.25321    0.02831   8.944   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.05 on 997 degrees of freedom
#> Multiple R-squared:  0.2969, Adjusted R-squared:  0.2955 
#> F-statistic: 210.5 on 2 and 997 DF,  p-value: < 2.2e-16
#> 

# Multiple regression - with more than one related PV variable
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 
#> $Reading
#> [1] "Reading1" "Reading2" "Reading3"
#> 

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       method = "ICILS") # the name of the method aka the study name
#> 153 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> replm(formula = Math ~ 1 + GENDER + Reading, pvs = pvs, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.6512     0.0613  10.625    0.000
#> GENDER       -1.2266     0.1024 -11.983    0.000
#> Reading1      0.3501     0.1122   3.119    0.002
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model for first plausible value combination:
#> 
#> Call:
#> replm(formula = Math ~ 1 + GENDER + Reading, pvs = pvs, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#>     NA     NA     NA     NA     NA 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.60574    0.03945  15.355   <2e-16 ***
#> GENDER      -1.14005    0.05638 -20.221   <2e-16 ***
#> Reading1     0.25321    0.02831   8.944   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.05 on 997 degrees of freedom
#> Multiple R-squared:  0.2969, Adjusted R-squared:  0.2955 
#> F-statistic: 210.5 on 2 and 997 DF,  p-value: < 2.2e-16
#> 

# Multiple regression - with more than UNrelated PV variables
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 
#> $Reading
#> [1] "Reading1" "Reading2" "Reading3"
#> 

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       relatedpvs = FALSE, # Unrelated PVs
       method = "ICILS") # the name of the method aka the study name
#> 459 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> replm(formula = Math ~ 1 + GENDER + Reading, pvs = pvs, relatedpvs = FALSE, 
#>     df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.6575     0.0583  11.278        0
#> GENDER       -1.2378     0.0963 -12.858        0
#> Reading1      0.3687     0.0723   5.101        0
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model for first plausible value combination:
#> 
#> Call:
#> replm(formula = Math ~ 1 + GENDER + Reading, pvs = pvs, relatedpvs = FALSE, 
#>     df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#>     NA     NA     NA     NA     NA 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.60574    0.03945  15.355   <2e-16 ***
#> GENDER      -1.14005    0.05638 -20.221   <2e-16 ***
#> Reading1     0.25321    0.02831   8.944   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.05 on 997 degrees of freedom
#> Multiple R-squared:  0.2969, Adjusted R-squared:  0.2955 
#> F-statistic: 210.5 on 2 and 997 DF,  p-value: < 2.2e-16
#> 


### Groups ----

# Simple regression - weights within df
replm(formula = Math1 ~ 1 + GENDER,
      wt = "wt", # Name of total weight column within df
      repwt = "REPWT", # Common names of replicate weights within df
      df = cbind(repdata2,RW), # Data frame
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name
#> 50 weights found.
#> Estimating pooled model.
#> Estimating group 1 of 3.
#> Estimating group 2 of 3.
#> Estimating group 3 of 3.
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.4976     0.0339  14.695
#> GENDER       -0.9731     0.0476 -20.462
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> replm(formula = Math1 ~ 1 + GENDER, df = cbind(repdata2, RW), 
#>     wt = "wt", repwt = "REPWT", group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   GENDER       
#> GR1 -0.076 (0.06) -0.976 (0.08)
#> GR3  1.048 (0.06) -1.002 (0.08)
#> GR2  0.520 (0.06) -0.941 (0.09)

# Simple regression - weights as a separate data frame
replm(formula = Math1 ~ 1 + GENDER,
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Estimating group 1 of 3.
#> Estimating group 2 of 3.
#> Estimating group 3 of 3.
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.4976     0.0339  14.695
#> GENDER       -0.9731     0.0476 -20.462
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> replm(formula = Math1 ~ 1 + GENDER, df = repdata2, wt = "wt", 
#>     repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   GENDER       
#> GR1 -0.076 (0.06) -0.976 (0.08)
#> GR3  1.048 (0.06) -1.002 (0.08)
#> GR2  0.520 (0.06) -0.941 (0.09)

# Multiple regression
replm(formula = Math1 ~ 1 + GENDER + Reading1,
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Estimating group 1 of 3.
#> Estimating group 2 of 3.
#> Estimating group 3 of 3.
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.5047     0.0379  13.322
#> GENDER       -0.9593     0.0483 -19.854
#> Reading1     -0.0201     0.0311  -0.646
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> replm(formula = Math1 ~ 1 + GENDER + Reading1, df = repdata2, 
#>     wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   GENDER        Reading1     
#> GR1 -0.010 (0.08) -1.013 (0.08)  0.066 (0.05)
#> GR3  1.048 (0.06) -1.002 (0.08)  0.000 (0.06)
#> GR2  0.475 (0.06) -0.863 (0.09) -0.127 (0.05)

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 

replm(formula = Math ~ 1 + GENDER + Reading1, # Math1 now is "Math"
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      pvs = pvs, # Named list
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Estimating group 1 of 3.
#> Estimating group 2 of 3.
#> Estimating group 3 of 3.
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.5157     0.0473  10.899
#> GENDER       -0.9789     0.0611 -16.018
#> Reading1      0.0594     0.0592   1.003
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> replm(formula = Math ~ 1 + GENDER + Reading1, pvs = pvs, df = repdata2, 
#>     wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  GENDER        Reading1     
#> GR1 0.027 (0.10) -0.974 (0.09)  0.126 (0.09)
#> GR3 1.019 (0.07) -1.022 (0.10)  0.058 (0.08)
#> GR2 0.501 (0.07) -0.941 (0.12) -0.005 (0.13)

# Multiple regression - with more than one related PV variable
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 
#> $Reading
#> [1] "Reading1" "Reading2" "Reading3"
#> 

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      pvs = pvs, # Named list
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Estimating group 1 of 3.
#> Estimating group 2 of 3.
#> Estimating group 3 of 3.
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.5584     0.0875   6.382
#> GENDER       -1.0350     0.0818 -12.654
#> Reading1      0.1110     0.1030   1.078
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> replm(formula = Math ~ 1 + GENDER + Reading, pvs = pvs, df = repdata2, 
#>     wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  GENDER        Reading1    
#> GR1 0.131 (0.23) -1.061 (0.11) 0.206 (0.19)
#> GR3 1.018 (0.08) -1.061 (0.11) 0.089 (0.13)
#> GR2 0.526 (0.10) -0.983 (0.19) 0.038 (0.21)

# Multiple regression - with UNrelated PV variables
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 
#> $Reading
#> [1] "Reading1" "Reading2" "Reading3"
#> 

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      pvs = pvs, # Named list
      relatedpvs = FALSE, # Unrelated PVs
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Estimating group 1 of 3.
#> Estimating group 2 of 3.
#> Estimating group 3 of 3.
#> 1836 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.5648     0.0654   8.638
#> GENDER       -1.0512     0.0715 -14.703
#> Reading1      0.1368     0.0669   2.044
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> replm(formula = Math ~ 1 + GENDER + Reading, pvs = pvs, relatedpvs = FALSE, 
#>     df = repdata2, wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  GENDER        Reading1    
#> GR1 0.146 (0.17) -1.070 (0.13) 0.221 (0.12)
#> GR3 1.011 (0.07) -1.079 (0.10) 0.118 (0.09)
#> GR2 0.537 (0.08) -1.004 (0.14) 0.071 (0.13)
```
