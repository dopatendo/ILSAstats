# Generalized Linear Models with Replicate Weights

Fits a generalized linear model using
[glm](https://rdrr.io/r/stats/glm.html) for replicate weights. For a
detailed explanation on how the standard errors are estimated see
[`repse`](https://dopatendo.github.io/ILSAstats/reference/repse.md).

## Usage

``` r
repglm(
  formula,
  family = stats::gaussian,
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
  method
)
```

## Arguments

- formula:

  an object of class
  `"`[`formula`](https://rdrr.io/r/stats/formula.html)`"` (or one that
  can be coerced to that class): a symbolic description of the model to
  be fitted. The details of model specification are given under
  ‘Details’.

- family:

  a description of the error distribution and link function to be used
  in the model. For `glm` this can be a character string naming a family
  function, a family function or the result of a call to a family
  function. For `glm.fit` only the third option is supported. (See
  [`family`](https://rdrr.io/r/stats/family.html) for details of family
  functions.)

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

## Value

a list with the standard errors and the total weights models.

## Examples

``` r
# Less data for shorter example
repdata2 <- repdata[1:1000,]

# Creation of replicate weights
RW <- repcreate(df = repdata2, # the data frame with all the information
                wt = "wt", # the total weights column name
                jkzone = "jkzones", # the jkzones column name
                jkrep = "jkrep", # the jkreps column name
                repwtname = "REPWT", # the desired name for the rep weights
                reps = 50, # the number of replications
                method = "ICILS") # the name of the method aka the study name

### No groups ----

# Simple regression - default family = gaussian
repglm(formula = GENDER ~ 1 + Math1,
        family = gaussian, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        method = "ICILS") # the name of the method aka the study name
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 51 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> repglm(formula = GENDER ~ 1 + Math1, family = gaussian, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.5092     0.0137  37.131        0
#> Math1        -0.2404     0.0100 -24.089        0
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model:
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = gaussian, df = repdata2, 
#>     wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.50917    0.01380   36.91   <2e-16 ***
#> Math1       -0.24035    0.01352  -17.78   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.2854853)
#> 
#>     Null deviance: 375.15  on 999  degrees of freedom
#> Residual deviance: 284.91  on 998  degrees of freedom
#> AIC: 1201.2
#> 
#> Number of Fisher Scoring iterations: 2
#> 


# Simple regression - change link function
repglm(formula = GENDER ~ 1 + Math1,
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        method = "ICILS") # the name of the method aka the study name
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 51 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> repglm(formula = GENDER ~ 1 + Math1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.0542     0.0743    0.73    0.465
#> Math1        -1.3100     0.0855  -15.33    0.000
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model:
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.05424    0.07225   0.751    0.453    
#> Math1       -1.31000    0.09499 -13.791   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.466464)
#> 
#>     Null deviance: 2080.3  on 999  degrees of freedom
#> Residual deviance: 1661.5  on 998  degrees of freedom
#> AIC: NA
#> 
#> Number of Fisher Scoring iterations: 4
#> 

# Multiple regression
repglm(formula = GENDER ~ 1 + Math1 + Reading1,
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        method = "ICILS") # the name of the method aka the study name
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 51 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> repglm(formula = GENDER ~ 1 + Math1 + Reading1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.1080     0.0859   1.258    0.209
#> Math1        -1.6827     0.1194 -14.098    0.000
#> Reading1      1.0634     0.0965  11.023    0.000
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model:
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1 + Reading1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.10804    0.07994   1.352    0.177    
#> Math1       -1.68273    0.11580 -14.531   <2e-16 ***
#> Reading1     1.06341    0.09644  11.027   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.487731)
#> 
#>     Null deviance: 2080.3  on 999  degrees of freedom
#> Residual deviance: 1426.5  on 997  degrees of freedom
#> AIC: NA
#> 
#> Number of Fisher Scoring iterations: 5
#> 

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 

repglm(formula = GENDER ~ 1 + Math + Reading1, # Math1 now is "Math"
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        pvs = pvs, # Named list
        method = "ICILS") # the name of the method aka the study name
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 153 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> repglm(formula = GENDER ~ 1 + Math + Reading1, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.1114     0.0886   1.257    0.209
#> Math1        -1.7798     0.1604 -11.098    0.000
#> Reading1      1.1877     0.1628   7.295    0.000
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model for first plausible value combination:
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading1, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.10804    0.07994   1.352    0.177    
#> Math1       -1.68273    0.11580 -14.531   <2e-16 ***
#> Reading1     1.06341    0.09644  11.027   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.487731)
#> 
#>     Null deviance: 2080.3  on 999  degrees of freedom
#> Residual deviance: 1426.5  on 997  degrees of freedom
#> AIC: NA
#> 
#> Number of Fisher Scoring iterations: 5
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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        pvs = pvs, # Named list
        method = "ICILS") # the name of the method aka the study name
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 153 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.1489     0.1124   1.325    0.186
#> Math1        -2.0623     0.4660  -4.425    0.000
#> Reading1      1.6379     0.6107   2.682    0.007
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model for first plausible value combination:
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.10804    0.07994   1.352    0.177    
#> Math1       -1.68273    0.11580 -14.531   <2e-16 ***
#> Reading1     1.06341    0.09644  11.027   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.487731)
#> 
#>     Null deviance: 2080.3  on 999  degrees of freedom
#> Residual deviance: 1426.5  on 997  degrees of freedom
#> AIC: NA
#> 
#> Number of Fisher Scoring iterations: 5
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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        pvs = pvs, # Named list
        relatedpvs = FALSE, # Unrelated PVs
        method = "ICILS") # the name of the method aka the study name
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 459 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Replicate weights' model:
#>  
#> Call: 
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, relatedpvs = FALSE, df = repdata2, wt = "wt", 
#>     repwt = RW, method = "ICILS")
#> 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   0.1464     0.1016   1.442     0.15
#> Math1        -2.1058     0.3378  -6.234     0.00
#> Reading1      1.6859     0.4452   3.787     0.00
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model for first plausible value combination:
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, relatedpvs = FALSE, df = repdata2, wt = "wt", 
#>     repwt = RW, method = "ICILS")
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.10804    0.07994   1.352    0.177    
#> Math1       -1.68273    0.11580 -14.531   <2e-16 ***
#> Reading1     1.06341    0.09644  11.027   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.487731)
#> 
#>     Null deviance: 2080.3  on 999  degrees of freedom
#> Residual deviance: 1426.5  on 997  degrees of freedom
#> AIC: NA
#> 
#> Number of Fisher Scoring iterations: 5
#> 


### Groups ----

# Simple regression - default family = gaussian
repglm(formula = GENDER ~ 1 + Math1,
       family = gaussian, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 1 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 2 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 3 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.5080     0.0146  34.858
#> Math1        -0.2954     0.0115 -25.780
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = gaussian, df = repdata2, 
#>     wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  Math1        
#> GR1 0.363 (0.03) -0.285 (0.02)
#> GR3 0.658 (0.03) -0.292 (0.02)
#> GR2 0.503 (0.02) -0.309 (0.03)


# Simple regression - change link function
repglm(formula = GENDER ~ 1 + Math1,
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 1 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 2 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 3 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.0439     0.0906   0.484
#> Math1        -1.7653     0.1224 -14.424
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1        
#> GR1 -0.819 (0.15) -1.734 (0.19)
#> GR3  0.925 (0.18) -1.719 (0.18)
#> GR2  0.025 (0.14) -1.843 (0.26)

# Multiple regression
repglm(formula = GENDER ~ 1 + Math1 + Reading1,
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 1 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 2 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 3 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.0996     0.1135   0.878
#> Math1        -1.7881     0.1312 -13.627
#> Reading1      0.9595     0.1185   8.097
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1 + Reading1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1 -0.137 (0.20) -1.832 (0.24) 1.021 (0.21)
#> GR3  0.361 (0.23) -1.739 (0.18) 0.954 (0.20)
#> GR2  0.075 (0.15) -1.793 (0.26) 0.903 (0.21)

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3))
pvs
#> $Math
#> [1] "Math1" "Math2" "Math3"
#> 

repglm(formula = GENDER ~ 1 + Math + Reading1, # Math1 now is "Math"
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 1 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 2 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 3 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.1334     0.1366   0.977
#> Math1        -1.8618     0.1619 -11.500
#> Reading1      1.1163     0.1548   7.211
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading1, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, group = "GROUP", 
#>     method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  Math1         Reading1    
#> GR1 0.014 (0.30) -1.733 (0.26) 1.126 (0.23)
#> GR3 0.348 (0.23) -1.847 (0.22) 1.089 (0.25)
#> GR2 0.039 (0.16) -2.006 (0.35) 1.134 (0.32)

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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 1 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 2 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 3 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.1632     0.2183   0.748
#> Math1        -2.0752     0.2580  -8.043
#> Reading1      1.6648     0.4650   3.581
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, group = "GROUP", 
#>     method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  Math1         Reading1    
#> GR1 0.241 (0.54) -2.001 (0.33) 1.631 (0.72)
#> GR3 0.147 (0.33) -2.000 (0.44) 1.574 (0.72)
#> GR2 0.102 (0.17) -2.225 (0.54) 1.790 (0.95)

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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       relatedpvs = FALSE, # Unrelated PVs
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name
#> Estimating pooled model.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 1 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 2 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Estimating group 3 of 3.
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> Warning: observations with zero weight not used for calculating dispersion
#> 1836 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.1554     0.1665   0.933
#> Math1        -2.1117     0.2188  -9.649
#> Reading1      1.7107     0.3408   5.020
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, relatedpvs = FALSE, df = repdata2, wt = "wt", 
#>     repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  Math1         Reading1    
#> GR1 0.223 (0.37) -2.024 (0.37) 1.655 (0.55)
#> GR3 0.140 (0.29) -2.033 (0.34) 1.631 (0.53)
#> GR2 0.102 (0.17) -2.277 (0.42) 1.846 (0.68)
```
