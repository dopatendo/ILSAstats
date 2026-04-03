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
repdata2 <- repdata[1:200,]

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
#> (Intercept)   0.4690     0.0314  14.956        0
#> Math1        -0.2325     0.0211 -11.036        0
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
#> (Intercept)  0.46905    0.03045  15.404  < 2e-16 ***
#> Math1       -0.23246    0.02760  -8.423 7.42e-15 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.2776452)
#> 
#>     Null deviance: 74.670  on 199  degrees of freedom
#> Residual deviance: 54.974  on 198  degrees of freedom
#> AIC: 238.36
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
#> (Intercept)  -0.1602     0.1749  -0.916    0.361
#> Math1        -1.3159     0.2167  -6.074    0.000
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
#> (Intercept)  -0.1602     0.1612  -0.994    0.322    
#> Math1        -1.3159     0.2030  -6.483 6.99e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.400462)
#> 
#>     Null deviance: 414.39  on 199  degrees of freedom
#> Residual deviance: 321.44  on 198  degrees of freedom
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
#> (Intercept)  -0.0600     0.1977  -0.303    0.762
#> Math1        -1.7319     0.2872  -6.031    0.000
#> Reading1      1.1746     0.2710   4.335    0.000
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
#> (Intercept) -0.05998    0.18856  -0.318    0.751    
#> Math1       -1.73192    0.26209  -6.608 3.56e-10 ***
#> Reading1     1.17458    0.24219   4.850 2.50e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.552154)
#> 
#>     Null deviance: 414.39  on 199  degrees of freedom
#> Residual deviance: 272.85  on 197  degrees of freedom
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
#> (Intercept)  -0.0343     0.2097  -0.163     0.87
#> Math1        -1.6702     0.2840  -5.882     0.00
#> Reading1      1.1384     0.2653   4.290     0.00
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
#> (Intercept) -0.05998    0.18856  -0.318    0.751    
#> Math1       -1.73192    0.26209  -6.608 3.56e-10 ***
#> Reading1     1.17458    0.24219   4.850 2.50e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.552154)
#> 
#>     Null deviance: 414.39  on 199  degrees of freedom
#> Residual deviance: 272.85  on 197  degrees of freedom
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
#> (Intercept)   0.1186     0.3635   0.326    0.745
#> Math1        -2.1348     0.6569  -3.250    0.001
#> Reading1      1.7444     0.7259   2.403    0.017
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
#> (Intercept) -0.05998    0.18856  -0.318    0.751    
#> Math1       -1.73192    0.26209  -6.608 3.56e-10 ***
#> Reading1     1.17458    0.24219   4.850 2.50e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.552154)
#> 
#>     Null deviance: 414.39  on 199  degrees of freedom
#> Residual deviance: 272.85  on 197  degrees of freedom
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
#> (Intercept)   0.1223     0.2917   0.419    0.675
#> Math1        -2.1572     0.5851  -3.687    0.000
#> Reading1      1.7873     0.6706   2.665    0.008
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
#> (Intercept) -0.05998    0.18856  -0.318    0.751    
#> Math1       -1.73192    0.26209  -6.608 3.56e-10 ***
#> Reading1     1.17458    0.24219   4.850 2.50e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.552154)
#> 
#>     Null deviance: 414.39  on 199  degrees of freedom
#> Residual deviance: 272.85  on 197  degrees of freedom
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
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.4776     0.0326  14.654
#> Math1        -0.2834     0.0230 -12.304
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = gaussian, df = repdata2, 
#>     wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  Math1        
#> GR1 0.354 (0.05) -0.275 (0.03)
#> GR3 0.652 (0.07) -0.321 (0.04)
#> GR2 0.427 (0.05) -0.255 (0.05)


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
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)  -0.1015     0.2125  -0.478
#> Math1        -1.7226     0.2717  -6.340
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1        
#> GR1 -0.832 (0.29) -1.618 (0.39)
#> GR3  0.893 (0.47) -2.080 (0.56)
#> GR2 -0.365 (0.31) -1.470 (0.45)

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
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)  -0.0509     0.2563  -0.199
#> Math1        -1.7695     0.2989  -5.920
#> Reading1      1.2161     0.3101   3.922
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1 + Reading1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.114 (0.37) -1.922 (0.50) 1.406 (0.46)
#> GR3  0.113 (0.59) -2.074 (0.57) 1.462 (0.60)
#> GR2 -0.379 (0.32) -1.313 (0.47) 0.780 (0.54)

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
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)  -0.0137     0.2776  -0.049
#> Math1        -1.6748     0.3647  -4.592
#> Reading1      1.2346     0.3224   3.829
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading1, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, group = "GROUP", 
#>     method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.309 (0.47) -1.776 (0.64) 1.416 (0.56)
#> GR3  0.021 (0.59) -1.812 (0.67) 1.352 (0.58)
#> GR2 -0.371 (0.36) -1.436 (0.57) 0.936 (0.54)

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
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.0906     0.4678   0.194
#> Math1        -2.1550     0.9821  -2.194
#> Reading1      2.1304     1.5459   1.378
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, group = "GROUP", 
#>     method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.596 (0.80) -2.126 (1.05) 1.872 (1.01)
#> GR3  0.006 (0.59) -2.278 (0.81) 2.001 (1.13)
#> GR2 -0.330 (0.98) -2.061 (2.63) 2.518 (4.38)

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
#> 1836 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.1112     0.3706   0.300
#> Math1        -2.1381     0.8096  -2.641
#> Reading1      2.1387     1.1308   1.891
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, relatedpvs = FALSE, df = repdata2, wt = "wt", 
#>     repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.589 (0.63) -2.121 (0.89) 1.914 (0.92)
#> GR3  0.034 (0.63) -2.334 (1.19) 2.034 (1.05)
#> GR2 -0.289 (0.66) -1.960 (1.92) 2.468 (3.09)
```
