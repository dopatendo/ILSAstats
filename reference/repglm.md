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
repdata2 <- repdata[1:500,]

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
#> (Intercept)   0.5000     0.0172  29.129        0
#> Math1        -0.2345     0.0145 -16.184        0
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
#> (Intercept)  0.50000    0.01942   25.74   <2e-16 ***
#> Math1       -0.23452    0.01822  -12.87   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.2834801)
#> 
#>     Null deviance: 188.16  on 499  degrees of freedom
#> Residual deviance: 141.17  on 498  degrees of freedom
#> AIC: 597.4
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
#> (Intercept)   0.0046     0.0922   0.050    0.961
#> Math1        -1.3091     0.1396  -9.378    0.000
#> 
#> -------------------------------------------------------------------------------- 
#> Total weights' model:
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, method = "ICILS")
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.004565   0.103681   0.044    0.965    
#> Math1       -1.309150   0.133845  -9.781   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.493772)
#> 
#>     Null deviance: 1043.53  on 499  degrees of freedom
#> Residual deviance:  823.07  on 498  degrees of freedom
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
#> (Intercept)   0.0338     0.1141   0.296    0.767
#> Math1        -1.6906     0.1949  -8.673    0.000
#> Reading1      1.0094     0.1528   6.605    0.000
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
#> (Intercept)  0.03381    0.11595   0.292    0.771    
#> Math1       -1.69060    0.16592 -10.189  < 2e-16 ***
#> Reading1     1.00937    0.13871   7.277 1.34e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.583703)
#> 
#>     Null deviance: 1043.5  on 499  degrees of freedom
#> Residual deviance:  716.8  on 497  degrees of freedom
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
#> (Intercept)   0.0384     0.1151   0.333    0.739
#> Math1        -1.7658     0.1969  -8.967    0.000
#> Reading1      1.0935     0.1854   5.896    0.000
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
#> (Intercept)  0.03381    0.11595   0.292    0.771    
#> Math1       -1.69060    0.16592 -10.189  < 2e-16 ***
#> Reading1     1.00937    0.13871   7.277 1.34e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.583703)
#> 
#>     Null deviance: 1043.5  on 499  degrees of freedom
#> Residual deviance:  716.8  on 497  degrees of freedom
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
#> (Intercept)   0.0848     0.1398   0.606    0.545
#> Math1        -2.1205     0.5870  -3.612    0.000
#> Reading1      1.6681     0.7211   2.313    0.021
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
#> (Intercept)  0.03381    0.11595   0.292    0.771    
#> Math1       -1.69060    0.16592 -10.189  < 2e-16 ***
#> Reading1     1.00937    0.13871   7.277 1.34e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.583703)
#> 
#>     Null deviance: 1043.5  on 499  degrees of freedom
#> Residual deviance:  716.8  on 497  degrees of freedom
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
#> (Intercept)   0.0888     0.1361   0.653    0.514
#> Math1        -2.1552     0.4362  -4.941    0.000
#> Reading1      1.7197     0.5725   3.004    0.003
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
#> (Intercept)  0.03381    0.11595   0.292    0.771    
#> Math1       -1.69060    0.16592 -10.189  < 2e-16 ***
#> Reading1     1.00937    0.13871   7.277 1.34e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for quasibinomial family taken to be 1.583703)
#> 
#>     Null deviance: 1043.5  on 499  degrees of freedom
#> Residual deviance:  716.8  on 497  degrees of freedom
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
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.5023     0.0209  23.993
#> Math1        -0.2783     0.0180 -15.427
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = gaussian, df = repdata2, 
#>     wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)  Math1        
#> GR1 0.410 (0.04) -0.257 (0.02)
#> GR3 0.633 (0.04) -0.281 (0.03)
#> GR2 0.464 (0.03) -0.297 (0.04)


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
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)  -0.0072     0.1214  -0.059
#> Math1        -1.6468     0.1959  -8.406
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1        
#> GR1 -0.558 (0.20) -1.542 (0.26)
#> GR3  0.729 (0.25) -1.594 (0.27)
#> GR2 -0.192 (0.17) -1.804 (0.46)

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
#> 204 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)  -0.0390     0.1469  -0.266
#> Math1        -1.6917     0.2041  -8.289
#> Reading1      1.0154     0.1666   6.095
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math1 + Reading1, family = quasibinomial, 
#>     df = repdata2, wt = "wt", repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.071 (0.27) -1.575 (0.27) 0.828 (0.28)
#> GR3  0.034 (0.30) -1.705 (0.29) 1.244 (0.30)
#> GR2 -0.222 (0.18) -1.795 (0.47) 0.974 (0.28)

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
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)  -0.0181     0.1715  -0.105
#> Math1        -1.7842     0.2329  -7.662
#> Reading1      1.1512     0.2110   5.455
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading1, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, group = "GROUP", 
#>     method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.173 (0.32) -1.541 (0.29) 0.877 (0.29)
#> GR3  0.103 (0.33) -1.866 (0.32) 1.320 (0.34)
#> GR2 -0.330 (0.23) -1.945 (0.55) 1.256 (0.45)

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
#> 612 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.0459     0.2183   0.210
#> Math1        -2.1022     0.4260  -4.935
#> Reading1      1.8638     0.6283   2.966
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, df = repdata2, wt = "wt", repwt = RW, group = "GROUP", 
#>     method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.403 (0.45) -1.822 (0.59) 1.436 (0.74)
#> GR3 -0.009 (0.39) -2.240 (0.90) 1.931 (1.04)
#> GR2 -0.256 (0.26) -2.244 (0.69) 2.224 (1.38)

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
#> 1836 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Composite coefficients:
#>             Estimate Std. Error t value
#> (Intercept)   0.0485     0.2006   0.242
#> Math1        -2.1172     0.3430  -6.172
#> Reading1      1.9053     0.4872   3.911
#> -------------------------------------------------------------------------------- 
#> 
#> Call:
#> repglm(formula = GENDER ~ 1 + Math + Reading, family = quasibinomial, 
#>     pvs = pvs, relatedpvs = FALSE, df = repdata2, wt = "wt", 
#>     repwt = RW, group = "GROUP", method = "ICILS")
#> 
#> Coefficients and standard errors by group:
#>     (Intercept)   Math1         Reading1    
#> GR1  0.410 (0.40) -1.889 (0.55) 1.528 (0.73)
#> GR3 -0.006 (0.37) -2.193 (0.57) 1.932 (0.71)
#> GR2 -0.258 (0.26) -2.270 (0.65) 2.256 (1.05)
```
