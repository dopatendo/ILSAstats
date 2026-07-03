# Linear Mixed-Models with Plausible Values

Fits a linear mixed-effects model using
[lmer](https://rdrr.io/pkg/lme4/man/lmer.html) and plausible values.

## Usage

``` r
lmerPV(
  formula,
  data = NULL,
  weights = NULL,
  pvs,
  relatedpvs = TRUE,
  grandmean = NULL,
  groupmean = NULL,
  group = NULL,
  nullmodel = FALSE,
  barnardrubin = TRUE,
  ...
)
```

## Arguments

- formula:

  a two-sided linear formula object describing both the fixed-effects
  and random-effects part of the model, with the response on the left of
  a `~` operator and the terms, separated by `+` operators, on the
  right. Random-effects terms are distinguished by vertical bars (`|`)
  separating expressions for design matrices from grouping factors. By
  default, non-scalar random effects (where the design matrix has more
  than one column, e.g. `(1+x|f)`) are fitted with *unstructured*
  (general positive semidefinite) covariance matrices.

  - Two vertical bars (`||`) can be used to specify multiple
    uncorrelated random effects for the same grouping variable. With
    default settings, the `||`-syntax *works only for design matrices
    containing numeric (continuous) predictors*; to fit models with
    independent categorical effects, use `diag(f|g)` or set
    `options(lme4.doublevert.default = "diag_special")` (see
    [`getDoublevertDefault`](https://rdrr.io/pkg/lme4/man/getDoublevertDefault.html)).

  - Tags preceding a random effect term specify covariance structure:

    - `us` (default: `us(f|g)` is equivalent to `(f|g)`): unstructured,
      positive semi-definite

    - `diag`: diagonal (all correlations set to zero). Specify
      `diag(f|g, hom = TRUE)` to fit a homogeneous diagonal covariance
      matrix

    - `cs`: compound symmetric (all pairwise correlations set
      identical). Specify `cs(f|g, hom = TRUE)` for homogeneous
      variances.

    - `ar1`: autoregressive order 1. Note that AR1 models are
      homogeneous by default; specify `ar1(f|g, hom = FALSE)` for
      heterogeneous variances.

- data:

  an optional data frame containing the variables named in `formula`. By
  default the variables are taken from the environment from which `lmer`
  is called. While `data` is optional, the package authors *strongly*
  recommend its use, especially when later applying methods such as
  `update` and `drop1` to the fitted model (*such methods are not
  guaranteed to work properly if `data` is omitted*). If `data` is
  omitted, variables will be taken from the environment of `formula` (if
  specified as a formula) or from the parent frame (if specified as a
  character vector).

- weights:

  an optional vector of ‘prior weights’ to be used in the fitting
  process. Should be `NULL` or a numeric vector. Prior `weights` are
  *not* normalized or standardized in any way. In particular, the
  diagonal of the residual covariance matrix is the squared residual
  standard deviation parameter
  [`sigma`](https://rdrr.io/pkg/lme4/man/sigma.html) times the vector of
  inverse `weights`. Therefore, if the `weights` have relatively large
  magnitudes, then in order to compensate, the
  [`sigma`](https://rdrr.io/pkg/lme4/man/sigma.html) parameter will also
  need to have a relatively large magnitude.

- pvs:

  a list indicating which variables from `formula` should be replaced by
  which plausible values variables. For more details check the examples.

- relatedpvs:

  a logical value indicating if `pvs` are drawn from the same model, and
  have the same number of plausible values. If `TRUE` (default), a total
  of \\n\\ estimations will be done, where \\n\\ is the number of
  plausible values for each plausible value variable. If `FALSE`, a
  total of \\n_1 \times n_2 \times n\_...\\ estimations will be done,
  where \\n_i\\ is the number of plausible values in each plausible
  value variable.

- grandmean:

  a character vector indicating the names of columns of `data` to which
  grand-mean should be applied.

- groupmean:

  a character vector indicating the names of columns of `data` to which
  group-mean should be applied.

- group:

  a string specifying the variable name (within `df`) to be used for
  grouping. Categories in `group` are treated as independent, e.g.,
  countries.

- nullmodel:

  a logical value indicating if the null model should also be estimated.

- barnardrubin:

  a logical value indicating if Barnard & Rubin adjustment should be
  used for estimating the degrees of freedom. Default is `TRUE`.

- ...:

  Arguments passed on to
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html)

  `REML`

  :   logical scalar - Should the estimates be chosen to optimize the
      REML criterion (as opposed to the log-likelihood)?

  `control`

  :   a list (of correct class, resulting from
      [`lmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html)
      or
      [`glmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html)
      respectively) containing control parameters, including the
      nonlinear optimizer to be used and parameters to be passed through
      to the nonlinear optimizer, see the `*lmerControl` documentation
      for details.

  `start`

  :   a numeric vector or a named
      [list](https://rdrr.io/r/base/list.html) with one optional
      component named `par` or `theta`, giving starting values for
      covariance parameters. Numeric `start` is equivalent to
      `list(par = start)`. Parameters corresponding to unstructured
      covariance matrices are on the scale of the Cholesky factor of the
      relative covariance matrix. By default, all relative covariance
      matrices are identity matrices.

  `verbose`

  :   integer scalar. If `> 0` verbose output is generated during the
      optimization of the parameter estimates. If `> 1` verbose output
      is generated during the individual penalized iteratively
      reweighted least squares (PIRLS) steps.

  `subset`

  :   an optional expression indicating the subset of the rows of `data`
      that should be used in the fit. This can be a logical vector, or a
      numeric vector indicating which observation numbers are to be
      included, or a character vector of the row names to be included.
      All observations are included by default.

  `na.action`

  :   a function that indicates what should happen when the data contain
      `NA`s. The default action (`na.omit`, inherited from the 'factory
      fresh' value of `getOption("na.action")`) strips any observations
      with any missing values in any variables.

  `offset`

  :   this can be used to specify an *a priori* known component to be
      included in the linear predictor during fitting. This should be
      `NULL` or a numeric vector of length equal to the number of cases.
      One or more [`offset`](https://rdrr.io/r/stats/offset.html) terms
      can be included in the formula instead or as well, and if more
      than one is specified their sum is used. See
      [`model.offset`](https://rdrr.io/r/stats/model.extract.html).

  `contrasts`

  :   an optional list. See the `contrasts.arg` of
      `model.matrix.default`.

  `devFunOnly`

  :   logical - return only the deviance evaluation function. Note that
      because the deviance function operates on variables stored in its
      environment, it may not return *exactly* the same values on
      subsequent calls (but the results should always be within machine
      tolerance).

## Value

a list.

## Examples

``` r
# Null model - with PVs
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m1 <- lmerPV(formula = MATH ~ 1 + (1|GROUP), # Intercept varies across GROUP
      pvs = pvs, # Named list
      data = repdata, # Data frame
      weights = repdata$wt) # Weights vector
m1
#> 5 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Multilevel results with PVs:
#>  
#> Call: 
#> lmerPV(formula = MATH ~ 1 + (1 | GROUP), data = repdata, weights = repdata$wt, 
#>     pvs = pvs)
#> 
#> Random effects:
#>                    Variance
#> GROUP.(Intercept) 0.3498329
#> Residual          1.2185299
#> 
#> Fixed effects:
#>                Estimate Std. Error   t value       df Pr(>|t|)
#> (Intercept) 0.004824544   0.341747 0.0141173 4982.124 0.988737
#> --------------------------------------------------------------------------------
#> Estimated models:
#> Math1 ~ 1 + (1 | GROUP)
#> Math2 ~ 1 + (1 | GROUP)
#> Math3 ~ 1 + (1 | GROUP)
#> Math4 ~ 1 + (1 | GROUP)
#> Math5 ~ 1 + (1 | GROUP)

## Fixed effects
m1$fixef
#>                Estimate Std. Error   t value       df Pr(>|t|)
#> (Intercept) 0.004824544   0.341747 0.0141173 4982.124 0.988737

## Random effects
m1$ranef
#>                    Variance
#> GROUP.(Intercept) 0.3498329
#> Residual          1.2185299

## Models for each PV
summary(m1$models)
#>       Length Class   Mode
#> Math1 1      lmerMod S4  
#> Math2 1      lmerMod S4  
#> Math3 1      lmerMod S4  
#> Math4 1      lmerMod S4  
#> Math5 1      lmerMod S4  

# Multiple regression
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m2 <- lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata, # Data frame
             weights = repdata$wt) # Weights vector
m2
#> 5 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Multilevel results with PVs:
#>  
#> Call: 
#> lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1 | GROUP), 
#>     data = repdata, weights = repdata$wt, pvs = pvs)
#> 
#> Random effects:
#>                       Variance
#> GROUP.(Intercept) 0.0002756404
#> Residual          0.8381211936
#> 
#> Fixed effects:
#>                Estimate Std. Error    t value       df     Pr(>|t|)
#> (Intercept) -64.9990035 3.90253673 -16.655578 7.010247 6.773595e-07
#> GENDER       -0.8990175 0.07564016 -11.885452 4.651592 1.173057e-04
#> SES           0.2458926 0.06342425   3.876950 4.152052 1.666901e-02
#> schoolSES     1.0641781 0.12298663   8.652795 4.979792 3.476117e-04
#> --------------------------------------------------------------------------------
#> Estimated models:
#> Math1 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math2 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math3 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math4 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math5 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)

# Multiple regression with grandmean centering
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m3 <- lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata, # Data frame
             weights = repdata$wt,
             grandmean = c("SES","schoolSES"))
m3
#> 5 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Multilevel results with PVs:
#>  
#> Call: 
#> lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1 | GROUP), 
#>     data = repdata, weights = repdata$wt, pvs = pvs, grandmean = c("SES", 
#>         "schoolSES"))
#> 
#> Random effects:
#>                       Variance
#> GROUP.(Intercept) 0.0002756404
#> Residual          0.8381211936
#> 
#> Fixed effects:
#>               Estimate Std. Error    t value       df     Pr(>|t|)
#> (Intercept)  0.4589106 0.03739158  12.273100 6.660211 8.094989e-06
#> GENDER      -0.8990175 0.07564016 -11.885452 4.651592 1.173057e-04
#> SES          0.2458926 0.06342425   3.876950 4.152052 1.666901e-02
#> schoolSES    1.0641781 0.12298663   8.652795 4.979792 3.476117e-04
#> --------------------------------------------------------------------------------
#> Estimated models:
#> Math1 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math2 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math3 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math4 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math5 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)

# Multiple regression with groupmean centering
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m4 <- lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata, # Data frame
             weights = repdata$wt,
             grandmean = "schoolSES",
             groupmean = "SES",
             group = repdata$GROUP)
m4
#> 5 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Multilevel results with PVs:
#>  
#> Call: 
#> lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1 | GROUP), 
#>     data = repdata, weights = repdata$wt, pvs = pvs, grandmean = "schoolSES", 
#>     groupmean = "SES", group = repdata$GROUP)
#> 
#> Random effects:
#>                     Variance
#> GROUP.(Intercept) 0.00028522
#> Residual          0.83812119
#> 
#> Fixed effects:
#>               Estimate Std. Error    t value       df     Pr(>|t|)
#> (Intercept)  0.4589091 0.03743452  12.258983 6.690742 7.869635e-06
#> GENDER      -0.8990154 0.07564056 -11.885362 4.651586 1.173108e-04
#> SES          0.2458884 0.06342472   3.876855 4.152046 1.667043e-02
#> schoolSES    1.3102965 0.07847520  16.696950 7.038119 6.399173e-07
#> --------------------------------------------------------------------------------
#> Estimated models:
#> Math1 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math2 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math3 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math4 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math5 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)


```
