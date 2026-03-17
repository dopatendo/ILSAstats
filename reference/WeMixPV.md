# Survey Weighted Mixed-Effects Models with Plausible Values

Fits a linear mixed-effects model using
[mix](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)
and plausible values.

## Usage

``` r
WeMixPV(formula, data = NULL, weights = NULL, pvs, relatedpvs = TRUE, ...)
```

## Arguments

- formula:

  a formula object in the style of `lme4` that creates the model.

- data:

  a data frame containing the raw data for the model.

- weights:

  a character vector of names of weight variables found in the data
  frame starts with units (level 1) and increasing (larger groups).

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

- ...:

  Arguments passed on to
  [`WeMix::mix`](https://american-institutes-for-research.github.io/WeMix/reference/mix.html)

  `cWeights`

  :   logical, set to `TRUE` to use conditional weights. Otherwise,
      `mix` expects unconditional weights.

  `center_group`

  :   a list where the name of each element is the name of the
      aggregation level, and the element is a formula of variable names
      to be group mean centered; for example to group mean center gender
      and age within the group student: `list("student"= ~gender+age)`,
      default value of NULL does not perform any group mean centering.

  `center_grand`

  :   a formula of variable names to be grand mean centered, for example
      to center the variable education by overall mean of education:
      `~education`. Default is NULL which does no centering.

  `max_iteration`

  :   a optional integer, for non-linear models fit by adaptive
      quadrature which limits number of iterations allowed before
      quitting. Defaults to 10. This is used because if the likelihood
      surface is flat, models may run for a very long time without
      converging.

  `nQuad`

  :   an optional integer number of quadrature points to evaluate models
      solved by adaptive quadrature. Only non-linear models are
      evaluated with adaptive quadrature. See notes for additional
      guidelines.

  `run`

  :   logical; `TRUE` runs the model while `FALSE` provides partial
      output for debugging or testing. Only applies to non-linear models
      evaluated by adaptive quadrature.

  `verbose`

  :   logical, default `FALSE`; set to `TRUE` to print results of
      intermediate steps of adaptive quadrature. Only applies to
      non-linear models.

  `acc0`

  :   deprecated; ignored.

  `keepAdapting`

  :   logical, set to `TRUE` when the adaptive quadrature should adapt
      after every Newton step. Defaults to `FALSE`. `FALSE` should be
      used for faster (but less accurate) results. Only applies to
      non-linear models.

  `start`

  :   optional numeric vector representing the point at which the model
      should start optimization; takes the shape of c(coef, vars) from
      results (see help).

  `fast`

  :   logical; deprecated

  `family`

  :   the family; optionally used to specify generalized linear mixed
      models. Currently only
      [`binomial()`](https://rdrr.io/r/stats/family.html) and
      [`poisson()`](https://rdrr.io/r/stats/family.html) are supported.

## Value

a list.

## Examples

``` r
# Prepare data weights
repdata2 <- repdata
repdata2$wt1 <- repdata2$wt # weight level 1
repdata2$wt2 <- 1 # weight level 2


# Null model - with PVs
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m1 <- WeMixPV(formula = MATH ~ 1 + (1|GROUP), # Intercept varies across GROUP
             pvs = pvs, # Named list
             data = repdata2, # Data frame
             weights = c("wt1","wt2")) # Weights vector
m1
#> 5 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Multilevel results with PVs:
#>  
#> Call: 
#> WeMixPV(formula = MATH ~ 1 + (1 | GROUP), data = repdata2, weights = c("wt1", 
#>     "wt2"), pvs = pvs)
#> 
#> Random effects:
#>                    Estimate Std. Error   t value          df     Pr(>|t|)
#> GROUP.(Intercept) 0.2332214 0.11973218  1.947859 2282.198284 5.155428e-02
#> Residual          0.8132685 0.02299559 35.366280    7.250355 2.186994e-09
#> 
#> Fixed effects:
#>                Estimate Std. Error   t value        df  Pr(>|t|)
#> (Intercept) 0.004824544  0.3417468 0.0141173 181384504 0.9887364
#> --------------------------------------------------------------------------------
#> Estimated models:
#> Math1 ~ 1 + (1 | GROUP)
#> Math2 ~ 1 + (1 | GROUP)
#> Math3 ~ 1 + (1 | GROUP)
#> Math4 ~ 1 + (1 | GROUP)
#> Math5 ~ 1 + (1 | GROUP)

## Fixed effects
m1$fixef
#>                Estimate Std. Error   t value        df  Pr(>|t|)
#> (Intercept) 0.004824544  0.3417468 0.0141173 181384504 0.9887364

## Random effects
m1$ranef
#>                    Estimate Std. Error   t value          df     Pr(>|t|)
#> GROUP.(Intercept) 0.2332214 0.11973218  1.947859 2282.198284 5.155428e-02
#> Residual          0.8132685 0.02299559 35.366280    7.250355 2.186994e-09

## Models for each PV
summary(m1$models)
#>       Length Class        Mode
#> Math1 22     WeMixResults list
#> Math2 22     WeMixResults list
#> Math3 22     WeMixResults list
#> Math4 22     WeMixResults list
#> Math5 22     WeMixResults list

# Multiple regression
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m2 <- WeMixPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata2, # Data frame
             weights = c("wt1","wt2")) # Weights vector
m2
#> 5 models were estimated.
#> -------------------------------------------------------------------------------- 
#> Multilevel results with PVs:
#>  
#> Call: 
#> WeMixPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1 | 
#>     GROUP), data = repdata2, weights = c("wt1", "wt2"), pvs = pvs)
#> 
#> Random effects:
#>                       Estimate   Std. Error    t value         df     Pr(>|t|)
#> GROUP.(Intercept) 1.666843e-05 8.968157e-05  0.1858623 413.549345 8.526439e-01
#> Residual          5.591163e-01 3.337543e-02 16.7523312   5.787741 3.995963e-06
#> 
#> Fixed effects:
#>                Estimate Std. Error    t value       df     Pr(>|t|)
#> (Intercept) -64.9984631 3.52733990 -18.427048 4.705471 1.440966e-05
#> GENDER       -0.8991677 0.08072093 -11.139215 6.108438 2.761341e-05
#> SES           0.2458968 0.06537683   3.761223 4.814884 1.408541e-02
#> schoolSES     1.0641649 0.11767005   9.043634 4.211655 6.498273e-04
#> --------------------------------------------------------------------------------
#> Estimated models:
#> Math1 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math2 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math3 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math4 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)
#> Math5 ~ 1 + GENDER + SES + schoolSES + (1 | GROUP)



```
