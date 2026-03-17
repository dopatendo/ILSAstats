# Intraclass Correlation Coefficient

Calculates the intraclass correlation coefficient (ICC) fitting a linear
mixed-effects model using
[lmer](https://rdrr.io/pkg/lme4/man/lmer.html).

## Usage

``` r
icc(x, PV = FALSE, group, data, weights = NULL, ...)
```

## Arguments

- x:

  a string vector specifying variable names (within `data`).

- PV:

  a logical value indicating if the variables in `x` are plausible
  values.

- group:

  a string specifying the variable name (within `data`) to be used for
  grouping.

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

- ...:

  Arguments passed on to
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html)

  `formula`

  :   a two-sided linear formula object describing both the
      fixed-effects and random-effects part of the model, with the
      response on the left of a `~` operator and the terms, separated by
      `+` operators, on the right. Random-effects terms are
      distinguished by vertical bars (`|`) separating expressions for
      design matrices from grouping factors. Two vertical bars (`||`)
      can be used to specify multiple uncorrelated random effects for
      the same grouping variable. (Because of the way it is implemented,
      the `||`-syntax *works only for design matrices containing numeric
      (continuous) predictors*; to fit models with independent
      categorical effects, see
      [`dummy`](https://rdrr.io/pkg/lme4/man/dummy.html) or the
      `lmer_alt` function from the
      [afex](https://CRAN.R-project.org/package=afex) package.)

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

  :   a named [`list`](https://rdrr.io/r/base/list.html) of starting
      values for the parameters in the model. For `lmer` this can be a
      numeric vector or a list with one component named `"theta"`.

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

a numeric value or a list.

## Examples

``` r
# ICC of one variable
icc(x = "Math1",group = "GROUP", weights = repdata$wt, data = repdata)
#> $Math1
#> [1] 0.2277918
#> 


# ICC of more than one variable
icc(x = c("Math1","Math2","Math3","Math4","Math5","SES"),
    group = "GROUP", weights = repdata$wt, data = repdata)
#> $Math1
#> [1] 0.2277918
#> 
#> $Math2
#> [1] 0.211539
#> 
#> $Math3
#> [1] 0.252869
#> 
#> $Math4
#> [1] 0.1997381
#> 
#> $Math5
#> [1] 0.2229097
#> 
#> $SES
#> [1] 0.1334063
#> 


# ICC of PVs
icc(x = c("Math1","Math2","Math3","Math4","Math5"), PV = TRUE,
    group = "GROUP", weights = repdata$wt, data = repdata)
#> $Average
#> [1] 0.2229695
#> 
#> $Math1
#> [1] 0.2277918
#> 
#> $Math2
#> [1] 0.211539
#> 
#> $Math3
#> [1] 0.252869
#> 
#> $Math4
#> [1] 0.1997381
#> 
#> $Math5
#> [1] 0.2229097
#> 
```
