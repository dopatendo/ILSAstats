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
  separating expressions for design matrices from grouping factors. Two
  vertical bars (`||`) can be used to specify multiple uncorrelated
  random effects for the same grouping variable. (Because of the way it
  is implemented, the `||`-syntax *works only for design matrices
  containing numeric (continuous) predictors*; to fit models with
  independent categorical effects, see
  [`dummy`](https://rdrr.io/pkg/lme4/man/dummy.html) or the `lmer_alt`
  function from the [afex](https://CRAN.R-project.org/package=afex)
  package.)

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

a list.
