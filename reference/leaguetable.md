# ILSA's league tables

Estimates the mean score for all countries within a cycle of an ILSA.
Arguments `method`, `reps`, and `var`, are extracted from
[`autoILSA`](https://dopatendo.github.io/ILSAstats/reference/autoILSA.md)
and can be overridden by the user.

## Usage

``` r
leaguetable(
  df,
  study = NULL,
  year,
  subject = NULL,
  specification = NULL,
  addCI = TRUE,
  alpha = 0.05,
  method = NULL,
  reps = NULL,
  fixN = TRUE
)
```

## Arguments

- df:

  a data frame.

- study:

  an optional character vector indicating the ILSA name, for a list of
  available ILSA, check
  [`autoILSA`](https://dopatendo.github.io/ILSAstats/reference/autoILSA.md).
  If `NULL`, the ILSA name will be determined by the column names in the
  data frame.

- year:

  a numeric vector indicating the ILSA name, for a list of available
  cycles, check
  [`autoILSA`](https://dopatendo.github.io/ILSAstats/reference/autoILSA.md).

- subject:

  an optional character vector indicating the subject for a list of
  available ILSA, check
  [`autoILSA`](https://dopatendo.github.io/ILSAstats/reference/autoILSA.md).

- specification:

  a character value indicating extra specification like grade (e.g.,
  `"G8"` for TIMSS) or subject (e.g., `"Math"` for TIMSSADVANCED).

- addCI:

  a logical value indicating if confidence intervals should be added.
  Defaults is `TRUE`.

- alpha:

  a numeric value indicating confidence level.

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

- reps:

  an integer indicating the number of replications to be created. If
  `NULL` the maximum number of zones will be used.

- fixN:

  a logical value indicating if data should be "fixed" to meet official
  criteria. For example, reducing the sample for certain countries in
  TIMSS 1995. Default is `TRUE`.

## Value

a data frame.

## Examples

``` r
data(timss99)
leaguetable(df = timss99, year = 1999)
#>   study study2 year subject  group    N     mean      se   CIdown     CIup
#> 1 TIMSS     G8 1999    math  Chile 1076 392.7611 5.45224 382.0749 403.4473
#> 2 TIMSS     G8 1999    math  Japan  885 578.4152 2.95280 572.6278 584.2026
#> 3 TIMSS     G8 1999    math Taiwan 1039 590.4357 4.94637 580.7410 600.1304
#> 4 TIMSS     G8 1999 science  Chile 1076 420.5514 4.78608 411.1708 429.9319
#> 5 TIMSS     G8 1999 science  Japan  885 550.6614 2.59819 545.5690 555.7538
#> 6 TIMSS     G8 1999 science Taiwan 1039 573.3537 5.09858 563.3606 583.3467
```
