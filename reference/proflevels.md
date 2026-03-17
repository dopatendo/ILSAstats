# ILSA's proficiency levels

Estimates the proficiency levels for all countries within a cycle of an
ILSA. Arguments `method`, and `reps`, are extracted from
[`autoILSA`](https://dopatendo.github.io/ILSAstats/reference/autoILSA.md)
and can be overridden by the user.

## Usage

``` r
proflevels(
  df,
  study = NULL,
  year,
  subject = NULL,
  method = NULL,
  reps = NULL,
  type = c("long", "wide1", "wide2"),
  separateSE = TRUE,
  fixN = TRUE,
  accumulated = FALSE
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

- type:

  a character value indicating the type of table to produce. Options
  include: `"long"`, for a long table with a column with the proportions
  and another one for the standard error; `"wide1"` for a wide table
  where groups are distributed in lines; `"wide2"` for a wide table
  where groups are distributed in columns.

- separateSE:

  a logical value indicating if standard errors should be separated from
  proportions, each as an element from a list. Only works for wide
  tables. Default is `TRUE`.

- fixN:

  a logical value indicating if data should be "fixed" to meet official
  criteria. For example, reducing the sample for certain countries in
  TIMSS 1995. Default is `TRUE`.

- acumulated:

  a logical value indicating if proficiency levels should be
  accumulated.

## Value

a data frame or a list.

## Examples

``` r
data(timss99)

proflevels(timss99,year = 1999,type = "long",subject = "math")
#>     group category                  level    prop      se
#> 1   Chile        0    Below Low Benchmark 0.53764 0.02425
#> 2   Japan        0    Below Low Benchmark 0.01731 0.00536
#> 3  Taiwan        0    Below Low Benchmark 0.04854 0.01000
#> 4   Chile        1          Low Benchmark 0.30316 0.01657
#> 5   Japan        1          Low Benchmark 0.08603 0.01251
#> 6  Taiwan        1          Low Benchmark 0.08282 0.01075
#> 7   Chile        2 Intermediate Benchmark 0.12673 0.02006
#> 8   Japan        2 Intermediate Benchmark 0.23168 0.01508
#> 9  Taiwan        2 Intermediate Benchmark 0.17295 0.01455
#> 10  Chile        3         High Benchmark 0.02846 0.00878
#> 11  Japan        3         High Benchmark 0.38279 0.01848
#> 12 Taiwan        3         High Benchmark 0.31132 0.01833
#> 13  Chile        4     Advanced Benchmark 0.00400 0.00356
#> 14  Japan        4     Advanced Benchmark 0.28218 0.01601
#> 15 Taiwan        4     Advanced Benchmark 0.38437 0.02309
```
