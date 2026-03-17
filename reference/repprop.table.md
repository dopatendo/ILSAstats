# Tables for Proportions with Replicate Weights

Creates tables for proportions using replicate weights for a variable or
a group of plausible values variables and for one or more populations.

## Usage

``` r
repprop.table(x, type = c("long", "wide1", "wide2"), separateSE = TRUE)
```

## Arguments

- x:

  a list produced by
  [`repprop`](https://dopatendo.github.io/ILSAstats/reference/repprop.md).

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

## Value

a adata frame or a list.

## Examples

``` r
# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                wt = "wt", # the total weights column name
                jkzone = "jkzones", # the jkzones column name
                jkrep = "jkrep", # the jkreps column name
                repwtname = "REPWT", # the desired name for the rep weights
                reps = 50, # the number of replications
                method = "ICILS") # the name of the method aka the study name




x = repprop(x = c("item01"),
            group = "GROUP",
            repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
            method = "ICILS")

repprop.table(x, type = "long")
#>        group category    prop      se
#> 1     Pooled        1 0.01532 0.00184
#> 2  Composite        1 0.01531 0.00188
#> 3        GR1        1 0.01336 0.00360
#> 4        GR2        1 0.02077 0.00341
#> 5        GR3        1 0.01180 0.00270
#> 6     Pooled        2 0.05747 0.00320
#> 7  Composite        2 0.05749 0.00335
#> 8        GR1        2 0.05008 0.00503
#> 9        GR2        2 0.06000 0.00662
#> 10       GR3        2 0.06240 0.00563
#> 11    Pooled        3 0.22237 0.00648
#> 12 Composite        3 0.22240 0.00668
#> 13       GR1        3 0.22040 0.01161
#> 14       GR2        3 0.21432 0.01117
#> 15       GR3        3 0.23248 0.01195
#> 16    Pooled        4 0.70484 0.00664
#> 17 Composite        4 0.70480 0.00701
#> 18       GR1        4 0.71616 0.01133
#> 19       GR2        4 0.70491 0.01276
#> 20       GR3        4 0.69332 0.01230


repprop.table(x, type = "wide1", separateSE = TRUE)
#> $prop
#>       group  prop.1  prop.2  prop.3  prop.4
#> 1    Pooled 0.01532 0.05747 0.22237 0.70484
#> 2 Composite 0.01531 0.05749 0.22240 0.70480
#> 3       GR1 0.01336 0.05008 0.22040 0.71616
#> 4       GR2 0.02077 0.06000 0.21432 0.70491
#> 5       GR3 0.01180 0.06240 0.23248 0.69332
#> 
#> $se
#>       group    se.1    se.2    se.3    se.4
#> 1    Pooled 0.00184 0.00320 0.00648 0.00664
#> 2 Composite 0.00188 0.00335 0.00668 0.00701
#> 3       GR1 0.00360 0.00503 0.01161 0.01133
#> 4       GR2 0.00341 0.00662 0.01117 0.01276
#> 5       GR3 0.00270 0.00563 0.01195 0.01230
#> 


repprop.table(x, type = "wide1", separateSE = FALSE)
#>       group  prop.1    se.1  prop.2    se.2  prop.3    se.3  prop.4    se.4
#> 1    Pooled 0.01532 0.00184 0.05747 0.00320 0.22237 0.00648 0.70484 0.00664
#> 2 Composite 0.01531 0.00188 0.05749 0.00335 0.22240 0.00668 0.70480 0.00701
#> 3       GR1 0.01336 0.00360 0.05008 0.00503 0.22040 0.01161 0.71616 0.01133
#> 4       GR2 0.02077 0.00341 0.06000 0.00662 0.21432 0.01117 0.70491 0.01276
#> 5       GR3 0.01180 0.00270 0.06240 0.00563 0.23248 0.01195 0.69332 0.01230


repprop.table(x, type = "wide2", separateSE = TRUE)
#> $prop
#>   category     Pooled  Composite        GR1        GR2        GR3
#> 1        1 0.01532088 0.01531238 0.01336472 0.02077053 0.01180189
#> 3        2 0.05747261 0.05749194 0.05007533 0.06000229 0.06239819
#> 5        3 0.22236927 0.22239926 0.22039974 0.21431643 0.23248162
#> 7        4 0.70483724 0.70479642 0.71616021 0.70491075 0.69331830
#> 
#> $se
#>   category      Pooled   Composite         GR1         GR2         GR3
#> 2        1 0.001836850 0.001882179 0.003604455 0.003406165 0.002699875
#> 4        2 0.003202504 0.003345585 0.005028686 0.006617893 0.005626032
#> 6        3 0.006475113 0.006684990 0.011608645 0.011167149 0.011947218
#> 8        4 0.006639949 0.007012546 0.011329773 0.012761530 0.012302917
#> 


repprop.table(x, type = "wide2", separateSE = FALSE)
#>   category statistic      Pooled   Composite         GR1         GR2
#> 1        1      prop 0.015320884 0.015312381 0.013364717 0.020770534
#> 2        1        se 0.001836850 0.001882179 0.003604455 0.003406165
#> 3        2      prop 0.057472611 0.057491938 0.050075329 0.060002290
#> 4        2        se 0.003202504 0.003345585 0.005028686 0.006617893
#> 5        3      prop 0.222369269 0.222399263 0.220399744 0.214316429
#> 6        3        se 0.006475113 0.006684990 0.011608645 0.011167149
#> 7        4      prop 0.704837237 0.704796418 0.716160210 0.704910747
#> 8        4        se 0.006639949 0.007012546 0.011329773 0.012761530
#>           GR3
#> 1 0.011801893
#> 2 0.002699875
#> 3 0.062398195
#> 4 0.005626032
#> 5 0.232481616
#> 6 0.011947218
#> 7 0.693318296
#> 8 0.012302917
```
