# Creation of Replicate Weights

Creates replicate weights given jackknife replicates and jackknife
zones.

## Usage

``` r
repcreate(df, wt, jkzone, jkrep, repwtname = "RWT", reps = NULL, method)

repcreateILSA(study, year, df, repwtname = "RWT")
```

## Arguments

- df:

  a data frame.

- wt:

  a string specifying the name of the column (within `df`) with the
  total weights.

- jkzone:

  a string specifying the name of the column in `df` that contains the
  jackknife zone information.

- jkrep:

  a string specifying the name of the column in `df` that contains the
  jackknife replicate information.

- repwtname:

  a string specifying the variable names for the replicate weights.

- reps:

  an integer indicating the number of replications to be created. If
  `NULL` the maximum number of zones will be used.

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

- study:

  a string indicating the study name. For checking available studies use
  `ILSAinfo$weights`.

- year:

  a numeric value indicating the study year. For checking available
  years use `ILSAinfo$weights`.

## Value

a data frame.

## Examples

``` r
head(repdata)
#>   GROUP ID GENDER      SES schoolSES item01 item02 item03 item04 item05 item06
#> 1   GR1  1      0 49.46451  49.51990      4      4     NA      4     NA      3
#> 2   GR3  2      0 49.81929  50.42130      3      2      4      4      1      1
#> 3   GR1  3      0 46.71804  49.51990     NA      3      4      4      2      2
#> 4   GR3  4      1 50.08840  50.42130      4      4      1     NA      3      1
#> 5   GR2  5      1 51.56532  49.95242      4     NA      4      4      2      2
#> 6   GR1  6      0 49.12905  49.51990      3     NA      4      4     NA      1
#>   item07 item08 item09 item10 item11 item12 item13 item14 item15 item16 item17
#> 1      2      1      1     NA      2      1      2      2      1     NA      1
#> 2      3      2      2      1      2      1      1      2      1      2     NA
#> 3      3      2      3      2      3      2      3     NA      1      1      1
#> 4      2     NA      1      1      1      1      1      1      1      2      1
#> 5      3     NA      3      2      2     NA      2      2      2      1      1
#> 6     NA      2      2      2     NA      1      1      1      1      1      1
#>   item18 item19 item20 item21 item22 item23 item24 item25   Math1   Math2
#> 1      2      1      2      1      3     NA      4      1  0.3169  0.3743
#> 2      2      1      1      1      2      2      4      4  0.2472  0.5552
#> 3      1      1      1     NA      3      3      4      4  0.3419 -0.7735
#> 4      1     NA      1      1      2      2      4      4  0.3367  0.2251
#> 5      1      1      1      1      1      2      3      4  0.7410  0.5050
#> 6      1      1      1     NA      3      2      4      4 -0.5120 -0.7480
#>     Math3   Math4   Math5 Reading1 Reading2 Reading3 Reading4 Reading5 CatMath1
#> 1 -0.5159  0.5385  0.2755  -1.5008  -1.0121  -1.0543  -1.1637  -1.6070        3
#> 2 -0.0666  0.4742  0.3283  -0.0367  -0.3127   0.1088   0.2674  -0.1438        3
#> 3 -0.0289 -0.0628 -0.7322  -1.7979  -2.4500  -2.3124  -3.1504  -2.3207        3
#> 4  1.4625  1.1424  1.2526   0.9898   0.5570   0.4939   0.8705   0.3729        3
#> 5  0.4783  0.3551  1.1254   0.5875   0.5454   0.8774   1.9997   1.0121        3
#> 6 -0.0752 -0.6820 -0.2089  -2.0346  -2.2791  -3.0703  -2.5179  -2.2149        2
#>   CatMath2 CatMath3 CatMath4 CatMath5 CatReading1 CatReading2 CatReading3
#> 1        3        2        3        3           1           1           1
#> 2        3        2        3        3           2           2           3
#> 3        2        2        2        2           1           1           1
#> 4        3        4        4        4           3           3           3
#> 5        3        3        3        4           3           3           3
#> 6        2        2        2        2           1           1           1
#>   CatReading4 CatReading5       wt jkzones jkrep
#> 1           1           1 1.491828       1     0
#> 2           3           2 1.183717       2     0
#> 3           1           1 1.351983       3     0
#> 4           3           3 1.219501       4     0
#> 5           4           4 1.298092       5     0
#> 6           1           1 1.767518       6     0

# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

head(RW)
#>     REPWT1   REPWT2   REPWT3   REPWT4   REPWT5   REPWT6   REPWT7   REPWT8
#> 1 0.000000 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828
#> 2 1.183717 0.000000 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717
#> 3 1.351983 1.351983 0.000000 1.351983 1.351983 1.351983 1.351983 1.351983
#> 4 1.219501 1.219501 1.219501 0.000000 1.219501 1.219501 1.219501 1.219501
#> 5 1.298092 1.298092 1.298092 1.298092 0.000000 1.298092 1.298092 1.298092
#> 6 1.767518 1.767518 1.767518 1.767518 1.767518 0.000000 1.767518 1.767518
#>     REPWT9  REPWT10  REPWT11  REPWT12  REPWT13  REPWT14  REPWT15  REPWT16
#> 1 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828
#> 2 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717
#> 3 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983
#> 4 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501
#> 5 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092
#> 6 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518
#>    REPWT17  REPWT18  REPWT19  REPWT20  REPWT21  REPWT22  REPWT23  REPWT24
#> 1 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828
#> 2 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717
#> 3 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983
#> 4 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501
#> 5 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092
#> 6 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518
#>    REPWT25  REPWT26  REPWT27  REPWT28  REPWT29  REPWT30  REPWT31  REPWT32
#> 1 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828
#> 2 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717
#> 3 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983
#> 4 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501
#> 5 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092
#> 6 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518
#>    REPWT33  REPWT34  REPWT35  REPWT36  REPWT37  REPWT38  REPWT39  REPWT40
#> 1 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828
#> 2 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717
#> 3 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983
#> 4 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501
#> 5 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092
#> 6 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518
#>    REPWT41  REPWT42  REPWT43  REPWT44  REPWT45  REPWT46  REPWT47  REPWT48
#> 1 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828 1.491828
#> 2 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717 1.183717
#> 3 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983 1.351983
#> 4 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501 1.219501
#> 5 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092 1.298092
#> 6 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518 1.767518
#>    REPWT49  REPWT50
#> 1 1.491828 1.491828
#> 2 1.183717 1.183717
#> 3 1.351983 1.351983
#> 4 1.219501 1.219501
#> 5 1.298092 1.298092
#> 6 1.767518 1.767518
```
