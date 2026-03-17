# Prepare ILSA Data

Modifies ILSA data to meet official participation cases, selects columns
and transforms data into simple data frames converting missing values to
NAs.

## Usage

``` r
prepILSA(
  df,
  study = NULL,
  year = NULL,
  specification = NULL,
  fixN = TRUE,
  columns = NULL
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

- specification:

  a character value indicating extra specification like grade (e.g.,
  `"G8"` for TIMSS) or subject (e.g., `"Math"` for TIMSSADVANCED).

- fixN:

  a logical value indicating if data should be "fixed" to meet official
  criteria. For example, reducing the sample for certain countries in
  TIMSS 1995. Default is `TRUE`.

- columns:

  a character vector indicating which columns should be selected. If
  `NULL`, all columns will be selected.

## Examples

``` r
data(timss99)

head(timss99)
#>      IDCNTRY CNTRY IDCNTRY_STR    TOTWGT JKZONE JKREP BSMMAT01 BSMMAT02
#> 5047     152   CHL       Chile  30.45785      6     0   386.66   364.03
#> 2642     152   CHL       Chile  27.59028     41     0   432.69   452.06
#> 1414     152   CHL       Chile  40.13670     21     1   342.83   366.80
#> 8586     392   JPN       Japan 279.79560     41     1   446.48   472.37
#> 2818     152   CHL       Chile  36.95128     45     1   247.62   250.38
#> 4707     152   CHL       Chile  32.72866     75     1   498.33   514.51
#>      BSMMAT03 BSMMAT04 BSMMAT05 BSSSCI01 BSSSCI02 BSSSCI03 BSSSCI04 BSSSCI05
#> 5047   383.67   412.24   353.79   442.62   446.64   343.96   385.70   426.69
#> 2642   451.47   456.83   440.25   516.88   469.45   486.77   568.88   501.23
#> 1414   419.76   410.11   352.05   464.16   245.48   423.77   390.25   240.32
#> 8586   429.80   443.71   463.09   497.88   473.06   350.56   452.29   544.69
#> 2818   284.34   176.84   224.88   395.63   341.74   271.16   351.62   262.55
#> 4707   497.91   481.09   480.25   521.61   531.39   542.79   511.49   500.97

newdata <- prepILSA(df = timss99, columns = paste0("BSMMAT0",1:5),fixN = FALSE)

head(newdata)
#>   BSMMAT01 BSMMAT02 BSMMAT03 BSMMAT04 BSMMAT05
#> 1   386.66   364.03   383.67   412.24   353.79
#> 2   432.69   452.06   451.47   456.83   440.25
#> 3   342.83   366.80   419.76   410.11   352.05
#> 4   446.48   472.37   429.80   443.71   463.09
#> 5   247.62   250.38   284.34   176.84   224.88
#> 6   498.33   514.51   497.91   481.09   480.25
```
