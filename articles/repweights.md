# Replicate weights

Depending on our data, it may contain or not replicate weights. If it
does not, we need to create them. For this we will need variables
containing the jackknife zones, the jackknife replications, and the
total weights. Normally, all this information is provided within the
data.

For example, using the data `timss99`, we can find the variables
`"JKZONE"`, `"JKREP"`, and `"TOTWGT"`:

``` r
colnames(repdata)
```

    ##  [1] "GROUP"       "ID"          "GENDER"      "SES"         "schoolSES"  
    ##  [6] "item01"      "item02"      "item03"      "item04"      "item05"     
    ## [11] "item06"      "item07"      "item08"      "item09"      "item10"     
    ## [16] "item11"      "item12"      "item13"      "item14"      "item15"     
    ## [21] "item16"      "item17"      "item18"      "item19"      "item20"     
    ## [26] "item21"      "item22"      "item23"      "item24"      "item25"     
    ## [31] "Math1"       "Math2"       "Math3"       "Math4"       "Math5"      
    ## [36] "Reading1"    "Reading2"    "Reading3"    "Reading4"    "Reading5"   
    ## [41] "CatMath1"    "CatMath2"    "CatMath3"    "CatMath4"    "CatMath5"   
    ## [46] "CatReading1" "CatReading2" "CatReading3" "CatReading4" "CatReading5"
    ## [51] "wt"          "jkzones"     "jkrep"

Now, we would need also to select some options for the creation of the
replicate weights: a) we need to establish the number of replications;
and b) we need to decide which method we are going to use.

For the number of replications, we use the argument `reps`, that is
`NULL` by default (we recommend using this option). If left `NULL`, the
number of replications will be determined using the `"jkrep"` variable.
Nevertheless, we can adjust this if needed.

## Method

More importantly, we need to decide on the method. All methods are
available for all ILSAs, so you are able to use the official one or
experiment with other. For creating replicate weights, 4 methods are
implemented, you can use the method name or the ILSA names:

- `"JK-full"`, corresponds to `"TIMSS"`, `"PIRLS"` and `"LANA"`.
- `"JK-half"`, corresponds to `"ICILS"`, `"ICCS"`, and `"CIVED"`.
- `"JK2-half-1PV"`, corresponds to `"oldTIMSS"`, `"oldPIRLS"`, and
  `"RLII"` (for TIMSS and PIRLS conducted before 2015).
- `"FAY-0.5"`, corresponds to `"PISA"` and `"TALIS"`.

Our data corresponds to TIMSS 1999, so we are advised to use either the
old method of TIMSS to replicate results (`"oldTIMSS"`), or the new
method of TIMSS to compare it to new results (`"TIMSS"`).

## Creation of weights

So, now knowing our options we can use `timss99` data and
[`repcreate()`](https://dopatendo.github.io/ILSAstats/reference/repcreate.md)
for creating the replicate weights:

``` r
RW1 <- repcreate(df = timss99,
                jkzone = "JKZONE",
                jkrep = "JKREP",
                wt = "TOTWGT",
                method = "oldTIMSS")
```

We can see that this new object is a data frame with 75 columns
(representing the 75 replications):

``` r
class(RW1)
```

    ## [1] "data.frame"

``` r
ncol(RW1)
```

    ## [1] 75

``` r
head(colnames(RW1))
```

    ## [1] "RWT1" "RWT2" "RWT3" "RWT4" "RWT5" "RWT6"

## Automatic creation of weights

Some ILSA information are already included in `ILSAstats`, so if our
data corresponds to any of the included ones, we can create the
replicate weights easier, using
[`repcreateILSA()`](https://dopatendo.github.io/ILSAstats/reference/repcreate.md).
To check which ILSA information is included we can use
`ILSAinfo$weights`:

``` r
ILSAinfo$weights
```

    ##            study  study2 year   method country jkzones jkreps reps totalweight
    ## 1          ICILS       - 2013    ICILS   CNTRY JKZONES JKREPS   75     TOTWGTS
    ## 2          ICILS       - 2018    ICILS   CNTRY JKZONES JKREPS   75     TOTWGTS
    ## 3          ICILS       - 2023    ICILS   CNTRY JKZONES JKREPS   75     TOTWGTS
    ## 4           ICCS      G8 2009    ICILS COUNTRY JKZONES JKREPS   75     TOTWGTS
    ## 5           ICCS      G9 2009    ICILS COUNTRY JKZONES JKREPS   75     TOTWGTS
    ## 6           ICCS       - 2016    ICILS COUNTRY JKZONES JKREPS   75     TOTWGTS
    ## 7           ICCS       - 2022    ICILS COUNTRY JKZONES JKREPS   75     TOTWGTS
    ## 8          CIVED      G8 1999    CIVED IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 9          CIVED     G12 1999    CIVED IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 10          RLII       - 1991     RLII IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 11          RLII       - 2001     RLII IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 12         TIMSS      G4 1995 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 13         TIMSS      G8 1995 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 14         TIMSS      G8 1999 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 15         TIMSS      G4 2003 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 16         TIMSS      G8 2003 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 17         TIMSS      G4 2007 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 18         TIMSS      G8 2007 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 19         TIMSS      G4 2011 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 20         TIMSS      G8 2011 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 21         TIMSS      G4 2015    TIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 22         TIMSS      G8 2015    TIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 23         TIMSS      G4 2019    TIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 24         TIMSS      G8 2019    TIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 25         TIMSS      G4 2023    TIMSS IDCNTRY  JKZONE  JKREP  125      TOTWGT
    ## 26         TIMSS      G8 2023    TIMSS IDCNTRY  JKZONE  JKREP  125      TOTWGT
    ## 27 TIMSSADVANCED    math 1995 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 28 TIMSSADVANCED    math 2008 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 29 TIMSSADVANCED    math 2015    TIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 30 TIMSSADVANCED physics 1995 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 31 TIMSSADVANCED physics 2008 oldTIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 32 TIMSSADVANCED physics 2015    TIMSS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 33         PIRLS       - 2001 oldPIRLS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 34         PIRLS       - 2006 oldPIRLS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 35         PIRLS       - 2011 oldPIRLS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 36         PIRLS       - 2016    PIRLS IDCNTRY  JKZONE  JKREP   75      TOTWGT
    ## 37         PIRLS       - 2021    PIRLS IDCNTRY  JKZONE  JKREP  125      TOTWGT
    ## 38          LANA       - 2023     LANA IDCNTRY  JKZONE  JKREP  125      TOTWGT
    ## 39     TIMSSLONG     G89 2023    TIMSS IDCNTRY  JKZONE  JKREP  125      TOTWGT
    ## 40     TIMSSLONG     G45 2023    TIMSS IDCNTRY  JKZONE  JKREP  125      TOTWGT

So this information already contains the name of the method, the
jackknife zones, the jackknife replications, and the total weights from
a total of 40 cycles.

So, for our TIMSS 1999 case, we would need to specify the name of the
study, the cycle and the data:

``` r
RW2 <- repcreateILSA(study = "TIMSS",year = 1999,df = timss99)
```

This will produce the same results as `RW1`:

``` r
identical(RW1,RW2)
```

    ## [1] TRUE
