# League tables

To calculate league tables of an ILSA we need to obtain the country mean
of all participant country of a cycle. Of course, we can do this using
[`repmean()`](https://dopatendo.github.io/ILSAstats/reference/repmean.md),
nevertheless, for certain ILSA all the relevant information is already
incorporated into `ILSAstats`. We can check all the available cycles
with
[`autoILSA()`](https://dopatendo.github.io/ILSAstats/reference/autoILSA.md):

``` r
autoILSA()
```

    ## Options for leaguetable(...):
    ## 
    ## year: 
    ## CIVED = 1999. 
    ## ICCS = 2009, 2016, 2022. 
    ## ICILS = 2013, 2018, 2023. 
    ## LANA = 2023. 
    ## PIRLS = 2001, 2006, 2011, 2016, 2021. 
    ## RLII = 1991, 2001. 
    ## TIMSS = 1995, 1999, 2003, 2007, 2011, 2015, 2019, 2023. 
    ## TIMSSADVANCED = 1995, 2008, 2015. 
    ## TIMSSLONG = 2023.
    ## 
    ## specification: 
    ## CIVED = G12, G8. 
    ## ICCS = NULL, G8, G9. 
    ## ICILS = NULL. 
    ## LANA = NULL. 
    ## PIRLS = NULL. 
    ## RLII = NULL. 
    ## TIMSS = G4, G8. 
    ## TIMSSADVANCED = math, physics. 
    ## TIMSSLONG = G45, G89.
    ## 
    ## subject: 
    ## CIVED = civic knowledge. 
    ## ICCS = civic knowledge. 
    ## ICILS = cil, ct. 
    ## LANA = math, reading. 
    ## PIRLS = reading. 
    ## RLII = reading. 
    ## TIMSS = math, science. 
    ## TIMSSADVANCED = math, physics. 
    ## TIMSSLONG = math23, math24, mathchange, science23, science24, sciencechange.

There is a total of 68 combinations.

For producing an automatic league table we can use
[`leaguetable()`](https://dopatendo.github.io/ILSAstats/reference/leaguetable.md),
which contains the following arguments:

- `df`: the data frame with the ILSA data.
- `study`: the name of the study as it appears in `ILSAinfo$pvs$study`.
- `year`: the year of the cycle.
- `subject`: the name of the subject as it appears in
  `ILSAinfo$pvs$year`. If left `NULL` all subjects will be estimated.
- `specification`: a string indicating extra specification for the
  identifying the study if it is necessary, as it appears in
  `ILSAinfo$pvs$study2`. If left `NULL` it is equivalent to `"-"`.
- `method`: the method to use. If left `NULL` the default method will be
  used.
- `fixN`: a logical value indicating if cases should be removed to match
  official results. This is necessary in some studies like TIMSS 1995,
  in which the public data includes cases that are not used for
  producing the league tables.

## Official league tables

We can use `timss99` data to calculate the official league table
published in the TIMSS 1999 report using:

``` r
# Only math
leaguetable(df = timss99, study = "TIMSS", year = 1999, subject = "math")
```

    ##   study study2 year subject  group    N     mean      se   CIdown     CIup
    ## 1 TIMSS     G8 1999    math  Chile 1076 392.7611 5.45224 382.0749 403.4473
    ## 2 TIMSS     G8 1999    math  Japan  885 578.4152 2.95280 572.6278 584.2026
    ## 3 TIMSS     G8 1999    math Taiwan 1039 590.4357 4.94637 580.7410 600.1304

``` r
# Only math
leaguetable(df = timss99, study = "TIMSS", year = 1999, subject = "science")
```

    ##   study study2 year subject  group    N     mean      se   CIdown     CIup
    ## 1 TIMSS     G8 1999 science  Chile 1076 420.5514 4.78608 411.1708 429.9319
    ## 2 TIMSS     G8 1999 science  Japan  885 550.6614 2.59819 545.5690 555.7538
    ## 3 TIMSS     G8 1999 science Taiwan 1039 573.3537 5.09858 563.3606 583.3467

``` r
# Both subjects
leaguetable(df = timss99, study = "TIMSS", year = 1999)
```

    ##   study study2 year subject  group    N     mean      se   CIdown     CIup
    ## 1 TIMSS     G8 1999    math  Chile 1076 392.7611 5.45224 382.0749 403.4473
    ## 2 TIMSS     G8 1999    math  Japan  885 578.4152 2.95280 572.6278 584.2026
    ## 3 TIMSS     G8 1999    math Taiwan 1039 590.4357 4.94637 580.7410 600.1304
    ## 4 TIMSS     G8 1999 science  Chile 1076 420.5514 4.78608 411.1708 429.9319
    ## 5 TIMSS     G8 1999 science  Japan  885 550.6614 2.59819 545.5690 555.7538
    ## 6 TIMSS     G8 1999 science Taiwan 1039 573.3537 5.09858 563.3606 583.3467

## Un-official league tables

Nevertheless, we can also change the method. For example, we can do this
to compare TIMSS 1999 to TIMSS 2023, therefore, we would need to use the
method `"TIMSS"` instead of `"oldTIMSS"`:

``` r
leaguetable(df = timss99, study = "TIMSS", year = 1999, method = "TIMSS")
```

    ##   study study2 year subject  group    N     mean      se   CIdown     CIup
    ## 1 TIMSS     G8 1999    math  Chile 1076 392.7611 5.46288 382.0540 403.4681
    ## 2 TIMSS     G8 1999    math  Japan  885 578.4152 2.80539 572.9167 583.9137
    ## 3 TIMSS     G8 1999    math Taiwan 1039 590.4357 5.00692 580.6223 600.2491
    ## 4 TIMSS     G8 1999 science  Chile 1076 420.5514 4.86369 411.0187 430.0840
    ## 5 TIMSS     G8 1999 science  Japan  885 550.6614 2.53919 545.6847 555.6381
    ## 6 TIMSS     G8 1999 science Taiwan 1039 573.3537 5.01030 563.5337 583.1737
