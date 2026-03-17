# Setup

The majority of `ILSAstats` functions will need us to provide
information about which data we are using (`df`), what are the total
weights (`wt`), what are the replicate weights `(repwt)`, which
jackknife method is used `(method)`, how to group the results `(group)`,
and which groups should be excluded from pooled and composite aggregates
`(exclude)`.

Since this is a lot of information, and most probably we will use it
several times, we can allocate all of it within a single object of class
`"repsetup"`.

## Creation of setup object

For creating the setup we will need both the replicate weights, and the
data. So, using the included `timss99` data and
[`repsetup()`](https://dopatendo.github.io/ILSAstats/reference/repsetup.md),
we would:

``` r
RW2 <- repcreateILSA(study = "TIMSS", year = 1999, df = timss99)
```

Then, we would need to specify:

``` r
ST1 <- repsetup(repwt = RW2,
                wt = "TOTWGT",
                df = timss99,
                method = "oldTIMSS",
                group = "IDCNTRY_STR")
```

We can print `ST1` to check what information it holds:

``` r
ST1
```

    ## repsetup:
    ## data = timss99: 16 columns, and 3000 rows.
    ## total weights = "TOTWGT".
    ## replicate weights = 75 weights.
    ## method = "oldTIMSS".
    ## groups = "IDCNTRY_STR".
    ## excluded groups = NULL.

## Automatic creation of setup object

As we did for the automatic creation of weights in `RW2`, we can also
create automatically the setup using
[`repsetupILSA()`](https://dopatendo.github.io/ILSAstats/reference/repsetup.md):

``` r
ST2 <- repsetupILSA(study = "TIMSS",
                    year = 1999,
                    repwt = RW2,
                    df = timss99,
                    group = "IDCNTRY_STR")
```

And, as we can see, we would obtain the same result as `ST1`:

``` r
identical(ST1,ST2)
```

    ## [1] TRUE
