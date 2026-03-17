# Options for automatic functions

Shows all options for the automatic functions:
[`leaguetable`](https://dopatendo.github.io/ILSAstats/reference/leaguetable.md),
and
[`proflevels`](https://dopatendo.github.io/ILSAstats/reference/proflevels.md).

## Usage

``` r
autoILSA(func = c("leaguetable", "proflevels"), study = NULL)
```

## Arguments

- func:

  a character value indicating the options of which function should be
  printed.

- study:

  an optional character vector indicating the ILSA name.

## Value

a list.

## Examples

``` r
autoILSA()
#> Options for leaguetable(...):
#> 
#> year: 
#> CIVED = 1999. 
#> ICCS = 2009, 2016, 2022. 
#> ICILS = 2013, 2018, 2023. 
#> LANA = 2023. 
#> PIRLS = 2001, 2006, 2011, 2016, 2021. 
#> RLII = 1991, 2001. 
#> TIMSS = 1995, 1999, 2003, 2007, 2011, 2015, 2019, 2023. 
#> TIMSSADVANCED = 1995, 2008, 2015. 
#> TIMSSLONG = 2023.
#> 
#> specification: 
#> CIVED = G12, G8. 
#> ICCS = NULL, G8, G9. 
#> ICILS = NULL. 
#> LANA = NULL. 
#> PIRLS = NULL. 
#> RLII = NULL. 
#> TIMSS = G4, G8. 
#> TIMSSADVANCED = math, physics. 
#> TIMSSLONG = G45, G89.
#> 
#> subject: 
#> CIVED = civic knowledge. 
#> ICCS = civic knowledge. 
#> ICILS = cil, ct. 
#> LANA = math, reading. 
#> PIRLS = reading. 
#> RLII = reading. 
#> TIMSS = math, science. 
#> TIMSSADVANCED = math, physics. 
#> TIMSSLONG = math23, math24, mathchange, science23, science24, sciencechange.
```
