# Changelog

## ILSAstats 0.4.5

- Added extensive documentation.
- Added support for TIMSS Longitudinal 2023.
- Added
  [`autoILSA()`](https://dopatendo.github.io/ILSAstats/reference/autoILSA.md).

## ILSAstats 0.4.4

CRAN release: 2026-03-13

- Added support for LaNA.
- Added
  [`repmeanCI()`](https://dopatendo.github.io/ILSAstats/reference/repmeanCI.md).
- Added `prepdata()`.
- `untidy()` converts missings to NAs.
- Argument `year` now accepts characters.
- Improved efficiency of
  [`repcreate()`](https://dopatendo.github.io/ILSAstats/reference/repcreate.md).
- Added
  [`proflevels()`](https://dopatendo.github.io/ILSAstats/reference/proflevels.md).
- Added
  [`proflevels.get()`](https://dopatendo.github.io/ILSAstats/reference/proflevels.get.md).
- Added
  [`repprop.table()`](https://dopatendo.github.io/ILSAstats/reference/repprop.table.md).
- Added argument `aggregates` to
  [`replm()`](https://dopatendo.github.io/ILSAstats/reference/replm.md),
  and
  [`repprop()`](https://dopatendo.github.io/ILSAstats/reference/repprop.md).
- Added `leaguetables()`.
- Fixed a bug in
  [`repquant()`](https://dopatendo.github.io/ILSAstats/reference/repquant.md)when
  using tibbles.
- Fixed a bug in
  [`repmean()`](https://dopatendo.github.io/ILSAstats/reference/repmean.md)
  when `by` and `aggregates = NULL`.
- Added `RLII` as method.
- Fixed `CIVED` as method.
- Function `myfunc()` deleted. This function was added my mistake in
  0.4.0.

## ILSAstats 0.4.0

CRAN release: 2025-07-15

- Added
  [`ILSAinfo()`](https://dopatendo.github.io/ILSAstats/reference/ILSAinfo.md)
  with ILSA information, and pre-setups with
  [`repcreateILSA()`](https://dopatendo.github.io/ILSAstats/reference/repcreate.md)
  and `setupILSA()`.
- Added support for CIVED, TIMSS and PIRLS pre 2015.
- Added aggregate options to
  [`repmean()`](https://dopatendo.github.io/ILSAstats/reference/repmean.md)
  and
  [`reprho()`](https://dopatendo.github.io/ILSAstats/reference/reprho.md).

## ILSAstats 0.3.8

CRAN release: 2025-04-16

- Fixed a bug in
  [`repglm()`](https://dopatendo.github.io/ILSAstats/reference/repglm.md).
- Fixed a bug in the method of
  [`repmean()`](https://dopatendo.github.io/ILSAstats/reference/repmean.md).

## ILSAstats 0.3.7

CRAN release: 2025-02-21

- Initial version
