# ILSAstats 0.4.5
- Added `repmeanfast()`. - Missing option for PVs and by.
- Added creation of replicate weights by index.
- Added extensive documentation.
- Added support for TIMSS Longitudinal 2023.
- Added `autoILSA()`.

# ILSAstats 0.4.4
- Added support for LaNA.
- Added `repmeanCI()`.
- Added `prepdata()`.
- `untidy()` converts missings to NAs.
- Argument `year` now accepts characters.
- Improved efficiency of `repcreate()`.
- Added `proflevels()`.
- Added `proflevels.get()`.
- Added `repprop.table()`.
- Added argument `aggregates` to `replm()`, and `repprop()`.
- Added `leaguetables()`.
- Fixed a bug in `repquant()`when using tibbles.
- Fixed a bug in `repmean()` when `by` and `aggregates = NULL`.
- Added `RLII` as method.
- Fixed `CIVED` as method.
- Function `myfunc()` deleted. This function was added my mistake in 0.4.0.


# ILSAstats 0.4.0
- Added `ILSAinfo()` with ILSA information, and pre-setups  with `repcreateILSA()` and `setupILSA()`.
- Added support for CIVED, TIMSS and PIRLS pre 2015.
- Added aggregate options to `repmean()` and `reprho()`.


# ILSAstats 0.3.8
- Fixed a bug in `repglm()`.
- Fixed a bug in the method of `repmean()`.


# ILSAstats 0.3.7
 - Initial version
