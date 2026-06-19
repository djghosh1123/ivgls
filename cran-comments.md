## Test environments
* macOS Tahoe 26.5.1, R 4.6.0 (local)
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Notes for CRAN reviewers
* `glmgraph` is listed in `Suggests` but is not currently on CRAN.
  It is available at https://github.com/cran/glmgraph. Both `ivgl()`
  and `ivgl_s()` check for its presence via `requireNamespace()` and
  throw an informative error if absent. All other functions work
  without it.

## Downstream dependencies
* None.
