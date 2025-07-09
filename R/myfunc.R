#' My Wrapper Function
#'
#' @param x Main input parameter
#' @param method,setup Parameters passed to [repse()]. Defaults:
#'   - method: `r eval(formals(repse)$method)`
#'   - setup: `r eval(formals(repse)$setup)`
#' @param ... Additional arguments passed to [repse()]
#'
#' @return Description of return value
#' @export
myfunc <- function(x,
                        method = formals(repse)$method,
                        setup = formals(repse)$setup,
                        ...) {
  repse(x, method = method, setup = setup, ...)
}
