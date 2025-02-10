#' Standard Error for Plausible Values
#'
#'
#' @param PVer to FILL
#' @param PVe0 to FILL
#' @param df description
#'
#'
#' @return the standard error and degrees of freedom.
#'
#'
#' @name pvse
#'
NULL

#' @rdname pvse
#' @export
#'
#'

pvse <- function(PVer,PVe0,df = FALSE){
  
  if(df)
    return(.pvsedf(PVer,PVe0))
  
  return(.pvse(PVer,PVe0))
  
}

.pvse <- function(PVer,PVe0){
  
  m <- length(PVe0)
  B <- stats::var(PVe0)
  Ubar <- mean(PVer**2,na.rm = TRUE)
  
  sqrt(Ubar + (1+1/m)*B)
  
}

.pvsedf <- function(PVer,PVe0){
  
  m <- length(PVe0)
  B <- stats::var(PVe0)
  U <- mean(PVer**2,na.rm = TRUE)
  
  eee = 1+m*U/((m+1)*B)
  eee = (m-1)*(eee)**2
  
  c(sqrt(U + (1+1/m)*B),(m-1)*(1+(m*U)/((m+1)*B))**2)
  
}





