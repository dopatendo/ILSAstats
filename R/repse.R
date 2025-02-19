#' Standard Error for Estimates with Replicate Weights and Plausible Values
#'
#'
#' Calculates the standard error given a vector or list of previous estimations.
#'
#' @param er a vector or a list containing any statistic of interest
#' (e.g., percent, mean, variance, regression coefficient).
#' If it is a vector or list of \code{length==1}, the function estimates
#' standard errors without plausible values. If it is a list with
#' \code{length>1}, it estimates standard errors with plausible values.
#' @param e0 a numeric vector or a vector containing any statistic of interest
#' (e.g., percent, mean, variance, regression coefficient), computed using
#' total weights. For scenarios without plausible values,
#' \code{e0} should be a single value. For scenarios with plausible values,
#' \code{e0} should be a vector of the same length as \code{er}.
#' @param method a string indicating the name of the large-scale assessment
#' to determine the replication method to use. Available options are:
#' \code{"TIMSS"}, \code{"PIRLS"}, \code{"ICILS"}, \code{"ICCS"},
#' \code{"PISA"}, and \code{"TALIS"}. Note that \code{"TIMSS"} and \code{"PIRLS"}
#' refer to the method used from 2016 onwards.
#' Their method has not yet been implemented for previous cycles.
#' @param se a numeric vector with standard errors,
#' used by \code{repsecomp()} to estimate a composite standard error.
#' @param setup an optional list produced by \code{\link{repsetup}}.
#' @param PVse a numeric vector containing the standard errors of the estimates of each
#' plausible value.
#' @param PVe0 a numeric vector containing the point estimates of each plausible value.
#' @param df a logical value indicating if degrees should be calculated.
#'
#'
#' @return the standard error.
#'
#' @details
#'
#' The standard errors are calculated using a modifier \eqn{m}, for TIMSS
#' and ICILS: \eqn{m = 0.5}; for ICILS and ICCS: \eqn{m = 1}; and for PISA and TALIS:
#' \eqn{\frac{1}{R(1-0.5)^2}}. Depending on the statistic, one of the following
#' formulas is used.
#'
#'
#'
#' The standard error not involving plausible values is calculated by:
#'
#' \deqn{\sqrt{m\times \sum_{r=1}^{R}(\varepsilon_r-\varepsilon_0)^2}.}
#'
#' The standard error involving plausibles values and replicate weights is calculated by:
#'
#' \deqn{\sqrt{\left[ \sum_{p=1}^{P} \left( m\times \sum_{r=1}^{R}(\varepsilon_{rp}-\varepsilon_{0p})^2 \right)  \dfrac{1}{P}\right]+  \left[ \left(1+ \dfrac{1}{P} \right) \dfrac{\sum_{p=1}^{P} (\varepsilon_{0p}-\overline{\varepsilon}_{0p})^{2}}{P-1} \right]}.}
#'
#' The standard error involving plausibles values without replicate weights is calculated by:
#'
#' \deqn{\sqrt{  \dfrac{\sum_{p=1}^{P} SE^2_{\varepsilon_{0P}}}{P}+  \left[ \left(1+ \dfrac{1}{P} \right) \dfrac{\sum_{p=1}^{P} (\varepsilon_{0p}-\overline{\varepsilon}_{0p})^{2}}{P-1} \right]}.}
#'
#'
#' The standard error of the difference of
#' two statistics (\eqn{a} and \eqn{b}) from independent samples is calculated by:
#'
#' \deqn{\sqrt{SE_a^{2}+SE_b^{2}}.}
#'
#'
#' The standard error of the difference of
#' two statistics (\eqn{a} and \eqn{b}) from dependent samples
#' not involving plausible values
#' is calculated by:
#'
#' \deqn{\sqrt{m\times \sum_{r=1}^R((a_r-b_r)-(a_0-b_0))^2}.}
#'
#' The standard error of the difference of
#' two statistics (\eqn{a} and \eqn{b}) from dependent samples
#' involving plausible values
#' is calculated by:
#'
#' \deqn{\sqrt{\left[ \sum_{p=1}^{P} \left( m\times \sum_{r=1}^{R}((a_{rp}-b_{rp})-(a_{0p}-b_{0p}))^2 \right)  \dfrac{1}{P}\right]+  \left[ \left(1+ \dfrac{1}{P} \right) \dfrac{\sum_{p=1}^{P} \left((a_{0p}-b_{0p})- ( \overline{a}_{0p}-\overline{b}_{0p}) \right)^{2}}{P-1} \right]}.}
#'
#' The standard error of a composite estimate is calculated by:
#'
#' \deqn{\sqrt{\dfrac{\sum_{c=1}^CSE^2_{\varepsilon_c}}{C^{2}}}.}
#'
#' The standard error of the difference between an element (\eqn{a}) of the composite
#' and the composite is calculated by:
#'
#' \deqn{\sqrt{\dfrac{\sum_{c=1}^CSE^2_{\varepsilon_c}}{C^{2}}+\left(\dfrac{(C-1)^2-1}{C^2}\right)SE^2_a}.}
#'
#'
#'
#' Where
#' \eqn{\varepsilon} represents a statistic of interest,
#' the subindex \eqn{0} indicates an estimate using the total weights,
#' \eqn{r} indicates a replicate from a total of \eqn{R},
#'  \eqn{p} indicates a plausible value from a total of \eqn{P},
#'  and \eqn{c} indicates an element in a composite estimate from value a total of \eqn{C}.
#'
#' @example inst/examples/repse_example.R
#'
#'
#' @name repse
#'
NULL

#' @rdname repse
#' @export
#'
repse <- function(er,e0, setup = NULL,
                  method = c("TIMSS","PIRLS","ICILS","ICCS","PISA","TALIS")){


  # Newchecks ----






  if(!is.null(setup)){
    # if(setup$repwt.type!="df"){repwt <- setup$repwt}else{repwt <- get(setup$repwt)}
    # wt <- setup$wt
    method <- setup$method
    # group <- setup$group
    # exclude <- setup$exclude
    # df <- get(setup$df)
  }

  frm <- formals(repse)
  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = frm$method)


  # Checks ----

  ## er
  if(!(is.vector(er)|is.list(er)))
      stop(c("\nInvalid input for 'er'.",
             "\nIt should be a vector or a list."),call. = FALSE)

  if(!is.vector(er)&!is.numeric(er))
    stop(c("\nInvalid input for 'er'.",
           "\nIt should be numeric."),call. = FALSE)


  if(!is.list(er)&min(sapply(er,is.numeric))!=1)
    stop(c("\nInvalid input for 'er'.",
           "\nIt should be numeric."),call. = FALSE)

  ## e0
  if(!is.vector(e0)&is.numeric(e0))
    stop(c("\nInvalid input for 'e0'.",
           "\nIt should be numeric vector."),call. = FALSE)

  ## er & e0
  if(!is.vector(er)&length(e0)!=1)
    stop(c("\nInvalid input for 'e0'.",
           "\nIf 'er' is a single value, 'e0' must have length 1."),call. = FALSE)

  ## er & e0
  if(is.list(er))
    if((length(e0)!=length(er)))
      stop(c("\nInvalid input for 'e0'.",
           "\nIf 'er' is a list, 'e0' must have length as 'er'."),call. = FALSE)

  ## method
  if(!is.vector(method)&!is.character(method))
    stop(c("\nInvalid input for 'method'.",
           "\nIt should be a character vector."),call. = FALSE)

  # if(min(method%in%c('TIMSS','PIRLS','ICILS','ICCS','PISA', 'TALIS'))!=1)
  #   stop(c("\nInvalid input for 'method'.",
  #          "\nIt should be a 'TIMSS','PIRLS','ICILS','ICCS', 'PISA', or 'TALIS'."),call. = FALSE)
  #
  # method <- match.arg(method,c('TIMSS','PIRLS','ICILS','ICCS','PISA','TALIS'))

  # is pv

  if(is.atomic(er)&is.vector(er)){
    lg <- length(er)
    pv <- FALSE
  }else{

    if(length(er)==1){
      pv = FALSE
    }else{
      pv <- TRUE
    }

    lg <- length(er[[1]])

  }


  # Process ----

  method <- tolower(method)

  if(method%in%c('timss','pirls')){
    mod <- 2
  }

  if(method%in%c('icils','iccs')){
    mod <- 1
  }

  if(method%in%c('pisa','talis')){
    mod <- length(er)*((1-0.5)**2)
  }

  # Process & Output ----

  if(!pv)
    return(sqrt(sum((unlist(er)-e0)**2)/mod))


  P <- length(er)

  ebar <- sum(e0)/P

  t1 <- sum(sapply(1:P,function(x){
    sum((er[[x]]-e0[x])**2)/mod
  }))/P
  t2 <- sum((e0-sum(e0)/P)**2)/(P-1)
  t2 <- (1+1/P)*t2

  return(sqrt(t1+t2))

  # # Output ----
  #
  #
  # return(sqrt(sum((er-e0)**2)/mod))
}

#' @rdname repse
#' @export
#'
repsecomp <- function(se){



  # Checks ----

  if(!is.vector(se)&&!is.numeric(se))
    stop(c("\nInvalid input for 'se'.",
           "\nIt should be a numeric vector."),call. = FALSE)

  # Process ----

  se <- se[!is.na(se)]




  # Output ----

  return(sqrt(sum(se**2)/length(se)**2))

}



.repse <- function(er,e0,
                  method = c('timss','pirls','icils','iccs','pisa','talis')){


  # Checks ----


  if(is.atomic(er)&is.vector(er)){
    lg <- length(er)
    pv <- FALSE
  }else{

    if(length(er)==1){
      pv = FALSE
    }else{
      pv <- TRUE
    }

    lg <- length(er[[1]])

  }

  method <- tolower(method)


  if(method%in%c('timss','pirls')){
    mod <- 2
  }

  if(method%in%c('icils','iccs')){
    mod <- 1
  }

  if(method%in%c('pisa','talis')){
    mod <- lg*((1-0.5)**2)
  }

  # Process & Output ----

  if(!pv)
    return(sqrt(sum((unlist(er)-e0)**2)/mod))


  P <- length(er)

  ebar <- sum(e0)/P

  t1 <- sum(sapply(1:P,function(x){
    sum((er[[x]]-e0[x])**2)/mod
  }))/P
  t2 <- sum((e0-sum(e0)/P)**2)/(P-1)
  t2 <- (1+1/P)*t2

  return(sqrt(t1+t2))



}


.repsecomp <- function(se){

se <- se[!is.na(se)]

sqrt(sum(se**2)/length(se)**2)

}

#' @rdname repse
#' @export
#'
pvse <- function(PVse,PVe0,df = FALSE){

  if(df)
    return(.pvsedf(PVse,PVe0))

  return(.pvse(PVse,PVe0))

}

.pvse <- function(PVse,PVe0){

  m <- length(PVe0)
  B <- stats::var(PVe0)
  Ubar <- mean(PVse**2,na.rm = TRUE)

  sqrt(Ubar + (1+1/m)*B)

}

.pvsedf <- function(PVse,PVe0){

  m <- length(PVe0)
  B <- stats::var(PVe0)
  Ubar <- mean(PVse**2,na.rm = TRUE)

  eee = 1+m*Ubar/((m+1)*B)
  eee = (m-1)*(eee)**2

  c(sqrt(Ubar + (1+1/m)*B),(m-1)*(1+(m*Ubar)/((m+1)*B))**2)

}




