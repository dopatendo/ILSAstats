#' Standard Error for Estimates with Replicate Weights
#'
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
#' \code{"PISA"}, and \code{"TALIS"}.
#' @param se a numeric vector with standard errors,
#' used by \code{repsecomp()} to estimate a composite standard error.
#' @param setup an optional list produced by \code{\link{repsetup}}.
#'
#'
#' @return the standard error.
#'
#' @example inst/examples/repse_example.R
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





