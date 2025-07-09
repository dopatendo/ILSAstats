#' Setup for Analysis with Replicate Weights
#'
#' Creates a list with common arguments used for analysis with replicate
#' weights.
#'
#' @inheritParams repmean
#' @param study a string indicating the study name. For checking available studies use
#' \code{ILSAinfo$weights}.
#' @param year a numeric value indicating the study year. For checking available
#' years use
#' \code{ILSAinfo$weights}.
#'
#' @return a list to be used in other functions.
#'
#' @example inst/examples/repsetup_example.R
#'
#' @export

repsetup <- function(repwt, wt, df,
                     method,
                     group = NULL,
                     exclude = NULL){



  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))


  if(!ischaval(repwt)){
    repwt <- deparse(substitute(repwt))
    repwt.type <- "df"
  }else{
    repwt.type <- "character"
  }



  out <- list(repwt = repwt, wt = wt, df = deparse(substitute(df)),
       method = method, group = group,
       exclude = exclude,
       repwt.type = repwt.type)

  structure(out,class = "repsetup")
}







#' @export
print.repsetup <- function(x, ...){

  dims <- dim(get(x$df))

  if(x$repwt.type!="df"){
    nc <- sum(grepl(x$repwt,colnames(get(x$df))))
  }else{nc <- ncol(get(x$repwt))}

  cat("repsetup:")
  cat("\ndata = ",x$df,": ",dims[2]," columns, and ",dims[1]," rows.",sep = "")
  cat("\ntotal weights = ",deparse(x$wt),".",sep = "")
  cat("\nreplicate weights = ",nc," weights.",sep = "")
  cat("\nmethod = ",deparse(x$method),".",sep = "")
  cat("\ngroups = ",deparse(x$group),".",sep = "")
  cat("\nexcluded groups = ",deparse(x$exclude),".",sep = "")
}



#' @rdname repsetup
#' @export
repsetupILSA <- function(study,
                         year,
                         repwt,
                         df,
                         group = NULL,
                         exclude = NULL){


  # Checks ----
  returnis(ischaval,study)
  returnis(isnumval,year)


  # Process ----


  x <- ILSAstats::ILSAinfo$weights
  x <- x[(tolower(x$study)%in%tolower(study))&((x$year)%in%(year)),]

  if(nrow(x)!=1){
    stop("\nCombination of study and year not available.\nCheck available studies using ILSAinfo$weights.")
  }



  if(!ischaval(repwt)){
    repwt <- deparse(substitute(repwt))
    repwt.type <- "df"
  }else{
    repwt.type <- "character"
  }



  out <- list(repwt = repwt, wt = x$totalweight,
              df = deparse(substitute(df)),
              method = x$method, group = group,
              exclude = exclude,
              repwt.type = repwt.type)




  # repsetup(repwt = repwt,
  #          wt = x$totalweight,
  #          df = df,
  #          method = x$method,
  #          group = group,
  #          exclude = exclude)

  structure(out,class = "repsetup")

}








