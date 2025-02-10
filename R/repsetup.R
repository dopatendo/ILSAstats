#' Setup for Analysis with Replicate Weights
#'
#' Creates a list with common arguments used for analysis with replicate
#' weights.
#'
#' @inheritParams repmean
#'
#'
#' @return a list to be used in other functions.
#'
#' @example inst/examples/repsetup_example.R
#'
#' @export

repsetup <- function(repwt, wt, df,
                     method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA","TALIS"),
                     group = NULL,
                     exclude = NULL){


  frm <- formals(repsetup)
  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = frm$method)


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









