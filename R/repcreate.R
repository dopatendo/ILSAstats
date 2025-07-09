#' Creation of Replicate Weights
#'
#' Creates replicate weights given jackknife replicates and jackknife zones.
#'
#' @param jkzone a string specifying the name of the column in \code{df}
#' that contains the jackknife zone information.
#' @param jkrep a string specifying the name of the column in \code{df}
#' that contains the jackknife replicate information.
#' @param repwtname a string specifying the variable names for the
#' replicate weights.
#' @param reps an integer indicating the number of replications to be created.
#' If \code{NULL} the maximum number of zones will be used.
#' @inheritParams repmean
#' @param study a string indicating the study name. For checking available studies use
#' \code{ILSAinfo$weights}.
#' @param year a numeric value indicating the study year. For checking available
#' years use
#' \code{ILSAinfo$weights}.
#' @param method a string indicating the name of the replication method.
#' Available options are:
#' \code{"JK2-full"}, \code{"JK2-half"},
# #' \code{"FAY-0.5"},
#' and \code{"JK2-half-1PV"}. \cr\cr
#' Additionally, ILSA names can be used, defaulting into:
#' \itemize{
#' \item \code{"TIMSS"} or \code{"PIRLS"} for \code{"JK2-full"};
#' \item \code{"ICILS"}, \code{"ICCS"}, or \code{"CIVED"} for \code{"JK2-half"};
# #' \item \code{"PISA"} or \code{"TALIS"} for \code{"FAY-0.5"};
#' and \code{"oldTIMSS"} or \code{"oldPIRLS"} for \code{"JK2-half-1PV"}.
#' }
#' Note that \code{"oldTIMSS"} and \code{"oldPIRLS"}
#' refer to the method used for TIMSS and PIRLS before 2015,
#' where within imputation variance
#' is estimated using only 1 plausible value.
#'
#' @return a data frame.
#'
#' @example inst/examples/repcreate_example.R
#'
#' @export
#'

repcreate <- function(df,
                      wt,
                      jkzone,
                      jkrep,
                      repwtname = "RWT",
                      reps = NULL,
                      method){
#
#   if(!is.null(setup)){
#     if(setup$repwt.type!="df"){repwt <- setup$repwt}else{repwt <- get(setup$repwt)}
#     wt <- setup$wt
#     method <- setup$method
#     # group <- setup$group
#     # exclude <- setup$exclude
#     df <- get(setup$df)
#   }


  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = ILSAmethods(repse = FALSE))


  # source("R/argchecks.R")
  # source("R/internal.R")

  # Check arguments ----

  ## df
  if(!is.data.frame(df))
    stop('df is not a data frame.')
  df <- untidy(df)


  ## jkzone
  if(!is.character(jkzone))
    stop('jkzone is not a character vector.')
  if(!(length(jkzone)==1&is.atomic(jkzone)))
    stop('Invalid input for jkzone')
  if(min(jkzone%in%names(df))==0)
    stop('jkzone not in data frame.')

  if(is.null(reps)){
    reps <- max(df[,jkzone],na.rm = TRUE)
  }

  ## reps
  if(!(is.numeric(reps)&is.atomic(reps)&length(reps)==1))
    stop('Invalid input for reps.')


  ## repwtname
  if(!is.character(repwtname))
    stop('repwtname is not a character vector.')
  if(!(length(repwtname)==1&is.atomic(repwtname)))
    stop('Invalid input for repwtname')

  ## jkrep
  if(!is.character(jkrep))
    stop('jkrep is not a character vector.')
  if(!(length(jkrep)==1&is.atomic(jkrep)))
    stop('Invalid input for jkrep')
  if(min(jkrep%in%names(df))==0)
    stop('jkrep not in data frame.')



  ## wt
  if(!is.character(wt))
    stop('wt is not a character vector.')
  if(!(length(wt)==1&is.atomic(wt)))
    stop('Invalid input for wt')
  if(min(wt%in%names(df))==0)
    stop('wt not in data frame.')

  ## method
  # if(!(is.character(method)&length(method)==1&is.atomic(method)))
  #   stop('Invalid input for method.')
  #
  # if(min(method%in%c('TIMSS','PIRLS','ICILS','ICCS',"oldTIMSS","oldPIRLS"))==0)
  #   stop('Invalid input for method')

  # methf <- method[1]

  # ## method
  # if(!(is.character(method)&length(method)==1&is.atomic(method)))
  #   stop('Invalid input for method.')
  #
  # if(!'JK2'%in%method)
  #   stop('Invalid input for method.')
  #
  # ## simple
  # if(!islova(simple))
  #   stop('Invalid input for simple')


  method <- tolower(method)

  # Process ----

  if(method%in%c("jk2-full",'timss','pirls')){
    simple <- FALSE
  }

  if(method%in%c("jk2-half","jk2-half-1pv",
                 'icils','iccs',"cived",
                 "oldtimss","oldpirls")){
    simple <- TRUE
  }



  RWT1 <- matrix(df[,wt],ncol = reps,nrow = nrow(df))
  for(i in 1:reps){
    RWT1[df[,jkzone]==i&df[,jkrep]%in%0,i] <- RWT1[df[,jkzone]==i&df[,jkrep]%in%0,i]*0
    RWT1[df[,jkzone]==i&df[,jkrep]%in%1,i] <- RWT1[df[,jkzone]==i&df[,jkrep]%in%1,i]*2
  }
  if(!simple){
    RWT2 <- matrix(df[,wt],ncol = reps,nrow = nrow(df))
    for(i in 1:reps){
      RWT2[df[,jkzone]==i&df[,jkrep]%in%0,i] <- RWT2[df[,jkzone]==i&df[,jkrep]%in%0,i]*2
      RWT2[df[,jkzone]==i&df[,jkrep]%in%1,i] <- RWT2[df[,jkzone]==i&df[,jkrep]%in%1,i]*0
    }
    RWT <- cbind.data.frame(RWT1,RWT2)
  }else{
    RWT <- as.data.frame(RWT1)
  }
  colnames(RWT) <- paste0(repwtname,1:ncol(RWT))

  RWT
}


#' @rdname repcreate
#' @export
repcreateILSA <- function(study,
                          year,
                          df,
                          repwtname = "RWT"){


  # Checks ----
  returnis(ischaval,study)
  returnis(isnumval,year)


  # Process ----


  x <- ILSAstats::ILSAinfo$weights
  x <- x[(tolower(x$study)%in%tolower(study))&((x$year)%in%(year)),]

  if(nrow(x)!=1){
    stop("\nCombination of study and year not available.\nCheck available studies using ILSAinfo$weights.")
  }

  repcreate(df = df,
            repwtname = repwtname,
            wt = x$totalweight,
            jkzone = x$jkzones,
            jkrep = x$jkreps,
            method = x$method,
            reps = x$reps)


}






