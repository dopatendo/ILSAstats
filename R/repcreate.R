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
#' @param index a logical value indicating if the result should be just an index
#' of zero and double weights instead of a matrix. Default is \code{FALSE}.
#'
#' @return a data frame or a list.
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
                       method,
                      index = FALSE){


  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = ILSAmethods(repse = FALSE))


  # source("R/argchecks.R")
  # source("R/internal.R")

  # Check arguments ----

  ## df
  if(!is.data.frame(df))
    stop('df is not a data frame.')

  if(!isdfonly(df)){


    df <- df[,c(wt,jkrep,jkzone)]
    df <- untidy(df)
  }



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



  method <- tolower(method)

  # Process ----

  if(method%in%c("jk2-full",'timss','pirls','lana')){
    simple <- FALSE
  }

  if(method%in%c("jk2-half","jk2-half-1pv",
                 'icils','iccs',"cived","rlii",
                 "oldtimss","oldpirls")){
    simple <- TRUE
  }




  if(index){
    out <- .repcreateIndex(reps = reps,
                           jzn = df[,jkzone],
                           jre = df[,jkrep],
                           method = method, simple = simple)
    return(out)
  }


  RWT <- matrix(df[,wt],ncol = reps,nrow = nrow(df))
  MM <- matrix(1,ncol = reps,nrow = nrow(df))
  # JKZ <- matrix(df[,jkzone],ncol = reps,nrow = nrow(df))
  # JKR <- matrix(df[,jkrep],ncol = reps,nrow = nrow(df))
  # JKi <- matrix(1:reps,nrow = nrow(df),ncol = reps,byrow = T)

  # isjkzi <- (JKZ==JKi)

  isjkzi = matrix(rep(df[,jkzone],each=reps)==rep(1:reps,nrow(df)),nrow = nrow(df),byrow = TRUE)

  # table(c(aa==isjkzi))

  # MM[isjkzi&(JKR==0)] <- 0
  # MM[isjkzi&(JKR==1)] <- 2
  # rm(JKZ,JKi,JKR)

  MM[(isjkzi*(df[,jkrep]==0))==1] <- 0
  MM[(isjkzi*(df[,jkrep]==1))==1] <- 2



  # for(i in 1:reps){
  #   RWT1[df[,jkzone]==i&df[,jkrep]%in%0,i] <- RWT1[df[,jkzone]==i&df[,jkrep]%in%0,i]*0
  #   RWT1[df[,jkzone]==i&df[,jkrep]%in%1,i] <- RWT1[df[,jkzone]==i&df[,jkrep]%in%1,i]*2
  # }
  if(!simple){

    # for(i in 1:reps){

    # MM[(JKZ==JKi)&(JKR==0)] <- 2
    # MM[(JKZ==JKi)&(JKR==1)] <- 0

    # RWT2[df[,jkzone]==i&df[,jkrep]%in%0,i] <- RWT2[df[,jkzone]==i&df[,jkrep]%in%0,i]*2
    # RWT2[df[,jkzone]==i&df[,jkrep]%in%1,i] <- RWT2[df[,jkzone]==i&df[,jkrep]%in%1,i]*0
    # }
    RWT <- cbind.data.frame(cbind(RWT,RWT)*cbind(MM,abs(MM-2)))
  }else{
    RWT <- as.data.frame(RWT*MM)
  }
  colnames(RWT) <- paste0(repwtname,1:ncol(RWT))

  RWT
}



#' @rdname repcreate
#' @export
repcreateILSA <- function(study,
                          year,
                          df,
                          repwtname = "RWT", index = FALSE){


  # Checks ----
  # returnis(isval,year);year <- as.numeric(year)
  returnis(ischaval,study)
  # returnis(isnumval,year)
  returnis(ischaeqnum,year)


  # Process ----


  x <- ILSAstats::ILSAinfo$weights
  x <- x[(tolower(x$study)%in%tolower(study))&((x$year)%in%(year)),]
  x <- unique(x[,!colnames(x)%in%"study2"])

  if(nrow(x)!=1){
    stop("\nCombination of study and year not available.\nCheck available studies using ILSAinfo$weights.")
  }

  repcreate(df = df,
            repwtname = repwtname,
            wt = x$totalweight,
            jkzone = x$jkzones,
            jkrep = x$jkreps,
            method = x$method,
            reps = x$reps,
            index = index)


}



repcreateOLD <- function(df,
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

  if(!isdfonly(df)){


    df <- df[,c(wt,jkrep,jkzone)]
    df <- untidy(df)
  }



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



  method <- tolower(method)

  # Process ----

  if(method%in%c("jk2-full",'timss','pirls','lana')){
    simple <- FALSE
  }

  if(method%in%c("jk2-half","jk2-half-1pv",
                 'icils','iccs',"cived","rlii",
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


.repcreateIndex <- function(reps, jzn, jre, method, simple){
  # jzn <- df[,jkzones]
  # jre <- df[,jkrep]



  ze = lapply(1:reps, function(i){
    which(jzn==i&jre==0)
  })
  do = lapply(1:reps, function(i){
    which(jzn==i&jre==1)
  })

  if(!simple){
    zef <- c(ze,do)
    dof <- c(do,ze)
  }else{
    zef <- ze
    dof <- do
  }



  out <- list(zef, dof)
  attributes(out)$multiplier <- c(0,2)
  attributes(out)$method <- method
  attributes(out)$reps <- reps

  class(out) <- c("repweights.index","repweights",class(out))

  return(out)
}

#' @export
print.repweights.index <- function(x, ...){

cat(paste0('Replicate weights indices for ', attributes(x)$reps, " replications."))
# cat(paste0('Weights for ',length(x[[1]])," replications."))
cat(paste0("\nMethod: ",attributes(x)$method,"."))

}
