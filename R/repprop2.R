#' Proportions with Replicate Weights
#'
#' Estimates proportions using replicate weights
#' for a variable or a group of plausible values variables and for one or more
#' populations.
#' For a detailed explanation on how the standard errors are estimated
#' see \code{\link{repse}}.
#'
#' @inheritParams repmean
#' @param categories a vector indicating all possible response categories.
#' If \code{NULL}, categories will be derived from the data.
#'
#'
#' @return a list.
#'
#' @example inst/examples/repprop_example.R
#' @export
#'

repprop2 <- function(x,categories = NULL,
                    setup = NULL,
                     repwt,wt,df,
                     method,
                     group = NULL, exclude = NULL,
                    aggregates = c("pooled", "composite")){



  if(!is.null(setup)){
    if(setup$repwt.type!="df"){repwt <- setup$repwt}else{repwt <- get(setup$repwt)}
    wt <- setup$wt
    method <- setup$method
    group <- setup$group
    exclude <- setup$exclude
    df <- get(setup$df)
  }


  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))

  returnis(isdf,df)
  if(!isdfonly(df)){

    if(ischavec(repwt)){
      rena = repwt
    }else{
      rena = NULL
    }

    df <- df[,c(x,group,wt,rena)]
    df <- untidy(df)
  }


  # Checks ----

  ## df
  if(!isdf(df))
    stop(c("\nInvalid input for 'df'.",
           "\nIt should be a data frame."))

  if(min(sapply(c(df),function(i) length(class(i))))!=1){
    df <- untidy(df)
  }




  ## Check if they are in df (x,pv,pv2,wt,group)

  indf <- c(x,wt,group)

  if(min(indf%in%colnames(df))==0)
    stop(c("\nInvalid input.",
           "\n",
           paste(paste0(indf[!indf%in%colnames(df)],collapse = ', '),"not found in 'df'")))


  if(length(x)>1)
    message("More than one variable provided. 'x' treated as PVs.")

  # Process ----

  ## Transformation of arguments ----

  X = df[,x,drop = FALSE]
  TW <- df[,wt]

  if(!is.null(group)){
    GR = df[,group]

    if('Pooled'%in%GR)
      stop(c("\nInvalid input for 'group'.",
             "\nNo group names should be 'Pooled'."))

    if('Composite'%in%GR)
      stop(c("\nInvalid input for 'group'.",
             "\nNo group names should be 'Composite'."))



  }else{
    GR = NULL
  }



  if(length(method)>1){
    method = method[1]
    message(paste0('More than one method provided. Only first method will be used ',
                   '(',method,').'))
  }



  if((is.atomic(repwt)&is.vector(repwt))){

    if(length(repwt==1)){
      RW = df[,grepl(repwt,colnames(df))]
    }else{
      RW = df[,repwt]
    }

  }else{
    RW <- repwt
  }

  if(!is.null(categories)){
    psb <- sort(categories)
  }else{
    psb <- sort(unique(omitna(unlist(c(X)))))
  }


  ## Proportions via means ----


  XX <- lapply(1:ncol(X),function(i){
    matrix(rep(X[,i],length(psb)),ncol = length(psb))
  })

  XX <- lapply(XX,function(j){
    sapply(X = 1:length(psb),
           FUN = function(i){
             outraw <- j[,i]
             out <- rep(0,length(outraw))
             out[outraw%in%psb[i]] <- 1
             out[is.na(outraw)] <- NA
             out


           })
  })

  xx <- XX[[1]]

  XX <- do.call(cbind,lapply(XX,function(i){as.vector(unlist(i))}))
  XN <- paste0('XX',1:ncol(XX))
  colnames(XX) <- XN



  XX <- cbind.data.frame(BY = rep(psb,each = nrow(X)),
                         XX,
                         GR = switch(is.null(GR)+1,rep(GR,length(psb)),NA),
                         # GR = switch(is.null(GR),NA,rep(GR,length(psb))),
                         # GR = rep(GR,length(psb)),
                         # GR = NA,
                         TW = rep(TW,length(psb)),
                         RRW = do.call(rbind,lapply(1:length(psb),function(i) RW)))
  colnames(XX)[1] <- ifelse(length(x)>1,'PVs',x)


  out <- sm(.repmean0(x = XN,
                    PV = ifelse(length(x)>1,TRUE,FALSE),
                    repwt = 'RRW',
                    wt = 'TW',
                    df = XX,
                    method = method,
                    # group = NULL,
                    group = switch(is.null(GR)+1,'GR',NULL),
                    by = ifelse(length(x)>1,'PVs',x),
                    exclude = exclude,
                    var = 0,
                    aggregates = aggregates))[-1]

  ## Fix of N and df remove ----

  out <- lapply(1:length(psb), function(i){


    outi <- out[[i]]
    outi <- outi[,!substr(colnames(outi),1,3)%in%c('df_','pva')]

    if(ifelse(length(x)>1,TRUE,FALSE)){
      return(outi[,!colnames(outi)=='N'])
    }


    fri <- .simplefreq(xx[,i],group = GR,exclude = exclude)

    if(is.null(GR)){
      outi[1,'N'] <- fri[fri$category%in%1,'frequency']
    }else{
      outi[c(1,3:nrow(outi)),'N'] <- fri[fri$category%in%1,'frequency']
    }



    outi
  })

  ## Names of list ----
  names(out) <- paste0(ifelse(length(x)>1,'PVs',x),'==',psb)

  # Output ----

  out

}
