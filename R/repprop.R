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

repprop <- function(x,categories = NULL,
                    setup = NULL,
                     repwt,wt,df,
                     method,
                     group = NULL, exclude = NULL){



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
                    var = 0))[-1]

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


.simplefreq <- function(x,
                        group = NULL, wt = NULL,
                        valid = NULL, invalid = NULL,
                        exclude = NULL){

  # Checks -----

  if(isdf(x))
    if(min(sapply(c(x),function(i) length(class(i))))!=1){
      x <- untidy(x)
    }


  # vector to data frame
  DF <- as.matrix(x)
  nDF <- ncol(DF)
  rDF <- nrow(DF)

  if(is.null(colnames(DF))){
    colnames(DF) <- paste0('V',1:nDF)
  }



  # group
  if(!is.null(group)){

    GR <- group



    if(!is.null(exclude)){


      exclude <- unique(exclude)

    }


  }else{
    GR <- rep("ALL",rDF)
  }


  # wt

    WT <- rep(1,rDF)




  # Process ----



  ## Possible and valid values ----

  # psb: possible value = data + valid + invalid
  # psv: valid value = if(valid){valid}{data}


  # Deciding which values are possible
  # if(is.null(valid)){
  # possible values by column
  psb <- lapply(1:nDF,function(x) {
    out <- unique(c(NA,DF[,x]))
    out[order(out)]

  })
  psv <- lapply(psb,function(x)
    x[!is.na(x)]
  )
  # }else{

  if(!is.null(valid)){

    if((is.atomic(valid)&is.vector(valid))){
      valid <- list(valid)
    }


    psb <- lapply(1:nDF,function(x){
      out <- unique(c(psb[[x]],valid[[x]],NA))
      out[order(out)]
    })

    psv <- lapply(1:nDF,function(x){

      if(is.null(valid[[x]])){
        out <- omitna(psb[[x]])
      }else{
        out <- unique(valid[[x]])
      }


      out[order(out)]
    })



  }

  # Deciding which values are valid
  if(!is.null(invalid)){

    if((is.atomic(invalid)&is.vector(invalid))){
      invalid <- list(invalid)
    }

    psb <- lapply(1:length(psb),function(x){
      out <- unique(c(psb[[x]],invalid[[x]]))
      out[order(out)]
    })

    psv <-  lapply(1:length(psv), function(i){
      psv[[i]][!psv[[i]]%in%invalid[[i]]]
    })



  }

  names(psb) <- colnames(DF)
  names(psv) <- colnames(DF)








  ## Frequencies ----




  UG <- unique(GR)
  UG <- UG[order(UG)]


  frq <-  lapply(1:nDF, function(k){

    indv <- psb[[k]]%in%psv[[k]]

    outj <- lapply(UG, function(j){

      DFj <- as.matrix(DF[GR%in%j,])
      WTj <- WT[GR%in%j]

      outi <- sapply(psb[[k]],function(i){

        # print(i)

        ind <- (DFj[,k]%in%i)
        cbind(N = sum(ind),W = sum(WTj[ind]))


      })

      outi <- cbind(t(outi),
                    NTOT = sum(outi[1,]),
                    WTOT = sum(outi[2,]),
                    NVAL = sum(outi[1,indv]),
                    WVAL = sum(outi[2,indv]))


    })

    if(length(outj)>1){
      outj <- do.call(rbind,c(list(Reduce("+",outj[!UG%in%exclude])),outj))
      rownames(outj) <- rep(c('ALL',UG),each = length(indv))
    }else{
      outj <- outj[[1]]
      rownames(outj) <- rep(c(UG),each = length(indv))
    }



    outj <- cbind.data.frame(item = colnames(DF)[k],
                             category = psb[[k]],
                             valid = indv,
                             group = rownames(outj),
                             frequency = outj[,1],
                             frequencyWT = outj[,2],
                             outj[,-(1:2)])
    rownames(outj) <- NULL
    outj

  })
  frq <- do.call(rbind,frq)

  frq <- cbind.data.frame(frq[,c('item','category','valid',
                                 'group','frequency')],
                          prop = frq$frequency/frq$NTOT,
                          propvalid = frq$frequency/frq$NVAL,
                          frequencyWT = frq$frequencyWT,
                          propWT = frq$frequencyWT/frq$WTOT,
                          propvalidWT = frq$frequencyWT/frq$WVAL)
  frq$propvalid[!frq$valid] <- NA
  frq$propvalidWT[!frq$valid] <- NA







  ## Merge ----

  out <- lapply(1:length(psb),function(x){
    cbind.data.frame(item = names(psb)[x],
                     category = psb[[x]])
  })
  out <- do.call(rbind.data.frame,out)
  out <- cbind.data.frame(out,group = rep(unique(frq$group),each = nrow(out)))
  out$ind <- 1:nrow(out)

  out <- merge(out,frq,sort = FALSE,all.x = TRUE)
  out <- out[order(out$ind),]
  rownames(out) <- NULL

  if(is.null(wt)){
    out <- out[,c('item','category','valid','group',
                  'frequency','prop','propvalid')]
  }else{
    out <- out[,c('item','category','valid','group',
                  'frequency','prop','propvalid',
                  'frequencyWT','propWT','propvalidWT')]
  }
  repl <- is.na(out)
  repl[,1:3] <- FALSE
  out[repl] <- 0
  out[!out$valid,grepl('propvalid',colnames(out))] <- NA

  # Output ----

  out

}


