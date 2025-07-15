#' Mean, Variance and Standard Deviation with Replicate Weights
#'
#' Estimates the mean, variance and standard deviation with replicate weights
#' for a variable or a group of variables and for one or more
#' populations. For a detailed explanation on how the standard errors are estimated
#' see \code{\link{repse}}.
#'
#' @param x a string vector specifying variable names (within \code{df}) for
#' analysis.
#' @param PV a logical value indicating if the variables in \code{x} are
#' plausible values.
#' @param repwt a string indicating the common names for the replicate weights
#' columns (within \code{df}), or a data frame with the replicate weights.
#' @param wt a string specifying the name of the column (within \code{df}) with the
#'  total weights.
#' @param df a data frame.
#' @param var a string indicating the method to use for the variance:
#' \code{"unbiased"} calculates the unbiased estimate (n-1); \code{"ML"}
#' calculates the maximum likelihood estimate.
#' @param group a string specifying the variable name (within \code{df}) to be
#' used for grouping. Categories in \code{group} are treated as independent, e.g.,
#' countries.
#' @param by a string specifying a second variable (within \code{df}) for grouping.
#' Categories used in \code{by} are not considered independent, e.g., gender
#' within a country. If used,
#' the output will be a list with the same length as the  unique values of
#' \code{by}. This can only be used for analyses with one variable or a group
#' of PVs.
#' @param exclude a vector indicating which groups
#' (in the same format as \code{group})
#' should be excluded from the pooled and composite estimates.
#' @param aggregates a string vector indicating which aggregates should be
#' included, options are \code{"pooled"} and \code{"composite"}, both options can be
#' used at the same time. If \code{NULL} no aggregate will be estimated.
#' @param zones a string specifying the name of the variable containing the
#' replicate zones.
#' Used for calculating the number of zones to be used by variable and group.
#' If \code{NULL}, zones are not be calculated.
#' @inheritParams repse
#'
#'
#' @return a data frame or a list.
#'
#'
#'
#'
#' @example inst/examples/repmean_example.R
#' @export
#'


repmean <- function(x, PV = FALSE, setup = NULL, repwt, wt, df,
                    method,
                    var = c("unbiased","ML"), group = NULL, by = NULL,
                    exclude = NULL,
                    aggregates = c("pooled","composite"),
                    zones = NULL){

  if(!is.null(setup)){
    if(setup$repwt.type!="df"){repwt <- setup$repwt}else{repwt <- get(setup$repwt)}
    wt <- setup$wt
    method <- setup$method
    group <- setup$group
    exclude <- setup$exclude
    # aggregates <- setup$aggregates
    df <- get(setup$df)
  }

  # returnis(ischavec, method)
  # method <- returnis(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))



  # Checks ----

  frm <- formals(repmean)



  ## class

  returnis(ischavec,x)
  returnis(islova,PV)
  returnis(is.chavec.or.dfonly,repwt)
  returnis(ischaval,wt)
  returnis(isdf,df)
  if(!isdfonly(df)){

   if(ischavec(repwt)){
     rena = repwt
   }else{
     rena = NULL
   }

    df <- df[,c(x,group,wt,by,zones,rena)]
    df <- untidy(df)
  }
  returnisNULL(ischaval, group)
  returnisNULL(ischaval, by)
  returnisNULL(ischavec, exclude)
  returnisNULL(ischaval, zones)
  returnis(ischavec, method)
  returnis(ischavec, var)



  ## match option
  method <- returnis(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))
  var <- returnis(isinvec,x = var[1L],choices = frm$var)

  aggregates <- returnisNULL(isinvecmul,x = aggregates, choices = frm$aggregates)

  ## Combinations

  ### x & PV
  if(length(x)==1&PV)
    stop(c("\nInvalid input for 'x'.",
           "\nOnly one PV provided."),call. = FALSE)

  ### PV are complete
  if(PV){
    if(!all(rowSums(is.na(df[,x]))%in%c(0,length(x))))
      stop(c("\nInvalid input for 'x'.",
             "\nThere are cases with incomplete PVs."),call. = FALSE)
  }

  ### by for 1 variable
  if(length(x)>1&!PV&!is.null(by))
    stop(paste0("Invalid input for 'by'.",
                "\nIt is only implemented for one variable."),call. = FALSE)


  ### x, wt, group, by, zones in df
  indf <- c(x, wt, group, by, zones)



  if(!all(indf%in%colnames(df)))
    stop(c("\nInvalid input.",
           "\n",
           paste(paste0(indf[!indf%in%colnames(df)],collapse = ', '),"not found in 'df'")),
         call. = FALSE)

  if(length(indf)!=length(unique(indf)))
    stop(paste0("\nInvalid input. Repeated arguments."),
         call. = FALSE)


  ### exclude + group
  if(!is.null(exclude)&&min(!exclude%in%df[,group])==1)
    warning("Check 'exclude', values not found in 'group'.",call. = FALSE)


  ### repwt(df) + df
  if(!is.vector(repwt)&&(nrow(repwt)!=nrow(df)))
    stop(c("\nInvalid input for 'repwt'.",
           "\nIf it is a data frame it should have the same number of rows as 'df'."),call. = FALSE)





  # Process ----

  # Transformation of arguments ----

    ### X - x
  X <- df[,x]
  TW <- df[,wt]

    ### GR - groups
    if(!is.null(group)){
      GR = df[,group]

      if('Pooled'%in%GR)
        stop(c("\nInvalid input for 'group'.",
               "\nNo group names should be 'Pooled'."),call. = FALSE)

      if('Composite'%in%GR)
        stop(c("\nInvalid input for 'group'.",
               "\nNo group names should be 'Composite'."),call. = FALSE)


    }else{
      GR = NULL
    }

    ### ZN - zones
    if(!is.null(zones)){
      ZN = df[,zones]
    }else{
      ZN = NULL
    }

    ### RW - replicate weights
    if(is.vector(repwt)){

      if(length(repwt==1)){
        RW <- df[,grepl(repwt,colnames(df)),drop = FALSE]

        if(ncol(RW)<2)
          stop(c("\nInvalid input for 'repwt'.",
                 "\nLess than 2 column names found in 'df'."),call. = FALSE)

        message(paste0(ncol(RW)," weights found."))

      }else{

        if(min(repwt%in%colnames(df))!=1)
          stop(c("\nInvalid input for 'repwt'.",
                 "\nColumns not found in 'df'."),call. = FALSE)

        RW = df[,repwt]

        if(ncol(RW)<2)
          stop(c("\nInvalid input for 'repwt'.",
                 "\nLess than 2 column names found in 'df'."),call. = FALSE)
      }

    }else{

      if(ncol(repwt)<2)
        stop(c("\nInvalid input for 'repwt'.",
               "\nLess than 2 columns."),call. = FALSE)


      if(nrow(repwt)!=nrow(df))
        stop(c("\nInvalid input for 'repwt'.",
               "\nIf a matrix or data frame is provided it should have the same number of rows as 'df'."),call. = FALSE)

      RW <- repwt
    }

  RW <- as.matrix(RW);colnames(RW) <- NULL



  # pooled ----

  if(var!=0){
    outp <- .repmean(X = X, RW = RW,TW = TW,method = method,PV = PV,var = var,
                     group = GR,exclude = exclude,zones = ZN,outrep = FALSE,
                     aggregates = aggregates)
  }else{
    outp <- NA
  }



  if(is.null(by)|(length(x)>1&!PV)){



    out <- outp[[1]]
    class(out) <- c("repmean", class(out))

    if(!"variable"%in%colnames(out)){
      class(out) <- c("repmean.single", class(out))
    }

    return(out)


  }

# by ----

  # message("dfs and pvalues are experimental.")
  bys <- sort(omitna(unique(df[,by])))

  outg <- vector(mode = "list", length = length(bys))

  ugr <- sort(unique(GR))

  for(i in 1:length(bys)){

    I <- bys[i]

    Xii <- as.data.frame(X)[df[,by]%in%I,]
    RWii <- RW[df[,by]%in%I,,drop = FALSE]
    TWii <- TW[df[,by]%in%I]
    GRii <- GR[df[,by]%in%I]
    ZNii <- ZN[df[,by]%in%I]

    # dim(Xii)
    # dim(RWii)
    # length(TWii)
    # length(GRii)
    # length(ZNii)

    # .repmean(X = Xii, RW = RWii,TW = TWii,
    #          method = method,
    #          PV = PV,
    #          group =GRi,
    #          var = var)

    # X = Xii
    # RW = RWii
    # TW = TWii
    # group = GRii
    # zones = ZNii

    outgi <- .repmean(X = Xii, RW = RWii,TW = TWii,
                          method = method,PV = PV,var = var,
                          group = GRii,
                      exclude = exclude,zones = ZNii,outrep = TRUE,
                      aggregates = aggregates)

    if(!is.null(ugr)){
      if(lu(GRii)!=length(ugr)){




        fixon <-  ugr[which(!ugr%in%names(outgi[[2]]))]


        fixo <-  lapply(1:length(fixon), function(i){
          lapply(outgi[[2]][[1]], `is.na<-`)
        })
        names(fixo) <- fixon




        outgi[[1]] <-  rbind(outgi[[1]][1:2,],merge(cbind.data.frame(group = c(ugr)),
                                                    outgi[[1]],sort = TRUE,all.x=TRUE))

        outgi[[1]][is.na(outgi[[1]][,'N']),'N'] <- 0



        outgi[[2]] <-  c(outgi[[2]],fixo)[c('Pooled',ugr)]

      }
    }


    outg[[i]] <- outgi

  }


  kom <- utils::combn(bys,2)

  outd <- lapply(1:ncol(kom),{
    function(k){
      G1 <- outg[[which(bys%in%kom[1,k])]][[-1]]
      G2 <- outg[[which(bys%in%kom[2,k])]][[-1]]


      if(is.null(GR)){


        Difk <- lapply(1:length(G1),function(i){
          G1[[i]]-G2[[i]]
        })


        Difk <- do.call(rbind,lapply(1:1, function(j){
          c(mean(sapply(Difk,function(i) i[1])),
            .repse(er = lapply(Difk,function(i) i[-1]),
                   e0 = sapply(Difk,function(i) i[1]),
                   method = method))

        })
        )



      }else{
        Difk <- lapply(1:length(G1),function(j){
          # Difk <- lapply(1:27,function(j){
          lapply(1:length(G1[[1]]),function(i){
            sw(G1[[j]][[i]]-G2[[j]][[i]])
          })
        })


        Difk <- do.call(rbind,lapply(1:length(Difk), function(j){
          c(mean(sapply(Difk[[j]],function(i) i[1])),
            .repse(er = lapply(Difk[[j]],function(i) i[-1]),
                   e0 = sapply(Difk[[j]],function(i) i[1]),
                   method = method))

        })
        )


        rbind(Difk[1,],
              c(mean(Difk[-1,1][!ugr%in%exclude],na.rm = T),
                .repsecomp(Difk[-1,2][!ugr%in%exclude])),
              Difk[-1,])

      }






    }
  })



  outg <- lapply(1:length(bys),function(k){



    wh <- which(kom[1,]%in%bys[k]|kom[2,]%in%bys[k])

    cbind.data.frame(outg[[k]][[1]],do.call(cbind,lapply(1:length(wh),function(j){


      whj <- wh[j]
      komj <- kom[,whj]

      difj <- outd[[whj]]
      if(bys[k]%in%komj[2]){difj[,1] = -(difj[,1])}

      difj <- cbind(outg[[which(bys==komj[!komj%in%bys[k]])]][[1]][,c('mean')],difj)
      difj <- cbind(difj,difj[,2]/difj[,3])
      dfs <- outg[[k]][[1]][,'N']+outg[[which(bys==komj[!komj%in%bys[k]])]][[1]][,'N']-2
      dfs[dfs<0] <- NA
      difj <- cbind(difj,dfs,2*(stats::pt(q = abs(difj[,4]), df = dfs, lower.tail = FALSE)))

      colnames(difj) <- paste0(c('mean_',"meandiff_","meandiffse_",'tvalue_','df_','pvalue_'),
                               komj[!komj%in%bys[k]])
      difj

    })
    ))


  })


  names(outg) <- paste0(by,'==',bys)


  # Output ----


  out <- c(ALL=outp,outg)

  out <- lapply(out,function(i){
    outi <- i
    class(outi) <- c("repmean", class(outi))
    if(!"variable"%in%colnames(outi)){
      class(outi) <- c("repmean.single", class(outi))
    }
    outi

  })

  class(out) <- c("repmean.list","repmean", class(out))





  # class(out) <- c(class(out),"repmean")
  #
  # if(length(x)==1|PV){
  #   class(out) <- c(class(out),"repmean.single")
  # }

  return(out)



}


# NEED TO BE MODIFIED WHEN repprop is modified for POOL

.repmean0 <- function(x, PV = FALSE, repwt, wt, df,
                    method,
                    var = 0, group = NULL, by = NULL,
                    exclude = NULL, zones = NULL){

  # if(!is.null(setup)){
  #   if(setup$repwt.type!="df"){repwt <- setup$repwt}else{repwt <- get(setup$repwt)}
  #   wt <- setup$wt
  #   method <- setup$method
  #   group <- setup$group
  #   exclude <- setup$exclude
  #   df <- get(setup$df)
  # }

  # Checks ----

  frm <- formals(repmean)




  ## class
  returnis(ischavec,x)
  returnis(islova,PV)
  returnis(is.chavec.or.dfonly,repwt)
  returnis(ischaval,wt)
  returnis(isdf,df)
  if(!isdfonly(df)){

    if(ischavec(repwt)){
      rena = repwt
    }else{
      rena = NULL
    }

    df <- df[,c(group,wt,by,zones,rena)]
    df <- untidy(df)
  }
  returnisNULL(ischaval, group)
  returnisNULL(ischaval, by)
  returnisNULL(ischaval, exclude)
  returnisNULL(ischaval, zones)
  returnis(ischavec, method)



  ## match option
  method <- returnis(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))


  ## Combinations

  ### x & PV
  if(length(x)==1&PV)
    stop(c("\nInvalid input for 'x'.",
           "\nOnly one PV provided."),call. = FALSE)

  ### PV are complete
  if(PV){
    if(!all(rowSums(is.na(df[,x]))%in%c(0,length(x))))
      stop(c("\nInvalid input for 'x'.",
             "\nThere are cases with incomplete PVs."),call. = FALSE)
  }

  ### by for 1 variable
  if(length(x)>1&!PV&!is.null(by))
    stop(paste0("Invalid input for 'by'.",
                "\nIt is only implemented for one variable."),call. = FALSE)


  ### x, wt, group, by, zones in df
  indf <- c(x, wt, group, by, zones)



  if(!all(indf%in%colnames(df)))
    stop(c("\nInvalid input.",
           "\n",
           paste(paste0(indf[!indf%in%colnames(df)],collapse = ', '),"not found in 'df'")),
         call. = FALSE)

  if(length(indf)!=length(unique(indf)))
    stop(paste0("\nInvalid input. Repeated arguments."),
         call. = FALSE)


  ### exclude + group
  if(!is.null(exclude)&&min(!exclude%in%df[,group])==1)
    warning("Check 'exclude', values not found in 'group'.",call. = FALSE)


  ### repwt(df) + df
  if(!is.vector(repwt)&&(nrow(repwt)!=nrow(df)))
    stop(c("\nInvalid input for 'repwt'.",
           "\nIf it is a data frame it should have the same number of rows as 'df'."),call. = FALSE)



  # Process ----

  # Transformation of arguments ----

  ### X - x
  X <- df[,x]
  TW <- df[,wt]

  ### GR - groups
  if(!is.null(group)){
    GR = df[,group]

    if('Pooled'%in%GR)
      stop(c("\nInvalid input for 'group'.",
             "\nNo group names should be 'Pooled'."),call. = FALSE)

    if('Composite'%in%GR)
      stop(c("\nInvalid input for 'group'.",
             "\nNo group names should be 'Composite'."),call. = FALSE)


  }else{
    GR = NULL
  }

  ### ZN - zones
  if(!is.null(zones)){
    ZN = df[,zones]
  }else{
    ZN = NULL
  }

  ### RW - replicate weights
  if(is.vector(repwt)){

    if(length(repwt==1)){
      RW <- df[,grepl(repwt,colnames(df)),drop = FALSE]

      if(ncol(RW)<2)
        stop(c("\nInvalid input for 'repwt'.",
               "\nLess than 2 column names found in 'df'."),call. = FALSE)

      message(paste0(ncol(RW)," weights found."))

    }else{

      if(min(repwt%in%colnames(df))!=1)
        stop(c("\nInvalid input for 'repwt'.",
               "\nColumns not found in 'df'."),call. = FALSE)

      RW = df[,repwt]

      if(ncol(RW)<2)
        stop(c("\nInvalid input for 'repwt'.",
               "\nLess than 2 column names found in 'df'."),call. = FALSE)
    }

  }else{

    if(ncol(repwt)<2)
      stop(c("\nInvalid input for 'repwt'.",
             "\nLess than 2 columns."),call. = FALSE)


    if(nrow(repwt)!=nrow(df))
      stop(c("\nInvalid input for 'repwt'.",
             "\nIf a matrix or data frame is provided it should have the same number of rows as 'df'."),call. = FALSE)

    RW <- repwt
  }

  RW <- as.matrix(RW);colnames(RW) <- NULL



  # pooled ----

  if(var!=0){
    outp <- .repmean(X = X, RW = RW,TW = TW,method = method,PV = PV,var = var,
                     group = GR,
                     exclude = exclude,zones = ZN,outrep = FALSE)
  }else{
    outp <- NA
  }



  if(is.null(by)|(length(x)>1&!PV)){



    out <- outp[[1]]
    class(out) <- c("repmean", class(out))

    if(!"variable"%in%colnames(out)){
      class(out) <- c("repmean.single", class(out))
    }

    return(out)


  }

  # by ----

  message("dfs and pvalues are experimental.")
  bys <- sort(omitna(unique(df[,by])))

  outg <- vector(mode = "list", length = length(bys))

  ugr <- sort(unique(GR))

  for(i in 1:length(bys)){

    I <- bys[i]

    Xii <- as.data.frame(X)[df[,by]%in%I,]
    RWii <- RW[df[,by]%in%I,]
    TWii <- TW[df[,by]%in%I]
    GRii <- GR[df[,by]%in%I]
    ZNii <- ZN[df[,by]%in%I]

    # dim(Xii)
    # dim(RWii)
    # length(TWii)
    # length(GRii)
    # length(ZNii)

    # .repmean(X = Xii, RW = RWii,TW = TWii,
    #          method = method,
    #          PV = PV,
    #          group =GRi,
    #          var = var)

    outgi <- .repmean(X = Xii, RW = RWii,TW = TWii,
                      method = method,PV = PV,var = var,
                      group = GRii,exclude = exclude,zones = ZNii,outrep = TRUE)

    if(!is.null(ugr)){
      if(lu(GRii)!=length(ugr)){




        fixon <-  ugr[which(!ugr%in%names(outgi[[2]]))]


        fixo <-  lapply(1:length(fixon), function(i){
          lapply(outgi[[2]][[1]], `is.na<-`)
        })
        names(fixo) <- fixon




        outgi[[1]] <-  rbind(outgi[[1]][1:2,],merge(cbind.data.frame(group = c(ugr)),
                                                    outgi[[1]],sort = TRUE,all.x=TRUE))

        outgi[[1]][is.na(outgi[[1]][,'N']),'N'] <- 0



        outgi[[2]] <-  c(outgi[[2]],fixo)[c('Pooled',ugr)]

      }
    }


    outg[[i]] <- outgi

  }


  kom <- utils::combn(bys,2)

  outd <- lapply(1:ncol(kom),{
    function(k){
      G1 <- outg[[which(bys%in%kom[1,k])]][[-1]]
      G2 <- outg[[which(bys%in%kom[2,k])]][[-1]]


      if(is.null(GR)){


        Difk <- lapply(1:length(G1),function(i){
          G1[[i]]-G2[[i]]
        })


        Difk <- do.call(rbind,lapply(1:1, function(j){
          c(mean(sapply(Difk,function(i) i[1])),
            .repse(er = lapply(Difk,function(i) i[-1]),
                   e0 = sapply(Difk,function(i) i[1]),
                   method = method))

        })
        )



      }else{
        Difk <- lapply(1:length(G1),function(j){
          lapply(1:length(G1[[1]]),function(i){
            G1[[j]][[i]]-G2[[j]][[i]]
          })
        })


        Difk <- do.call(rbind,lapply(1:length(Difk), function(j){
          c(mean(sapply(Difk[[j]],function(i) i[1])),
            .repse(er = lapply(Difk[[j]],function(i) i[-1]),
                   e0 = sapply(Difk[[j]],function(i) i[1]),
                   method = method))

        })
        )


        rbind(Difk[1,],
              c(mean(Difk[-1,1][!ugr%in%exclude],na.rm = T),
                .repsecomp(Difk[-1,2][!ugr%in%exclude])),
              Difk[-1,])

      }






    }
  })



  outg <- lapply(1:length(bys),function(k){



    wh <- which(kom[1,]%in%bys[k]|kom[2,]%in%bys[k])

    cbind.data.frame(outg[[k]][[1]],do.call(cbind,lapply(1:length(wh),function(j){


      whj <- wh[j]
      komj <- kom[,whj]

      difj <- outd[[whj]]
      if(bys[k]%in%komj[2]){difj[,1] = -(difj[,1])}

      difj <- cbind(outg[[which(bys==komj[!komj%in%bys[k]])]][[1]][,c('mean')],difj)
      difj <- cbind(difj,difj[,2]/difj[,3])
      dfs <- outg[[k]][[1]][,'N']+outg[[which(bys==komj[!komj%in%bys[k]])]][[1]][,'N']-2
      dfs[dfs<0] <- NA
      difj <- cbind(difj,dfs,2*(stats::pt(q = abs(difj[,4]), df = dfs, lower.tail = FALSE)))

      colnames(difj) <- paste0(c('mean_',"meandiff_","meandiffse_",'tvalue_','df_','pvalue_'),
                               komj[!komj%in%bys[k]])
      difj

    })
    ))


  })


  names(outg) <- paste0(by,'==',bys)


  # Output ----


  out <- c(ALL=outp,outg)

  out <- lapply(out,function(i){
    outi <- i
    class(outi) <- c("repmean", class(outi))
    if(!"variable"%in%colnames(outi)){
      class(outi) <- c("repmean.single",class(outi))
    }
    outi

  })

  class(out) <- c("repmean","repmean.list", class(out))





  # class(out) <- c(class(out),"repmean")
  #
  # if(length(x)==1|PV){
  #   class(out) <- c(class(out),"repmean.single")
  # }

  return(out)



}

.repmeanX <- function(X, RW, TW,method,var = 'unbiased',zones = NULL,
                      outrep = FALSE){

  TRW <- cbind(TW,RW)
  RE <- ncol(TRW)

  if(is.null(ncol(RW))){
    RE <- 1+length(RW)
  }



  mod <- ifelse(var=='unbiased',1,0)

  if(is.atomic(X)&&is.vector(X)){
    X <- list(X)
    nzones <- length(unique(zones[!is.na(X)]))
  }else{
    nzones <- length(unique(zones[!is.na(X[[1]])]))
  }

  if(is.null(zones)){nzones = NULL}



  ER <-  lapply(1:length(X),function(j){

    if(length(X[[j]])==1){
      wm <- rep(X[[j]],RE)
    }else{

      wm = as.vector(colSums((X[[j]]*TRW),na.rm = TRUE)/colSums(TRW[!is.na(X[[j]]),,drop = FALSE],na.rm = TRUE))
      # wv = colSums(TRW*(X[[j]]%*%t(rep(1,RE))-t(t(rep(1,length(X[[j]]))))%*%wm)**2,na.rm = TRUE)/(colSums(TRW[!is.na(X[[j]]),,drop = FALSE],na.rm = TRUE)-mod)
    }



    if(var%in%c('ML','unbiased')){

      # if(length(X[[j]])==1&&is.na(X[[j]])){
        if(length(X[[j]])==1){
        wv <- rep(NA,RE)
      }else{
        wv = colSums(TRW*(X[[j]]%*%t(rep(1,RE))-t(t(rep(1,length(X[[j]]))))%*%wm)**2,na.rm = TRUE)/(colSums(TRW[!is.na(X[[j]]),,drop = FALSE],na.rm = TRUE)-mod)
      }


      return(cbind(wm,wv))
    }


    as.matrix(wm)



  })




  me <- mean(sapply(ER,function(y) y[1,1]))

  if(var%in%c('ML','unbiased')){
    va <- mean(sapply(ER,function(y) y[1,2]))
    sd <- mean(sapply(ER,function(y) sw(sqrt(y[1,2]))))

    out <- cbind(N = sum(!is.na(X[[1]])),
                 nzones = nzones,
                 mean = mean(sapply(ER,function(y) y[1,1])),
                 se = .repse(er = lapply(ER,function(y) y[-1,1]),
                             e0 = sapply(ER,function(y) y[1,1]),
                             method = method),
                 sd = mean(sapply(ER,function(y) sw(sqrt(y[1,2])))),
                 sdse = .repse(er = lapply(ER,function(y) sw(sqrt(y[-1,2]))),
                               e0 = (sapply(ER,function(y) sw(sqrt(y[1,2])))),
                               method = method),
                 var = mean(sapply(ER,function(y) y[1,2])),
                 varse = .repse(er = lapply(ER,function(y) (y[-1,2])),
                                e0 = (sapply(ER,function(y) y[1,2])),
                                method = method))
  }else{

    out <- cbind(N = sum(!is.na(X[[1]])),
                 nzones = nzones,
                 mean = mean(sapply(ER,function(y) y[1,1])),
                 se = .repse(er = lapply(ER,function(y) y[-1,1]),
                             e0 = sapply(ER,function(y) y[1,1]),
                             method = method))

  }







  if(!outrep){
    return(list(out))
  }

  c(list(out),list(lapply(ER,function(i) i[,1])))



}




.repmeanG <- function(X, RW, TW,method,var = 'unbiased',
                      group = NULL,exclude = NULL,
                      zones = NULL,outrep = FALSE,
                      aggregates = c("pooled","composite")){

  gopo = "pooled"%in%aggregates
  goco = "composite"%in%aggregates




  if(is.null(group)){

    out <- .repmeanX(X = X,
                     RW = RW,TW = TW,method = method,
                     var = var,zones=zones,outrep = outrep)
    out[[1]] <- as.data.frame(out[[1]])

    return(out)
  }


  ugr <- unique(group)
  ugr <- ugr[order(ugr)]

  if(is.atomic(X)&&is.vector(X)){
    XX <- list(X)
  }else{
    XX <- X
  }







outr <- lapply((abs(gopo-1)):length(ugr),function(i){

    # print(i)

    if(i==0){
      # Xi <- X[!group%in%exclude]

      # Xi <- lapply(1:ncol(XX),function(k) k[!group%in%exclude])
      Xi <- lapply(XX,function(k) k[!group%in%exclude])

      RWi <- RW[!group%in%exclude,]
      TWi <- TW[!group%in%exclude]
      Zi <- zones[!group%in%exclude]
    }else{
      # Xi <- X[group%in%ugr[i]]
      # Xi <- lapply(1:ncol(XX),function(k) k[group%in%ugr[i]])
      Xi <- lapply(XX,function(k) k[group%in%ugr[i]])
      RWi <- RW[group%in%ugr[i],]
      TWi <- TW[group%in%ugr[i]]
      Zi <- zones[group%in%ugr[i]]
    }

    .repmeanX(X = Xi,RW = RWi,TW = TWi,
              method = method,var = var,zones=Zi,
              outrep = outrep)



  })

oute <- do.call(rbind,lapply(outr,function(i) i[[1]]))






  # composite

if(goco){

  if(gopo){
    notexc <- oute[-1,][!ugr%in%exclude,]
  }else{
    notexc <- oute[!ugr%in%exclude,]
  }



  if(var%in%c('ML','unbiased')){
    comp <- c(colMeans(notexc[,c('mean','sd','var')],na.rm = TRUE),
              apply(notexc[,c('se','sdse','varse')],2,function(x){

                .repsecomp(se = x)


              }))[c(1,4,2,5,3,6)]

    comp <- c(rep(NA,ncol(oute)-6),comp)

  }else{

    comp <- c(mean(notexc[,c('mean')],na.rm = TRUE),
              .repsecomp(se = notexc[,c('se')]))

    comp <- c(rep(NA,ncol(oute)-2),comp)

  }


  if(gopo){
    oute <- rbind(oute[1,],comp,oute[-1,])
    rownames(oute) <- NULL

    oute <- cbind.data.frame(group = c('Pooled','Composite',ugr),oute)
  }else{
    cln <- colnames(oute)
    oute <- rbind(comp,oute)
    colnames(oute) <- cln
    rownames(oute) <- NULL

    oute <- cbind.data.frame(group = c('Composite',ugr),oute)
  }


}else{

  if(gopo){

    rownames(oute) <- NULL

    oute <- cbind.data.frame(group = c('Pooled',ugr),oute)
  }else{

    rownames(oute) <- NULL

    oute <- cbind.data.frame(group = c(ugr),oute)
  }

}







  # oute <- rbind(oute[1,],comp,oute[-1,])
  # rownames(oute) <- NULL
  #
  # oute <- cbind.data.frame(group = c('Pooled','Composite',ugr),oute)
  attributes(oute)$excluded = exclude

  if(!outrep){
    return(list(oute))
  }

  outr <- lapply(outr,function(i) i[[-1]])

  if(gopo){
    names(outr) <- c('Pooled',ugr)
  }else{
    names(outr) <- c(ugr)
  }

  c(list(oute),list(outr))


}




.repmean <- function(X,RW,TW,method,PV = FALSE,var = 'unbiased',
                     group = NULL,exclude = NULL,zones = NULL,
                     outrep = FALSE,aggregates = c("pooled","composite")){

  if(PV|is.null(ncol(X))){
    outr <- .repmeanG(X = X,
                     RW = RW,
                     TW = TW,
                     method=method,
                     var=var,
                     zones = zones,
                     group=group,
                     exclude=exclude,
                     outrep = outrep,
                     aggregates = aggregates)
    return(outr)

  }



  nms <- colnames(X)

list(  do.call(rbind,lapply(1:ncol(X),function(y){

  cbind.data.frame(variable = nms[y],  .repmeanG(X = X[,y],
                                                 RW = RW,
                                                 TW = TW,
                                                 method=method,
                                                 var=var,
                                                 zones = zones,
                                                 group=group,
                                                 exclude=exclude,
                                                 outrep = FALSE,
                                                 aggregates = aggregates))

})
))


}
