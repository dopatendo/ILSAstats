#' Mean with Replicate Weights
#'
#' Estimates only the mean with replicate weights
#' for a variable or a group of variables and for one or more
#' populations. Please be aware, this function is under construction and does not
#' have argument checks yet.
#'
#' For a detailed explanation on how the standard errors are estimated
#' see \code{\link{repse}}.
#'
#' @param repindex a \code{repweights.index} object generate
#' with \code{\link{repcreate}(..., index = TRUE)}.
#' @param simplify a logical value indicating if only the summary
#' statistics should be printed. If \code{FALSE} estimations for
#' all replicated will be provided and no aggregates will
#' be estimated. Default is \code{TRUE}.
#' @inheritParams repmean
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


repmeanfast <- function(x,
                         PV = FALSE,
                         setup = NULL,
                         repindex,
                         wt,
                         df,
                         groups = NULL,
                         by = NULL,
                         exclude = NULL,
                         aggregates = c("pooled", "composite"),
                         simplify = TRUE){



  if(!is.null(groups)){
    GR = df[,groups, drop = FALSE]
    GR <- do.call(paste,c(GR,sep = "_"))
  }else{
    GR <- NULL
    aggregates <- NULL
    exclude <- NULL
  }


  Xs = df[,x, drop = FALSE]
  W = df[,wt]



  if(is.null(by)){


    estm <- .repmeanfastXsG(PV = PV,
                            GR = GR,
                            uGR = NULL,
                            exclude = exclude,
                            aggregates = aggregates,
                            Xs = Xs,
                            XsN = x,
                            W = W,
                            RWL = repindex,
                            simplify = simplify)

    if(is.null(GR)){
      attributes(estm)$groups <- FALSE
    }else{
      attributes(estm)$groups <- TRUE
    }

    return(estm)
  }


  BY <- df[,by]

  if(!PV){
    estm <- .repmeanfastBYX(GR = GR,
                    Xs = Xs,
                    W = W,
                    BY = BY,
                    exclude = exclude,
                    aggregates = aggregates,
                    x = x,
                    repindex = repindex)


    if(is.null(GR)){
      attributes(estm)$groups <- FALSE
    }else{
      attributes(estm)$groups <- TRUE
    }

    return(estm)
  }


}




.repmeanfastX <- function(X, W, RWL, add0 = NULL,
                          simplify = TRUE){

  mi = ncol(X)

  if(nrow(X)==0){

    if(simplify)
      return(cbind.data.frame(N = NA,
                              mean = NA,
                              se = NA))


    out <- list(N = NA,
                mean = NA,
                meanpvs = rep(NA,ncol(X)),
                se = NA,
                meanreps = NA)

    if(mi==1){
      out <- out[-3]
    }

    return(out)

  }

  # X = de$BSMMAT01
  # W = de$TOTWGT
  # RWL = rw2
  # add0 = NULL

  # if(is.null(ncol(X))){
  #   mi <- 1
  # }else{

  # }

  ER <- vector("list",mi)
  E0 <- vector("numeric",mi)
  for(i in 1:mi){

    if(is.null(ncol(X))){
      Xi <- X
    }else{
      Xi <- X[,i]
    }


    RWL0 <- RWL[[1]]
    RWL2 <- RWL[[2]]
    me <- attributes(RWL)$method

    isn <- union(which(is.na(Xi)|is.na(W)),add0)
    # isn <- which(is.na(Xi)|is.na(W))

    LX <- length(Xi)-length(isn)

    Xi[isn] <- 0
    W[isn] <- 0

    XW <- Xi*W
    TS <- sum(XW)
    TW <- sum(W)

    # names(XW) <- ids
    # names(W) <- ids


    er <- sapply(1:length(RWL0),function(i){
      (TS-sum(XW[RWL0[[i]]])+sum(XW[RWL2[[i]]]))/(TW-sum(W[RWL0[[i]]])+sum(W[RWL2[[i]]]))
    })
    e0 <- TS/TW

    E0[i] <- e0
    ER[[i]] <- er
  }


  if(length(ER)==1){ER = ER[[1]]}

  SE <- repse(er = ER, e0 = E0, method = me)

  if(simplify)
    return(cbind.data.frame(N = LX,
                            mean = mean(E0,na.rm = TRUE),
                            se = SE))


  out <- list(N = LX,
              mean = mean(E0,na.rm = TRUE),
              meanpvs = E0,
              se = SE,
              meanreps = ER)

  if(mi==1){
    out <- out[-3]
  }

  return(out)



}


.repmeanfastXG <- function(GR = NULL, uGR = NULL,
                           exclude = NULL,
                           aggregates = NULL,
                              X, W, RWL,
                           simplify = TRUE){

  if(is.null(GR)){

    rpf <- .repmeanfastX(X,
                         W = W,
                         RWL = RWL,
                         add0 = NULL,
                         simplify = simplify)

    cls <- ifelse(simplify,"repmean.single","repmean.complex")

    class(rpf) <- c(cls,"repmean",class(rpf))

    return(rpf)

  }

  if(is.null(uGR)){
    uGR <- sort(unique(GR))
  }

  spl <- splitindex(repindex = RWL, group = GR)
  rpf <- vector("list",length(uGR))
  for(i in 1:length(uGR)){

    coui <- which(GR%in%uGR[i])

    rpf[[i]] <- .repmeanfastX(X = X[coui,,drop = FALSE],
                              W = W[coui],
                              RWL = spl[[as.character(uGR[i])]],
                              simplify = simplify)



  }

  if(simplify){
    rpf <- do.call(rbind,rpf)

    if("composite"%in%aggregates){
      rpfi <- rpf[!uGR%in%exclude,]
      rpfi <- c(N = NA,
                mean = mean(rpfi[,2],na.rm = TRUE),
                se = repsecomp(rpfi[,3]))
      rpf <- rbind(rpfi,rpf)
    }
    if("pooled"%in%aggregates){
      add0 <- which(GR%in%exclude)

      rpfi <- .repmeanfastX(X,
                            W = W,
                            RWL = RWL,
                            add0 = add0,
                            simplify = simplify)
      rpf <- rbind(rpfi,rpf)
    }
    rpf <- cbind.data.frame(group = c(upperfirst(aggregates),uGR),rpf)
    class(rpf) <- c("repmean.single","repmean",class(rpf))

    return(rpf)
  }

  names(rpf) <- uGR

  if("pooled"%in%aggregates){

    NGR <- (!GR%in%exclude)*1

    spl <- splitindex(repindex = RWL, group = NGR)


    coui <- which(NGR%in%1)

   pool <- .repmeanfastX(X = X[coui,,drop = FALSE],
                              W = W[coui],
                              RWL = spl[[as.character(1)]],
                              simplify = FALSE)
    rpf <- c(list(Pooled = pool),rpf)
  }


  if("composite"%in%aggregates){

    cli <- class(rpf)
    start <- ifelse("pooled"%in%aggregates,2,1)

    me <- mean(sapply(rpf[start:length(rpf)],function(i) i$mean),na.rm = TRUE)
    se <- repsecomp(sapply(rpf[start:length(rpf)],function(i) i$se))



    rpf <- c(rpf[start-1],
             list(Composite = list(N = NA, mean = me, se = se, meanreps = NA)),
             rpf[start:length(rpf)])
    class(rpf) <- cli

  }


  class(rpf) <- c("repmean.complex","repmean",class(rpf))
  return(rpf)


}


.repmeanfastXsG <- function(PV = FALSE,
                            GR = NULL, exclude = NULL, aggregates = NULL,
                            Xs,XsN, W, RWL,uGR = NULL,
                            simplify = TRUE){

  Xs <- as.matrix(Xs)


  if(ncol(Xs)==1|PV==TRUE){
    rpf <- .repmeanfastXG(GR = GR, uGR = uGR,
                          exclude = exclude, aggregates = aggregates,
                   X = Xs, W = W, RWL = RWL,
                   simplify = simplify)
    return(rpf)
  }


  rpf <- vector("list",ncol(Xs))
  for(i in 1:ncol(Xs)){
    rpfi <- .repmeanfastXG(GR = GR, uGR = uGR,
                           exclude = exclude, aggregates = aggregates,
                           X = Xs[,i,drop = FALSE], W = W, RWL = RWL,
                           simplify = simplify)
    if(simplify){
      rpfi <- cbind.data.frame(variable = XsN[i],rpfi)
    }

    rpf[[i]] <- rpfi
  }

  if(simplify){
    rpf <- do.call(rbind,rpf)
    class(rpf) <- c("repmean",class(rpf))

    return(rpf)
  }

  names(rpf) <- XsN
  class(rpf) <- c("repmean.complex","repmean",class(rpf))

  return(rpf)



}

.repmeanfastBYX <- function(GR,Xs,W,BY,exclude,aggregates,x,repindex){

  ugr <- sort(unique(GR))
  uby <- sort(unique(BY))
  simplify <- FALSE
  PV <- FALSE

  estm <- .repmeanfastXsG(PV = PV,
                          GR = GR,
                          uGR = ugr,
                          exclude = exclude,
                          aggregates = aggregates,
                          Xs = Xs,
                          XsN = x,
                          W = W,
                          RWL = repindex,
                          simplify = simplify)

  indi <- splitindex(repindex = repindex, group = BY)
  estmi <- vector("list",length(uby))
  names(estmi) <- uby

  for(i in 1:length(uby)){
    estmi[[i]] <- .repmeanfastXsG(PV = PV,
                                  GR = GR[BY%in%uby[i]],
                                  uGR = ugr,
                                  exclude = exclude,
                                  aggregates = aggregates,
                                  Xs = Xs[BY%in%uby[i],,drop = FALSE],
                                  XsN = x,
                                  W = W[BY%in%uby[i]],
                                  RWL = indi[[i]],
                                  simplify = simplify)
  }

  if(is.list(estmi[[1]][[1]])){
    estmii <- lapply(estmi, function(i){
      do.call(cbind,lapply(i, function(j){
        c(j$N,j$mean,j$meanreps)
      }))
    } )
  }else{
    estmii <- lapply(estmi, function(i){
      # do.call(cbind,lapply(i, function(j){
        cbind(c(i$N,i$mean,i$meanreps))
      # }))
    } )
  }



  i=1
  out <- vector("list",(length(estmii)+1))
  for(i in 1:length(estmii)){

    yest <- matrix(rep(estmii[[i]],(length(uby)-1)),nrow = nrow(estmii[[i]]))
    nest <- do.call(cbind,estmii[-i])
    dest <- yest-nest

    sei <- sapply(1:ncol(dest),function(j){
      .repse(er = dest[-(1:2),j], e0 = dest[2,j],method = method)
    })

    ids <- rep(uby[-i],each = ncol(estmii[[1]]))

    mei <- nest[2,]
    mdi <- dest[2,]
    dfi <- yest[1,]+nest[1,]-2
    tvi <- mdi/sei
    pvi <- 2*(stats::pt(q = abs(tvi),
                        df = dfi,
                        lower.tail = FALSE))

    mei <- split(mei,ids)
    mdi <- split(mdi,ids)
    sei <- split(sei,ids)
    tvi <- split(tvi,ids)
    dfi <- split(dfi,ids)
    pvi <- split(pvi,ids)

    outi <- do.call(cbind,lapply(1:length(mei), function(j){

      xo <- t(rbind(mei[[j]],
                    mdi[[j]],
                    sei[[j]],
                    tvi[[j]],
                    dfi[[j]],
                    pvi[[j]]))

      colnames(xo) <- paste0(c("mean",
                               "meandiff",
                               "meandiffse",
                               "tvalue",
                               "df",
                               "pvalue"),"_",uby[-i][j])
      rownames(xo) <- NULL
      xo

    })
    )

    outi <- cbind.data.frame(summary(estmi[[i]]),outi)
    class(outi) <- c("repmean.single","repmean",class(outi))

    out[[i+1]] <- outi

  }


  out[[1]] <- summary(estm)
  names(out) <- c("ALL",paste0(by,"==",uby))

  class(out) <- c("repmean.list","repmean",class(out))


  if("composite"%in%aggregates){
    start <- ifelse("pooled"%in%aggregates,2,1)
    cli <- class(out)
    out <- lapply(out,function(i){
      ci <- rep(NA,ncol(i))
      wi <- which(substr(colnames(i),1,5)%in%c("mean","mean_"))
      ci[wi] <- colMeans(i[start:nrow(i),wi,drop = FALSE],na.rm = TRUE)

      wi <- which(substr(colnames(i),1,9)%in%c("mean","meandiff_"))
      di <- colMeans(i[start:nrow(i),wi,drop = FALSE],na.rm = TRUE)
      si <- apply(i[start:nrow(i),wi+1,drop = FALSE],2,repsecomp)
      ci[wi+0] <- di
      ci[wi+1] <- si
      ci[wi[-1]+2] <- (di/si)[-1]

      outi <- rbind(i[start-1,],
                    ci,
                    i[start:nrow(i),])
      outi[start,1] <- "Composite"
      rownames(outi) <- NULL

      outi
    })

    class(out) <- cli

  }


  return(out)
}



#' @export
summary.repmean.complex <- function(x){

  de <- depth(x)
  gr <- attr(x,"groups")

  if(is.null(gr)){gr <- TRUE}


  if(de==2&&gr){
    xo <- do.call(rbind,lapply(x,function(i){

      unlist(i[setdiff(names(i),c("meanpvs","meanreps"))])


    }))

    xo <- cbind.data.frame(group = rownames(xo),xo)
    class(xo) <- c("repmean.single","repmean",class(xo))
    rownames(xo) <- NULL

  }

  if(de==1&&gr){
    xo <- unlist(x[setdiff(names(x),c("meanpvs","meanreps"))])
    xo <- as.data.frame(t(xo))
    class(xo) <- c("repmean.single","repmean",class(xo))
    rownames(xo) <- NULL
  }

  if(de==3){

    de <- depth(x)
    gr <- attr(x,"groups")
    de
    gr

    xo <- do.call(rbind,lapply(1:length(x),function(j){

      xoi <- do.call(rbind,lapply(x[[j]],function(i){

        unlist(i[setdiff(names(i),c("meanpvs","meanreps"))])


      }))

      xoi <- cbind.data.frame(group = rownames(xoi),xoi)

      # group = rownames(xo)

      xoi <- cbind.data.frame(variable = names(x[j]),xoi)

      xoi
    }))

    class(xo) <- c("repmean",class(xo))
    rownames(xo) <- NULL


  }

  if(de==2&&!gr){

    xo <- do.call(rbind,lapply(x,function(i){

      unlist(i[setdiff(names(i),c("meanpvs","meanreps"))])


    }))

    xo <- cbind.data.frame(variable = rownames(xo),xo)
    class(xo) <- c("repmean",class(xo))
    rownames(xo) <- NULL

  }



  xo
}

#' @export
print.repmean.complex <- function(x,...){
  attr(x,"groups") <- NULL
  print(unclass(x))
}

