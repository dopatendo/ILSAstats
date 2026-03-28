


repmeanfast <- function(x,
                        PV = FALSE,
                        setup = NULL,
                        repindex,
                        wt,
                        df,
                        method,
                        var = c("unbiased","ML","none"),
                        group = NULL,
                        by = NULL,
                        exclude = NULL,
                        aggregates = c("pooled", "composite"),
                        simplify = TRUE){

  assignsetup(repmeanfast,setup = setup,mc = match.call())


  var <- var[1L]
  if(var%in%"none"){mod <- -1}
  if(var%in%"unbiased"){mod <- 1}
  if(var%in%"ML"){mod <- 0}



  if(!is.null(group)){
    GR = df[,group, drop = FALSE]
    GR <- do.call(paste,c(GR,sep = "_"))
  }else{
    GR <- NULL
    aggregates <- NULL
    exclude <- NULL
  }

  if(is.null(aggregates)){
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
                                   mod = mod, method = method,
                                   simplify = simplify)

    if(is.null(GR)){
      attributes(estm)$groups <- FALSE
    }else{
      attributes(estm)$groups <- TRUE
      attributes(estm)$excluded <- exclude
    }

    return(estm)
  }


  BY <- df[,by]

  if(!PV){
    estm <- .repmeanfastBYX(GR = GR,
                                   Xs = Xs,
                                   W = W,
                                   BY = BY,
                                   by = by,
                                   exclude = exclude,
                                   aggregates = aggregates,
                                   x = x,
                                   mod = mod, method = method,
                                   repindex = repindex)


    if(is.null(GR)){
      attributes(estm)$groups <- FALSE
    }else{
      attributes(estm)$groups <- TRUE
      attributes(estm)$excluded <- exclude
    }

    return(estm)
  }



  estm <- .repmeanfastBYPV(GR = GR,
                          Xs = Xs,
                          W = W,
                          BY = BY,
                          by = by,
                          exclude = exclude,
                          aggregates = aggregates,
                          x = x,
                          mod = mod, method = method,
                          repindex = repindex)


  if(is.null(GR)){
    attributes(estm)$groups <- FALSE
  }else{
    attributes(estm)$groups <- TRUE
    attributes(estm)$excluded <- exclude
  }

  return(estm)



}









#' @export
summary.repmean.complex <- function(x){

  de <- depth(x)
  gr <- attr(x,"groups")

  if(is.null(gr)){gr <- TRUE}


  if(de==2&&gr){
    xo <- do.call(rbind,lapply(x,function(i){

      unlist(i[setdiff(names(i),c("meanpvs","meanreps","varpvs","varreps"))])


    }))

    xo <- cbind.data.frame(group = rownames(xo),xo)
    class(xo) <- c("repmean.single","repmean",class(xo))
    rownames(xo) <- NULL

  }

  if(de==1){
    xo <- unlist(x[setdiff(names(x),c("meanpvs","meanreps","varpvs","varreps"))])
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

        unlist(i[setdiff(names(i),c("meanpvs","meanreps","varpvs","varreps"))])


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

      unlist(i[setdiff(names(i),c("meanpvs","meanreps","varpvs","varreps"))])


    }))

    xo <- cbind.data.frame(variable = rownames(xo),xo)
    class(xo) <- c("repmean",class(xo))
    rownames(xo) <- NULL

  }


  attr(xo,"groups") <- attr(x,"groups")
  attr(xo,"excluded") <- attr(x,"excluded")

  xo
}


.repmeanfastX <- function(X, W, RWL, add0 = NULL,method,
                                 simplify = TRUE, mod = -1){

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


  ER <- vector("list",mi)
  E0 <- vector("numeric",mi)
  ER2 <- vector("list",mi)
  E02 <- vector("numeric",mi)
  for(i in 1:mi){

    if(is.null(ncol(X))){
      Xi <- X
    }else{
      Xi <- X[,i]
    }


    RWL0 <- RWL[[1]]
    RWL2 <- RWL[[2]]
    me <-  method

    isn <- union(which(is.na(Xi)|is.na(W)),add0)
    # isn <- which(is.na(Xi)|is.na(W))

    LX <- length(Xi)-length(isn)

    Xi[isn] <- 0
    W[isn] <- 0

    XW <- Xi*W
    TS <- sum(XW)
    TW <- sum(W)

    wemean <- TS/TW



    er <- sapply(1:length(RWL0),function(i){
      (TS-sum(XW[RWL0[[i]]])+sum(XW[RWL2[[i]]]))/(TW-sum(W[RWL0[[i]]])+sum(W[RWL2[[i]]]))
    })
    e0 <- wemean


    if(mod>=0){

      mod <- 1

      XWV <- W*(Xi-wemean)**2
      TXWV <- sum(XWV)

      er2 <- sapply(1:length(RWL0),function(i){
        # er[i]
        XWVi <- W*(Xi-er[i])**2
        TXWVi <- sum(XWVi)
        (TXWVi-sum(XWVi[RWL0[[i]]])+sum(XWVi[RWL2[[i]]]))/(TW-sum(W[RWL0[[i]]])+sum(W[RWL2[[i]]])-mod)
      })
      e02 <- TXWV/(TW-1)

      E02[i] <- e02
      ER2[[i]] <- er2

    }



    E0[i] <- e0
    ER[[i]] <- er


  }


  if(length(ER)==1){ER = ER[[1]]}
  if(length(ER2)==1){ER2 = ER2[[1]]}

  SE <- repse(er = ER, e0 = E0, method = me)


  if(mod<0){
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

  if(is.list(ER2)){
    SD2 <- lapply(ER2,sqrt)
  }else{
    SD2 <- sqrt(ER2)
  }

  SEV <- repse(er = ER2, e0 = E02, method = me)
  SED <- repse(er = SD2, e0 = sqrt(E02), method = me)

  if(simplify)
    return(cbind.data.frame(N = LX,
                            mean = mean(E0,na.rm = TRUE),
                            se = SE,
                            sd = mean(sqrt(E02),na.rm = TRUE),
                            sdse = SED,
                            var = mean(E02,na.rm = TRUE),
                            varse = SEV))



  out <- list(N = LX,
              mean = mean(E0,na.rm = TRUE),
              meanpvs = E0,
              se = SE,
              meanreps = ER,
              sd = mean(sqrt(E02),na.rm = TRUE),
              sdse = SED,
              var = mean(E02,na.rm = TRUE),
              varse = SEV,
              varpvs = E02,
              varreps = ER2)

  if(mi==1){
    out <- out[-c(3,10)]
  }

  return(out)




}

.repmeanfastXG <- function(GR = NULL, uGR = NULL,
                                  exclude = NULL, method,
                                  aggregates = NULL,
                                  X, W, RWL,
                                  simplify = TRUE,
                                  mod = -1){

  if(is.null(GR)){

    rpf <- .repmeanfastX(X,
                                W = W,
                                RWL = RWL,
                                add0 = NULL,
                         method = method,
                                simplify = simplify,
                                mod = mod)

    cls <- ifelse(simplify,"repmean.single","repmean.complex")

    class(rpf) <- c(cls,"repmean",class(rpf))

    if(is.null(GR)){
      attributes(rpf)$groups <- FALSE
    }else{
      attributes(rpf)$groups <- TRUE
      attributes(rpf)$excluded <- exclude
    }

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
                                     W = W[coui], method = method,
                                     RWL = spl[[as.character(uGR[i])]],
                                     simplify = simplify, mod = mod)



  }

  if(simplify){
    rpf <- do.call(rbind,rpf)

    if("composite"%in%aggregates){
      rpfi <- rpf[!uGR%in%exclude,]

      if(mod<0){
        rpfi <- c(N = NA,
                  mean = mean(rpfi[,2],na.rm = TRUE),
                  se = repsecomp(rpfi[,3]))
      }else{

        rpfi <- c(N = NA,
                  mean = mean(rpfi[,2],na.rm = TRUE),
                  se = repsecomp(rpfi[,3]),
                  sd = mean(rpfi[,4],na.rm = TRUE),
                  sdse = repsecomp(rpfi[,5]),
                  var = mean(rpfi[,6],na.rm = TRUE),
                  varse = repsecomp(rpfi[,7]))

      }


      rpf <- rbind(rpfi,rpf)
    }
    if("pooled"%in%aggregates){
      add0 <- which(GR%in%exclude)

      rpfi <- .repmeanfastX(X,method = method,
                                   W = W,
                                   RWL = RWL,
                                   add0 = add0,mod = mod,
                                   simplify = simplify)
      rpf <- rbind(rpfi,rpf)
    }
    rpf <- cbind.data.frame(group = c(upperfirst(aggregates),uGR),rpf)
    class(rpf) <- c("repmean.single","repmean",class(rpf))

    if(is.null(GR)){
      attributes(rpf)$groups <- FALSE
    }else{
      attributes(rpf)$groups <- TRUE
      attributes(rpf)$excluded <- exclude
    }

    return(rpf)
  }

  names(rpf) <- uGR

  if("pooled"%in%aggregates){

    NGR <- (!GR%in%exclude)*1

    spl <- splitindex(repindex = RWL, group = NGR)


    coui <- which(NGR%in%1)

    pool <- .repmeanfastX(X = X[coui,,drop = FALSE],
                                 W = W[coui],
                          method = method,
                                 RWL = spl[[as.character(1)]],
                                 simplify = FALSE, mod = mod)
    rpf <- c(list(Pooled = pool),rpf)
  }


  if("composite"%in%aggregates){

    cli <- class(rpf)
    start <- ifelse("pooled"%in%aggregates,2,1)

    rpfex <- rpf[start:length(rpf)]
    rpfex <- rpfex[!names(rpfex)%in%exclude]

    me <- mean(sapply(rpfex,function(i) i$mean),na.rm = TRUE)
    se <- repsecomp(sapply(rpfex,function(i) i$se))


    if(mod<0){
      rpf <- c(rpf[start-1],
               list(Composite = list(N = NA, mean = me, se = se, meanreps = NA)),
               rpf[start:length(rpf)])
      class(rpf) <- cli
    }else{

      sd <- mean(sapply(rpfex,function(i) i$sd),na.rm = TRUE)
      sdse <- repsecomp(sapply(rpfex,function(i) i$sdse))
      var <- mean(sapply(rpfex,function(i) i$var),na.rm = TRUE)
      varse <- repsecomp(sapply(rpfex,function(i) i$varse))

      rpf <- c(rpf[start-1],
               list(Composite = list(N = NA, mean = me, se = se, meanreps = NA,
                                     sd = sd, sdse = sdse, var = var,
                                     varse=varse, varreps = NA)),
               rpf[start:length(rpf)])
      class(rpf) <- cli
    }




  }


  class(rpf) <- c("repmean.complex","repmean",class(rpf))
    if(is.null(GR)){
    attributes(rpf)$groups <- FALSE
  }else{
    attributes(rpf)$groups <- TRUE
    attributes(rpf)$excluded <- exclude
  }
  return(rpf)


}

.repmeanfastXsG <- function(PV = FALSE, method,
                                   GR = NULL, exclude = NULL, aggregates = NULL,
                                   Xs,XsN, W, RWL,uGR = NULL,
                                   simplify = TRUE, mod = -1){

  Xs <- as.matrix(Xs)


  if(ncol(Xs)==1|PV==TRUE){
    rpf <- .repmeanfastXG(GR = GR, uGR = uGR, method = method,
                                 exclude = exclude, aggregates = aggregates,
                                 X = Xs, W = W, RWL = RWL,
                                 simplify = simplify, mod = mod)

    if(is.null(GR)){
      attributes(rpf)$groups <- FALSE
    }else{
      attributes(rpf)$groups <- TRUE
      attributes(rpf)$excluded <- exclude
    }

    return(rpf)
  }


  rpf <- vector("list",ncol(Xs))
  for(i in 1:ncol(Xs)){
    rpfi <- .repmeanfastXG(GR = GR, uGR = uGR, method = method,
                                  exclude = exclude, aggregates = aggregates,
                                  X = Xs[,i,drop = FALSE], W = W, RWL = RWL,
                                  simplify = simplify, mod = mod)
    if(simplify){
      rpfi <- cbind.data.frame(variable = XsN[i],rpfi)
    }

    rpf[[i]] <- rpfi
  }

  if(simplify){
    rpf <- do.call(rbind,rpf)
    class(rpf) <- c("repmean",class(rpf))

    if(is.null(GR)){
      attributes(rpf)$groups <- FALSE
    }else{
      attributes(rpf)$groups <- TRUE
      attributes(rpf)$excluded <- exclude
    }

    return(rpf)
  }

  names(rpf) <- XsN
  class(rpf) <- c("repmean.complex","repmean",class(rpf))

  if(is.null(GR)){
    attributes(rpf)$groups <- FALSE
  }else{
    attributes(rpf)$groups <- TRUE
    attributes(rpf)$excluded <- exclude
  }

  return(rpf)



}

.repmeanfastBYX <- function(GR,Xs,W,by,BY,exclude = NULL, method,
                                   aggregates = NULL,x,repindex,mod = -1){

  ugr <- sort(unique(GR))
  uby <- sort(unique(BY))
  simplify <- FALSE
  PV <- FALSE
  method <- method

  estm <- .repmeanfastXsG(PV = PV,
                                 GR = GR,
                                 uGR = ugr, method = method,
                                 exclude = exclude,
                                 aggregates = setdiff(aggregates,"composite"),
                                 Xs = Xs,
                                 XsN = x,
                                 W = W,
                                 RWL = repindex,
                                 simplify = simplify, mod = mod)

  indi <- splitindex(repindex = repindex, group = BY)
  estmi <- vector("list",length(uby))
  names(estmi) <- uby

  for(i in 1:length(uby)){
    estmi[[i]] <- .repmeanfastXsG(PV = PV, method = method,
                                         GR = GR[BY%in%uby[i]],
                                         uGR = ugr,
                                         exclude = exclude,
                                         aggregates = setdiff(aggregates,"composite"),
                                         Xs = Xs[BY%in%uby[i],,drop = FALSE],
                                         XsN = x,
                                         W = W[BY%in%uby[i]],
                                         RWL = indi[[i]],
                                         simplify = simplify, mod = mod)
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

    if(is.null(GR)){
      attributes(outi)$groups <- FALSE
    }else{
      attributes(outi)$groups <- TRUE
      attributes(outi)$excluded <- exclude
    }


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


      iex <- i[start:nrow(i),]
      iex <- iex[!iex$group%in%exclude,]

      wi <- which(substr(colnames(iex),1,5)%in%c("mean","mean_"))
      ci[wi] <- colMeans(iex[,wi,drop = FALSE],na.rm = TRUE)

      wi <- which(substr(colnames(iex),1,9)%in%c("mean","meandiff_"))
      di <- colMeans(iex[,wi,drop = FALSE],na.rm = TRUE)
      si <- apply(iex[,wi+1,drop = FALSE],2,repsecomp)
      ci[wi+0] <- di
      ci[wi+1] <- si
      ci[wi[-1]+2] <- (di/si)[-1]


      if(mod>=0){
        wi <- which(colnames(iex)%in%c("sd","var"))
        ci[wi] <- colMeans(iex[,wi,drop = FALSE],na.rm = TRUE)
        wi <- which(colnames(iex)%in%c("sdse","varse"))
        ci[wi] <- apply(iex[,wi,drop = FALSE],2,repsecomp)
      }

      outi <- rbind(i[start-1,],
                    ci,
                    i[(start):nrow(i),])
      outi[start,1] <- "Composite"
      rownames(outi) <- NULL

      outi
    })

    class(out) <- cli

  }

  if(is.null(GR)){
    attributes(out)$groups <- FALSE
  }else{
    attributes(out)$groups <- TRUE
    attributes(out)$excluded <- exclude
  }

  return(out)
}

.repmeanfastBYPV <- function(GR,Xs,W,by,BY,exclude = NULL, method,
                             aggregates = NULL,x,repindex,mod = -1){

  ugr <- sort(unique(GR))
  uby <- sort(unique(BY))
  simplify <- FALSE
  PV <- TRUE
  method <- method

  estm <- .repmeanfastXsG(PV = PV,
                          GR = GR,
                          uGR = ugr, method = method,
                          exclude = exclude,
                          aggregates = setdiff(aggregates,"composite"),
                          Xs = Xs,
                          XsN = x,
                          W = W,
                          RWL = repindex,
                          simplify = simplify, mod = mod)

  indi <- splitindex(repindex = repindex, group = BY)
  estmi <- vector("list",length(uby))
  names(estmi) <- uby

  for(i in 1:length(uby)){
    estmi[[i]] <- .repmeanfastXsG(PV = PV, method = method,
                                  GR = GR[BY%in%uby[i]],
                                  uGR = ugr,
                                  exclude = exclude,
                                  aggregates = setdiff(aggregates,"composite"),
                                  Xs = Xs[BY%in%uby[i],,drop = FALSE],
                                  XsN = x,
                                  W = W[BY%in%uby[i]],
                                  RWL = indi[[i]],
                                  simplify = simplify, mod = mod)
  }

  if(is.list(estmi[[1]][[1]])){
    estmii <- lapply(estmi, function(i){
      lapply(i, function(j){
        rbind(j$N,j$meanpvs,do.call(cbind,j$meanreps))
      })
    } )
  }else{
    estmii <- lapply(estmi, function(i){
      # do.call(cbind,lapply(i, function(j){
      rbind(i$N,i$meanpvs,do.call(cbind,i$meanreps))
      # }))
    } )
  }



  i=1
  out <- vector("list",(length(estmii)+1))
  for(i in 1:length(estmii)){


    yesti <- estmii[[i]]
    nesti <- estmii[-i]

    if(!is.list(yesti)){
      yesti <- list(yesti)
    }

    j=1
    # for(j in 1:length(nesti)){}

    nesti <- lapply(1:length(nesti),function(j){

      nestj <- nesti[[j]]

      if(!is.list(nestj)){
        nestj <- list(nestj)
      }

      nestj <- do.call(rbind,    lapply(1:length(yesti),function(k){

        dk <- yesti[[k]]-nestj[[k]]


        outk <- c(mean = mean(nestj[[k]][2,]),
                  meandiff = mean(dk[2,]),
                  meandiffse = .repse(er = c(as.data.frame(dk[-(1:2),])),e0 = dk[2,],method = method),
                  tvalue = NA,
                  df = yesti[[k]][1,1]+nestj[[k]][1,1]-2,
                  pvalue = NA)
        outk

      })
      )

      nestj[,4] <- nestj[,2]/nestj[,3]
      nestj[,6] <- 2*(stats::pt(q = abs(nestj[,4]), df = nestj[,5],lower.tail = FALSE))

      colnames(nestj) <- paste0(colnames(nestj),"_",names(nesti)[j])
      nestj

    })
    sesti <- summary(estmi[[i]])
    nesti <- cbind(sesti,nesti)
    class(nesti) <- class(sesti)

    if(is.null(GR)){
      attributes(nesti)$groups <- FALSE
    }else{
      attributes(nesti)$groups <- TRUE
      attributes(nesti)$excluded <- exclude
    }



    out[[i+1]] <- nesti

  }


  out[[1]] <- summary(estm)
  names(out) <- c("ALL",paste0(by,"==",uby))

  class(out) <- c("repmean.list","repmean",class(out))


  if("composite"%in%aggregates){
    start <- ifelse("pooled"%in%aggregates,2,1)
    cli <- class(out)
    out <- lapply(out,function(i){
      ci <- rep(NA,ncol(i))


      iex <- i[start:nrow(i),]
      iex <- iex[!iex$group%in%exclude,]

      wi <- which(substr(colnames(iex),1,5)%in%c("mean","mean_"))
      ci[wi] <- colMeans(iex[,wi,drop = FALSE],na.rm = TRUE)

      wi <- which(substr(colnames(iex),1,9)%in%c("mean","meandiff_"))
      di <- colMeans(iex[,wi,drop = FALSE],na.rm = TRUE)
      si <- apply(iex[,wi+1,drop = FALSE],2,repsecomp)
      ci[wi+0] <- di
      ci[wi+1] <- si
      ci[wi[-1]+2] <- (di/si)[-1]


      if(mod>=0){
        wi <- which(colnames(iex)%in%c("sd","var"))
        ci[wi] <- colMeans(iex[,wi,drop = FALSE],na.rm = TRUE)
        wi <- which(colnames(iex)%in%c("sdse","varse"))
        ci[wi] <- apply(iex[,wi,drop = FALSE],2,repsecomp)
      }

      outi <- rbind(i[start-1,],
                    ci,
                    i[(start):nrow(i),])
      outi[start,1] <- "Composite"
      rownames(outi) <- NULL

      outi
    })

    class(out) <- cli

  }


  if(is.null(GR)){
    attributes(out)$groups <- FALSE
  }else{
    attributes(out)$groups <- TRUE
    attributes(out)$excluded <- exclude
  }

  return(out)
}

#' @export
print.repmean.complex <- function(x,...){
  attr(x,"groups") <- NULL
  attr(x,"excluded") <- NULL
  print(unclass(x))
}

