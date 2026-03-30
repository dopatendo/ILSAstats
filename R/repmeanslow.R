repmeanslow <- function(x,
                        PV = FALSE,
                        setup = NULL,
                        repwt,
                        wt,
                        df,
                        method,
                        var = c("unbiased","ML","none"),
                        group = NULL,
                        by = NULL,
                        exclude = NULL,
                        aggregates = c("pooled", "composite"),
                        simplify = TRUE){

  assignsetup(repmeanslow,setup = setup,mc = match.call())


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


    estm <- .repmeanslowXsG(PV = PV,
                            GR = GR,
                            uGR = NULL,
                            exclude = exclude,
                            aggregates = aggregates,
                            Xs = Xs,
                            XsN = x,
                            W = W,
                            RW = repwt,
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
    estm <- .repmeanslowBYX(GR = GR,
                            Xs = Xs,
                            W = W,
                            BY = BY,
                            by = by,
                            exclude = exclude,
                            aggregates = aggregates,
                            x = x,
                            mod = mod, method = method,
                            repwt = repwt)


    if(is.null(GR)){
      attributes(estm)$groups <- FALSE
    }else{
      attributes(estm)$groups <- TRUE
      attributes(estm)$excluded <- exclude
    }

    return(estm)
  }



  estm <- .repmeanslowBYPV(GR = GR,
                           Xs = Xs,
                           W = W,
                           BY = BY,
                           by = by,
                           exclude = exclude,
                           aggregates = aggregates,
                           x = x,
                           mod = mod, method = method,
                           repwt = repwt)


  if(is.null(GR)){
    attributes(estm)$groups <- FALSE
  }else{
    attributes(estm)$groups <- TRUE
    attributes(estm)$excluded <- exclude
  }

  return(estm)



}



.repmeanslowX <- function(X, W, RW, XRW, XW, add0 = NULL,method,
                          simplify = TRUE, mod = -1){

  mi = ncol(X)

  if(nrow(X)==0){

    if(simplify)
      return(cbind.data.frame(N = NA,
                              mean = NA,
                              se = NA))


    if(mod<0){
      out <- list(N = NA,
                  mean = NA,
                  meanpvs = rep(NA,ncol(X)),
                  se = NA,
                  meanreps = NA)
      if(mi==1){
        out <- out[-3]
      }

    }else{

      out <- list(N = NA,
                  mean = NA,
                  meanpvs = rep(NA,ncol(X)),
                  se = NA,
                  meanreps = NA,
                  sd = NA,
                  sdse = NA,
                  var = NA,
                  varse = NA,
                  varpvs = rep(NA,ncol(X)),
                  varreps = NA)
      if(mi==1){
        out <- out[-c(3,10)]
      }



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



    me <-  method

    # isn <- union(which(is.na(Xi)|is.na(W)),add0)
    #
    #
    #
    #
    # Xi[isn] <- 0
    # RWi <- RW
    # RWi[isn,] <- 0
    # Wi <- W
    # Wi[isn] <- 0
    #
    # XRW <- Xi*RWi

    er <- as.vector(colSums(XRW[[i]])/colSums(RW))
    e0 <- sum(XW[[i]])/sum(W)



    if(mod>=0){

      XRW2 <- RW*(outer(Xi,er,FUN = "-"))**2
      er2 <- colSums(XRW2)/(colSums(RW)-mod)
      e02 <- sum(W*(Xi-e0)**2)/(sum(W)-mod)

      E02[i] <- e02
      ER2[[i]] <- as.vector(er2)

    }



    E0[i] <- e0
    ER[[i]] <- er


  }


  if(length(ER)==1){ER = ER[[1]]}
  if(length(ER2)==1){ER2 = ER2[[1]]}

  SE <- repse(er = ER, e0 = E0, method = me)


  if(mod<0){
    if(simplify)
      return(cbind.data.frame(N = nrow(X),
                              mean = mean(E0,na.rm = TRUE),
                              se = SE))


    out <- list(N = nrow(X),
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
    return(cbind.data.frame(N = nrow(X),
                            mean = mean(E0,na.rm = TRUE),
                            se = SE,
                            sd = mean(sqrt(E02),na.rm = TRUE),
                            sdse = SED,
                            var = mean(E02,na.rm = TRUE),
                            varse = SEV))



  out <- list(N = nrow(X),
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



.repmeanslowXG <- function(GR = NULL, uGR = NULL,
                           exclude = NULL, method,
                           aggregates = NULL,
                           X, W, RW,
                           # XRW, XW,
                           simplify = TRUE,
                           mod = -1){


  isn <- which(is.na(X[,1])|is.na(W))

  X[isn,] <- 0
  RW[isn,] <- 0
  W[isn] <- 0

  XRW <- vector("list",ncol(X))
  XW <-  vector("list",ncol(X))
  for(i in 1:ncol(X)){
    XRW[[i]] <- X[,i]*RW
    XW[[i]] <- X[,i]*W
  }





  if(is.null(GR)){

    rpf <- .repmeanslowX(X,
                         W = W,
                         RW = RW,
                         add0 = NULL,
                         XRW = XRW, XW = XW,
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


  rpf <- vector("list",length(uGR))
  for(i in 1:length(uGR)){

    coui <- which(GR%in%uGR[i])
    XRWi <- lapply(XRW,function(k) k[coui,])
    XWi <- lapply(XW,function(k) k[coui])

    rpf[[i]] <- .repmeanslowX(X = X[coui,,drop = FALSE],
                              W = W[coui], method = method,
                              RW = RW[coui,],
                              XRW = XRWi,XW = XWi,
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
      NGR <- (!GR%in%exclude)*1




      coui <- which(NGR%in%1)
      XRWi <- lapply(XRW,function(k) k[coui,])
      XWi <- lapply(XW,function(k) k[coui])

      rpfi <- .repmeanslowX(X[coui,,drop = FALSE],method = method,
                            W = W[coui],
                            RW = RW[coui,],
                            XRW = XRWi,XW = XWi,
                            add0 = NULL,mod = mod,
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




    coui <- which(NGR%in%1)
    XRWi <- lapply(XRW,function(k) k[coui,])
    XWi <- lapply(XW,function(k) k[coui])

    pool <- .repmeanslowX(X = X[coui,,drop = FALSE],
                          W = W[coui],
                          method = method,
                          RW = RW[coui,],
                          XRW = XRWi,XW = XWi,
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


.repmeanslowXsG <- function(PV = FALSE, method,
                            GR = NULL, exclude = NULL, aggregates = NULL,
                            Xs,XsN, W, RW,uGR = NULL,
                            simplify = TRUE, mod = -1){

  Xs <- as.matrix(Xs)


  if(ncol(Xs)==1|PV==TRUE){



    # isn <- union(which(is.na(Xs[,1])|is.na(W)),add0)
    #
    # Xs[isn,] <- 0
    # Wi <- W
    # Wi[isn] <- 0
    # RWi <- RW
    # RWi[isn,] <- 0
    # XRWi <- Xs*RWi
    # XWi <- Xs*Wi


    rpf <- .repmeanslowXG(GR = GR, uGR = uGR, method = method,
                          exclude = exclude, aggregates = aggregates,
                          X = Xs, W = W, RW = RW,
                          # XRW = XRWi, XW = Xwi,
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

    # isn <- union(which(is.na(Xs[,i])|is.na(W)),add0)
    #
    # Xsi <- Xs[,i,drop = FALSE]
    #
    # Xsi[isn,] <- 0
    # Wi <- W
    # Wi[isn] <- 0
    # RWi <- RW
    # RWi[isn,] <- 0
    # XRWi <- Xsi*RWi
    # XWi <- Xsi*Wi
    #
    rpfi <- .repmeanslowXG(GR = GR, uGR = uGR, method = method,
                           exclude = exclude, aggregates = aggregates,
                           X = Xs[,i,drop = FALSE], W = W, RW = RW,
                           # XRW = XRWi, XW = Xwi,
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


.repmeanslowBYX <- function(GR,Xs,W,by,BY,exclude = NULL, method,
                            aggregates = NULL,x,repwt,mod = -1){

  ugr <- sort(unique(GR))
  uby <- sort(unique(BY))
  simplify <- FALSE
  PV <- FALSE
  method <- method

  estm <- .repmeanslowXsG(PV = PV,
                          GR = GR,
                          uGR = ugr, method = method,
                          exclude = exclude,
                          aggregates = setdiff(aggregates,"composite"),
                          Xs = Xs,
                          XsN = x,
                          W = W,
                          RW = repwt,
                          simplify = simplify, mod = mod)


  estmi <- vector("list",length(uby))
  names(estmi) <- uby

  for(i in 1:length(uby)){
    estmi[[i]] <- .repmeanslowXsG(PV = PV, method = method,
                                  GR = GR[BY%in%uby[i]],
                                  uGR = ugr,
                                  exclude = exclude,
                                  aggregates = setdiff(aggregates,"composite"),
                                  Xs = Xs[BY%in%uby[i],,drop = FALSE],
                                  XsN = x,
                                  W = W[BY%in%uby[i]],
                                  RW = repwt[BY%in%uby[i],],
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

    outi <- cbind.data.frame(summary.repmean.complex(estmi[[i]]),outi)
    class(outi) <- c("repmean.single","repmean",class(outi))

    if(is.null(GR)){
      attributes(outi)$groups <- FALSE
    }else{
      attributes(outi)$groups <- TRUE
      attributes(outi)$excluded <- exclude
    }


    out[[i+1]] <- outi

  }


  out[[1]] <- summary.repmean.complex(estm)
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


.repmeanslowBYPV <- function(GR,Xs,W,by,BY,exclude = NULL, method,
                             aggregates = NULL,x,repwt,mod = -1){

  ugr <- sort(unique(GR))
  uby <- sort(unique(BY))
  simplify <- FALSE
  PV <- TRUE
  method <- method

  estm <- .repmeanslowXsG(PV = PV,
                          GR = GR,
                          uGR = ugr, method = method,
                          exclude = exclude,
                          aggregates = setdiff(aggregates,"composite"),
                          Xs = Xs,
                          XsN = x,
                          W = W,
                          RW = repwt,
                          simplify = simplify, mod = mod)


  estmi <- vector("list",length(uby))
  names(estmi) <- uby

  for(i in 1:length(uby)){
    estmi[[i]] <- .repmeanslowXsG(PV = PV, method = method,
                                  GR = GR[BY%in%uby[i]],
                                  uGR = ugr,
                                  exclude = exclude,
                                  aggregates = setdiff(aggregates,"composite"),
                                  Xs = Xs[BY%in%uby[i],,drop = FALSE],
                                  XsN = x,
                                  W = W[BY%in%uby[i]],
                                  RW = repwt[BY%in%uby[i],],
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
    sesti <- summary.repmean.complex(estmi[[i]])
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


  out[[1]] <- summary.repmean.complex(estm)
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

