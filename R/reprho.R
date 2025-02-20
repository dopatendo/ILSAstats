#' Correlations with Replicate Weights
#'
#' Estimates correlation coefficients using replicate weights.
#' For a detailed explanation on how the standard errors are estimated
#' see \code{\link{repse}}.
#'
#' @param x a string vector specifying variable names (within \code{df}) for analysis.
#' If \code{pv} is \code{NULL}, this function estimates correlations between all variables in the vector.
#' If \code{pv2} is NOT \code{NULL}, then \code{x} should be set to \code{NULL}.
#' @param pv a string vector indicating the variable names for all plausible values
#' of a construct. If not \code{NULL}, this function estimates correlations only between
#' \code{x} and the plausible values construct.
#' @param pv2 a string vector indicating the variable names for all plausible values
#' of a second construct (distinct from \code{pv}).
#' @param relatedpvs a logical value indicating if \code{pv} and \code{pv2} are drawn
#' from the same model, and have the same number of plausible values.
#' If \code{TRUE} (default), a total of \eqn{n} estimations will be done,
#' where \eqn{n} is the number of plausible values of each.
#' If \code{FALSE}, a total of \eqn{n_1 \times n_2}
#' estimations will be done, where \eqn{n_1} is the number of plausible values in \code{pv}
#' and \eqn{n_2} is the number of plausible values in \code{pv2}.
#' @param repwt a string indicating the common names for the replicate weights
#' columns within \code{df}, or a data frame with the replicate weights.
#' @param wt a string specifying the name of the column in \code{df} that contains the total weights.
#' @param df a data frame.
#' @param rho a string indicating the correlation coefficient to be computed:
#' \code{"pearson"}, \code{"polychoric"}, or \code{"spearman"} (lower or uppercase).
#' @param method a string indicating the name of the large-scale assessment to
#' determine the replication method to use. Available options are:
#' \code{"TIMSS"}, \code{"PIRLS"}, \code{"ICILS"}, \code{"ICCS"}, and \code{"PISA"}.
#' @param group a string specifying the variable name (within \code{df}) to be used for grouping.
#' @param exclude a vector indicating which groups (in the same format as \code{group})
#' should be excluded from the estimation of pooled and composite estimates.
#' @inheritParams repmean
#'
#' @return a data frame.
#'
#' @example inst/examples/reprho_example.R
#'
#' @export
#'


reprho <- function(x = NULL,pv = NULL, pv2 = NULL,relatedpvs = TRUE,
                   setup = NULL,
                   repwt, wt, df,
                   rho = c("pearson", "spearman", "polychoric"),
                   method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA","TALIS"),
                   group = NULL,exclude = NULL){




  if(!is.null(setup)){
    if(setup$repwt.type!="df"){repwt <- setup$repwt}else{repwt <- get(setup$repwt)}
    wt <- setup$wt
    method <- setup$method
    group <- setup$group
    exclude <- setup$exclude
    df <- get(setup$df)
  }

  frm <- formals(reprho)
  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = frm$method)


  # Checks -----

  ## df
  if(!isdf(df))
    stop(c("\nInvalid input for 'df'.",
           "\nIt should be a data frame."))

  if(isdf(df))
    if(length(class(df))>1){
      df <- untidy(df)
    }


  if(!is.null(pv2)&is.null(pv))
    stop(c("\nInvalid input for 'pv'.",
           "\nIf 'pv2' is used, also 'pv' should be used."))

  # if(relatedpvs&(is.null(pv)|is.null(pv2)))
  #   stop(c("\nInvalid input for 'relatedpvs'.",
  #          "\nIt should only be used with 'pv' and 'pv2'."))

  totalv = length(x)+as.numeric(!is.null(pv))+as.numeric(!is.null(pv2))
  if(totalv<2)
    stop(c("\nNot enough variables."))



  ## Check if they are in df (x,pv,pv2,wt,group)

  indf <- c(x,pv,pv2,wt,group)

  if(min(indf%in%colnames(df))==0)
    stop(c("\nInvalid input.",
           "\n",
           paste(paste0(indf[!indf%in%colnames(df)],collapse = ', '),"not found in 'df'")))

  # Process----


  ## Transformation of arguments ----

  X <- as.matrix(df[,x])
  PV <- as.matrix(df[,pv])
  PV2 <- as.matrix(df[,pv2])
  TW <- as.matrix(df[,wt])

  # if(!is.null(group)){
  #   group = df[,group]
  # }

  if(!is.null(group)){
    GR = df[,group]

    # GR[GR%in%'Pooled'] <- "Pooled (name of group)"
    # GR[GR%in%'Composite'] <- "Composite (name of group)"

    if('Pooled'%in%GR)
      stop(c("\nInvalid input for 'group'.",
             "\nNo group names should be 'Pooled'."))

    if('Composite'%in%GR)
      stop(c("\nInvalid input for 'group'.",
             "\nNo group names should be 'Composite'."))


  }else{
    GR = NULL
  }

  rho <- tolower(rho)

  if(length(method)>1){
    method = method[1]
    message(paste0('More than one method provided. Only first method will be used ',
                   '(',method,').'))
  }

  if(length(rho)>1){
    rho = rho[1]
    message(paste0('More than one correlation method provided. Only first method will be used ',
                   '(',rho,').'))
  }

  # SUGGESTS
  if(!rho%in%"pearson")
    if(length(find.package("wCorr",quiet = T))==0)
      stop(paste0("\nFor running non-Pearson correlations, package 'wCorr' is required.",
                  "\nPlease, install it."),call. = FALSE)






  if((is.atomic(repwt)&is.vector(repwt))){

    if(length(repwt==1)){
      RW = df[,grepl(repwt,colnames(df))]
    }else{
      RW = df[,repwt]
    }

  }else{
    RW <- repwt
  }

  RW <- as.matrix(RW)




  ## Decide which case ----
  # 1 - 2X no PV
  # 2 - XX no PV
  # 3 - 1X w/ PV
  # 4 - XX w/ PV
  # 5 - 2 PVs related
  # 6 - 2 PVs unrelated is case 5



  if(is.null(pv)){

    lx <- length(x)

    if(lx==1)
      stop(c("\nInvalid input for 'x'.",
             "\nOnly 1 variable provided."))

    case <- as.numeric(!lx%in%2)+1

  }else{

    if(!is.null(pv2)){

      if(!is.null(x))
        stop(c("\nInvalid input for 'x'.",
               "\nFor a correlation between PVs, 'x' should be NULL."))

      case <- 5
    }else{
      case <- as.numeric(!length(x)%in%1)+3
    }


  }


  ## 1 - Correlation of two variables, no PV ----
  if(case==1){
    out <- .reprhoXYG(X = X[,1],Y = X[,2],RW = RW,TW = TW, method = method,
                      rho = rho,group = GR,exclude = exclude)

    # w/ groups
    if(is.data.frame(out)){
      out <- cbind(variable1 = colnames(X)[1],variable2 = colnames(X)[2],(out))
    }else{

      # wo/ groups
      out <- cbind.data.frame(colnames(X)[1],colnames(X)[2],
                              t(out))
      colnames(out) <- c('variable1','variable2','rho','se','n')
    }

  }


  ## 2 - Correlation of X variables, no PV ----

  if(case == 2){
    # for non PVs

    kom <- (utils::combn(1:ncol(X),2))


    kor <- lapply(1:ncol(kom),function(i){


      .reprhoXYG(X = X[,kom[1,i]],Y = X[,kom[2,i]],
                 RW = RW,TW = TW, method = method,
                 rho = rho,group = GR,exclude = exclude)

    })


    out <- kor
    # w/ groups
    if(is.data.frame(out[[1]])){
      out <- do.call(rbind,lapply(1:length(out),function(i){

        cbind.data.frame(variable1 = x[kom[1,i]],
                         variable2 = x[kom[2,i]],
                         out[[i]])


      }))

    }else{

      # wo/ groups
      out <- do.call(rbind,lapply(1:length(out),function(i){

        cbind.data.frame(x[kom[1,i]],
                         x[kom[2,i]],
                         t(out[[i]]))


      }))
      colnames(out) <- c('variable1','variable2','rho','se','n')
    }
  }

  ## 3 - Correlation of 1 variable against PVs ----

  if(case==3){

    XX <- as.matrix(X)

    out <- .reprhoXYG(X = XX[,1],Y = PV,RW = RW,TW = TW, method = method,
                      rho = rho,group = GR,exclude = exclude)

    # w/ groups
    if(is.data.frame(out)){
      out <- cbind(variable1 = x,variable2 = 'PVs',(out))
    }else{

      # wo/ groups
      out <- cbind.data.frame(x,'PVs',
                              t(out))
      colnames(out) <- c('variable1','variable2','rho','se','n')
    }
  }

  ## 4 - Correlation of X variables against PVs ----

  if(case==4){
    out <- lapply(1:ncol(X),function(i){

      .reprhoXYG(X = X[,i],Y = PV,RW = RW,TW = TW, method = method,
                 rho = rho,group = GR,exclude = exclude)


    })


    # w/ groups
    if(is.data.frame(out[[1]])){
      out <- do.call(rbind,lapply(1:length(out),function(i){

        cbind.data.frame(variable1 = x[i],
                         variable2 = 'PVs',
                         out[[i]])


      }))

    }else{

      # wo/ groups
      out <- do.call(rbind,lapply(1:length(out),function(i){

        cbind.data.frame(x[i],
                         'PVs',
                         t(out[[i]]))


      }))
      colnames(out) <- c('variable1','variable2','rho','se','n')
    }
  }


  ## 5- Correlation between 2 PVs ----


  if(case==5){
    out <- .reprhoPVG(PV1 = PV, PV2 = PV2,related = relatedpvs, RW = RW,TW = TW, method = method,
                      group = GR,exclude = exclude)

    # w/ groups
    if(is.null(GR)){
      out <- as.data.frame(t(out))
      colnames(out) <- c('rho','se','n')
    }


  }

  # Output ----

  out$tvalue <- c(out$rho/out$se)
  out$pvalue <- sapply(1:nrow(out),function(x){
    2*(stats::pt(q = abs(out$tvalue[x]), df = out$n[x], lower.tail = FALSE))

  })

  out



}


.reprhoXY <- function(X,Y,RW,TW,method,rho = 'pearson'){

  # Y can be a matrix of PV

  TRW <- cbind(TW,RW)
  RE <- ncol(TRW)



  XYW <- stats::na.omit(cbind(X,Y,TRW))
  XY <- XYW[,(1:(ncol(XYW)-RE))]
  WW <- XYW[,-(1:(ncol(XYW)-RE))]
  N <- nrow(XY)

  if(N==0){
    return(rep(NA,RE))
  }


  if(rho=='pearson'){

    SW <- colSums(WW)

    DX <- tcrossprod(XY[,1],rep(1,RE))-tcrossprod(rep(1,N),colSums(XY[,1]*WW)/SW)
    DY <- tcrossprod(XY[,2],rep(1,RE))-tcrossprod(rep(1,N),colSums(XY[,2]*WW)/SW)


    ER <- colSums(DX*DY*WW)/(sqrt(colSums(WW*DX**2))*sqrt(colSums(WW*DY**2)))


    return(c(ER[1],.repse(er = ER[-1],e0 = ER[1],method = method),N))



  }else{
    ER <- do.call(rbind,lapply(1:RE,function(i){



      # # print(i)
      # if(rho=='pearson'){
      #   return(stats::cov.wt(x = XY, wt = WW[,i],cor = TRUE)$cor[-1,1])
      # }


      # sapply(2:ncol(XY),function(j){
      wCorr::weightedCorr(x = XY[,1],
                          y = XY[,2],
                          method = rho,
                          weights = WW[,i])
      # })



    }))


    c(mean(ER[1,]),
      .repse(er = c(as.data.frame(ER[-1,])),
             e0 = as.vector(ER[1,]),
             method = method),nrow(XYW))
  }





}






.reprhoXYG <- function(X,Y,RW,TW,method,rho = 'pearson',
                       group = NULL,exclude = NULL,
                       PV1 = NULL,PV2 = NULL){



  if(is.null(group)){



    return((.reprhoXY(X = X,Y = Y,
                      RW = RW,TW = TW,method = method,
                      rho = rho)))
  }

  Y = as.matrix(Y)

  ugr <- unique(group)
  ugr <- ugr[order(ugr)]



  out <- do.call(rbind,lapply(0:length(ugr),function(i){

    if(i==0){


      # Xi <- lapply(XX,function(x) X[!group%in%exclude])
      Xi <- X[!group%in%exclude]
      Yi <- Y[!group%in%exclude,]
      RWi <- RW[!group%in%exclude,]
      TWi <- TW[!group%in%exclude]

    }else{

      # Xi <- lapply(XX,function(x) X[group%in%ugr[i]])
      Xi <- X[group%in%ugr[i]]
      Yi <- Y[group%in%ugr[i],]
      RWi <- RW[group%in%ugr[i],]
      TWi <- TW[group%in%ugr[i]]

    }

    .reprhoXY(X = Xi,Y = Yi,
              RW = RWi,TW = TWi,method = method,rho = rho)

  })
  )


  notexc <- out[-1,][!ugr%in%exclude,]

  comp <- c(mean(notexc[,1],na.rm = TRUE),.repsecomp(se = notexc[,2]),NA)


  out <- rbind(out[1,],comp,out[-1,])
  rownames(out) <- NULL
  colnames(out) <- c('rho','se','n')

  cbind.data.frame(group = c('Pooled','Composite',ugr),out)
}





.reprhoPV <- function(PV1, PV2, related = TRUE, RW, TW, method,rho = 'pearson'){




  P1 = ncol(PV1)
  P2 = ncol(PV2)

  TRW <- cbind(TW,RW)
  RE <- ncol(TRW)

  if(related){
    grd <- expand.grid(1:P1,1,stringsAsFactors = FALSE)
    grd[,2] <- grd[,1]
  }else{
    grd <- expand.grid(1:P1,1:P2,stringsAsFactors = FALSE)
  }


  ER <- lapply(1:nrow(grd),function(j){



    # xx <- PV1[,grd[j,1]]
    # yy <- PV2[,grd[j,2]]


    XYW <- omitna(cbind(PV1[,grd[j,1]],PV2[,grd[j,2]],TRW))
    XY <- XYW[,(1:(ncol(XYW)-RE))]
    WW <- XYW[,-(1:(ncol(XYW)-RE))]
    N <- nrow(XY)


    if(N==0){
      return(rep(NA,RE))
    }

    # print(i)
    if(rho=='pearson'){

      SW <- colSums(WW)

      DX <- tcrossprod(XY[,1],rep(1,RE))-tcrossprod(rep(1,N),colSums(XY[,1]*WW)/SW)
      DY <- tcrossprod(XY[,2],rep(1,RE))-tcrossprod(rep(1,N),colSums(XY[,2]*WW)/SW)


      return(colSums(DX*DY*WW)/(sqrt(colSums(WW*DX**2))*sqrt(colSums(WW*DY**2))))


      # return(stats::cov.wt(x = XY, wt = WW[,i],cor = TRUE)$cor[-1,1])
    }


    sapply(1:RE,function(i){

      XYW <- stats::na.omit(cbind(XY[,1],XY[,2],WW[,i]))






      # sapply(2:ncol(XY),function(j){
      wCorr::weightedCorr(x = XY[,1],
                          y = XY[,2],
                          method = rho,
                          weights = WW[,i])
      # })



    })


  })

  # ER = e1

  c(mean((sapply(ER,function(i) i[1]))),
    .repse(er = lapply(ER,function(i) i[-1]),
           e0 = (sapply(ER,function(i) i[1])),
           method = method),nrow(stats::na.omit(cbind(PV1,PV2))))
  # ER=e2
  # c(mean((sapply(ER,function(i) i[1]))),
  #   .repse(er = lapply(ER,function(i) i[-1]),
  #          e0 = (sapply(ER,function(i) i[1])),
  #          method = method),nrow(stats::na.omit(cbind(PV1,PV2))))
}



.reprhoPVG <- function(PV1, PV2, related = TRUE, RW,TW,method,rho = 'pearson',
                       group = NULL,exclude = NULL){

  if(is.null(group)){
    return((.reprhoPV(PV1 = PV1,PV2 = PV2, related = related,
                      RW = RW,TW = TW,method = method,
                      rho = rho)))


  }




  ugr <- unique(group)
  ugr <- ugr[order(ugr)]




  out <- do.call(rbind,lapply(0:length(ugr),function(i){

    if(i==0){


      # Xi <- lapply(XX,function(x) X[!group%in%exclude])
      PV1i <- PV1[!group%in%exclude,]
      PV2i <- PV2[!group%in%exclude,]
      RWi <- RW[!group%in%exclude,]
      TWi <- TW[!group%in%exclude]

    }else{

      # Xi <- lapply(XX,function(x) X[group%in%ugr[i]])
      PV1i <- PV1[group%in%ugr[i],]
      PV2i <- PV2[group%in%ugr[i],]
      RWi <- RW[group%in%ugr[i],]
      TWi <- TW[group%in%ugr[i]]

    }

    .reprhoPV(PV1 = PV1i,PV2 = PV2i,
              RW = RWi,TW = TWi,method = method,
              rho = rho, related = related)

  })
  )


  notexc <- out[-1,][!ugr%in%exclude,]

  comp <- c(mean(notexc[,1],na.rm = TRUE),.repsecomp(se = notexc[,2]),NA)


  out <- rbind(out[1,],comp,out[-1,])
  rownames(out) <- NULL
  colnames(out) <- c('rho','se','n')

  cbind.data.frame(group = c('Pooled','Composite',ugr),out)
}


