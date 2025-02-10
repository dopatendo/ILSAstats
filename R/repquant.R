#' Quantiles with Replicate Weights
#'
#' Estimates quantiles with replicate weights
#' for a variable or a group of variables and for one or more
#' populations.
#'
#' @param qtl a numeric vector indicating the desired quantiles (between 0 and 1).
#' @inheritParams repmean
#'
#' @return a data frame or a list.
#'
#' @example inst/examples/repquant_example.R
#' @export
#'

repquant <- function(x,qtl = c(0.05, 0.25, 0.75, 0.95),PV = FALSE,
                     setup = NULL,
                     repwt, wt, df,
                     method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA","TALIS"),
                     group = NULL,by = NULL, exclude = NULL){




  if(!is.null(setup)){
    if(setup$repwt.type!="df"){repwt <- setup$repwt}else{repwt <- get(setup$repwt)}
    wt <- setup$wt
    method <- setup$method
    group <- setup$group
    exclude <- setup$exclude
    df <- get(setup$df)
  }

  frm <- formals(repquant)
  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = frm$method)


  # Checks ----



  ## class
  returnis(ischavec,x)
  returnis(islova,PV)
  returnis(is.chavec.or.dfonly,repwt)
  returnis(ischaval,wt)
  returnis(isdf,df)
  if(!isdfonly(df)){
    df <- untidy(df)
  }
  returnisNULL(ischaval, group)
  returnisNULL(ischaval, by)
  returnisNULL(ischavec, exclude)
  returnis(ischavec, method)

  returnis(isnumvec,qtl)
  returnis(isnumbet,qtl,from = 0, to = 1)

  ## match option
  method <- returnis(isinvec,x = method[1L],choices = frm$method)

  ## Combinations

  ### PV are complete
  if(PV){
    if(!all(rowSums(is.na(df[,x]))%in%c(0,length(x))))
      stop(c("\nInvalid input for 'x'.",
             "\nThere are cases with incomplete PVs."),call. = FALSE)
  }

  ### x, wt, group, by in df
  indf <- c(x, wt, group, by)

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


  ### x & PV
  if(length(x)==1&PV)
    stop(c("\nInvalid input for 'x'.",
           "\nOnly one PV provided."),call. = FALSE)

  if(length(x)>1&!PV)
    stop("\nThis function can only handle one non-PV variable.",call. = FALSE)


  # Process ----

  ## Transformation of arguments ----

  # X <- as.matrix(df[,x])
  X <- df[,x,drop = FALSE]
  TW <- df[,wt]

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




  if(is.vector(repwt)){

    if(length(repwt==1)){
      RW = df[,grepl(repwt,colnames(df))]
    }else{
      RW = df[,repwt]
    }

  }else{
    RW <- repwt
  }

  ## Quantiles ----


  outt <- .repquant(X = X, PV = PV, RW = RW, TW = TW,method = method,
                    GR = GR, exclude = exclude,qtl = qtl)

  if(is.null(by)){
    return(outt)
  }


  bys <- sort(omitna(unique(df[,by])))

  out <- vector(mode = "list", length = length(bys))

  for(i in 1:length(bys)){

    I <- bys[i]


    Xii <- X[df[,by]%in%I,,drop = FALSE]
    RWii <- RW[df[,by]%in%I,]
    TWii <- TW[df[,by]%in%I]
    GRii <- GR[df[,by]%in%I]

    # .repquant(X = Xii, qtl = qtl, RW = RWii,TW = TWii,
    #           method = method,PV = PV,
    #           GR = GRii,exclude = exclude)



    out[[i]] <- .repquant(X = Xii, qtl = qtl, RW = RWii,TW = TWii,
                          method = method,PV = PV,
                          GR = GRii,exclude = exclude)



  }

  # Output ----

  names(out) <- paste0(by,'==',bys)

  c(list(ALL=outt),out)
}


# Main function (before by)

.repquant <- function(X,qtl = c(0.05, 0.25, 0.75, 0.95),PV = FALSE, RW,TW,
                      method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA"),
                      GR = NULL,
                      exclude = NULL){

  # Example
# X = Xii
# RW = RWii
# TW = TWii
# GR = GRii

  # Process


  if(PV){
    xna <- 'PVs'
  }else{
    xna <- colnames(X)
  }


  if(is.null(GR)){

    return(cbind.data.frame(variable = xna,
                            .repquantX(X = X, RW = RW, TW = TW,PV = PV,
                                       qtl = qtl,method = method)))
  }



  UGR <- sort(unique(GR))


  # i = 1

  raw <- lapply(0:length(UGR),function(i){

    # print(i)



    if(i==0){

      Xi <- (X[!GR%in%exclude,,drop = FALSE])

      RWi <- RW[!GR%in%exclude,]
      TWi <- TW[!GR%in%exclude]

    }else{

      Xi <- (X[GR%in%UGR[i],,drop = FALSE])
      RWi <- RW[GR%in%UGR[i],]
      TWi <- TW[GR%in%UGR[i]]

    }

    .repquantX(X = Xi, RW = RWi,TW = TWi, qtl = qtl,method = method,PV = PV)

  })

  # raw <- do.call(rbind,raw)

  comp <- do.call(rbind,lapply(1L:nrow(raw[[1L]]),function(j){

    # print(j)

    outi <- do.call(rbind,  lapply(2L:length(raw),function(i){
      raw[[i]][j,]
    }))[!UGR%in%exclude,]

    if(ncol(outi)==2L){
      c(mean(outi[,1L],na.rm = TRUE),
        repsecomp(outi[,2L]))
    }else{
      sei <- apply(outi[,1L:length(qtl)*2L],2,function(i) repsecomp(i))
      outi <- colMeans(outi[,1L:length(qtl)*2L-1L],na.rm = TRUE)

      c(rbind(outi,sei))
    }




  })
  )
  colnames(comp) <- colnames(raw[[1L]])

  out <- c(raw[1L],list(comp),raw[-1L])



  out <- do.call(rbind,lapply(1:length(xna),function(j){

    do.call(rbind,lapply(out,function(i) i[j,]))
  }))


  cbind.data.frame(variable = rep(xna,each = (length(UGR)+2)),
                   group = rep(c('Pooled','Composite',c(UGR)),length(xna)),
                   out)


}


# Basic repquant function (no groups)
.repquantX <- function(X,PV = FALSE,RW,TW,qtl = c(0.05, 0.25, 0.75, 0.95),
                       method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA")){
  lqt <- length(qtl)
  TRW <- cbind(TW,RW)


  ER <- lapply(1:ncol(X),function(j){

    do.call(rbind,lapply(1:ncol(TRW),function(i){

      # i = 1; j=1

      .stquantwX(x = X[,j], wt = TRW[,i],qtl = qtl)


    })
    )
  })

  if(!PV){
    SE <- do.call(rbind,lapply(ER,function(j){
      unlist(lapply(1:lqt,function(i){
        .repse(er = j[-1,i],e0 = j[1,i],method = method)
      }))
    })
    )

    E0 <- do.call(rbind,lapply(ER,function(i) i[1,]))

  }else{
    ####### pvs

    ER <- lapply(1:lqt,function(j){
      lapply(ER,function(i){
        i[,j]
      })
    })

    E0 <- sapply(ER,function(j){lapply(j,function(i) i[1])})
    E0 <- colMeans(apply(E0,2,function(i) as.numeric(i)))

    SE <- sapply(ER,function(j){
      .repse(er = lapply(j,function(i) i[-1]),
                     e0 = sapply(j,function(i) i[1]),
                     method = method)
    })

    SE <- as.matrix(t(SE))
    E0 <- as.matrix(t(E0))

  }



  out <- cbind(E0,SE)[,c(t(matrix(1:(lqt*2),ncol=2)))]
  if(!is.matrix(out)){
    out <- t(out)
  }

  colnames(out) <- paste0(rep(paste0("P",
                                     gsub(" ", 0, format(c(0.99, qtl) * 100))[-1]),
                              each = 2),
                          rep(c('','se'),lqt))


  out
}

# Basic quantile function (no repwt)
.stquantwX <- function(x,wt,qtl = c(0.05, 0.25, 0.75, 0.95)){
  xx <- omitna(cbind(x,wt))

  if(nrow(xx)<2){
    return(rep(NA,length(qtl)))
  }

  xx <- xx[order(xx[,1L]),]
  cs <- cumsum(xx[,2L])/sum(xx[,2L])

  unlist(lapply(qtl, function(y) xx[which(cs>=y)[1],1L]),use.names = FALSE)
}
