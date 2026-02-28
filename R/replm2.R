#' Linear Models with Replicate Weights
#'
#' Fits a linear model using \link[stats]{lm} for replicate weights.
#' For a detailed explanation on how the standard errors are estimated
#' see \code{\link{repse}}.
#'
#' @inheritParams repmean
#' @inheritParams stats::lm
#' @param formula an object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param pvs if plausible values are not used, this should be \code{NULL}.
#' Otherwise it is a list indicating which variables from \code{formula}
#' should be replaced by which plausible values variables. For more details
#' check the examples.
#' @param relatedpvs a logical value indicating if \code{pvs} are drawn
#' from the same model. If \code{TRUE} (default), a total of \eqn{n} estimations will be done,
#' where \eqn{n} is the number of plausible values for each plausible value variable.
#' If \code{FALSE}, a total of \eqn{n_1 \times n_2 \times n_...}
#' estimations will be done, where \eqn{n_i} is the number of plausible values in each plausible value variable.
#' @param quiet a logical value indicating if progress status should be shown
#' while estimating models by group. Default is \code{FALSE}.
#' @param summarize a logical value indicating if \code{lm} objects should be
#' converted to \code{summary.lm} or \code{summary.glm} objects and stripped from certain elements
#' to reduce the size of the output object. Default is \code{TRUE}.
#'
#' @return a list.
#'
#' @example inst/examples/replm_example.R
#' @export
#'


replm2 <- function(formula, pvs = NULL, relatedpvs = TRUE, quiet = FALSE,
                  summarize = TRUE,
                  setup = NULL, df, wt, repwt,
                  group = NULL, exclude = NULL,
                  na.action = getOption("na.action"),
                  method,
                  aggregates = c("pooled", "composite")){


  frm <- formals(replm2)

  # source("R/argchecks.R")
  # source("R/internal.R")

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

  aggregates <- returnisNULL(isinvecmul,x = aggregates, choices = frm$aggregates)

  # Checks ----
  returnis(isformula,formula)
  returnis(isdf,df)
  if(!isdfonly(df)){
    df <- untidy(df)
  }
  returnis(ischaval,wt)
  returnisNULL(islist, pvs)
  # returnis(ischavec, method)
  returnis(islova,relatedpvs)
  returnis(is.chavec.or.dfonly,repwt)

  # ## match option


  ## Consistency of pvs
  if(!is.null(pvs)){
    vfo <- setdiff(all.vars(formula),names(pvs))
    if(relatedpvs)
      if(lu(sapply(pvs,length))>1)
        stop(paste0("\nInvalid input for 'pvs'.",
                    "\nAll elements from list should have the same length."),call. = FALSE)

    if(length(names(pvs))!=length(pvs))
      stop(paste0("\nInvalid input for 'pvs'.",
                  "\nAll elements should be named."),call. = FALSE)

    if(!all(names(pvs)%in%all.vars(formula)))
      stop(paste0("\nInvalid input for 'pvs'.",
                  "\nVariables not found in 'formula'."),call. = FALSE)

    if(!all(unlist(pvs)%in%colnames(df)))
      stop(paste0("\nInvalid input for 'pvs'.",
                  "\nVariables not found in 'df'."),call. = FALSE)
  }else{
    vfo <- all.vars(formula)
  }

  ## formula, wt %in% df
  indf <- c(vfo, wt)

  if(!all(indf%in%colnames(df)))
    stop(c("\nInvalid input.",
           "\n",
           paste(paste0(indf[!indf%in%colnames(df)],collapse = ', '),"not found in 'df'")),
         call. = FALSE)

  if(length(indf)!=length(unique(indf)))
    stop(paste0("\nInvalid input. Repeated arguments."),
         call. = FALSE)

  ### repwt(df) + df
  if(!is.vector(repwt)&&(nrow(repwt)!=nrow(df)))
    stop(c("\nInvalid input for 'repwt'.",
           "\nIf it is a data frame it should have the same number of rows as 'df'."),call. = FALSE)



  # New objects -------------------------------------------------------------

  ## RW - replicate weights ----
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

  ## The rest ----

  WT = df[,wt]
  TRW <- cbind(WT,RW)


  if(is.null(group)){
    GR <- NULL
  }else{
    GR <- df[,group]
  }


  UGR <- unique(GR)
  UGRp <- UGR[!UGR%in%exclude]



  # Formulas with PVs -------------------------------------------------------


  if(is.null(pvs)){
    nfo <- list(formula)
  }else{

    if(relatedpvs){
      pvs <- do.call(cbind.data.frame,pvs)

      nfo <- lapply(1:nrow(pvs),function(i){
        out <- (do.call(substitute,list(formula,as.list(pvs[i,,drop = FALSE]))))
        stats::as.formula(gsub("\"","",deparse(out)))
      })
    }else{
      pvs <- expand.grid(pvs,stringsAsFactors = FALSE)
      nfo <- lapply(1:nrow(pvs),function(i){
        out <- (do.call(substitute,list(formula,c(pvs[i,]))))
        stats::as.formula(gsub("\"","",deparse(out)))
      })
    }


  }


  # Estimation --------------------------------------------------------------


  # estmG <- vector("list",length(UGR)+1)

  if("pooled"%in%aggregates|is.null(UGR)){
    ming <- 0
    estmG <- vector("list",length(UGR)+1)
  }else{
    ming <- 1
    estmG <- vector("list",length(UGR))
  }



  for(j in ming:length(UGR)){

    if(j==0){
      if(is.null(group)){
        dfj <- df
        TRWj <- TRW
        GRj <- GR
      }else{
        dfj <- df[GR%in%UGRp,]
        TRWj <- TRW[GR%in%UGRp,]
        GRj <- GR[GR%in%UGRp]
      }

      if(!quiet&&!is.null(group))
        cat("Estimating pooled model.\n")

    }else{
      dfj <- df[GR%in%UGR[j],]
      TRWj <- TRW[GR%in%UGR[j],]
      GRj <- GR[GR%in%UGR[j]]

      if(!quiet)
        cat(paste0("Estimating group ",j," of ",length(UGR),
                   ".\n"))

    }

    estm <- vector("list",length(nfo))
    for(i in 1:length(nfo)){

      if(summarize){strip <- TRUE}else{strip <- !j==0}

      trlm <- try(.replmX(nfo[[i]],dfj,TRWj,na.action = na.action, strip = strip),silent = TRUE)
      if(!inherits(trlm,"try-error")){
        estm[[i]] <-trlm
        estm[[i]][[1]]$call <- match.call()
      }

    }

    if(ming==0){
      estmG[[j+1]] <- estm
    }else{
      estmG[[j]] <- estm
    }


  }

  # estm <- vector("list",length(nfo))
  # for(i in 1:length(nfo)){
  #   estm[[i]] <- .replmX(nfo[[i]],df,TRW,na.action = na.action)
  #   estm[[i]][[1]]$call <- match.call()
  # }


  # incrn <- sum(c("pooled","composite")%in%aggregates)+1
  # # normally I would need coef+errors, + pooled and composite

  incrN1 <- "composite"%in%aggregates+2

  outj <- vector("list",length(estmG)+incrN1)

  for(k in 1:length(estmG)){
    estm <- estmG[[k]]

    if(!is.null(estm[[1]])){
      # Coefficients ------------------------------------------------------------

      # Total coef
      tc <- do.call(rbind,lapply(estm,function(i){

        if(inherits(i[[1]],"lm")){
          stats::coef(i[[1]])
        }else{
          stats::coef(i[[1]])[,1]
        }


      }))
      if(ncol(tc)==1){colnames(tc) <- "(Intercept)"}

      # Rep coef
      rc <- lapply(estm,function(i) do.call(rbind,i[-1]))


      # Standard errors ---------------------------------------------------------



      ster <- sapply(1:ncol(tc), function(i){
        .repse(er = lapply(rc,function(j) j[,i]), e0 = tc[,i], method = method)
      })
      ster <- cbind(Estimate = (colMeans(tc,na.rm = TRUE)),"Std. Error" = ster,
                    "t value" = NA, "Pr(>|t|)" = NA)
      ster[,3] <- ster[,1]/ster[,2]
      ster[,4] <- 2 * stats::pt(abs(ster[,3]), estm[[1]][[1]]$df.residual,
                                lower.tail = FALSE)


      # Output ------------------------------------------------------------------


      # Output
      nmodels <- sum(sapply(estm,length))
      if(length(estm)==1){
        outj[[k+incrN1]] <- (structure(list(repmodel = ster, totalmodel = estm[[1]][[1]],
                                       nmodels = nmodels),class = "replm"))
      }else{
        estm <- lapply(estm,function(i) i[[1]])
        names(estm) <- paste0("totalmodel",1:length(estm))

        outj[[k+incrN1]] <- structure(c(list(repmodel = ster), estm,
                                   nmodels = nmodels),class = "replm.pv")
      }
    }


  }


  # Output --------------------------------------------------------------

  if(incrN1==3){
    incrN2 <- 4
  }else{
    incrN2 <- 3
  }

  if(length(outj)==incrN2)
    return(outj[[incrN2]])


  noutj <- c("Composite","Coefficients","StdErrors","Pooled")
  noutj <- noutj[tolower(noutj)%in%c(aggregates,c("stderrors","coefficients"))]

  # aa = outj
  # names(outj) <- c("Composite","Coefficients","StdErrors","Pooled",UGR)
  names(outj) <- c(noutj,UGR)


  coes <- lapply(outj[setdiff(names(outj[-(1:length(noutj))]),exclude)],function(i) i[[1]][,1:2])
  coena <- unique(unlist(lapply(coes,rownames)))
  coes <- lapply(coes,function(i){

    if(is.null(i))
      return(i)

    stdi <- setdiff(coena,rownames(i))
    if(length(stdi)==0)
      return(i)

    out <- rbind(i,matrix(NA,
                          nrow = length(stdi),
                          ncol = 2,dimnames = list(stdi,colnames(i))))
    out[coena,]

  })

  ster <- do.call(rbind,lapply(coes,function(i) i[,2]))
  coes <- do.call(rbind,lapply(coes,function(i) i[,1]))
  outj[["composite"%in%aggregates+1]] <- coes
  outj[["composite"%in%aggregates+2]] <- ster

  if("composite"%in%aggregates){
    ster <- cbind(Estimate = colMeans(coes,na.rm = TRUE),
                  "Std. Error" = apply(ster,2,.repsecomp),
                  "t value" = NA)
    ster[,3] <- ster[,1]/ster[,2]
    outj[[1]] <- ster
  }



  return(structure(outj,class = "replm.group2"))

}




#' @export
print.replm.group2 <- function(x, ...){

  if(names(x)[1]=="Composite"){
    cat(paste0(sum(unlist(lapply(x[-(1:3)],function(i) i$nmodels))),
               " models were estimated.\n"))
    cat(paste0(rep("-",getOption("width")),collapse = ""),
        "\nComposite coefficients:\n")
    print(cbind(round(x[[1]][,1:2,drop = FALSE],4),round(x[[1]][,3,drop = FALSE],3)))
  }else{
    cat(paste0(sum(unlist(lapply(x[-(1:2)],function(i) i$nmodels))),
               " models were estimated.\n"))
  }

  wierr <- which(names(x)%in%"StdErrors"[1])
  wicoe <- which(names(x)%in%"Coefficients"[1])



  cat(paste0(rep("-",getOption("width")),collapse = ""),
      "\n\nCall:\n")

  print(x[[wierr+1]][[2]]$call)
  cat("\nCoefficients and standard errors by group:\n")
  coepr <- do.call(cbind,lapply(1:ncol(x[[wicoe]]),function(i){
    out <- paste0(format(round(x[[wicoe]][,i],3)),
                  " (",format(round(x[[wierr]][,i],2)),")")
    out[is.na(x[[wicoe]][,i])] <- NA
    out
  }))
  colnames(coepr) <- colnames(x[[wicoe]])
  rownames(coepr) <- rownames(x[[wicoe]])
  print(as.table(coepr))

}

