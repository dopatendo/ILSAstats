#' Generalized Linear Models with Replicate Weights
#'
#' Fits a generalized linear model using \link[stats]{glm} for replicate weights.
#' For a detailed explanation on how the standard errors are estimated
#' see \code{\link{repse}}.
#'
#' @inheritParams repmean
#' @inheritParams stats::glm
#' @inheritParams replm
#'
#' @return a list with the standard errors and the total weights models.
#'
#' @example inst/examples/replm_example.R
#' @export
#'


repglm <- function(formula, family = stats::gaussian,
                   pvs = NULL, relatedpvs = TRUE, quiet = FALSE,summarize = TRUE,
                   setup = NULL,
                  df, wt, repwt,
                  group = NULL, exclude = NULL,
                  na.action = getOption("na.action"),
                  method = c("TIMSS", "PIRLS", "ICILS", "ICCS", "PISA","TALIS")){


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

  frm <- formals(repglm)
  returnis(ischavec, method)
  method <- returnis(isinvec,x = method[1L],choices = frm$method)


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
  # frm <- formals(replm)
  # method <- returnis(isinvec,x = method[1L],choices = frm$method)

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

  estmG <- vector("list",length(UGR)+1)

  for(j in 0:length(UGR)){

    if(j==0){
      if(is.null(group)){
        dfj <- df
        TRWj <- TRW
      }else{
        dfj <- df[GR%in%UGRp,]
        TRWj <- TRW[GR%in%UGRp,]
      }

      if(!quiet&&!is.null(group))
        cat("Estimating pooled model.\n")

    }else{
      dfj <- df[GR%in%UGR[j],]
      TRWj <- TRW[GR%in%UGR[j],]

      if(!quiet)
        cat(paste0("Estimating group ",j," of ",length(UGR),
                   ".\n"))

    }

    estm <- vector("list",length(nfo))
    for(i in 1:length(nfo)){

      if(summarize){strip <- TRUE}else{strip <- !j==0}

      trlm <- try(.repglmX(nfo[[i]],dfj,TRWj,na.action = na.action, strip = strip),silent = TRUE)
      if(!inherits(trlm,"try-error")){
        estm[[i]] <-trlm
        estm[[i]][[1]]$call <- match.call()
      }

    }
    estmG[[j+1]] <- estm

  }

  # estm <- vector("list",length(nfo))
  # for(i in 1:length(nfo)){
  #   estm[[i]] <- .replmX(nfo[[i]],df,TRW,na.action = na.action)
  #   estm[[i]][[1]]$call <- match.call()
  # }



  outj <- vector("list",length(estmG)+3)

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
        outj[[k+3]] <- (structure(list(repmodel = ster, totalmodel = estm[[1]][[1]],
                                       nmodels = nmodels),class = "repglm"))
      }else{
        estm <- lapply(estm,function(i) i[[1]])
        names(estm) <- paste0("totalmodel",1:length(estm))

        outj[[k+3]] <- structure(c(list(repmodel = ster), estm,
                                   nmodels = nmodels),class = "repglm.pv")
      }
    }


  }


  if(length(outj)==4)
    return(outj[[4]])


  # aa = outj
  names(outj) <- c("Composite","Coefficients","StdErrors","Pooled",UGR)


  coes <- lapply(outj[setdiff(names(outj[-(1:4)]),exclude)],function(i) i[[1]][,1:2])
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
  outj[[2]] <- coes
  outj[[3]] <- ster
  ster <- cbind(Estimate = colMeans(coes,na.rm = TRUE),
                "Std. Error" = apply(ster,2,.repsecomp),
                "t value" = NA)
  ster[,3] <- ster[,1]/ster[,2]
  outj[[1]] <- ster


  return(structure(outj,class = "repglm.group"))

}



tr.glm <- function(FO, FM, DF, WT = NULL,NAC = getOption("na.action")){
  FO <- stats::as.formula(FO)
  environment(FO) <- environment()

  if(is.null(WT)){
    WT <- rep(1, nrow(DF))
  }
  environment(WT) <- environment()

  stats::glm(formula = FO, family = FM, data = DF, weights = WT, na.action = NAC)
}



.repglmX <- function(formula, DF, WT,na.action = getOption("na.action"), strip = FALSE){
  out <- vector("list",ncol(WT))
  for(i in 1:ncol(WT)){
    mod <- tr.lm(FO = formula, DF = DF, WT = WT[,i],NAC = na.action)

    if(strip){
      dfr <- mod$df.residual
      mod <- summary(mod)
      mod$residuals <- NULL
      mod$weights <- NULL
      mod$na.action <- NULL
      attributes(mod$terms)$.Environment <- NULL
      coe <- stats::coef(mod)[,1]
      mod$df.residual <- dfr
    }else{
      coe <- stats::coef(mod)
    }

    if(i==1){
      out[[i]] <- (mod)
    }else{
      out[[i]] <- coe
    }

  }
  out

}


#' @export
print.repglm <- function(x, ...){
  cat(paste0(x$nmodels," models were estimated.\n"))
  cat(paste0(rep("-",getOption("width")),collapse = ""),
      "\nReplicate weights' model:\n",
      "\nCall:","\n")
  print(x[[2]]$call)
  cat("\n\nCoefficients:\n")
  print(cbind(round(x[[1]][,1:2,drop = FALSE],4),round(x[[1]][,3:4,drop = FALSE],3)))
  cat("\n")
  cat(paste(rep("-",getOption("width")),collapse = ""),
      "\nTotal weights' model:\n")
  if(inherits(x[[2]],"lm")){
    print(summary(x[[2]]))
  }else{
    print((x[[2]]))
  }


}

#' @export
print.repglm.pv <- function(x, ...){
  cat(paste0(x$nmodels," models were estimated.\n"))
  cat(paste0(rep("-",getOption("width")),collapse = ""),
      "\nReplicate weights' model:\n",
      "\nCall:","\n")
  print(x[[2]]$call)
  cat("\n\nCoefficients:\n")
  print(cbind(round(x[[1]][,1:2,drop = FALSE],4),round(x[[1]][,3:4,drop = FALSE],3)))
  cat("\n")
  cat(paste(rep("-",getOption("width")),collapse = ""),
      "\nTotal weights' model for first plausible value combination:\n")
  if(inherits(x[[2]],"lm")){
    print(summary(x[[2]]))
  }else{
    print((x[[2]]))
  }


}

#' @export
print.repglm.group <- function(x, ...){
  cat(paste0(sum(unlist(lapply(x[-(1:3)],function(i) i$nmodels))),
             " models were estimated.\n"))
  cat(paste0(rep("-",getOption("width")),collapse = ""),
      "\nComposite coefficients:\n")
  # print(x[[2]]$call)
  # cat("\nCoefficients:\n")
  print(cbind(round(x[[1]][,1:2,drop = FALSE],4),round(x[[1]][,3,drop = FALSE],3)))
  cat(paste0(rep("-",getOption("width")),collapse = ""),
      "\n\nCall:\n")

  print(x[[4]][[2]]$call)
  cat("\nCoefficients and standard errors by group:\n")
  coepr <- do.call(cbind,lapply(1:ncol(x[[2]]),function(i){
    out <- paste0(format(round(x[[2]][,i],3)),
                  " (",format(round(x[[3]][,i],2)),")")
    out[is.na(x[[2]][,i])] <- NA
    out
  }))
  colnames(coepr) <- colnames(x[[2]])
  rownames(coepr) <- rownames(x[[2]])
  print(as.table(coepr))

}


