#' Survey Weighted Mixed-Effects Models with Plausible Values
#'
#' Fits a linear mixed-effects model using \link[WeMix]{mix} and plausible values.
#'
#' @inheritParams WeMix::mix
#' @inheritDotParams WeMix::mix
#' @param pvs a list indicating which variables from \code{formula}
#' should be replaced by which plausible values variables. For more details
#' check the examples.
#' @param relatedpvs a logical value indicating if \code{pvs} are drawn
#' from the same model, and have the same number of plausible values.
#' If \code{TRUE} (default), a total of \eqn{n} estimations will be done,
#' where \eqn{n} is the number of plausible values for each plausible value variable.
#' If \code{FALSE}, a total of \eqn{n_1 \times n_2 \times n_...}
#' estimations will be done, where \eqn{n_i} is the number of plausible values
#' in each plausible value variable.
#'
#' @return a list.
#'
#' @example inst/examples/WeMixPV_example.R
#' @export
#'



WeMixPV <- function(formula, data = NULL, weights = NULL,
                     pvs, relatedpvs = TRUE,...){


  # SUGGESTS
  if(!"WeMix"%in%rownames(utils::installed.packages()))
    stop(paste0("\nFor running WeMixPV, package 'WeMix' is required.",
                "\nPlease, install it."),call. = FALSE)



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



  modi <- vector("list",length(nfo))


    ndf <- untidy(data)



  for(i in 1:length(modi)){

    modi[[i]] <- (WeMix::mix(formula = nfo[[i]],data = ndf, weights = weights,...))

    # modi[[i]] <- lme4::lmer(formula = nfo[[i]], data = ndf, weights = WT,...)

    modi[[i]]$call$formula <- nfo[[i]]

  }
  names(modi) <- do.call(paste,c(pvs,sep = "_"))
  # sumi <- lapply(modi,summary)



  coei <- vector("list",length(modi))
  rani <- vector("list",length(modi))



  for(i in 1:length(modi)){
    sumi <- summary(modi[[i]])

    coei[[i]] <- as.matrix(sumi$coef)
    rani[[i]] <- as.matrix(sumi$vars)
    # reni <- as.data.frame(sumi$varcor)
    # rm(sumi)
    #
    # rani[[i]] <- reni[is.na(reni$var2),];rm(reni)
  }


  coei <- do.call(rbind,coei)
  coei <- split.data.frame(coei,f = rownames(coei),drop = FALSE)
  coe <- do.call(rbind,lapply(coei,function(i){
    sedf <- pvse(i[,2],i[,1],df = TRUE)
    esti <- mean(i[,1],na.rm = TRUE)
    cbind.data.frame("Estimate" = esti,
                     "Std. Error" = sedf[1],
                     "t value" = esti/sedf[1],
                     "df" = sedf[2],
                     "Pr(>|t|)" = 2 * stats::pt(abs(esti/sedf[1]),
                                                df = sedf[2],
                                                lower.tail = FALSE))

  })
  )

  rani0 <- rani
  names(rani0) <- names(modi)
  rani <- do.call(rbind,rani)
  rani <- split.data.frame(rani,f = rownames(rani),drop = FALSE)
  ran <- do.call(rbind,lapply(rani,function(i){
    sedf <- pvse(i[,2],i[,1],df = TRUE)
    esti <- mean(i[,1],na.rm = TRUE)
    cbind.data.frame("Estimate" = esti,
                     "Std. Error" = sedf[1],
                     "t value" = esti/sedf[1],
                     "df" = sedf[2],
                     "Pr(>|t|)" = 2 * stats::pt(abs(esti/sedf[1]),
                                                df = sedf[2],
                                                lower.tail = FALSE))

  })
  )


  CALL <- match.call()
  # CALL <- 1

  structure(list(fixef = coe,
                        ranef = ran,
                        models = modi,
                        fixefbyX = coei,
                        ranefbymodel = rani0,
                        # pseudoR2 = r2,
                        call = CALL),

                   class = "WeMixPV")


}


#' @export
print.WeMixPV <- function(x, ...){
  cat(paste0(length(x$models)," models were estimated.\n"))
  cat(paste0(rep("-",getOption("width")),collapse = ""),
      "\nMultilevel results with PVs:\n",
      "\nCall:","\n")
  print(x$call)
  cat("\nRandom effects:\n")
  print(x$ranef)
  cat("\nFixed effects:\n")
  print(x$fixef)

  cat(paste0(rep("-",getOption("width")),collapse = ""),
      "Estimated models:",sapply(x$models,function(i) deparse(stats::formula(i))),sep = "\n")

}




