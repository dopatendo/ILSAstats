#' Linear Mixed-Models with Plausible Values
#'
#' Fits a linear mixed-effects model using \link[lme4]{lmer} and plausible values.
#'
#' @inheritParams lme4::lmer
#' @inheritParams replm
#' @inheritParams center
#' @inheritDotParams lme4::lmer
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
#' @param grandmean a character vector indicating the names of columns of
#' \code{data} to which grand-mean should be applied.
#' @param groupmean a character vector indicating the names of columns of
#' \code{data} to which group-mean should be applied.
#'
#' @return a list.
#'
#' @example inst/examples/lmerPV_example.R
#' @export
#'


lmerPV <- function(formula, data = NULL, weights = NULL,
                   pvs, relatedpvs = TRUE,
                   grandmean = NULL, groupmean = NULL, group = NULL, ...){


  # SUGGESTS
  if(!"lme4"%in%rownames(utils::installed.packages()))
    stop(paste0("\nFor running lmerPV, package 'lme4' is required.",
                "\nPlease, install it."),call. = FALSE)

  # SUGGESTS
  if(!"MuMIn"%in%rownames(utils::installed.packages()))
    stop(paste0("\nFor running lmerPV, package 'MuMIn' is required.",
                "\nPlease, install it."),call. = FALSE)




  WT <- weights





  .lmerPV(formula = formula,
          pvs = pvs,
          relatedpvs = relatedpvs,
          df = data,
          grandmean = grandmean,
          groupmean = groupmean,
          group = group,
          WT = WT, CALL = match.call(), ...)
}

.lmerPV <- function(formula, pvs, relatedpvs = TRUE,
                    df = NULL, WT = NULL,
                    grandmean = NULL, groupmean = NULL, group = NULL,
                    CALL = NULL,
                    # TEST = FALSE,
                    # TESTtype = "CR1",
                    ...){

  # # SUGGESTS
  # if(TEST)
  #   if(!"clubSandwich"%in%rownames(utils::installed.packages()))
  #     stop(paste0("\nFor running lmerPV(TEST = TRUE), package 'clubSandwich' is required.",
  #                 "\nPlease, install it."),call. = FALSE)


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

  if(is.null(df)){
    ndf <- df
  }else{
    ndf <- center(df,grandmean = grandmean,groupmean = groupmean,group=group,wt = WT)
  }


  for(i in 1:length(modi)){

    # if(TEST){
      # modi[[i]] <- lme4::lmer(formula = nfo[[i]], data = df, REML = FALSE)
    # }else{
      modi[[i]] <- lme4::lmer(formula = nfo[[i]], data = ndf, weights = WT,...)
      # modi[[i]] <- lmer(formula = nfo[[i]], data = df, weights = wt)
    # }


    modi[[i]]@call$formula <- nfo[[i]]
  }
  names(modi) <- do.call(paste,c(pvs,sep = "_"))
  # modi


  coei <- vector("list",length(modi))
  rani <- vector("list",length(modi))

  for(i in 1:length(modi)){
    sumi <- summary(modi[[i]])

    coei[[i]] <- sumi$coefficients
    reni <- as.data.frame(sumi$varcor)
    rm(sumi)
    # coei[[i]] <- cbind(fixef(modi[[i]]),sqrt(diag(stats::vcov(modi[[i]]))))
    # reni <- as.data.frame(VarCorr(modi[[i]]))
    rani[[i]] <- reni[is.na(reni$var2),];rm(reni)
  }


  coei <- do.call(rbind,coei)

  # for TEST
  # if(TEST){
  #   CR1 <- lapply(modi,function(i) sqrt(diag(clubSandwich::vcovCR(i,type = TESTtype))))
  #   coei[,2] <- do.call(c,CR1)
  # }

  coei <- split.data.frame(coei,f = rownames(coei),drop = FALSE)


  ran <- cbind.data.frame("Variance" = sapply(1:nrow(rani[[1]]),function(i){
    mean(sapply(rani,function(j) j$vcov[i]),na.rm = TRUE)
  }))
  rownames(ran) <- apply(rani[[1]][,1:3],1,function(i) paste(stats::na.omit(i),collapse = "."))


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



  # colnames(coe) <- colnames(coei[[1]])
  if(is.null(CALL)){CALL <- match.call()}

  names(rani) <- names(modi)

  # MM <- m5$models
  r2 <- do.call(rbind,lapply(modi,MuMIn::r.squaredGLMM))
  r2 <- rbind.data.frame(colMeans(r2),r2)
  r2 <- cbind(c("Average",names(modi)),r2)
  colnames(r2) <- c("Model","MarginalR2","ConditionalR2")
  r2

  return(structure(list(fixef = coe,
                        ranef = ran,
                        models = modi,
                        fixefbyX = coei,
                        ranefbymodel = rani,
                        pseudoR2 = r2,
                        call = CALL),
                   class = "lmerPV"))

}

#' @export
print.lmerPV <- function(x, ...){
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
