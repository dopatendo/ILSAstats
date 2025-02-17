#' Intraclass Correlation Coefficient
#'
#' Calculates the intraclass correlation coefficient (ICC) fitting a
#' linear mixed-effects model using \link[lme4]{lmer}.
#'
#' @inheritParams repmean
#' @inheritParams lme4::lmer
#' @inheritDotParams lme4::lmer
#' @param x a string vector specifying variable names (within \code{data}).
#' @param group a string specifying the variable name (within \code{data}) to be used for grouping.
#'
#' @return a numeric value or a list.
#'
#' @example inst/examples/icc_example.R
#' @export
#'




icc <- function(x, PV = FALSE, group, data, weights = NULL,...){

  # SUGGESTS
  if(!"lme4"%in%rownames(utils::installed.packages()))
    stop(paste0("\nFor running icc, package 'lme4' is required.",
                "\nPlease, install it."),call. = FALSE)

  icci <- vector("list",length(x))
  for(i in 1:length(x)){

    ficc <- stats::as.formula(paste0(x[i]," ~ 1 + (1|",group,")"))

    WT <- weights

    mo <- lme4::lmer(formula = ficc,data = data,weights = WT,...)
    va <- as.data.frame(lme4::VarCorr(mo))


    icci[[i]] <- va$vcov[1]/sum(va$vcov)
  }
  names(icci) <- x

  if(!PV){
    return(icci)
  }

  return(c(list("Average"=mean(unlist(icci))),icci))


}





