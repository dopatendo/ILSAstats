#' Confidence Intervals for Replicated Means
#'
#'
#' Calculates the confidence intervals for a \code{\link{repmean}} object.
#'
#'
#' @param x an object produced by \code{\link{repmean}}.
#' @param alpha a numeric value indicating confidence level.
#' @param add a logical value indicating if the confidence intervals should be added
#' to the object or not. Defaults is \code{TRUE}.
#'
#'
#' @return a data frame or a list.
#'
#' @example inst/examples/repmeanCI_example.R
#' @export
#'


repmeanCI <- function(x, alpha = 0.05, add = TRUE){


  returnis(isrep.mean,x)

  returnis(isnumbet,alpha,from = 0, to = 1)
  returnis(islova, add)


  if(inherits(x,"repmean.list")){
    X <- lapply(x,function(i){
      .repmeanCI(x = i, alpha = alpha, add = add)
    })

    class(X) <- c("repmeanCI.list","repmeanCI","list")
  }else{
    X <- .repmeanCI(x = x, alpha = alpha, add = add)
    class(X) <- c("repmeanCI","data.frame")
  }




  return(X)

}

.repmeanCI <- function(x, alpha = 0.05, add = TRUE){

  M <- qnorm(1-alpha/2)
  up <- x$mean+M*x$se
  do <- x$mean-M*x$se


  if(add){
    br <- which(colnames(x)%in%"se")

    if(br==ncol(x)){
      X <- cbind.data.frame(x[,c(1:br)],
                            CIdown = do,
                            CIup = up)
    }else{
      X <- cbind.data.frame(x[,c(1:br)],
                            CIdown = do,
                            CIup = up,
                            x[,c((br+1):ncol(x))])
    }




  }else{
    X <- cbind.data.frame(x[,colnames(x)%in%"group",drop = FALSE],
                          CIdown = do,
                          CIup = up)
  }

  class(X) <- c("repmeanCI","data.frame")

  return(X)

}


#' @export
print.repmeanCI <- function(x, ...){

  class(x) <- setdiff(class(x),c("repmeanCI","repmeanCI.list"))

  print(x)

}


