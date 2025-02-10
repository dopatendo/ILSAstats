#' Centering
#'
#' Centers a vector, a matrix or a data frame to the grand mean or the group mean.
#'
#' @param x a vector, a matrix or a data frame.
#' @param X a matrix or a data frame.
#' @param wt a numeric vector of weights.
#' @param group a vector indicating the group for centering.
#' @param grandmean a numeric or character vector indicating the number or the
#' the name of columns of \code{X} to which grand-mean should be applied.
#' @param groupmean a numeric or character vector indicating the number or the
#' the name of columns of \code{X} to which group-mean should be applied.
#'
#' @return a data frame or a list.
#'
#' @example inst/examples/center_example.R
#' @export
#'


#' @rdname center
#' @export
center <- function(X, group = NULL, grandmean = NULL, groupmean = NULL, wt = NULL){
  # Checks ------------------------------------------------------------------

  ## 1 argument
  returnis(isdf.or.mat,X)
  returnisNULL(isvec,group)
  returnisNULL(isvec,grandmean)
  returnisNULL(isvec,groupmean)
  returnisNULL(isnumvec,wt)

  ## >2 arguments
  returnis2NULL(same.nrow.length,x = X, y = wt)
  returnis2NULL(same.nrow.length,x = X, y = group)

  if(!is.numeric(grandmean)){
    returnis2NULL(isindf,x = X, y = grandmean)
    grandmean <- which(colnames(X)%in%grandmean)
  }else{
    returnisNULL(isnumbet,grandmean,from = 1, to = ncol(X))
  }

  if(!is.numeric(groupmean)){
    returnis2NULL(isindf,x = X, y = groupmean)
    groupmean <- which(colnames(X)%in%groupmean)
  }else{
    returnisNULL(isnumbet,groupmean,from = 1, to = ncol(X))
  }

  if(length(intersect(grandmean,groupmean))>0)
    stop(paste0("\nInvalid input for 'grandmean' or 'groupmean'.",
                "\nOne or more columns are indicated in both arguments."),
         call. = FALSE)


  # Process -----------------------------------------------------------------

  if(length(grandmean)>0){
    X[,grandmean] <- grand.mean(x = X[,grandmean], wt = wt)
  }
  if(length(groupmean)>0){
    X[,groupmean] <- group.mean(x = X[,groupmean], group = group, wt = wt)
  }




  # Output ------------------------------------------------------------------

  return(X)
}

#' @rdname center
#' @export
grand.mean <- function(x, wt = NULL){

  ncx <- ncol(x)
  x <- cbind(x)


# Checks ------------------------------------------------------------------

  returnisNULL(isnumvec,wt)
  returnis2NULL(same.nrow.length,x = x, y = wt)





  # Process -----------------------------------------------------------------


  out <- .wtmean(x = x, wt = wt)

  out <- sweep(x = cbind(x), MARGIN = 2, STATS = out, FUN = "-")


  # Output ------------------------------------------------------------------


  if(is.null(ncx))
    return(as.vector(out))


  return(out)
}

#' @rdname center
#' @export
group.mean <- function(x, group, wt = NULL){


  # Checks ------------------------------------------------------------------
  ncx <- ncol(x)
  x <- cbind(x)

  returnis(isvec,group)
  returnisNULL(isnumvec,wt)


  returnis2NULL(same.nrow.length,x = x, y = wt)
  returnis2(same.nrow.length,x = x, y = group)

  # Process -----------------------------------------------------------------


  ugr <- unique(group)


  i=1
  for(i in 1:length(ugr)){
    coord <- group%in%ugr[i]
    xi <- x[coord,]
    wti <- wt[coord]
    x[coord,] <- grand.mean(x = xi, wt = wti)
  }


  # Output ------------------------------------------------------------------


  if(is.null(ncx))
    return(as.vector(x))


  return(x)
}



.wtmean <- function(x, wt = NULL){
  x <- cbind(x)
  if(is.null(wt)){
    return(colMeans(x,na.rm = TRUE))
  }
  colSums(x*wt,na.rm = TRUE)/colSums((!is.na(x))*wt,na.rm = TRUE)
}


