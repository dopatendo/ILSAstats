#' Internal functions
#'
#' @param x an R object.
#'
#' @noRd




lu <- function(x){
  length(unique(x[!is.na(x)]))
}

sm <- function(x){

  suppressMessages(x)

}

sw <- function(x){


  suppressWarnings(x)}


omitna <- function(x){

  if(is.vector(x)|is.factor(x))
    return(x[!is.na(x)])


  x[rowSums(is.na(x))<1,,drop = FALSE]

}


untidy <- function(x, mistoNAs = TRUE){

  if(is.data.frame(x)|is.matrix(x))
    return(.untidy2(x, mistoNAs = mistoNAs))

  lapply(x,function(i){
    .untidy2(i,mistoNAs = mistoNAs)
  })

}

.untidy <- function(x, mistoNAs = TRUE){
  out <- x
  out <- lapply(1:ncol(x),function(X){

    as.vector(out[,X,drop = TRUE])

  })

  out <- do.call(cbind.data.frame,out)
  colnames(out) <- colnames(x)
  if(mistoNAs){
    out[is.na(x)] <- NA
  }
  out

}


