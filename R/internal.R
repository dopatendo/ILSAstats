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
    return(.untidy(x, mistoNAs = mistoNAs))

  lapply(x,function(i){
    .untidy(i,mistoNAs = mistoNAs)
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


maxdec <- function(x, dec = 5){
  out <- do.call(cbind.data.frame,lapply(x,function(i){

    if(!is.numeric(i))
      return(i)

    round(i, digits = dec)


  }))

  class(out) <- class(x)

  return(out)

}

addcolumn <- function(df, xname, after, x = NA){

  cdf <- colnames(df)

  af <- which(cdf%in%after)

  out <- cbind.data.frame(df[,1:af,drop = FALSE],
                          x,
                          df[,((af+1):ncol(df)),drop = FALSE])
  colnames(out)[af+1] <- xname

  class(out) <- class(df)

  out


}
