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

upperfirst <- function(x){



  paste0(toupper(substr(x,1,1)),tolower(substring(x,2)))
}


splitindexold <- function(repindex, group){
  spl <- split(1:length(group),f = group)


  out <- lapply(spl,function(i){
    outi <- list(lapply(repindex[[1]],function(j) omitna(match(j,i))),
                 lapply(repindex[[2]],function(j) omitna(match(j,i))))

    attributes(outi) <- attributes(repindex)
    outi

  })

  out

}




splitindex <- function(repindex,group){
  att <- attributes(repindex)

  groupR <- rep(NA,length(group))
  ugr <- sort(unique(group))
  for(i in 1:length(ugr)){
    indi <- group%in%ugr[i]
    groupR[indi] <- 1:sum(indi)
  }

  out <- vector("list",length(ugr))
  names(out) <- ugr
  for(i in 1:length(ugr)){

    indi <- which(group%in%ugr[i])

    ugri1 <- vector("list",att$reps)
    ugri2 <- vector("list",att$reps)
    for(j in 1:att$reps){
      ugri1[[j]] <- groupR[repindex[[1]][[j]][group[repindex[[1]][[j]]]%in%ugr[i]]]
      ugri2[[j]] <- groupR[repindex[[2]][[j]][group[repindex[[2]][[j]]]%in%ugr[i]]]
    }

    out[[i]] <- list(ugri1,ugri2)
    attributes(out[[i]]) <- att

  }

  return(out)
}


depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)
