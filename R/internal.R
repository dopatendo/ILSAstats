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

  repsi <- lengths(repindex)[[1]]

  groupR <- rep(NA,length(group))
  ugr <- sort(unique(group))
  indix <- vector("list",length(ugr))
  for(i in 1:length(ugr)){
    indix[[i]] <- group%in%ugr[i]
    groupR[indix[[i]]] <- 1:sum(indix[[i]])
  }

  out <- vector("list",length(ugr))
  names(out) <- ugr
  for(i in 1:length(ugr)){

    indi <- which(indix[[i]])

    ugri1 <- vector("list",repsi)
    ugri2 <- vector("list",repsi)
    for(j in 1:repsi){
      ugri1[[j]] <- groupR[repindex[[1]][[j]][group[repindex[[1]][[j]]]%in%ugr[i]]]
      ugri2[[j]] <- groupR[repindex[[2]][[j]][group[repindex[[2]][[j]]]%in%ugr[i]]]
    }

    out[[i]] <- list(ugri1,ugri2)
    attributes(out[[i]]) <- att

  }

  return(out)
}


depth <- function(x){

  if(!is.list(x))
    return(1L)

  d <- 1L
  w <- TRUE
  y <- x
  while(w){
    if(is.list(y[[1]])){
      d <- d+1L
      y <- y[[1]]
    }else{
      w <- FALSE
    }
  }

  return(d)

}


assignsetup <- function(func,setup = NULL,mc){

  if(is.null(setup))
    return(NULL)

  # mc <- (match.call())
  args <- intersect(names(formals(func)),names(formals(repsetup)))
  # empt <- lapply(formals(repsetup),function(x){
  #
  #   x <- nzchar(x)
  #   if(length(x)==0){
  #     TRUE
  #   }else{FALSE}
  #
  #
  # })



  for(i in 1:length(args)){
    argi <- args[i]

    used <- argi%in%names(mc)


    if(!used&&!argi%in%c("df","repwt","repindex")){
      assign(argi,setup[[argi]], envir = parent.frame(n = 1))
    }


    if(!used&&argi=="df"){
      assign(argi,get(setup[[argi]],envir = parent.frame(n = 2)), envir = parent.frame(n = 1))
    }

    if(!used&&argi=="repindex"&!is.null(setup[[argi]])){
      assign(argi,get(setup[[argi]],envir = parent.frame(n = 2)), envir = parent.frame(n = 1))
    }


    if(!used&&argi=="repwt"&!is.null(setup[[argi]])){
      if(setup$repwt.type!="df"){
        assign(argi,setup[[argi]], envir = parent.frame(n = 1))
      }else{
        assign(argi,get(setup[[argi]],envir = parent.frame(n = 2)), envir = parent.frame(n = 1))
      }
    }


  }

}


repweitoindex <- function(RW,W,method){

  reps <- ncol(RW)

  if(tolower(method)%in%c("timss","pirls","lana","jk2-full")){
    reps <- reps/2
  }

  atri <- list(multiplier = c(0,2),
               method = tolower(method),
               reps = reps,
               class = c("repweights.index","repweights","list"))


  out <- list(vector("list",ncol(RW)),vector("list",ncol(RW)))
  W2 <- 2*W
  ptm = proc.time()
  # This finds the '2x' cases column by column
  for(j in 1:ncol(RW)){
    out[[1]][[j]] <- which(RW[, j] == 0)
    out[[2]][[j]] <- which(RW[, j] == (W2))
  }

  attributes(out) <- atri

  return(out)

}
