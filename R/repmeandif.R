#' Mean Difference of Independent Samples with Replicate Weights
#'
#'
#' Estimates the mean difference for a single variable with replicate weights.
#' For a detailed explanation on how the standard errors are estimated
#' see \code{\link{repse}}.
#'
#'
#' @param x a data frame produced by \code{\link{repmean}} for a single variable
#' or an object produced by \code{\link{leaguetable}}.
#'
#'
#' @return a data frame or a list.
#'
#' @example inst/examples/repmeandif_example.R
#' @export
#'


repmeandif <- function(x){

  if(inherits(x,"leaguetable")){

    return(.repmeandifLT(x))

  }


  returnis(isrep.mean,x)


  if(is.data.frame(x)){
    # message("dfs and pvalues are experimental.")
    return(.repmeandif(x))
  }


  # message("dfs and pvalues are experimental.")
  out <- lapply(x,function(i){
    .repmeandif(i)
  })

  class(out) <- c("repmeandif.list","repmeandif",class(out))

  out



}

.repmeandif <- function(x){

  returnis(isrep.meansingle,x)

  goco = "Composite"%in%unique(c(x$group1,x$group2))

  # returnis(isrepmean,x)

  # if(min(c('N','mean','se')%in%colnames(x))==0)
  #   stop('Invalid input for x.')
  #
  # if(!is.data.frame(x))
  #   stop('Invalid input for x.')
  #
  # if('variable'%in%colnames(x))
  #   if(lu(unlist(x$variable))>1)
  #     stop('Invalid input for x. Please use only one variable.')


  # zou <- x[!x[,1]%in%'ALL',]
  zou <- x
  zou2 <- x[!x$group%in%"Composite",]
  zou3 <- x[!x$group%in%c("Pooled","Composite"),]



  if(goco){
    # composite
    exc <- attr(x,'excluded')
    # grs <- zou2[-1,]
    grs <- zou3
    grs <- grs[!grs$group%in%exc,]
    grs <- grs[!is.na(grs$mean),]
    mult <- ((nrow(grs))**2-1)/(nrow(grs)**2)
    comp <- sqrt(x[2,'se']**2+mult*grs[,'se']**2)
  }



  # non composite

  if('variable'%in%colnames(x)){
    grn <- grep('variable',colnames(x))-1
  }else{
    grn <- grep('N',colnames(x))-1
  }

  if(grn==0)
    stop('Invalid input for x. No groups found.',call. = FALSE)


  mdif <- round(sapply(zou$mean,function(x) x-unlist(zou$mean)),5)

  mser <- round(sapply(zou$se,function(x){
    sqrt(unlist(zou$se)**2+x**2)
  }),5)

  ddff <- sapply(zou$N,function(x){
    (unlist(zou$N)-1)+(x-1)
  })

  mdif <- cbind.data.frame(dif = c(mdif),
                           se = c(mser),
                           tvalue = c(mdif/mser),
                           df = c(ddff))




  if(grn==1){
    nam <- unlist(zou$group)
  }else{
    nam <- apply(zou[,names(zou)[1:grn]],1,function(x) paste0(x,collapse = '_'))
    nam <- gsub(' ','',nam)
  }

  nam <- cbind.data.frame(rep(nam,each=nrow(zou)),rep(nam,nrow(zou)))
  colnames(nam) <- c('group1','group2')




  pv <- sapply(1:nrow(mdif),function(x){
    2*(stats::pt(q = abs(mdif$tvalue[x]), df = mdif$df[x], lower.tail = FALSE))

  })

  mdif <- cbind.data.frame(nam,mdif,pvalue = round(pv,5))


  if(goco){
    mdiftot <- mdif[mdif$group1!='Composite',]
    mdifcom <- mdif[mdif$group1=='Composite',]

    # group 1

    mcom <- mdifcom[mdifcom$group2%in%grs$group,]

    mcom$se <- comp
    mcom$tvalue <- mcom$dif/comp

    out <- rbind.data.frame(mcom,mdiftot)

    # group 2

    ## remove non comparable

    out = out[!((!out$group1%in%grs$group)&out$group2=='Composite'),]

    ## add group2
    out[out$group2=='Composite',"se"] <- comp
    out[out$group2=='Composite',"tvalue"] <- mcom$dif/comp
  }else{
    out = mdif
  }






  rownames(out) <- NULL
  class(out) <- c("repmeandif", class(out))
  out


}


.repmeandifLT <- function(x){

  xx <- split(x,f=x$subject)

  xx <- do.call(rbind,lapply(xx,function(i){
    x1 <- unique(i[,setdiff(colnames(x),c("group","mean","se","N","CIup","CIdown"))])
    x2 <- i[,c("group","N","mean","se")]

    class(x2) <- c("repmean.single","repmean","data.frame")
    sw(cbind(x1,repmeandif(x2)))
  })
  )
  rownames(xx) <- NULL
  class(xx) <- c("repmeandif",class(xx))

  xx

}



# .LTdif <- function(x){
#   spl <- split(x,f = x$subject)
#
#   out <- do.call(rbind,lapply(spl,function(i){
#
#     i2 <- i[,-(1:4)]
#     class(i2) <- c("repmean.single","repmean","data.frame")
#
#     sw(cbind.data.frame(i[,1:4],.repmeandif(i2)))
#
#
#   }))
#   rownames(out) <- NULL
#
#   class(out) <- c("repmeandif", class(out))
#
#   out
# }


#' @export
print.repmeandif <- function(x, ...){

  dec = 5

  class(x) <- setdiff(class(x),c("repmeandif","repmeandif.list"))

  if(inherits(x,"list")){


    print(    lapply(x,function(i){

      maxdec(i, dec = dec)

    }))

  }else{
    print(maxdec(x, dec = dec))
  }



}
