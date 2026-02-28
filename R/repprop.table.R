#' Tables for Proportions with Replicate Weights
#'
#' Creates tables for proportions using replicate weights
#' for a variable or a group of plausible values variables and for one or more
#' populations.
#'
#'
#' @param x a list produced by \code{\link{repprop}}.
#' @param type a charaacter value indicating the type of table to produce.
#' Options include: \code{"long"}, for a long table with a column with the proportions
#' and another one for the standard error; \code{"wide1"} for a wide table where groups
#' are distributed in lines; \code{"wide2"} for a wide table where groups are distributed in columns.
#' @param separateSE a logical value indicating if standard errors should be separated from
#' proportions, each as an element from a list. Only works for wide tables. Default is \code{TRUE}.
#'
#'
#' @return a adata frame or a list.
#'
#' @example inst/examples/repprop.table_example.R
#' @export
#'



repprop.table <- function(x,
                          type = c("long","wide1","wide2"),
                          separateSE = TRUE){


  # Argument checks ----

  returnis(isrep.prop,x)
  returnis(islova,separateSE)



  frm <- formals(repprop.table)
  returnis(ischavec, type)
  type <- returnis(isinvec,x = type[1L],choices = frm$type)

  # Process -----

  xxi <- x[-length(x)]
  xx <- do.call(rbind,lapply(xxi,function(i){
    i[,intersect(colnames(i),c("group","prop","se"))]
  }))

  if(!"group"%in%colnames(xx)){
    xx <- cbind.data.frame(group = "Pooled",xx)
  }

  xx$category <- rep(x$categories,each = nrow(x[[1]]))
  xx <- xx[,c("group","category","prop","se")]
  rownames(xx) <- NULL




  if(type%in%c("long"))
    return(xx)




  if(type%in%c("wide1")){


    w1 <- stats::reshape(xx,direction = "wide",idvar = "group",timevar = "category")

    if(separateSE){
      w1 <-  list(prop = w1[,c(1,seq(2,ncol(w1),by=2))],
                  se = w1[,c(1,seq(3,ncol(w1),by=2))])
    }


    return(w1)

  }



  if(type%in%"wide2"){

    w1 <- stats::reshape(xx,direction = "wide",idvar = "group",timevar = "category")

    w2 <- cbind.data.frame(category = rep(x$categories,each = 2),
                           statistic = c("prop","se"),
                           t(w1[,-1]))
    colnames(w2)[-(1:2)] <- w1$group
    rownames(w2) <- NULL
    out <- w2

    if(separateSE){
      out <- split(w2, ~ statistic)
      out <- lapply(out,function(i){
        i[,!colnames(w2)%in%"statistic"]
      })
    }

    return(out)

  }




}



