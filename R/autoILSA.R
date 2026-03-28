#' Options for automatic functions
#'
#' Shows all options for the automatic functions: \code{\link{leaguetable}}, and \code{\link{proflevels}}.
#'
#' @param func a character value indicating the options of which function should be printed.
#' @param study an optional character vector indicating the ILSA name.
#'
#' @return a list.
#'
#' @examples
#' autoILSA()
#'
#' @export
#'


autoILSA <- function(func = c("leaguetable","proflevels"),
                     study = NULL){


  ILSAinfo <- ILSAstats::ILSAinfo

# Argument checks ---------------------------------------------------------

  sch <- sort(unique(ILSAinfo$pvs$study))
  returnisNULL(isinvec,study,choices = sch)


  if(!is.null(study))
    study <- isinvec(study,sch)

  frm <- formals(autoILSA)
  returnis(ischavec, func)
  funci <- returnis(isinvec,func[1L],choices = frm$func)


# Select function ---------------------------------------------------------



  # funci <- "proflevels"

  # leaguetable
  if(funci=="leaguetable"){
    ili <- merge(ILSAstats::ILSAinfo$pvs,ILSAstats::ILSAinfo$weights,all.x = TRUE)
    ili$extravars[ili$extravars%in%"-"] <- NA
  }


  # proflevels
  if(funci=="proflevels"){
    ili <- merge(ILSAstats::ILSAinfo$pvs,ILSAstats::ILSAinfo$weights,all.x = TRUE)
    ili <- merge(ili,ILSAstats::ILSAinfo$levels,all.x = TRUE)
    ili <- stats::na.omit(ili)
  }



# Process -----------------------------------------------------------------

  us <- sort(unique(ili$study))

  i=1
  li <- lapply(1:length(us),function(i){
    ilii <- ili[ili$study%in%us[i],]
    list(year = (sort(unique(ilii$year))),
         specification = (gsub("-","NULL",sort(unique(ilii$study2)))),
         subject = (sort(unique(ilii$subject))))
  })
  names(li) <- us



# Output ------------------------------------------------------------------




  if(!is.null(study)){
    li <- (li[study])
  }

  class(li) <- c("autoILSA",class(li))
  attributes(li)$func <- funci

  return(li)

}


#' @export
print.autoILSA <- function(x, ...){



  if(length(x)==1){

    cat(paste0("Options for ",attributes(x)$func,"(...) in ",names(x),":"))

    cat("\nyear:",paste0(paste0(x[[1]][[1]],collapse = ", "),"."))

    cat("\nspecification:",paste0(paste0(x[[1]][[2]],collapse = ", "),"."))

    cat("\nsubject:",paste0(paste0(x[[1]][[3]],collapse = ", "),"."))
  }else{

    cat(paste0("Options for ",attributes(x)$func,"(...):"))

    cat("\n\nyear:",paste0("\n",names(x)," = ",sapply(sapply(x, `[[`, 1),function(i){
      paste0(i,collapse = ", ")
    }),"."))

    cat("\n\nspecification:",paste0("\n",names(x)," = ",sapply(sapply(x, `[[`, 2),function(i){
      paste0(i,collapse = ", ")
    }),"."))

    cat("\n\nsubject:",paste0("\n",names(x)," = ",sapply(sapply(x, `[[`, 3),function(i){
      paste0(i,collapse = ", ")
    }),"."))
  }




}



