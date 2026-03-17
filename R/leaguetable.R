#' ILSA's league tables
#'
#' Estimates the mean score for all countries within a cycle of an ILSA.
#' Arguments \code{method}, \code{reps}, and \code{var}, are extracted from
#' \code{\link{autoILSA}} and can be overridden by the user.
#'
#' @param study an optional character vector indicating the ILSA name, for a list of available
#'  ILSA, check \code{\link{autoILSA}}. If \code{NULL}, the ILSA name will be determined
#'  by the column names in the data frame.
#' @param year a numeric vector indicating the ILSA name, for a list of available
#'  cycles, check \code{\link{autoILSA}}.
#' @param subject an optional character vector indicating the subject for a list of available
#'  ILSA, check \code{\link{autoILSA}}.
#' @param fixN a logical value indicating if data should be "fixed" to meet official criteria.
#' For example, reducing the sample for certain countries in TIMSS 1995. Default is \code{TRUE}.
#' @param addCI a logical value indicating if confidence intervals should be added.
#' Defaults is \code{TRUE}.
#' @param specification a character value indicating extra specification like grade
#' (e.g., \code{"G8"} for TIMSS) or subject (e.g., \code{"Math"} for TIMSSADVANCED).
#' @inheritParams repmean
#' @inheritParams repcreate
#' @inheritParams repmeanCI
#'
#' @return a data frame.
#'
#' @examples
#' data(timss99)
#' leaguetable(df = timss99, year = 1999)
#'
#' @export
#'


leaguetable <- function(df,
                        study = NULL,
                        year,
                        subject = NULL,
                        specification = NULL,
                        addCI = TRUE,
                        alpha = 0.05,
                        method = NULL,
                        reps = NULL,
                        fixN = TRUE){

  # Argument checks ----

  # 1 - df
  # 2 - study
  # 3 - year
  # 4 - subject
  # 5 - method
  # 6 - reps
  # 7 - var - passes through repmean
  #
  # df = aa
  # study = "cived"
  # year = 1999
  # var = "ML"
  # subject = NULL
  # group = NULL
  # reps = NULL
  # method = NULL



  ili <- merge(ILSAstats::ILSAinfo$pvs,ILSAstats::ILSAinfo$weights,all.x = TRUE)
  ili$extravars[ili$extravars%in%"-"] <- NA
  cdf <- colnames(df)

  ## 3 - year, numeric value and within ILSAinfo ----
  # returnis(isval,year);year <- as.numeric(year)
  # returnis(isnumval,year)
  returnis(ischaeqnum,year)
  returnis(isinvec,x = year,choices = sort(unique(ili$year)))

  ili <- ili[ili$year%in%year,]


  ## 1 - df - check variables within df ----

  ilic <- lapply(1:nrow(ili), function(i){
    omitna(as.vector(unlist(lapply(ili[i,c("country","pvs","jkzones","jkreps",
                                           "totalweight","extravars")],
                                           strsplit,split = ";"))))
  })

  ili <- ili[sapply(ilic,function(i){all(i%in%cdf)}),]

  if(nrow(ili)==0)
    stop(paste0("\nInvalid input for 'df'.",
                "\nVariables in do not match conditions of any study.",
                "\nCheck needed variables in ILSAinfo$weights, and ILSAinfo$pvs"),
         call. = FALSE)





  ## 2 - study, character value and within ILSAinfo ----
  returnisNULL(ischaval,study)
  returnisNULL(isinvec,x = study,choices = sort(unique(ili$study)))

  if(!is.null(study)){
    study <- toupper(study)
    ili <- ili[ili$study%in%study,]

  }else{

    if(length(unique(ili$study))>1)
      stop(paste0("\nStudy can not be determined just using this 'df'.",
                  "\nSpecify the study name."),
           call. = FALSE)

  }



  ## 4 - subject, character value and within ILSAinfo ----
  returnisNULL(ischavec,subject)
  returnisNULL(isinvecmul,x = subject, choices = sort(unique(ili$subject)))

  if(!is.null(subject)){
    ili <- ili[ili$subject%in%subject,]
  }

  ## 5 - method ----
  returnisNULL(ischavec, method)
  returnisNULL(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))
  if(is.null(method)){method <- unique(ili$method)}

  ## 6 - reps ----
  returnisNULL(isnumval, reps)
  if(is.null(reps)){reps <- unique(ili$reps)}

  ## 7 - var - passes through repmean ----


  ## 8 - group ----

  if("IDCNTRY_STR"%in%colnames(df)){
    cou <- "IDCNTRY_STR"
  }else{
    cou <- unique(ili$country)
  }

  ## 6 - specification ----
  study2 <- specification





  if(!is.null(study2)){

    # specification = "-"
    returnisNULL(ischaval,specification)
    returnisNULL(isinvec,x = specification,choices = sort(unique(ili$study2)))

    study2 <- toupper(study2)
    ili <- ili[toupper(ili$study2)%in%study2,]

  }else{

    if(length(unique(ili$study2))>1){
      # specification = NULL
      returnisNULL(ischaval,specification)
      returnis(isinvec,specification,c(ili$study2))
      # specification <- unique(ili$study2)
      # returnis(ischaval,specification)
    }

  }





# Fixdata -----------------------------------------------------------------



  # if(!isdfonly(df)){
  #
  #
  #
  #   df <- df[,c(unlist(strsplit(ili$pvs,";")),
  #               ili$jkzones[1],
  #               ili$jkreps[1],
  #               ili$totalweight[1],
  #               cou)]
  #   df <- untidy(df)
  # }
  #

  evars <- strsplit(ili$extravars[1],";")[[1]]
  if(evars[1]%in%"-"){
    evars <- NULL
  }
  if(is.na(evars[1])){
    evars <- NULL
  }

  kolumns <- c(unlist(strsplit(ili$pvs,";")),
               ili$jkzones[1],
               ili$jkreps[1],
               ili$totalweight[1],
               evars,
               cou)

  # if(fixdata){
    df <- .fixdata(df = df,
                   study = ili$study[1],
                   year = ili$year[1],
                   specification = ili$study2[1],
                   columns = kolumns,
                   fixN = fixN)
  # }





  # Process -----------------------------------------------------------------

  rwi <- repcreate(df = df,
                   jkzone = ili$jkzones[1],
                   jkrep = ili$jkreps[1],
                   wt = ili$totalweight[1],
                   repwtname = "rwi",
                   reps = reps,
                   method = method)


  xx <- strsplit(ili$pvs,";")


  out <- vector("list",length(xx))
  for(i in 1:length(xx)){
    meai <- .repmean0(df = df,
                    x = xx[[i]],
                    PV = (length(xx[[i]])>1),
                    # setup = NULL,
                    repwt = rwi,
                    wt = ili$totalweight[i],
                    method = method,
                    var = -1,
                    group = cou,
                    by = NULL,
                    exclude = NULL,
                    aggregates = NULL,
                    zones = NULL)
    if(addCI){
      meai <- repmeanCI(x = meai, alpha = alpha, add = TRUE)
    }

    # if(includeid){
    meai <- cbind(study = ili$study[1],
                  study2 = ili$study2[1],
                  year = ili$year[1],
                  subject = ili$subject[i],
                  meai)
    # }

    out[[i]] <- meai
    rm(meai)

  }




  # Output ------------------------------------------------------------------

  out <- do.call(rbind,out)
  class(out) <- c("leaguetable",class(out))

  return(out)

}



#' @export
print.leaguetable <- function(x, ...){

  dec = 5

  class(x) <- setdiff(class(x),c("leaguetable"))

  if(inherits(x,"list")){


    print(    lapply(x,function(i){

      maxdec(i, dec = dec)

    }))

  }else{
    print(maxdec(x, dec = dec))
  }



}
