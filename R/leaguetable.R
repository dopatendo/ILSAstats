#' ILSA's league tables
#'
#' Estimates the mean score for all countries within a cycle of an ILSA.
#' Arguments \code{method}, \code{reps}, and \code{var}, are extracted from
#' \code{\link{ILSAinfo}} and can be overridden by the user.
#'
#' @param study an optional character vector indicating the ILSA name, for a list of available
#'  ILSA, check \code{\link{ILSAinfo}}. If \code{NULL}, the ILSA name will be determined
#'  by the column names in the data frame.
#' @param year a numeric vector indicating the ILSA name, for a list of available
#'  cycles, check \code{\link{ILSAinfo}}.
#' @param study an optional character vector indicating the subjects to be analyzed, for a list of available
#'  subjects, check \code{\link{ILSAinfo}}.
#' @param subject an optional character vector indicating the subject for a list of available
#'  ILSA, check \code{\link{ILSAinfo}}.
#' @inheritParams repmean
#' @inheritParams repcreate
#'
#' @return a data frame.
#'
#' @examples
#'
#' leaguetable(df = timss99, year = 1999)
#'
#' @export
#'


leaguetable <- function(df,
                        study = NULL,
                        year,
                        subject = NULL,
                        method = NULL,
                        reps = NULL,
                        var = c("unbiased", "ML")){

  # Argument checks ----

  # 1 - df
  # 2 - study
  # 3 - year
  # 4 - subject
  # 5 - method
  # 6 - reps
  # 7 - var - passes through repmean

  # df = df
  # study = "TIMSS"
  # year = 1999
  # var = "ML"
  # subject = "math"
  # group = NULL
  # reps = NULL
  # method = NULL



  ili <- merge(ILSAstats::ILSAinfo$pvs,ILSAstats::ILSAinfo$weights,all.x = TRUE)
  # ili$jkreps <- "caca"
  cdf <- colnames(df)

  ## 1 - df - check variables within df ----

  ilic <- lapply(1:nrow(ili), function(i){
    as.vector(unlist(lapply(ili[i,c("country","pvs","jkzones","jkreps","totalweight")],
                            strsplit,split = ";")))
  })

  ili <- ili[sapply(ilic,function(i){all(i%in%cdf)}),]

  if(nrow(ili)==0)
    stop(paste0("\nInvalid input for 'df'.",
                "\nVariables in do not match conditions of any study.",
                "\nCheck needed variables in ILSAinfo$weights and ILSAinfo$pvs"),
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

  ## 3 - year, numeric value and within ILSAinfo ----
  returnis(isnumval,year)
  returnis(isinvec,x = year,choices = sort(unique(ili$year)))

  ili <- ili[ili$year%in%year,]


  ## 4 - subject, numeric value and within ILSAinfo ----
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
    meai <- repmean(df = df,
                    x = xx[[i]],
                    PV = (length(xx[[i]])>1),
                    setup = NULL,
                    repwt = rwi,
                    wt = ili$totalweight[i],
                    method = method,
                    var = var,
                    group = cou,
                    by = NULL,
                    exclude = NULL,
                    aggregates = NULL,
                    zones = NULL)
    # if(includeid){
    meai <- cbind(study = ili$study[1],
                  year = ili$year[1],
                  subject = ili$subject[i],
                  meai)
    # }

    out[[i]] <- meai
    rm(meai)

  }




  # Output ------------------------------------------------------------------

  do.call(rbind,out)

}
