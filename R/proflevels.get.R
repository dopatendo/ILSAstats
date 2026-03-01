#' Get ILSA's proficiency levels
#'
#' Converts ILSA scores into proficiency levels.
#'
#'
#' @inheritParams leaguetable
#' @param combine a logical value indicating if subjects should be combined in
#' a single data frame.
#'
#' @return a data frame or a list.
#'
#' @examples
#'
#' proflevels.get(timss99, subject = "math")
#'
#' @export
#'


proflevels.get <- function(df, study = NULL, subject = NULL, combine = TRUE){
  # Argument checks ----

  returnis(islova,combine)

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



  # ili <- merge(ILSAstats::ILSAinfo$pvs,ILSAstats::ILSAinfo$weights,all.x = TRUE)
  ili <- merge(ILSAstats::ILSAinfo$pvs,ILSAstats::ILSAinfo$levels,all.x = TRUE)
  ili <- unique(ili[,!colnames(ili)%in%c("reps","year")])
  ili <- stats::na.omit(ili)
  cdf <- colnames(df)

  ## 1 - df - check variables within df ----

  ilic <- lapply(1:nrow(ili), function(i){
    as.vector(unlist(lapply(ili[i,c("country","pvs")],
                            strsplit,split = ";")))
  })

  ili <- ili[sapply(ilic,function(i){all(i%in%cdf)}),]

  if(nrow(ili)==0)
    stop(paste0("\nInvalid input for 'df'.",
                "\nVariables in do not match conditions of any study.",
                "\nCheck needed variables in ILSAinfo$levels, and ILSAinfo$pvs"),
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

  # ## 3 - year, numeric value and within ILSAinfo ----
  # returnis(isnumval,year)
  # returnis(isinvec,x = year,choices = sort(unique(ili$year)))
  #
  # ili <- ili[ili$year%in%year,]


  ## 4 - subject, character value and within ILSAinfo ----
  returnisNULL(ischavec,subject)
  returnisNULL(isinvecmul,x = subject, choices = sort(unique(ili$subject)))

  if(!is.null(subject)){
    ili <- ili[ili$subject%in%subject,]
  }

  # ## 5 - method ----
  # returnisNULL(ischavec, method)
  # returnisNULL(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))
  # if(is.null(method)){method <- unique(ili$method)}

  # ## 6 - reps ----
  # returnisNULL(isnumval, reps)
  # if(is.null(reps)){reps <- unique(ili$reps)}

  ## 7 - var - passes through repmean ----


  # ## 8 - group ----
  #
  # if("IDCNTRY_STR"%in%colnames(df)){
  #   cou <- "IDCNTRY_STR"
  # }else{
  #   cou <- unique(ili$country)
  # }
  #


  # Process -----------------------------------------------------------------



  lev <- vector("list",nrow(ili))
  for(i in 1:nrow(ili)){
    pvsi = strsplit(ili$pvs[i],";")[[1]]
    bchsi = as.numeric(strsplit(ili$cutoffs[i],";")[[1]])
    moreorequal = ili$moreorequal[i]==1

    lev[[i]] <- .proflevels.get(df = df,
                                pvs = pvsi,
                                bchs = bchsi,
                                moreorequal = moreorequal,
                                acumulated = FALSE)
  }



  # Output ------------------------------------------------------------------


  if(combine){
    lev <- do.call(cbind,lev)
    return(lev)
  }else{
    names(lev) <- ili$subject
    return(lev)
  }

}


.proflevels.get <- function(df,pvs,bchs,moreorequal = TRUE, acumulated = FALSE){

  sco <- untidy(df[,pvs])
  bch <- matrix(0,nrow = nrow(sco),ncol = ncol(sco))

  i=1
  if(moreorequal){
    for(i in 1:length(bchs)){
      bch[sco>=bchs[i]] <- i
    }
  }else{
    for(i in 1:length(bchs)){
      bch[sco>bchs[i]] <- i
    }
  }
  bch[is.na(sco)] <- NA
  bch <- as.data.frame(bch)
  rm(sco)
  colnames(bch) <- paste0(pvs,"_level")

  if(!acumulated)
    return(bch)



  bch <- lapply(1:length(bchs),function(i){
    (bch>=i)*1L
  })

  return(bch)


}



