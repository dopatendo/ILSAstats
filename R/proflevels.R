#' ILSA's proficiency levels
#'
#' Estimates the proficiency levels for all countries within a cycle of an ILSA.
#' Arguments \code{method}, \code{reps}, and \code{var}, are extracted from
#' \code{\link{ILSAinfo}} and can be overridden by the user.
#'
#'
#' @inheritParams leaguetable
#' @inheritParams repprop.table
#'
#' @return a data frame or a list.
#'
#' @examples
#'
#' proflevels(timss99,year = 1999,type = "long",subject = "math")
#'
#' @export
#'

proflevels <- function(df,
                       study = NULL, year, subject = NULL,
                       method = NULL, reps = NULL,
                       type = c("long","wide1","wide2"),
                       separateSE = TRUE){

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
  ili <- merge(ili,ILSAstats::ILSAinfo$levels,all.x = TRUE)
  # ili <- unique(ili[,!colnames(ili)%in%c("reps","year")])
  ili <- stats::na.omit(ili)
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
                "\nCheck needed variables in ILSAinfo$weights, ILSAinfo$levels, and ILSAinfo$pvs"),
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
  returnis(isnumval,year)
  returnis(isinvec,x = year,choices = sort(unique(ili$year)))

  ili <- ili[ili$year%in%year,]


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



  # Process -----------------------------------------------------------------

  levs <- df[,c(cou,unlist(c(ili[1,c("jkzones","jkreps","totalweight")]),
                           use.names = FALSE))]
  levs <- untidy(levs)

  levs <- cbind.data.frame(levs,
                           proflevels.get(df = df,study = study,combine = TRUE))

  rwi <- repcreate(df = levs,
                   jkzone = ili$jkzones[1],
                   jkrep = ili$jkreps[1],
                   wt = ili$totalweight[1],
                   repwtname = "rwi",
                   reps = reps,
                   method = method)

  out <- vector("list",nrow(ili))

  for(i in 1:length(out)){


    xi <- paste0(strsplit(ili[i,"pvs"],";")[[1]],"_level")
    ci <- 0:length(strsplit(ili[i,"cutoffs"],";")[[1]])

    pri <- sm(repprop(x = xi,
                      categories = ci,
                      setup = NULL,
                      repwt = rwi,
                      wt = ili$totalweight[1],
                      df = levs,
                      method = method,
                      group = cou,
                      exclude = NULL,
                      aggregates = NULL))
    pri <- repprop.table(x = pri, type = type, separateSE = separateSE)

    if(!type%in%"wide1"){
      if(islist(pri)){

        pri <- lapply(pri,function(j){

          cana <- rep(strsplit(ili[i,"names"],";")[[1]],each = nrow(j)/length(ci))

          nuca <- which(colnames(j)%in%"category")
          jj <- cbind(j[,c(1:nuca)],
                      level = cana,
                      j[,c((nuca+1):ncol(j))])
          colnames(jj)[nuca] <- c("category")
          jj

        })



      }else{

        cana <- rep(strsplit(ili[i,"names"],";")[[1]],each = nrow(pri)/length(ci))

        nuca <- which(colnames(pri)%in%"category")
        pri <- cbind(pri[,c(1:nuca)],
                     level = cana,
                     pri[,c((nuca+1):ncol(pri))])
        colnames(pri)[nuca] <- c("category")


      }
    }

    out[[i]] <- pri


  }

  if(length(out)==1)
    return(out[[1]])

  names(out) <- ili$subject
  return(out)

}
