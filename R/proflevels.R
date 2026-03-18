#' ILSA's proficiency levels
#'
#' Estimates the proficiency levels for all countries within a cycle of an ILSA.
#' Arguments \code{method}, and \code{reps}, are extracted from
#' \code{\link{autoILSA}} and can be overridden by the user.
#'
#'
#' @inheritParams leaguetable
#' @inheritParams repprop.table
#' @param acumulated a logical value indicating if proficiency levels should be accumulated.
#'
#' @return a data frame or a list.
#'
#' @examples
#' data(timss99)
#'
#' proflevels(timss99,year = 1999,type = "long",subject = "math")
#'
#' @export
#'



proflevels <- function(df,
                       study = NULL, year, subject = NULL,
                       method = NULL, reps = NULL,
                       type = c("long","wide1","wide2"),
                       separateSE = TRUE,
                       fixN = TRUE,
                       accumulated = FALSE){


  # rm(list = ls())
  # library(ILSAstats)
  # source("R/internal.R")
  # source("R/argchecks.R")
  # source("R/repmean.R")
  # source("R/repse.R")
  # source("R/prepILSA.R")
  # source("R/proflevels.get.R")
  # df = ILSAstats::timss99
  # study = NULL
  # year = 1999
  # subject = NULL
  # method = NULL
  # reps = NULL
  # type = c("wide2")
  # separateSE = FALSE
  # fixN = TRUE
  # accumulated = FALSE

  # Argument checks ----

  frm <- formals(repprop.table)
  returnis(ischavec, type)
  type <- returnis(isinvec,x = type[1L],choices = frm$type)

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

  # ## 3 - year, numeric value and within ILSAinfo --
  # returnis(isval,year);year <- as.numeric(year)
  # returnis(isnumval,year)
  returnis(ischaeqnum,year)
  returnis(isinvec,x = year,choices = sort(unique(ili$year)))

  ili <- ili[ili$year%in%year,]


  ## 1 - df - check variables within df --

  ilic <- lapply(1:nrow(ili), function(i){
    x <- omitna(as.vector(unlist(lapply(ili[i,c("country","pvs","jkzones","jkreps","totalweight","extravars")],
                                   strsplit,split = ";"))))
    setdiff(x,"-")
  })

  ili <- ili[sapply(ilic,function(i){all(i%in%cdf)}),]

  if(nrow(ili)==0)
    stop(paste0("\nInvalid input for 'df'.",
                "\nVariables in do not match conditions of any study.",
                "\nCheck needed variables in ILSAinfo$weights, ILSAinfo$levels, and ILSAinfo$pvs"),
         call. = FALSE)



  ## 2 - study, character value and within ILSAinfo --
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



  ## 4 - subject, character value and within ILSAinfo --
  returnisNULL(ischavec,subject)
  returnisNULL(isinvecmul,x = subject, choices = sort(unique(ili$subject)))

  if(!is.null(subject)){
    ili <- ili[ili$subject%in%subject,]
  }

  ## 5 - method --
  returnisNULL(ischavec, method)
  returnisNULL(isinvec,x = method[1L],choices = ILSAmethods(repse = TRUE))
  if(is.null(method)){method <- unique(ili$method)}

  ## 6 - reps --
  returnisNULL(isnumval, reps)
  if(is.null(reps)){reps <- unique(ili$reps)}

  ## 7 - var - passes through repmean --


  ## 8 - group --

  if("IDCNTRY_STR"%in%colnames(df)){
    cou <- "IDCNTRY_STR"
  }else{
    cou <- unique(ili$country)
  }


  # Fixdata -----------------------------------------------------------------



  if(fixN){
    df <- .fixdata(df = df, study = ili$study[1],
                   year = ili$year[1],
                   specification = ili$study2[1],fixN = fixN)
  }



  # Process -----------------------------------------------------------------

  # levs <- df[,c(cou,unlist(c(ili[1,c("jkzones","jkreps","totalweight")]),
  #                          use.names = FALSE))]
  # levs <- untidy(levs)
  #
  #
  # rwi <- repcreate(df = levs,
  #                  jkzone = ili$jkzones[1],
  #                  jkrep = ili$jkreps[1],
  #                  wt = ili$totalweight[1],
  #                  repwtname = "rwi",
  #                  reps = reps,
  #                  method = method)
  #
  #
  # if(acumulated)
  #
  # i=1
  # # for(i 1:nrow(ili)){}
  #
  # aci <- .proflevels.get(df = df,
  #                        pvs = strsplit(ili$pvs[i],";")[[1]],
  #                        bchs = as.numeric(strsplit(ili$cutoffs[i],";")[[1]]),
  #                        moreorequal = (ili$moreorequal[i]==1),
  #                        acumulated = TRUE)
  # ci <- 0:length(strsplit(ili[i,"cutoffs"],";")[[1]])
  #
  # j=1
  # pri <- vector("list",length(aci))
  # for(j in 1:length(aci)){
  #   dfj <- cbind(levs,aci[[j]])
  #
  #   xi <- colnames(aci[[j]])
  #
  #   prj <- sm(repprop(x = xi,
  #                     categories = 0:1,
  #                     setup = NULL,
  #                     repwt = rwi,
  #                     wt = ili$totalweight[i],
  #                     df = dfj,
  #                     method = method,
  #                     group = cou,
  #                     exclude = NULL,
  #                     aggregates = NULL))
  #   pri[[j]] <- prj[["PVs==1"]]
  #   names(pri)[j] <- paste0("PVs==",ci[j+1])
  # }
  # class(pri) <- c("repprop",class(pri))
  # attributes(pri)$categories <- ci[-1]



# Acumulated --------------------------------------------------------------



  if(accumulated){

    levs <- df[,c(cou,unlist(c(ili[1,c("jkzones","jkreps","totalweight")]),
                             use.names = FALSE))]
    levs <- untidy(levs)


    rwi <- repcreate(df = levs,
                     jkzone = ili$jkzones[1],
                     jkrep = ili$jkreps[1],
                     wt = ili$totalweight[1],
                     repwtname = "rwi",
                     reps = reps,
                     method = method)



    out <- vector("list",nrow(ili))

    i=1
    for(i in 1:length(out)){


      aci <- .proflevels.get(df = df,
                             pvs = strsplit(ili$pvs[i],";")[[1]],
                             bchs = as.numeric(strsplit(ili$cutoffs[i],";")[[1]]),
                             moreorequal = (ili$moreorequal[i]==1),
                             acumulated = TRUE)
      ci <- 0:length(strsplit(ili[i,"cutoffs"],";")[[1]])

      j=1
      pri <- vector("list",length(aci))
      for(j in 1:length(aci)){
        dfj <- cbind(levs,aci[[j]])

        xi <- colnames(aci[[j]])

        prj <- sm(repprop(x = xi,
                          categories = 0:1,
                          setup = NULL,
                          repwt = rwi,
                          wt = ili$totalweight[i],
                          df = dfj,
                          method = method,
                          group = cou,
                          exclude = NULL,
                          aggregates = NULL))
        pri[[j]] <- prj[["PVs==1"]]
        names(pri)[j] <- paste0("PVs==",ci[j+1])
      }
      class(pri) <- c("repprop",class(pri))
      attributes(pri)$categories <- ci[-1]

      pri <- repprop.table(x = pri, type = type, separateSE = separateSE)

      caN <- strsplit(ili[i,"names"],";")[[1]][-1]

      if(!type%in%"wide1"){
        if(islist(pri)){

          pri <- lapply(pri,function(j){

            cana <- rep(caN,each = nrow(j)/length(caN))

            nuca <- which(colnames(j)%in%"category")
            jj <- cbind(j[,c(1:nuca)],
                        level = cana,
                        j[,c((nuca+1):ncol(j))])
            colnames(jj)[nuca] <- c("category")
            jj

          })



        }else{

          cana <- rep(caN,each = nrow(pri)/length(caN))

          nuca <- which(colnames(pri)%in%"category")
          pri <- cbind(pri[,c(1:nuca)],
                       level = cana,
                       pri[,c((nuca+1):ncol(pri))])
          colnames(pri)[nuca] <- c("category")


        }
      }

      out[[i]] <- pri


    }
  }


# Not acumulated ----------------------------------------------------------



if(!accumulated){
  levs <- df[,c(cou,unlist(c(ili[1,c("jkzones","jkreps","totalweight")]),
                           use.names = FALSE))]
  levs <- untidy(levs)


  rwi <- repcreate(df = levs,
                   jkzone = ili$jkzones[1],
                   jkrep = ili$jkreps[1],
                   wt = ili$totalweight[1],
                   repwtname = "rwi",
                   reps = reps,
                   method = method)


  levs <- cbind.data.frame(levs,
                           proflevels.get(df = df,study = study,combine = TRUE))

  out <- vector("list",nrow(ili))

  i=1
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
}





# Output ------------------------------------------------------------------



  if(length(out)==1){
    out <- out[[1]]
    class(out) <- c("proflevels",class(out))
    return(out)
  }


  names(out) <- ili$subject
  class(out) <- c("proflevels",class(out))
  return(out)

}


#' @export
print.proflevels <- function(x, ...){

  dec = 5

  class(x) <- setdiff(class(x),c("proflevels"))

  if(inherits(x,"list")&&inherits(x[[1]],"list")){


    print( lapply(x,function(i){
      lapply(i,function(j){
        maxdec(j, dec = dec)
      })
    }))


  }else{
    if(inherits(x,"list")){


      print(    lapply(x,function(i){

        maxdec(i, dec = dec)

      }))

    }else{
      print(maxdec(x, dec = dec))
    }
  }






}
