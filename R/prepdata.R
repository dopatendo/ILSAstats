#' Prepare ILSA Data
#'
#'
#' Modifies ILSA data to meet official participation cases, selects columns and
#' transforms data into simple data frames converting missing values to NAs.
#'
#' @inheritParams leaguetable
#' @param specification a character value indicating extra specification like grade
#' (e.g., \code{"G8"} for TIMSS) or subject (e.g., \code{"Math"} for TIMSSADVANCED).
#' @param columns a character vector indicating which columns should be selected.
#' If \code{NULL}, all columns will be selected.
#'
#'
#' @examples
#'
#' data(timss99)
#'
#' head(timss99)
#'
#' newdata <- prepdata(df = timss99, columns = paste0("BSMMAT0",1:5),fixN = FALSE)
#'
#' head(newdata)
#'
#' @export
#'
#'





prepdata <- function(df,
                     study = NULL,
                     year = NULL,
                     specification = NULL,
                     fixN = TRUE,
                     columns = NULL){

  # Just select columns if fixN = FALSE
  ## 5 - fixN
  returnis(islova,fixN)

  ## 4 - columns

  cdf <- colnames(df)
  returnisNULL(isinvecmulExact, x = columns, choices = cdf)
  if(is.null(columns)){columns <- cdf}

  if(!fixN){
    xx <- df[,columns,drop = FALSE]
    xx <- untidy(xx,mistoNAs = TRUE)

    return(xx)
  }


  # Argument checks ----

  # 1 - year
  # 2 - df
  # 3 - study
  # 4 - columns


  ili <- merge(ILSAstats::ILSAinfo$pvs,ILSAstats::ILSAinfo$weights,all.x = TRUE)
  ili$extravars[ili$extravars%in%"-"] <- NA


  ## 1 - year, numeric value and within ILSAinfo ----
  returnis(ischaeqnum,year)
  returnis(isinvec,x = year,choices = sort(unique(ili$year)))
  ili <- ili[ili$year%in%year,]


  ## 2 - df - check variables within df ----
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


  ## 3 - study, character value and within ILSAinfo ----
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





  ## 6 - specification ----
  study2 <- specification





  if(!is.null(study2)){

    # specification = "-"
    returnisNULL(ischaval,specification)
    returnisNULL(isinvec,x = specification,choices = sort(unique(ili$study2)))

    study2 <- toupper(study2)
    ili <- ili[ili$study2%in%study2,]

  }else{

    # specification = NULL
    returnisNULL(ischaval,specification)
    specification <- unique(ili$study2)
    returnis(ischaval,specification)
  }


  # Process -----------------------------------------------------------------


  .fixdata(df = df,
           study = study,
           year = year,
           specification = study2[1],
           fixN = fixN,
           columns = columns)
}



.fixdata <- function(df, study, year, specification,fixN,columns){

  study <- toupper(study)
  specification <- toupper(specification)

  tofix <- cbind.data.frame(study = "TIMSS", year = c(1995), specification = c("G4","G8"))

  tofix <- tofix[tofix$study%in%study&tofix$year%in%year&tofix$specification%in%specification,]

  if((!fixN)|nrow(tofix)==0){

    if(!isdfonly(df)){
      df <- df[,columns,drop = FALSE]
      df <- untidy(df)
    }else{
      df <- df[,columns]
    }


    return(df)

  }




  if(study=="TIMSS"&year==1995&specification=="G4"){

    dati <- df[,unique(c(columns,"IDCNTRY","IDGRADER","IDGRADE")), drop = FALSE]
    # dati <- aa[,c(kolumns,"IDCNTRY","IDGRADER","IDGRADE")]
    if(!isdfonly(dati)){
      # dati <- dati[,kolumns]
      dati <- untidy(dati)
    }


    # Slovenia (ID 705) G4, IDGRADER==1, 3rd
    # Australia (ID 36) G4, IDGRADE==4, 4th

    excou <- c("Slovenia","Australia")
    excou <- c(705,36)

    v1 <- c("IDGRADER","IDGRADE")
    g1 <- c(1,4)


    for(w in 1:length(excou)){
      iscou <- dati$IDCNTRY%in%excou[w]
      iscas <- !dati[,v1[w],drop = TRUE]%in%g1[w]
      dati <- dati[-which(iscou&iscas),]
    }

    dati <- dati[,columns]
    return(dati)

  }

  if(study=="TIMSS"&year==1995&specification=="G8"){

    dati <- df[,unique(c(columns,"IDCNTRY","IDGRADER","IDGRADE")), drop = FALSE]
    if(!isdfonly(dati)){
      # dati <- dati[,kolumns]
      dati <- untidy(dati)
    }


    # Sweden (ID 752) G8, IDGRADER==3, 8th
    # Slovenia (ID 705) G8, IDGRADER==1, 7th
    # Colombia (ID 170) G8, IDGRADER==1, 7th
    # Australia (ID 36) G8, IDGRADE==8, 8th

    excou <- c("Sweden","Slovenia","Colombia","Australia")
    excou <- c(752,705,170,36)

    v1 <- c("IDGRADER","IDGRADER","IDGRADER","IDGRADE")
    g1 <- c(3,1,1,8)


    for(w in 1:length(excou)){
      iscou <- dati$IDCNTRY%in%excou[w]
      iscas <- !dati[,v1[w],drop = TRUE]%in%g1[w]
      dati <- dati[-which(iscou&iscas),]
    }

    dati <- dati[,columns]
    return(dati)



  }

}




