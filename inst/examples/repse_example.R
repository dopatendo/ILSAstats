# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                wt = "wt", # the total weights column name
                jkzone = "jkzones", # the jkzones column name
                jkrep = "jkrep", # the jkreps column name
                repwtname = "REPWT", # the desired name for the rep weights
                reps = 50, # the number of replications
                method = "ICILS") # the name of the method aka the study name

# Non-PVs ----

## Mean with total weights
E0 <- stats::weighted.mean(x = repdata$item01, w = repdata$wt, na.rm = TRUE)
E0

## Means by replication
ER <- as.vector(apply(RW,2,function(i){
  stats::weighted.mean(x = repdata$item01, w = i, na.rm = TRUE)
}))
ER

## Standard error by hand
repse(er = ER, e0 = E0, method = "ICILS")

## Standard error with repmean()
repmean(x = "item01",wt = "wt",repwt = RW,df = repdata, method = "ICILS")


# PVs ----

## Mean with total weights
E0 <- sapply(1:5,function(i){
  stats::weighted.mean(x = repdata[,paste0("Math",i)], w = repdata$wt,
                       na.rm = TRUE)
})
E0

## Means by replication
ER <- lapply(1:5, function(j){
  as.vector(apply(RW,2,function(i){
    stats::weighted.mean(x = repdata[,paste0("Math",j)], w = i, na.rm = TRUE)
  }))
})
ER

## Standard error by hand
repse(er = ER, e0 = E0, method = "ICILS")

## Standard error with repmean()
repmean(x = paste0("Math",1:5),wt = "wt",repwt = RW,df = repdata, method = "ICILS",PV = TRUE)


