# ICC of one variable
icc(x = "Math1",group = "GROUP", weights = repdata$wt, data = repdata)


# ICC of more than one variable
icc(x = c("Math1","Math2","Math3","Math4","Math5","SES"),
    group = "GROUP", weights = repdata$wt, data = repdata)


# ICC of PVs
icc(x = c("Math1","Math2","Math3","Math4","Math5"), PV = TRUE,
    group = "GROUP", weights = repdata$wt, data = repdata)
