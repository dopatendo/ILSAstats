# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----
stp1 <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "ICILS")
stp1

### Groups ----
stp2 <- repsetup(repwt = RW, wt = "wt", df = repdata, method = "ICILS",
                 group = "GROUP", exclude = "GR2")
stp2


### repmean ----

repmean(x = "Math1",setup = stp1)

repmean(x = "Math1",setup = stp2)
