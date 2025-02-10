# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name


### Groups ----

# One variable
reme <- repmean(x = c("item01"),
                PV = FALSE,
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)


# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)

### Groups and By ----

# One variable
reme <- repmean(x = c("item01"),
                PV = FALSE,
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)

# One PV variable
reme <- repmean(x = paste0("Math",1:5),
                PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
                repwt = RW, wt = "wt", df = repdata,
                method = "ICILS",var = "ML",zones = "jkzones",
                group = "GROUP",
                by = "GENDER", # results will be separated by GENDER
                exclude = "GR2") # GR2 will not be used for Pooled nor Composite

repmeandif(reme)
