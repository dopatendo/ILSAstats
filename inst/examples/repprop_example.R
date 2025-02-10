# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----

# One variable - weights within df
repprop(x = c("item01"),
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
        method = "ICILS")

# One variable - weights weights as a separate data frame
repprop(x = c("item01"),
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

# Multiple variables - PVs are assumed
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

### Groups ----

# One variable - weights within df
repprop(x = c("item01"),
        group = "GROUP",
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
        method = "ICILS")

# One variable - weights weights as a separate data frame
repprop(x = c("item01"),
        group = "GROUP",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

# Multiple variables - PVs are assumed
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        group = "GROUP",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

# Multiple variables - excluding one group
repprop(x = c("CatMath1","CatMath2","CatMath3"),
        group = "GROUP",
        exclude = "GR2",
        repwt = RW, wt = "wt", df = repdata,
        method = "ICILS")

