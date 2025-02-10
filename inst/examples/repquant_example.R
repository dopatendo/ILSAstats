
RWT <- repcreate(df = repdata, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----

# One variable - weights within df
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = FALSE,
        repwt = "REPWT", wt = "wt", df = cbind(repdata,RWT),
        method = "ICILS")

# One variable - weights as a separate data frame
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = FALSE,
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS")

# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS")

### Groups ----

# One variable
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = FALSE,
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite


# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite

### Groups and By ----

# One variable
repquant(x = c("item01"),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = FALSE,
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite

# One PV variable
repquant(x = paste0("Math",1:5),
        qtl = c(0.05, 0.25, 0.75, 0.95),
        PV = TRUE, # if set to TRUE, PVs will be treated as separate variables
        repwt = RWT, wt = "wt", df = repdata,
        method = "ICILS",
        group = "GROUP",
        by = "GENDER", # results will be separated by GENDER
        exclude = "GR2") # GR2 will not be used for Pooled nor Composite
