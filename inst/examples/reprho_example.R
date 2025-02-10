# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                wt = "wt", # the total weights column name
                jkzone = "jkzones", # the jkzones column name
                jkrep = "jkrep", # the jkreps column name
                repwtname = "REPWT", # the desired name for the rep weights
                reps = 50, # the number of replications
                method = "ICILS") # the name of the method aka the study name


### No groups ----

# Non PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = NULL,
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")

# X var and PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = paste0("Reading",1:5),
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")

# PVs and PVs (related)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")

# PVs and PVs (UNrelated)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       relatedpvs = FALSE,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       method = "ICILS")


### Groups ----

# Non PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = NULL,
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")

# X var and PVs
reprho(x = c("GENDER",paste0("Math",1:3)),
       pv = paste0("Reading",1:5),
       pv2 = NULL,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")

# PVs and PVs (related)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")

# PVs and PVs (UNrelated)
reprho(x = NULL,
       pv = paste0("Math",1:5),
       pv2 = paste0("Reading",1:5),
       relatedpvs = FALSE,
       rho = "pearson",
       repwt = RW,
       wt = "wt",
       df = repdata,
       group = "GROUP",
       method = "ICILS")