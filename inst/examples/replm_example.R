# Less data for shorter example
repdata2 <- repdata[1:1000,]

RW <- repcreate(df = repdata2, # the data frame with all the information
                 wt = "wt", # the total weights column name
                 jkzone = "jkzones", # the jkzones column name
                 jkrep = "jkrep", # the jkreps column name
                 repwtname = "REPWT", # the desired name for the rep weights
                 reps = 50, # the number of replications
                 method = "ICILS") # the name of the method aka the study name

### No groups ----

# Simple regression - weights within df
replm(formula = Math1 ~ 1 + GENDER,
      wt = "wt", # Name of total weight column within df
      repwt = "REPWT", # Common names of replicate weights within df
      df = cbind(repdata2,RW), # Data frame
      method = "ICILS") # the name of the method aka the study name

# Simple regression - weights as a separate data frame
replm(formula = Math1 ~ 1 + GENDER,
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       method = "ICILS") # the name of the method aka the study name

# Multiple regression
replm(formula = Math1 ~ 1 + GENDER + Reading1,
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       method = "ICILS") # the name of the method aka the study name

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3))
pvs

replm(formula = Math ~ 1 + GENDER + Reading1, # Math1 now is "Math"
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       method = "ICILS") # the name of the method aka the study name

# Multiple regression - with more than one related PV variable
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3))
pvs

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       method = "ICILS") # the name of the method aka the study name

# Multiple regression - with more than UNrelated PV variables
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3)) 
pvs

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       relatedpvs = FALSE, # Unrealted PVs
       method = "ICILS") # the name of the method aka the study name


### Groups ----

# Simple regression - weights within df
replm(formula = Math1 ~ 1 + GENDER,
      wt = "wt", # Name of total weight column within df
      repwt = "REPWT", # Common names of replicate weights within df
      df = cbind(repdata2,RW), # Data frame
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name

# Simple regression - weights as a separate data frame
replm(formula = Math1 ~ 1 + GENDER,
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name

# Multiple regression
replm(formula = Math1 ~ 1 + GENDER + Reading1,
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3))
pvs

replm(formula = Math ~ 1 + GENDER + Reading1, # Math1 now is "Math"
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      pvs = pvs, # Named list
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name

# Multiple regression - with more than one related PV variable
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3))
pvs

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      pvs = pvs, # Named list
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name

# Multiple regression - with more than UNrelated PV variables
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3),
           Reading = paste0("Reading",1:3)) 
pvs

replm(formula = Math ~ 1 + GENDER + Reading, # Reading1 now is "Reading"
      wt = "wt", # Name of total weight column within df
      repwt = RW, # Data frame of weights
      df = repdata2, # Data frame
      pvs = pvs, # Named list
      relatedpvs = FALSE, # Unrealted PVs
      group = "GROUP",
      method = "ICILS") # the name of the method aka the study name