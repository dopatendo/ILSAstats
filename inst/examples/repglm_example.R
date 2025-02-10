# Less data for shorter example
repdata2 <- repdata[1:1000,]

# Creation of replicate weights
RW <- repcreate(df = repdata2, # the data frame with all the information
                wt = "wt", # the total weights column name
                jkzone = "jkzones", # the jkzones column name
                jkrep = "jkrep", # the jkreps column name
                repwtname = "REPWT", # the desired name for the rep weights
                reps = 50, # the number of replications
                method = "ICILS") # the name of the method aka the study name

### No groups ----

# Simple regression - default family = gaussian
repglm(formula = GENDER ~ 1 + Math1,
        family = gaussian, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        method = "ICILS") # the name of the method aka the study name


# Simple regression - change link function
repglm(formula = GENDER ~ 1 + Math1,
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        method = "ICILS") # the name of the method aka the study name

# Multiple regression
repglm(formula = GENDER ~ 1 + Math1 + Reading1,
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        method = "ICILS") # the name of the method aka the study name

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3)) 
pvs

repglm(formula = GENDER ~ 1 + Math + Reading1, # Math1 now is "Math"
        family = quasibinomial, # Link function
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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
        family = quasibinomial, # Link function
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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
        family = quasibinomial, # Link function
        wt = "wt", # Name of total weight column within df
        repwt = RW, # Data frame of weights
        df = repdata2, # Data frame
        pvs = pvs, # Named list
        relatedpvs = FALSE, # Unrealted PVs
        method = "ICILS") # the name of the method aka the study name


### Groups ----

# Simple regression - default family = gaussian
repglm(formula = GENDER ~ 1 + Math1,
       family = gaussian, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name


# Simple regression - change link function
repglm(formula = GENDER ~ 1 + Math1,
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name

# Multiple regression
repglm(formula = GENDER ~ 1 + Math1 + Reading1,
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name

# Multiple regression - with PVs
## Named list, with element names matching formula variables
pvs = list(Math = paste0("Math",1:3)) 
pvs

repglm(formula = GENDER ~ 1 + Math + Reading1, # Math1 now is "Math"
       family = quasibinomial, # Link function
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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
       family = quasibinomial, # Link function
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

repglm(formula = GENDER ~ 1 + Math + Reading, # Reading1 now is "Reading"
       family = quasibinomial, # Link function
       wt = "wt", # Name of total weight column within df
       repwt = RW, # Data frame of weights
       df = repdata2, # Data frame
       pvs = pvs, # Named list
       relatedpvs = FALSE, # Unrealted PVs
       group = "GROUP",
       method = "ICILS") # the name of the method aka the study name