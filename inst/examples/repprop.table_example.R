# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                wt = "wt", # the total weights column name
                jkzone = "jkzones", # the jkzones column name
                jkrep = "jkrep", # the jkreps column name
                repwtname = "REPWT", # the desired name for the rep weights
                reps = 50, # the number of replications
                method = "ICILS") # the name of the method aka the study name




x = repprop(x = c("item01"),
            group = "GROUP",
            repwt = "REPWT", wt = "wt", df = cbind(repdata,RW),
            method = "ICILS")

repprop.table(x, type = "long")


repprop.table(x, type = "wide1", separateSE = TRUE)


repprop.table(x, type = "wide1", separateSE = FALSE)


repprop.table(x, type = "wide2", separateSE = TRUE)


repprop.table(x, type = "wide2", separateSE = FALSE)
