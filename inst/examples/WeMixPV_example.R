# Prepare data weights
repdata2 <- repdata
repdata2$wt1 <- repdata2$wt # weight level 1
repdata2$wt2 <- 1 # weight level 2


# Null model - with PVs
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m1 <- WeMixPV(formula = MATH ~ 1 + (1|GROUP), # Intercept varies across GROUP
             pvs = pvs, # Named list
             data = repdata2, # Data frame
             weights = c("wt1","wt2")) # Weights vector
m1

## Fixed effects
m1$fixef

## Random effects
m1$ranef

## Models for each PV
summary(m1$models)

# Multiple regression
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m2 <- WeMixPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata2, # Data frame
             weights = c("wt1","wt2")) # Weights vector
m2




