# Null model - with PVs
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m1 <- lmerPV(formula = MATH ~ 1 + (1|GROUP), # Intercept varies across GROUP
      pvs = pvs, # Named list
      data = repdata, # Data frame
      weights = repdata$wt) # Weights vector
m1

## Fixed effects
m1$fixef

## Random effects
m1$ranef

## Model for each PV
summary(m1$models)

# Multiple regression
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m2 <- lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata, # Data frame
             weights = repdata$wt) # Weights vector
m2

# Multiple regression with grandmean centering
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m3 <- lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata, # Data frame
             weights = repdata$wt,
             grandmean = c("SES","schoolSES"))
m3

# Multiple regression with groupmean centering
## Named list, with element names matching formula variables
pvs = list(MATH = paste0("Math",1:5))


m4 <- lmerPV(formula = MATH ~ 1 + GENDER + SES + schoolSES + (1|GROUP),
             pvs = pvs, # Named list
             data = repdata, # Data frame
             weights = repdata$wt,
             grandmean = "schoolSES",
             groupmean = "SES",
             group = repdata$GROUP)
m4



