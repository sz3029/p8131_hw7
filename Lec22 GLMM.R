# This code performs GLMM on epilepsy example
###############################################################
# load Epilepsy data
data("epilepsy", package = "HSAUR2")


per <- rep(log(2),nrow(epilepsy)) # 2 week per period, so log2 is the offset (for seisure rate per week)
epilepsy$period <- as.numeric(epilepsy$period)
names(epilepsy)[names(epilepsy) == "treatment"] <- "trt"


# fit GLMM 
library(lme4)
library(nlme)
# random intercept model
ep.GLMM1 <- glmer(seizure.rate  ~ base + age + trt*period + offset(per) + (1 | subject), 
                  family = 'poisson', data = epilepsy)
summary(ep.GLMM1) # correlation of fixed effects is related to Fisher information of estimates
random.effects(ep.GLMM1)
fixed.effects(ep.GLMM1)

# random intercept and slope (for period only)
ep.GLMM2 <- glmer(seizure.rate  ~ base + age + trt*period + offset(per) + (1+period | subject), 
                  family = 'poisson', data = epilepsy)
summary(ep.GLMM2) 
# explain exp(beta)-1 as %change in seizure rate (lam1-lam0)/lam0

