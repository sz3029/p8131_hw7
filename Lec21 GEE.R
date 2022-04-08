# this file contains code for GEE 
# Two discrete-response LDA examples: 
# Respiratory Infection and Epilepsy
# By Gen Li, 1/2/2018
# also check geepack

library(gee) # Note: data need to be sorted!!! make sure measures from the same subject are together!

# load RI data
data("respiratory", package = "HSAUR2")
head(respiratory)
resp <- subset(respiratory, month > "0")
dim(resp) # 444*7
resp$baseline <- rep(subset(respiratory, month == "0")$status,rep(4, 111)) # make month 0 as another covariate: baseline
resp$nstat <- as.numeric(resp$status == "good") # 1=good, 0=poor

# fit logistic GEE, with different association structure
resp_glm <- glm(status ~ centre + treatment + gender + baseline + age, data = resp, family = "binomial")
resp_gee1 <- gee(nstat ~ centre + treatment + gender + baseline + age, data = resp, family = "binomial", id = subject,
                 corstr = "independence", scale.fix = TRUE, scale.value = 1) # no correlation btw responses, exactly equiv to GLM
resp_gee2 <- gee(nstat ~ centre + treatment + gender + baseline + age, data = resp, family = "binomial", id = subject, 
                 corstr = "exchangeable", scale.fix = TRUE, scale.value = 1)
resp_gee3 <- gee(nstat ~ centre + treatment + gender + baseline + age, data = resp, family = "binomial", id = subject,
                 corstr = "unstructured", scale.fix = FALSE) # scale parameter is phi (over dispersion)

summary(resp_glm)
summary(resp_gee1)
summary(resp_gee2)
summary(resp_gee3)




# load Epilepsy data
data("epilepsy", package = "HSAUR2")
head(epilepsy)
itp <- interaction(epilepsy$treatment, epilepsy$period)
tapply(epilepsy$seizure.rate, itp, mean) # crude check of group mean by period*treatment

# fit log linear GEE
per <- rep(log(2),nrow(epilepsy)) # 2 week per period, so log2 is the offset (for seisure rate per week)
epilepsy$period <- as.numeric(epilepsy$period)
names(epilepsy)[names(epilepsy) == "treatment"] <- "trt" 
fm <- seizure.rate ~ base + age + period + trt + offset(per)
#
summary(glm(fm, data = epilepsy, family = "poisson"))
summary(glm(fm, data = epilepsy, family = quasi(link=log,variance=mu)))
epilepsy_gee1 <- gee(fm, data = epilepsy, family = "poisson", id = subject, corstr = "independence", scale.fix = TRUE, scale.value = 1)
summary(epilepsy_gee1) # equivalent to GLM (except for robust est of variance)
#
epilepsy_gee2 <- gee(fm, data = epilepsy, family = "poisson", id = subject, corstr = "exchangeable", scale.fix = TRUE, scale.value = 1)
summary(epilepsy_gee2)
epilepsy_gee3 <- gee(fm, data = epilepsy, family = "poisson", id = subject, corstr = "exchangeable", scale.fix = FALSE)
summary(epilepsy_gee3)


# fit GEE with interaction
epilepsy_gee4 <- gee(seizure.rate ~ offset(per)+ base + age + period*trt, data = epilepsy, family = "poisson", id = subject, corstr = "exchangeable", scale.fix = TRUE, scale.value = 1)
summary(epilepsy_gee4)




