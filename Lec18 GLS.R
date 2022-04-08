# this file contains code for linear marginal models for longitudinal data
# We test different covariance patterns and show how to fit model with WLS and REML
# also check lme4


library(nlme)
library(ggplot2)
opposites <- read.table("https://stats.idre.ucla.edu/stat/r/examples/alda/data/opposites_pp.txt",header=TRUE,sep=",")
head(opposites)

# spaghetti plot
p = ggplot(opposites, aes(time, opp, group=id)) + geom_line()
print(p)


# Fit different cov model with REML
###################################################
# unstructured covariance
unstruct <- gls(opp~time*ccog, opposites, correlation=corSymm(form = ~ 1 |id),  weights=varIdent(form = ~ 1| wave),method="REML")
# check ?gls ?corClasses ?corSymm
summary(unstruct) # focus on corr, var, (weight)
unstruct$modelStruct$corStruct # corr
unstruct$modelStruct$varStruct # variance:weight (variance ratios, with the first group being the reference; check ?varIdent)
unstruct$sigma # standard deviation
cov2cor(unstruct$varBeta) # fisher information -> correlation (distinguish from $modelStruct$corStruct)


# compound symmetry
comsym <- gls(opp~time,opposites, correlation=corCompSymm(form = ~ 1|id),   weights=varIdent(form = ~ 1| wave), method="REML")
summary(comsym)
corMatrix(comsym$modelStruct$corStruct)[[1]]



# AR(1)
auto1 <- gls(opp~time ,opposites, correlation=corAR1(form = ~ 1 |id), method="REML")
summary(auto1)
corMatrix(auto1$modelStruct$corStruct)[[1]]

