# this file contains code for linear mixed effects models of the orthodontic example
# Note: we use the R pacakge nlme here. An alternative is the lme4 package.
#       The lmer function in lme4 is more suitable for modeling multiple non-nested random effects

library(nlme)
library(ggplot2)
library (lattice)
Orthodont<- read.table ("orthodontic.dat",header=TRUE)
head(Orthodont)

Orth.new <- groupedData (distance ~ age | child,data = as.data.frame (Orthodont)) # not necessary for lme
head(Orth.new)
OrthFem <- subset(Orth.new, male==0)
OrthFem[1:5,]
plot(OrthFem) # default: ordered by max resp of each child
ggplot(Orth.new, aes(age, distance, group=child)) + geom_line()


######################################
# fit a random intercept model
LMM1 <- lme (distance ~ age, random = ~1 | child,  data = OrthFem, method='REML') 
summary (LMM1) # pay attention to: random effects, fixed effects, 
#
VarCorr(LMM1) # covariance estimates for random effects and variance for residuals
LMM1$sigma # std for residuals
vcov(LMM1) # covariance for fixed effects estimates (inv fisher info)
#
fixed.effects(LMM1) # fixed effects coeff 
random.effects(LMM1) # ordered random effects, BLUP (in this case, just b_i)
fitted(LMM1) # fixed+random for each subj in each visit
OrthFem$distance-fitted(LMM1) # residuals
LMM1$residuals
#
plot(OrthFem$age[1:4],OrthFem$distance[1:4],type='b',ylim=c(19,26),col=26) # original data for subj1
points(OrthFem$age[1:4],fixed.effects(LMM1)[1]+fixed.effects(LMM1)[2]*OrthFem$age[1:4], pch=4,col=26)
lines(OrthFem$age[1:4],fitted(LMM1)[1:4], lty=5,type='b',pch=4,col=26)
random.effects(LMM1) # check random effects b_i for the first subj: -1.229
#
lines(OrthFem$age[5:8],OrthFem$distance[5:8],type='b',col=20) # original data for subj2
points(OrthFem$age[5:8],fixed.effects(LMM1)[1]+fixed.effects(LMM1)[2]*OrthFem$age[5:8], pch=4,col=20)
lines(OrthFem$age[5:8],fitted(LMM1)[5:8], lty=5,type='b',pch=4,col=20)
random.effects(LMM1) #check random effects b_i for the second subj: 0.340


# check equivalence to marginal model with compound symmetry correlation structure
summary(gls(distance~age, OrthFem, correlation=corCompSymm(form = ~ 1 |child), method="REML"))
# check rho==sigma_b^2/(sigma_b^2+sigma^2)


###########################################
# compare models (likelihood ratio test)
LMM.1 <- lme (distance ~ age, random = ~ 1 | child,  data = OrthFem, method='ML') # do NOT use REML for likelihood ratio
LMM.2 <- lme (distance ~ 1, random = ~ 1 | child,  data = OrthFem, method='ML')
anova(LMM.2,LMM.1) 





#############################################
# fit a random intercept and slope model
LMM2 <- lme (distance ~ age, random = ~ 1+ age | child, data = OrthFem)
summary (LMM2) 
#
plot(OrthFem$age[1:4],OrthFem$distance[1:4],type='b',ylim=c(19,26),col=26) # original data for child1
points(OrthFem$age[1:4],fixed.effects(LMM2)[1]+fixed.effects(LMM2)[2]*OrthFem$age[1:4], pch=4,col=26)
lines(OrthFem$age[1:4],fitted(LMM2)[1:4], lty=5,type='b',pch=4,col=26)
random.effects(LMM2) # check BLUP 
lines(OrthFem$age[5:8],OrthFem$distance[5:8],type='b',col=20) 
points(OrthFem$age[5:8],fixed.effects(LMM2)[1]+fixed.effects(LMM2)[2]*OrthFem$age[5:8], pch=4,col=20)
lines(OrthFem$age[5:8],fitted(LMM2)[5:8], lty=5,type='b',pch=4,col=20)


/* generate sas data set */

data OrthFem;
input id distance age Subject $ Sex $;
datalines;
65      21.0   8     F01 Female
66      20.0  10     F01 Female
67      21.5  12     F01 Female
68      23.0  14     F01 Female
69      21.0   8     F02 Female
70      21.5  10     F02 Female
71      24.0  12     F02 Female
72      25.5  14     F02 Female
73      20.5   8     F03 Female
74      24.0  10     F03 Female
75      24.5  12     F03 Female
76      26.0  14     F03 Female
77      23.5   8     F04 Female
78      24.5  10     F04 Female
79      25.0  12     F04 Female
80      26.5  14     F04 Female
81      21.5   8     F05 Female
82      23.0  10     F05 Female
83      22.5  12     F05 Female
84      23.5  14     F05 Female
85      20.0   8     F06 Female
86      21.0  10     F06 Female
87      21.0  12     F06 Female
88      22.5  14     F06 Female
89      21.5   8     F07 Female
90      22.5  10     F07 Female
91      23.0  12     F07 Female
92      25.0  14     F07 Female
93      23.0   8     F08 Female
94      23.0  10     F08 Female
95      23.5  12     F08 Female
96      24.0  14     F08 Female
97      20.0   8     F09 Female
98      21.0  10     F09 Female
99      22.0  12     F09 Female
100     21.5  14     F09 Female
101     16.5   8     F10 Female
102     19.0  10     F10 Female
103     19.0  12     F10 Female
104     19.5  14     F10 Female
105     24.5   8     F11 Female
106     25.0  10     F11 Female
107     28.0  12     F11 Female
108     28.0  14     F11 Female
;
run;

/* sas code for intercept only random effects */

proc mixed data=OrthFem;
class Subject;
model distance = age/solution;
random int/ subject=Subject type=un;
run;


/* sas code for both intercept and slope random effects */

proc mixed data=OrthFem;
class Subject;
model distance = age/solution;
random int age/ subject=Subject type=un;
run;

/* sas code for independent intercept and slope random effects */

proc mixed data=OrthFem;
class Subject;
model distance = age/solution;
random int / subject=Subject type=un;
random age / subject=Subject type=un;
run;


