# this file contains two case studies for LMM: TURBT and GCase
library(nlme)
library(ggplot2)



# Ex1: TURBT
########################################
# load data
TURBT.data=read.table("TURBT.csv",header=TRUE,sep=',')
dim(TURBT.data) # 285*4
colnames(TURBT.data)

# process data
start='12/1/2017'
days = as.Date(as.character(TURBT.data$Date.Of.Procedure), format="%m/%d/%Y")-
  as.Date(as.character(start), format="%m/%d/%Y")
days = as.numeric(days) # days since the intervention
doc=as.character(TURBT.data$Performing.Provider)
unique(doc) # 18 physicians
score=TURBT.data$Count_OfElements
after=as.numeric(days>0)

# Q1: compare pre- and post-intervention scores
boxplot(score~after) # not quite right, need to account for subject-level change 
# fit GLS with compound symmetry cov
comsym <- gls(score~after, correlation=corCompSymm(form = ~ 1|doc),  method="REML")
summary(comsym) # interpret beta
# fit LMM with random intercept
LMM1 <- lme (score ~ after, random = ~1 | doc, method='REML') 
summary (LMM1) # Note: we can obtain the corr in gls with the two standard deviation estimates in LMM
#
fixed.effects(LMM1) # fixed effects coeff 
random.effects(LMM1) # ordered random effects, BLUP (in this case, just b_i)


# Q2: trend of post-intervention scores
data1=data.frame(days,score,doc)
data2=subset(data1,days>0)
# speghetti plot
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  ggtitle("Score change after intervention")
# fit LMM with random intercept
LMM2 = lme(score ~ days,random = ~ 1 | doc,  data = data2)
summary(LMM2)
# fit LMM with random intercept and slope
LMM3 <- lme(score ~ days, random = ~ 1+ days | doc, data = data2)
summary (LMM3) 
vcov(LMM3) # fisher information inverse
random.effects(LMM3) 
# translate this into plain language:  
# who has overall better performance? 
subj1=rownames(random.effects(LMM3))[which.max(random.effects(LMM3)[,1])]
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  geom_line(data=data2[data2$doc==subj1,],aes(x=days, y=score,color=subj1))+
  ggtitle("Score change after intervention")
# who improve over time
subj2=rownames(random.effects(LMM3))[fixed.effects(LMM3)[2]+random.effects(LMM3)[,2]>0]
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  geom_line(data=data2[data2$doc==subj2,],aes(x=days, y=score,color=subj2))+
  ggtitle("Score change after intervention")

# Caution: goodness of fit check is also important, but not covered here
# some outlying curves are not well represented in this example:
ggplot(data2) + 
  geom_path(aes(x = days, y = score, group = doc)) +
  geom_line(data=data2[data2$doc=='John Naitoh MD',],aes(x=days, y=score,color="John Naitoh MD"))+
  ggtitle("Score change after intervention")
# apparently this guy has increasing trend, but not reflected in the estimated parameters (due to model restriction)





# Ex2: GCase
######################################
# load data
GCase.data=read.table("GCase.csv",header=TRUE,sep=',')
dim(GCase.data) # 1562*6
colnames(GCase.data)
table(GCase.data$DIAGNOSIS)

# process data
GCase.data$visit=as.numeric(GCase.data$CLINICAL_EVENT)-1 # BL=0, V04=1, V06=2, V08=3
GCase.data$mutation=as.numeric(GCase.data$GBA.Mutation!='N') # 1=mutation, 0=no


# Q1: GCase ~ time + sex+PD+mutation
# fit a random intercept model
LMM1 <- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation, random = ~1 | PATNO,  data = GCase.data, method='REML') 
summary (LMM1) 



# Q2: PD rate of change =? Control rate of change
# fit subgroup analysis, separate LMM for PD and control
LMM1.PD <- lme (TESTVALUE ~ visit+GENDER+mutation, random = ~1 | PATNO,  data = GCase.data, subset=DIAGNOSIS=="PD", method='REML') 
summary (LMM1.PD)# rate of change: 0.24
LMM1.Control <- lme (TESTVALUE ~ visit+GENDER+mutation, random = ~1 | PATNO,  data = GCase.data, subset=DIAGNOSIS=="Control", method='REML') 
summary (LMM1.Control)# rate of change: 0.13 (not significant)

# fit one analysis, with interaction
LMM2.1 <- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation+DIAGNOSIS*visit, random = ~1 | PATNO,  data = GCase.data, method='ML') # do NOT use REML for likelihood ratio
LMM2.2 <- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation, random = ~1 | PATNO,  data = GCase.data, method='ML')
anova(LMM2.2,LMM2.1)  # LRT of interaction (not significant)
LMM2<- lme (TESTVALUE ~ visit+GENDER+DIAGNOSIS+mutation+DIAGNOSIS*visit, random = ~1 | PATNO,  data = GCase.data, method='REML') 
summary (LMM2) # Wald test of interaction ( not significant)

# discrepency? Why? 
# -- subgroup analysis is different from unified analysis (where effect sizes in other variables are enforced to be the same)
# -- (relevant to this case) the time effect for PD group is still significant in the unified model
beta_PD=fixed.effects(LMM2)[2]+fixed.effects(LMM2)[6] # beta_visit+beta_visit*diagPD
beta_PD_std=sqrt(vcov(LMM2)[2,2]+vcov(LMM2)[6,6]+2*vcov(LMM2)[2,6]) # std(beta_visit+beta_visit*diagPD)
1-pnorm(beta_PD/beta_PD_std) # approx wald p-value
