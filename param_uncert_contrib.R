
rm(list = ls())

library(ggplot2)
library(lme4)

#load from: 10.6084/m9.figshare.23587563

df=read.csv("pu.abund.giraffe.csv")

#total abundance

#per group
tot.group=aggregate(ID~run+par.samp+TIME.SIM+group,sum,data=df)

#all groups summed
tot=aggregate(ID~run+par.samp+TIME.SIM,sum,data=df)

tot.group$group=factor(tot.group$group)
tot.group$run=factor(tot.group$run)
tot.group$par.samp=factor(tot.group$par.samp)

tot$run=factor(tot$run)
tot$par.samp=factor(tot$par.samp)

# Mixed linear model

# GROUP-SPECIFIC ABUNDANCE (SEASON 121-150)
m1=lmer(log(ID)~1+(1|par.samp),data=tot.group)

summary(m1)

# % contribution of parameter uncertainty (par.samp) to total variance
summary(m1)[13]$varcor$par.samp[1]/(var(residuals(m1))+summary(m1)[13]$varcor$par.samp[1])


m2=lmer(log(ID)~group+(1|par.samp),data=tot.group) # with group as fixed effect

summary(m2)


## INCLUDING SCENARIOS

# final abundance
df=read.csv("tot.end.scen.pu.csv")


head(df)

df$par.samp=factor(df$par.samp)
df$run=factor(df$run)

# ABUNDANCE (SEASON 99)
m1=lmer(log(ID)~1+(1|par.samp),data=df)

summary(m1)

# % contribution of parameter uncertainty (par.samp) to total variance
summary(m1)[13]$varcor$par.samp[1]/(var(residuals(m1))+summary(m1)[13]$varcor$par.samp[1])


m2=lmer(log(ID)~scen+(1|par.samp),data=df) # with group as fixed effect

summary(m2)

# % contribution of parameter uncertainty (par.samp) to random variance one group effects are accounted for
summary(m2)[13]$varcor$par.samp[1]/(var(residuals(m2))+summary(m2)[13]$varcor$par.samp[1])

### Abundance through time
df=read.csv("ab.time.scen.pu.csv")

head(df)

df$par.samp=factor(df$par.samp)
df$run=factor(df$run)


m1=lmer(log(ID)~1+(1|par.samp),data=df)

summary(m1)

# % contribution of parameter uncertainty (par.samp) to total variance
summary(m1)[13]$varcor$par.samp[1]/(var(residuals(m1))+summary(m1)[13]$varcor$par.samp[1])


m2=lmer(log(ID)~scen+(1|par.samp),data=df) # with group as fixed effect

summary(m2)

summary(m2)[13]$varcor$par.samp[1]/(var(residuals(m2))+summary(m2)[13]$varcor$par.samp[1])


