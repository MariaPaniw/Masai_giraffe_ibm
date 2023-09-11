rm(list = ls())

### Plot the results of sensitivity analysis

library(ggplot2)

#load from GitHub

read.csv("sens.giraffe.base.increase.csv")

#total abundance

tot=aggregate(ID~run+TIME.SIM+group+scen,sum,data=df)

# mean size from season 121-150

mean.size=aggregate(ID~run+group+scen,mean,data=tot)

# CV size from season 121-150

CV <- function(x){
  (sd(x,na.rm = T)/mean(x,na.rm = T))*100
}

CV.size=aggregate(ID~run+group+scen,CV,data=tot)

mean.size$group=factor(mean.size$group)

mean.size$scen=factor(mean.size$scen,levels=c("control","survNeo","survJuv","survAd","Rep","Disp"))

ggplot(data=mean.size,aes(scen,ID))+
  geom_boxplot(fill="grey",alpha=0.7)+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Mean total abundance")+
  xlab("Perturbation")+
  facet_grid(group~.,scales = "free")+
  theme_bw(base_size = 20)+
  theme(axis.text.y=element_text(size=7))


CV.size$group=factor(CV.size$group)

CV.size$scen=factor(CV.size$scen,levels=c("control","survNeo","survJuv","survAd","Rep","Disp"))

ggplot(data=CV.size,aes(scen,ID))+
  geom_boxplot(fill="grey",alpha=0.7)+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("CV total abundance")+
  xlab("Perturbation")+
  facet_grid(group~.,scales = "free")+
  theme_bw(base_size = 20)+
  theme(axis.text.y=element_text(size=10))