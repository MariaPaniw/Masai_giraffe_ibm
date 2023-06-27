
rm(list = ls())

### Plot the results of base model for 1000 simulations

library(ggplot2)
library(dplyr)

#load from: 10.6084/m9.figshare.23587563

df=read.csv("abund.giraffe.base.csv")

# Plot mean abundance over time by age class

mean.size=aggregate(ID~TIME.SIM+group+class,mean,data=df)
mean.size$UB=aggregate(ID~TIME.SIM+group+class,quantile,0.975,data=df)$ID
mean.size$LB=aggregate(ID~TIME.SIM+group+class,quantile,0.025,data=df)$ID

mean.size$group=factor(mean.size$group)

mean.size$class=factor(mean.size$class,levels=c("C","SA1","SA2","SA3","SA4","A"))

abund.class <- ggplot(data=mean.size,aes(TIME.SIM,ID,col=group,group=group,fill=group))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Season")+
  facet_grid(class~.,scales = "free")+
  theme_bw(base_size = 20)+ 
  labs(col='Social Community') + labs(fill='Social Community')

#to match Fig 1b which is a color for each social community
abund.class + scale_color_manual(values=c("#E69F00", "#663300", "#009E73", "#D55E00",
                                          "#F0E442", "#9966CC", "#00ffff", "#0072B2",
                                          "#CC79A7"))+ scale_fill_manual(values=c("#E69F00", "#663300", "#009E73", "#D55E00",
                                                                                  "#F0E442", "#9966CC", "#00ffff", "#0072B2",
                                                                                  "#CC79A7"))

#  Plot mean total abundance over time

abundance.core=aggregate(ID~run+TIME.SIM+group, sum, data=df)
mean.size.all=aggregate(ID~TIME.SIM+group,sum,data=mean.size)
mean.size.all$UB=aggregate(ID~TIME.SIM+group,quantile,0.975,data=abundance.core)$ID
mean.size.all$LB=aggregate(ID~TIME.SIM+group,quantile,0.025,data=abundance.core)$ID

mean.size.all$TIME.SIM=as.numeric(mean.size.all$TIME.SIM)

tot.abundance <- ggplot(data=mean.size.all,aes(TIME.SIM,ID,col=group,group=group,fill=group))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line(lwd=1.3)+
  ylab("Abundance")+
  xlab("Season")+
  labs(col='Social Community',fill='Social Community')+ 
  theme_bw(base_size = 20)

tot.abundance + scale_color_manual(values=c("#E69F00", "#663300", "#009E73", "#D55E00",
                                            "#F0E442", "#9966CC", "#00ffff", "#0072B2",
                                            "#CC79A7"))
