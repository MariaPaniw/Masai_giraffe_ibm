
rm(list = ls())

library(ggplot2)
library(dplyr)

#load from: 10.6084/m9.figshare.23587563

# mean size 
mean.size=read.csv("mean.size.scen.csv")

CV.size=read.csv("CV.size.scen.csv")

mean.size$group=factor(mean.size$group)

mean.size$scen=factor(mean.size$scen,levels=c("1_control","2_town","3_disp","4_decr_antiTNP","5_decr_antiMR","6_decr_antiALL",
                                              "7_incr_antiTNP", "8_incr_antiMR", "9_incr_antiALL", "10_lose_migr", "11_lose_predator",
                                              "12_rain10","13_rain25", "14_town_disp", "15_town_decr_antiTNP","16_town_decr_antiMR",
                                              "17_town_decr_antiALL","18_town_disp_decr_antiALL", "19_town_incr_antiTNP","20_town_incr_antiMR",
                                              "21_town_incr_antiALL", "22_rain10_predmigr","23_rain25_predmigr","24_rain10_incr_antiTNP",
                                              "25_rain10_incr_antiMR","26_rain10_incr_antiALL","27_rain25_incr_antiTNP","28_rain25_incr_antiMR",
                                              "29_rain25_incr_antiALL","30_rain25_incr_antiTNP_predmigr","31_rain25_incr_antiMR_predmigr",
                                              "32_rain25_incr_antiALL_predmigr"))

levels(mean.size$scen)=c(1:32)

ggplot(data=mean.size,aes(scen,ID))+
  geom_boxplot(fill="grey",alpha=0.7)+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Mean total abundance")+
  xlab("Scenario")+
  facet_grid(group~.,scales = "free")+
  theme_bw(base_size = 20)+
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=8))


CV.size$group=factor(CV.size$group)

CV.size$scen=factor(CV.size$scen,levels=c("1_control","2_town","3_disp","4_decr_antiTNP","5_decr_antiMR","6_decr_antiALL",
                                          "7_incr_antiTNP", "8_incr_antiMR", "9_incr_antiALL", "10_lose_migr", "11_lose_predator",
                                          "12_rain10","13_rain25", "14_town_disp", "15_town_decr_antiTNP","16_town_decr_antiMR",
                                          "17_town_decr_antiALL","18_town_disp_decr_antiALL", "19_town_incr_antiTNP","20_town_incr_antiMR",
                                          "21_town_incr_antiALL", "22_rain10_predmigr","23_rain25_predmigr","24_rain10_incr_antiTNP",
                                          "25_rain10_incr_antiMR","26_rain10_incr_antiALL","27_rain25_incr_antiTNP","28_rain25_incr_antiMR",
                                          "29_rain25_incr_antiALL","30_rain25_incr_antiTNP_predmigr","31_rain25_incr_antiMR_predmigr",
                                          "32_rain25_incr_antiALL_predmigr"))

levels(CV.size$scen)=c(1:32)

ggplot(data=CV.size,aes(scen,ID))+
  geom_boxplot(fill="grey",alpha=0.7)+
  geom_point(size=0.01)+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("CV total abundance")+
  xlab("Scenario")+
  facet_grid(group~.,scales = "free")+
  theme_bw(base_size = 20)+
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=8))

### Extinction

ext.df=read.csv("ext.df.scen.csv")

ext.df.mu=aggregate(ext~group+scen,mean,data=ext.df)

ext.df.mu$scen=factor(ext.df.mu$scen,levels=c("1_control","2_town","3_disp","4_decr_antiTNP","5_decr_antiMR","6_decr_antiALL",
                                          "7_incr_antiTNP", "8_incr_antiMR", "9_incr_antiALL", "10_lose_migr", "11_lose_predator",
                                          "12_rain10","13_rain25", "14_town_disp", "15_town_decr_antiTNP","16_town_decr_antiMR",
                                          "17_town_decr_antiALL","18_town_disp_decr_antiALL", "19_town_incr_antiTNP","20_town_incr_antiMR",
                                          "21_town_incr_antiALL", "22_rain10_predmigr","23_rain25_predmigr","24_rain10_incr_antiTNP",
                                          "25_rain10_incr_antiMR","26_rain10_incr_antiALL","27_rain25_incr_antiTNP","28_rain25_incr_antiMR",
                                          "29_rain25_incr_antiALL","30_rain25_incr_antiTNP_predmigr","31_rain25_incr_antiMR_predmigr",
                                          "32_rain25_incr_antiALL_predmigr"))

levels(ext.df.mu$scen)=c(1:32)

# The variation in the boxplot is across the 9 groups 
ggplot(data=ext.df.mu,aes(scen,ext))+
  geom_boxplot(fill="grey",alpha=0.7)+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Extinction probability")+
  xlab("Scenario")+
  theme_bw(base_size = 20)+
  theme(axis.text.x=element_text(size=10))

### Total abundance 

end.abund=read.csv("tot.end.scen.csv")

end.abund$scen=factor(end.abund$scen,levels=c("1_control","2_town","3_disp","4_decr_antiTNP","5_decr_antiMR","6_decr_antiALL",
                                              "7_incr_antiTNP", "8_incr_antiMR", "9_incr_antiALL", "10_lose_migr", "11_lose_predator",
                                              "12_rain10","13_rain25", "14_town_disp", "15_town_decr_antiTNP","16_town_decr_antiMR",
                                              "17_town_decr_antiALL","18_town_disp_decr_antiALL", "19_town_incr_antiTNP","20_town_incr_antiMR",
                                              "21_town_incr_antiALL", "22_rain10_predmigr","23_rain25_predmigr","24_rain10_incr_antiTNP",
                                              "25_rain10_incr_antiMR","26_rain10_incr_antiALL","27_rain25_incr_antiTNP","28_rain25_incr_antiMR",
                                              "29_rain25_incr_antiALL","30_rain25_incr_antiTNP_predmigr","31_rain25_incr_antiMR_predmigr",
                                              "32_rain25_incr_antiALL_predmigr"))

levels(end.abund$scen)=c(1:32)



ggplot(data=end.abund,aes(scen,ID))+
  geom_boxplot(fill="grey",alpha=0.7)+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ylab("Metapopulation abundance after 50 yrs")+
  # ylim(110,160)+
  xlab("Scenario")+
  theme_bw(base_size = 20)+
  theme(axis.text.x=element_text(size=10),
  axis.title.y=element_text(size=10))

# Abundance through time

abund=read.csv("ab.time.scen.csv")


# plot the mitigating effects of increasing law enforcement on towns

sub=abund[abund$scen%in%c("1_control","2_town",
                      "19_town_incr_antiTNP","20_town_incr_antiMR","21_town_incr_antiALL"),]

sub.mu=aggregate(ID~TIME.SIM+scen,mean,data=sub)
sub.mu$UB=aggregate(ID~TIME.SIM+scen,quantile,0.975,data=sub)$ID
sub.mu$LB=aggregate(ID~TIME.SIM+scen,quantile,0.025,data=sub)$ID

sub.mu$scen=factor(sub.mu$scen,levels = c("1_control","2_town",
                                          "19_town_incr_antiTNP","20_town_incr_antiMR","21_town_incr_antiALL"))

ggplot(data=sub.mu,aes(TIME.SIM,ID,col=scen,group=scen,fill=scen))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Time")+
  theme_bw(base_size = 20)


# Climate

sub=abund[abund$scen%in%c("1_control","12_rain10","13_rain25"),]

sub.mu=aggregate(ID~TIME.SIM+scen,mean,data=sub)
sub.mu$UB=aggregate(ID~TIME.SIM+scen,quantile,0.975,data=sub)$ID
sub.mu$LB=aggregate(ID~TIME.SIM+scen,quantile,0.025,data=sub)$ID

sub.mu$scen=factor(sub.mu$scen,levels = c("1_control","12_rain10","13_rain25"))

ggplot(data=sub.mu,aes(TIME.SIM,ID,col=scen,group=scen,fill=scen))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Time")+
  theme_bw(base_size = 20)

# plotting the effects of increased or decreased law enforcement

sub=abund[abund$scen%in%c("1_control","4_decr_antiTNP","5_decr_antiMR","6_decr_antiALL",
                          "7_incr_antiTNP", "8_incr_antiMR", "9_incr_antiALL"),]

sub.mu=aggregate(ID~TIME.SIM+scen,mean,data=sub)
sub.mu$UB=aggregate(ID~TIME.SIM+scen,quantile,0.975,data=sub)$ID
sub.mu$LB=aggregate(ID~TIME.SIM+scen,quantile,0.025,data=sub)$ID

sub.mu$scen=factor(sub.mu$scen,levels = c("1_control","4_decr_antiTNP","5_decr_antiMR","6_decr_antiALL",
                                          "7_incr_antiTNP", "8_incr_antiMR", "9_incr_antiALL"))

ggplot(data=sub.mu,aes(TIME.SIM,ID,col=scen,group=scen,fill=scen))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Time")+
  theme_bw(base_size = 20) 


ggplot(data=sub.mu,aes(TIME.SIM,ID,col=scen,group=scen,fill=scen))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Time")+
  scale_fill_discrete(labels=c("Control","Less Law Enforcement TNP","Less Law Enforcement MR","Less Law Enforcement Both",
                              "More Law Enforcement TNP", "More Law Enforcement MR", "More Law Enforcement Both"))+
  theme_bw(base_size = 20)+ 
  theme(legend.title=element_blank())


# time to extinction for law enforcement scenarios; how much worse is the worst case?

sub=abund[abund$scen%in%c("6_decr_antiALL",
                          "17_town_decr_antiALL","18_town_disp_decr_antiALL"),]

sub.mu=aggregate(ID~TIME.SIM+scen,mean,data=sub)
sub.mu$UB=aggregate(ID~TIME.SIM+scen,quantile,0.975,data=sub)$ID
sub.mu$LB=aggregate(ID~TIME.SIM+scen,quantile,0.025,data=sub)$ID

sub.mu$scen=factor(sub.mu$scen,levels = c("6_decr_antiALL",
                                          "17_town_decr_antiALL","18_town_disp_decr_antiALL"))

ggplot(data=sub.mu,aes(TIME.SIM,ID,col=scen,group=scen,fill=scen))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Time")+
  theme_bw(base_size = 20)

# comparing towns and blocking movements 

sub=abund[abund$scen%in%c("1_control","2_town","3_disp","14_town_disp"),]

sub.mu=aggregate(ID~TIME.SIM+scen,mean,data=sub)
sub.mu$UB=aggregate(ID~TIME.SIM+scen,quantile,0.975,data=sub)$ID
sub.mu$LB=aggregate(ID~TIME.SIM+scen,quantile,0.025,data=sub)$ID

sub.mu$scen=factor(sub.mu$scen,levels = c("1_control","2_town","3_disp","14_town_disp"))

ggplot(data=sub.mu,aes(TIME.SIM,ID,col=scen,group=scen,fill=scen))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Time")+
  theme_bw(base_size = 20)


# losing predators 

sub=abund[abund$scen%in%c("1_control","11_lose_predator"),]

sub.mu=aggregate(ID~TIME.SIM+scen,mean,data=sub)
sub.mu$UB=aggregate(ID~TIME.SIM+scen,quantile,0.975,data=sub)$ID
sub.mu$LB=aggregate(ID~TIME.SIM+scen,quantile,0.025,data=sub)$ID

sub.mu$scen=factor(sub.mu$scen,levels = c("1_control","11_lose_predator"))

ggplot(data=sub.mu,aes(TIME.SIM,ID,col=scen,group=scen,fill=scen))+
  geom_ribbon(aes(ymin = LB, ymax = UB),col=NA,alpha=0.2)+
  geom_line()+
  ylab("Abundance")+
  xlab("Time")+
  theme_bw(base_size = 20)


