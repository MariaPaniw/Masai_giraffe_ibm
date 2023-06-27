
rm(list = ls())

### Plot the results of validation

valid=read.csv("validation.csv") #get from https://github.com/MariaPaniw/Masai_giraffe_ibm
valid$group <- as.character(valid$group)

df=read.csv("abund.giraffe.base.csv")  # get from 10.6084/m9.figshare.23587563

sub=df[df$TIME.SIM%in%c(1:15),]

sub$group <- as.character(sub$group)

aggregate_sub=aggregate(ID~run+TIME.SIM+group, sum, data=sub)

summary(aggregate_sub)

library(ggplot2)
library(tidyverse)
library(dplyr)


ggplot(data=valid,aes(x = TIME.SIM,y = ID, group = group, col=group))+
  geom_line()+
  theme_bw(base_size = 20) +
  labs(col='Social Community', x='Season', y='Abundance') 


validation <- ggplot(data=valid,aes(x = TIME.SIM,y = ID, group = group, col=group))+
  geom_line(lwd=1.3)+
  theme_bw(base_size = 20) +
  labs(col='Social Community', x='Season', y='Abundance') 
 

# Again, match colors to Fig 1b
validation + scale_color_manual(values=c("#E69F00", "#663300", "#009E73", "#D55E00",
                                         "#F0E442", "#9966CC", "#00ffff", "#0072B2",
                                         "#CC79A7"))
#Simulated abundance for 15 seasons

sim=read.csv("sim_valid.csv") #get from https://github.com/MariaPaniw/Masai_giraffe_ibm

sim$group <- as.character(sim$group)

der=read.csv("der_valid.csv") #get from https://github.com/MariaPaniw/Masai_giraffe_ibm

der$group <- as.character(der$group)


valid<-ggplot(data=der,aes(x = TIME.SIM,y = ID, group = group, col=group))+
  geom_line(lwd=1.3)+
  theme_bw(base_size = 20) +
  labs(col='Social Community', x='Season', y='Abundance') 


# Again, match colors to Fig 1b
valid + scale_color_manual(values=c("#E69F00", "#663300", "#009E73", "#D55E00",
                                         "#F0E442", "#9966CC", "#00ffff", "#0072B2",
                                         "#CC79A7"))

valid + geom_ribbon(data=sim,aes(x=TIME.SIM, ymin=MIN, ymax=MAX, group=group, col=group))

sim_valid <- aggregate_sub %>%
  group_by(run) %>%
  summarise(min_y=min(ID), max_y=max(ID))

p <- ggplot(aggregate_sub, aes(TIME.SIM, y=ID, group=group, col=group))
p + geom_ribbon(data=aggregate_sub,aes(ymin=min(ID), ymax=max(ID), group=group, col=group))
p + geom_area(aes(y=ID))

ggplot(data=aggregate_sub,aes(x = TIME.SIM, ymin=min(ID), ymax=max(ID), group = group, col=group))+
  geom_ribbon(aes(x = TIME.SIM, ymin=min(ID), ymax=max(ID), group = group, col=group))+
  theme_bw(base_size = 20) +
  labs(col='Social Community', x='Season', y='Abundance')


ggplot(data=aggregate_sub,aes(x = TIME.SIM,y = ID, group = group, col=group))+
  geom_line()+
  theme_bw(base_size = 20) +
  labs(col='Social Community', x='Season', y='Abundance')+ 
  facet_wrap(facets = vars(Type))


