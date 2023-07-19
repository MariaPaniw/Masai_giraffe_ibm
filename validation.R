
### R script to validate giraffe IBM core model with empirical data
## R version 4.2.0
# Code developed by Monica Bond, Maria Paniw, and Derek Lee 
# validation: Validate core model using empirical data from the Tarangire Ecosystem, Tanzania

# set path to location where all files necessary to run the script are and load the data

# load necessary libraries
library(tidyverse)
library(ggpattern)

# load data
dat<-read_csv("validation.csv")
dat

# clean up data
dat <-dat %>%
  rename(Season=Time,
         Abundance=DerAbund,
         `Social Community`=Community) %>%
  mutate(`Social Community`=as.factor(`Social Community`))

# define colour scheme
community_cols=c("Social Community 1"="#ef476fff",
                 "Social Community 2"="#83d483ff",
                 "Social Community 3"="#118ab2ff",
                 "Social Community 4"="#f78c6bff",
                 "Social Community 5"="#06d6a0ff",
                 "Social Community 6"="#0c637fff",
                 "Social Community 7"="#ffd166ff",
                 "Social Community 8"="#0cb0a9ff",
                 "Social Community 9"="#073b4cff")


# plot communities
dat %>%
  mutate(`Social Community`=paste("Social Community", `Social Community`)) %>%
  ggplot(aes(x=Season, color=`Social Community`, fill=`Social Community`)) +
  geom_ribbon(aes(ymin=`SimLo95%CI`, ymax=`SimHi95%CI`), alpha=0.4, color=NA) +
  geom_ribbon_pattern(aes(ymin=`DerLo95%CI`, ymax=`DerHi95%CI`,
                          pattern_fill=`Social Community`), 
                      pattern_color = NA,
                      pattern = "stripe", 
                      fill = NA,
                      linewidth=0.25) +
  geom_line(aes(y=Abundance), linewidth=1) +
  scale_x_continuous(n.breaks=8) +
  scale_y_continuous(n.breaks=10) +
  scale_color_manual(values = community_cols) +
  scale_fill_manual(values = community_cols) +
  scale_pattern_fill_manual(values = community_cols) +
  facet_wrap(~`Social Community`, ncol=3) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        #panel.grid.major = element_line(colour = "#FDF6D6"),
        #panel.border = element_rect(colour = "#FDF6D6", fill=NA, size=2),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(color="black"))
