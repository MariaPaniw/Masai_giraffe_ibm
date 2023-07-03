
### R script to construct and project giraffe IBM under parameter uncertainty
## R version 4.2.0
# Code developed by Monica Bond, Maria Paniw, and Derek Lee 

# options(warn=-1)

library(truncnorm)

##------Simple Core Model with reproduction, Aging, and Mortality loops---------
InitPopGiraffe_DF=read.csv("InitPopGiraffe.csv")


inds <- InitPopGiraffe_DF


inds$AFSC <- as.character(inds$AFSC)
inds$age <- as.numeric(inds$age)
colnames(inds)[3] <- "group"


#########################################################################################
# Model parameters based on Bonenfant et al.(2009), Lee et al.(2016a,b), Lee et al.(2017), Bond et al.(2021a,b) 
#   Lee & Bond(2022), Bond et al.(2023), 
# Seasonal survival estimates based on season of birth = S higher when born in dry, then long, then short [Lee et al.(2017)]
# Calf survival varies by age, season of birth, and AFC [Lee et al.(2017), Bond et al.(2021b)]

#########################################################################################

####################
### BASELINE IBM ###
####################

### INITIALIZE ###

pu.df=NULL

pu.sim=1:100

parameters=1:100 # # number of simulations 

# sub=as.numeric(cut(pu.sim,totalFold))
# 
# par.sub=pu.sim[sub==foldNumber]

par.sub=parameters

groups=c(1:9) # AF social community IDs  

seasons <- 150 # how many seasons of simulation

## COVARIATES

curr.seas.char=rep(c("short","long","dry"),seasons)

rain.obs=read.csv("rainfall_now.csv")  # load observed rainfall anomalies


### BEGIN #####
for(ex in 1:length(par.sub)){
  
  print(paste("Parameter draw", par.sub[ex]))
  
  AB.DATA=NULL  # create database to hold abundances
  
  ### SURVIVAL FUNCTIONS ###
  
  b_short=rtruncnorm(1, a=0.01, b=0.99, mean = 0.13, sd = 0.04)
  b_dry=rtruncnorm(1, a=0.01, b=0.99, mean = 0.26, sd = 0.04)
  # Calf: apply to individuals aged 1:3 
  
  survNeo <- function(age, age2, density, beta3, rain, com.effjuv, pred.eff){  
    
    beta0 = 1  # intercept, changes with density dependence: intercept-1 is the reduction in S with density dependence
    beta1 <- rnorm(1,mean = 0.3, sd = 0.072)  # effect of age [Lee et al. 2017]
    beta2 <- rnorm(1,mean = -0.005, sd = 0.004) # effect of age2 [Lee et al. 2017]
    rain.beta <- rnorm(1,mean = 0.061, sd = 0.03)
    # beta3 is season of birth [Lee et al. 2017]
    # com.effjuv is community effect [Bond et al. 2021b]
    # pred.eff is higher predation in short and long and lower predation in dry in TNP-no effect in MRC [Lee et al. 2016b]
    
    eq.density=150
    if(density<eq.density){
      mu.surv <- exp(beta0 + beta1*age + beta2*age2 + beta3 - rain.beta*rain) # rain [Bond et al. 2023]
    }else{
      
      mu.surv <- exp(beta0 -1 + beta1*age + beta2*age2 + beta3 - rain.beta*rain) # -1 reduces S at higher density [Bonenfant et al. 2009]
    }
    
    prob.surv=mu.surv/(1+mu.surv)
    
    prob.surv=prob.surv + com.effjuv + pred.eff
    prob.surv[prob.surv>1]=1
    prob.surv[prob.surv<0]=0
    
    return(prob.surv)
    
  }   
  
  # Juvenile/Subadult: apply to individuals aged 4:12 seasons
  
  survJuv <- function(age, age2, density, beta3, rain, com.effjuv){  
    
    beta0 = 1  # intercept, changes with density dependence: intercept-0.7 is the reduction in S with density dependence
    beta1 <- rnorm(1,mean = 0.3, sd = 0.072)  # effect of age 
    beta2 <- rnorm(1,mean = -0.005, sd = 0.004) # effect of age2 
    rain.beta <- rnorm(1,mean = 0.061, sd = 0.03)
    # beta3 is season of birth
    # com.effjuv is community effect [Bond et al. 2021b]
    
    eq.density=150
    if(density<eq.density){
      mu.surv <- exp(beta0 + beta1*age + beta2*age2 + beta3 - rain.beta*rain) # rain [Bond et al. 2023]
    }else{
      
      mu.surv <- exp(beta0 -0.7 + beta1*age + beta2*age2 + beta3 - rain.beta*rain)  # -0.7 reduces S at higher density
    }
    
    prob.surv=mu.surv/(1+mu.surv)
    
    prob.surv=prob.surv + com.effjuv
    prob.surv[prob.surv>1]=1
    prob.surv[prob.surv<0]=0
    
    return(prob.surv)
    
  }    
  
  # Adult: apply to individuals aged 13:87 seasons
  
  survAd <- function(season, rain, density, com.effaf){
    
    if(season%in%"short"){beta0<-4}else if(season%in%"long"){beta0<-4}else{beta0<-4} # season has no effect on AF survival
    
    rain.beta <- rnorm(1,mean = 0.466, sd = 0.107)
    eq.density=150  # density dependent survival
    if(density<eq.density){
      mu.surv <- exp(beta0 - rain.beta*rain)  # rain [Bond et al. 2023]
    }else{
      
      mu.surv <- exp(beta0 -0.5 - rain.beta*rain)   # -0.5 reduces S at higher density
    }
    
    prob.surv=mu.surv/(1+mu.surv)
    
    # com.effaf is community effect [Bond et al. 2021b]
    
    prob.surv=prob.surv + com.effaf
    prob.surv[prob.surv>1]=1
    prob.surv[prob.surv<0]=0
    
    return(prob.surv)
    
  }
  
  
  com.effjuv.vec=c(rtruncnorm(1, a=0.001, b=0.99, mean = 0.046, sd = 0.002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.036, sd = 0.002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.046, sd = 0.002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.036, sd = 0.002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.006, sd = 0.002),
                   rtruncnorm(1, a=(-0.99), b=0, mean = -0.074, sd = 0.002),rtruncnorm(1, a=(-0.99), b=0, mean = -0.034, sd = 0.002),rtruncnorm(1, a=(-0.99), b=0, mean = -0.044, sd = 0.002),rtruncnorm(1, a=(-0.99), b=0, mean = -0.015, sd = 0.002)) 
 
  # AF differs among communities [S from MARK results; Bond et al. 2021b]
  com.effaf.vec=c(rtruncnorm(1, a=0.001, b=0.99, mean = 0.005, sd = 0.0002),0,0,rtruncnorm(1, a=0.001, b=0.99, mean = 0.0075, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.005, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.0075, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.01, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.005, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.01, sd = 0.0002)) 
  
  # Vector of community and seasonal specific modifications to juvenile survival due to predation
  
  pred.eff.pu=rtruncnorm(1, a=(-0.99), b=0, mean = -0.1, sd = 0.06)
  pred.eff.vec=c(pred.eff.pu,0,0,0,pred.eff.pu,pred.eff.pu,pred.eff.pu,pred.eff.pu,pred.eff.pu)  # predation higher in TNP in long and short rains [Lee et al. 2016b]
  
  pred.eff.puV2=rtruncnorm(1, a=0.001, b=0.99, mean = 0.03, sd = 0.007)
  pred.eff.vecV2=c(pred.eff.puV2,0,0,0,pred.eff.puV2,pred.eff.puV2,pred.eff.puV2,pred.eff.puV2,pred.eff.puV2)  # predation lower in TNP in dry 
  
  
  for(pu in 1:length(parameters)){ # first loop: SIMULATIONS 
    
    SIM.DATA=NULL  # Empty files to hold results for each simulation run in
    
    DATA.G=inds
    
    DATA.G=DATA.G[,c("ID","age","repCount","season","SOB","group")] 
    
    
    for(i in 1:(seasons-1)){ # second loop: SEASONS
      
      
      DISP=NULL
      DATA.GROUP.TEMP=NULL
      
      if(curr.seas.char[i]=="short"){
        
        rain = sample(rain.obs$rainfall[rain.obs$season=="short"],1)
      }else if(curr.seas.char[i]=="long"){
        
        rain = sample(rain.obs$rainfall[rain.obs$season=="long"],1)
      }else if(curr.seas.char[i]=="dry"){
        
        rain = sample(rain.obs$rainfall[rain.obs$season=="dry"],1)
      }
      
      
      for(gg in 1:length(groups)){ # third loop: GROUP, which is social community
        
        g.name=groups[gg]
        
        ##### DEMOGRAPHIC RATES #####
        
        if(i==1){DATA<-DATA.G[DATA.G$group%in%g.name,] }else{DATA<-DATA.GROUP[DATA.GROUP$group%in%g.name,] }
        
        # SURVIVAL #
        
        # Survivors are drawn as random binomial variables for survival (0 or 1), and probability of getting a 1 is 
        # governed by the probability for age-class-specific survival with covariates
        
        DATA$surv=NA  # makes a column for survival, either 0 or 1
        
        DATA$beta3=NA # makes a column for effect of season of birth for ages below 13 seasons [Lee et al. 2017]
        
        DATA$beta3[DATA$age%in%c(1:12)&DATA$SOB%in%"short"] <- b_short
        DATA$beta3[DATA$age%in%c(1:12)&DATA$SOB%in%"dry"] <- b_dry
        DATA$beta3[DATA$age%in%c(1:12)&is.na(DATA$beta3)] <- 0
        
        # first season survival
        
        if(nrow(DATA[DATA$age%in%c(1:3),])>0){
          
          DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=survNeo(DATA$age[DATA$age%in%c(1:3)],DATA$age[DATA$age%in%c(1:3)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(1:3)],rain,com.effjuv.vec[gg], ifelse(any(i==c(3,6,9,12,15,18,21,24,27,30,
                                                                                                                                                                                                                                                                       33,36,39,42,45,48,51,54,57,
                                                                                                                                                                                                                                                                       60,63,66,69,72,75,78,81,84,
                                                                                                                                                                                                                                                                       87,90,93,96,99,102,105,108,
                                                                                                                                                                                                                                                                       111,114,117,120,123,126,129,
                                                                                                                                                                                                                                                                       132,135,138,141,144,147)),pred.eff.vecV2[gg],pred.eff.vec[gg])),size=1)
          
        } 
        
        # juvenile survival
        
        if(nrow(DATA[DATA$age%in%c(4:12),])>0){
          
          DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=survJuv(DATA$age[DATA$age%in%c(4:12)],DATA$age[DATA$age%in%c(4:12)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(4:12)],rain,com.effjuv.vec[gg]),size=1)
          
        }
        
        # adult survival
        
        if(nrow(DATA[DATA$age%in%c(13:87),])>0){
          
          DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=survAd(curr.seas.char[i],rain,nrow(DATA[DATA$age%in%c(1:87),]),com.effaf.vec[gg]),size=1)
          
        }
        
        DATA$surv[DATA$age>87] <-0  # remove individuals older than 87 seasons
        
        # REPRODUCTION #
        
        DATA$calf=NA
        
        # Make new calves
        
        DATA$calf[DATA$age%in%c(19:85)&DATA$repCount%in%0&DATA$surv%in%1] = 1
        
        if(length(DATA$calf[DATA$age%in%c(19:85)&DATA$repCount%in%0&DATA$surv%in%1])>0){
          
          DF=DATA[DATA$age%in%c(19:85)&DATA$surv%in%1&DATA$repCount%in%0,]
          
          CALF.DATA=DF[rep(row.names(DF), DF$calf),]  # adding new rows for each calf
          
          # Assign calf sex, so we can remove males
          
          CALF.DATA$Sex=NA
          CALF.DATA$Sex=sample(c("M","F"),nrow(CALF.DATA),prob=c(0.5,0.5),replace = T)
          CALF.DATA=CALF.DATA[CALF.DATA$Sex%in%"F",]  # immediately delete the male calves
          
        }  
        
        #### Save old #####
        
        TIME.SIM=i
        
        SIM.DATA=rbind(SIM.DATA,cbind(DATA,TIME.SIM))
        
        ##### New data for T + 1 #####
        
        DATA=DATA[DATA$surv%in%1,]  # next time step only includes the ones that survived
        
        DATA$age=DATA$age+1
        DATA$season=DATA$season+1
        
        ## Update RepCount: repCount goes down by one each time step until 0, then goes to 6 and repeat
        
        DATA$repCount[DATA$repCount%in%0] <-6
        DATA$repCount <- DATA$repCount-1
        
        ### DISPERSAL AFTER SURVIVAL ###
        
        #how many:
        
        DATA$disp=NA
        
        # Because dispersal is at the end of a time step, after the individual aged,
        # we set age of dispersal to 14 (they start as 14 seasons old in a new group)
        
        if(length(DATA$disp[DATA$age%in%14])>0){
          
          
          DATA$disp[DATA$age%in%14]=rbinom(length(DATA$disp[DATA$age%in%14]),prob=0.13,size=1)
          
          if (any(DATA$disp[DATA$age%in%14]%in%1)) DATA$ID[DATA$disp%in%1]<-paste(DATA$ID[DATA$disp%in%1],"Disperser")
          
          DISP=rbind(DISP,DATA[DATA$disp%in%1,c(1:6)]) # save the ones that dispersed in dispersal data
          
          DATA=DATA[!DATA$disp%in%1,] # keep the ones that didn't disperse 
          
        }
        
        
        DATA=DATA[,c("ID","age","repCount","season","SOB","group")]
        
        DISP=DISP[,c("ID","age","repCount","season","SOB","group")]
        
        
        ## Add calves to the data
        
        if(exists("CALF.DATA")){
          
          DATA=rbind(DATA,data.frame(ID=paste(c(1:length(CALF.DATA$ID)),g.name,i,parameters[pu]), # ID keeps name of natal group, season, and simulation
                                     age=1,
                                     repCount=19,
                                     season=i+1,
                                     SOB=curr.seas.char[i+1],
                                     group=g.name))
          
          rm(CALF.DATA)
          
          
        }
        
        if(nrow(DATA)<1) next  
        
        DATA.GROUP.TEMP = rbind(DATA.GROUP.TEMP,DATA)
        
      }
      
      DATA.GROUP=DATA.GROUP.TEMP # Assign new dataframe to loop through the groups
      
      # Now deal with dispersers. 
      # We randomly assign them to a different group name, done inside the group loop 
      
      if(nrow(DISP)>0){
        
        DISP$group.old=DISP$group
        for(gg in 1:length(groups)){
          
          DISP$group[DISP$group.old%in%gg]=sample(groups[groups!=gg],length(DISP$group[DISP$group.old%in%gg]),replace = T)
        }
        
        DATA.GROUP=rbind(DATA.GROUP,DISP[,-which(colnames(DISP)%in%"group.old")])
      }
      
      
    }
    
    SIM.DATA$run=parameters[pu]
    
    DATA.sub=SIM.DATA
    DATA.sub$class=NA
    DATA.sub$class[DATA.sub$age%in%c(1:3)]="C"
    DATA.sub$class[DATA.sub$age%in%c(4:6)]="SA1"
    DATA.sub$class[DATA.sub$age%in%c(7:9)]="SA2"
    DATA.sub$class[DATA.sub$age%in%c(10:12)]="SA3"
    DATA.sub$class[DATA.sub$age%in%c(13:15)]="SA4"
    DATA.sub$class[DATA.sub$age>15]="A"
    
    if(nrow(DATA.sub[DATA.sub$TIME.SIM>120,])>0){
      
      data.agg=aggregate(ID~group+class+TIME.SIM+run, data=DATA.sub[DATA.sub$TIME.SIM>120,], function(x) length(unique(x)))
      
      # IBM.DATA=rbind(IBM.DATA,SIM.DATA)
      AB.DATA=rbind(AB.DATA,data.agg)
    }
    
    
  }
  
  AB.DATA$par.samp=par.sub[ex]
  pu.df=rbind(pu.df,AB.DATA)
}


# The output of this simulation is saved as "pu.abund.giraffe.csv" and is available here: 10.6084/m9.figshare.23587563
# and can be directly accessed on FigShare for plotting


