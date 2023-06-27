
### R script to construct and project giraffe IBM: This is baseline
## R version 4.2.0
# Code developed by Monica Bond, Maria Paniw, and Derek Lee 
# giraffe_IBM_core_model_sensitivity: Test sensitivity of population dynamics to each demographic rate 

# We use females only because they contribute most to population dynamics and omitting males does not affect our 
#   inference because there are always males around (polygnous, nonterritorial, etc).

# We use a seasonal (4-month) time step with loops for reproduction, aging, survival/mortality.

# We run for 50 years (seasons = 150) which is 5 giraffe generations (Muller et al. 2018, Suraud et al. 2012)

# We start with an initial population of 1000 female giraffes. Traits include the following:
#  (1) unique ID
#  (2) age in seasons (=4-month steps) 
#  (3) social community (=9 AFCs) with natal dispersal movement among them for age 13 seasons and different demographic rates
#  (4) repCount = rep countdown for mothers to ensure those with calves cannot produce another for 5 seasons

# Females can begin giving birth at age 19 seasons (= age 6 yr)

# We randomly selected females in each AFC to have calves ages 1:5 in the initial population.

# Mothers were assigned repCount according to the age of their calves 1-5 (after she has a calf at repCount=0 mom gets repCount=6  
#   and she cannot birth again until her repCount returns to 0). This is to ensure moms with little calves are not available to 
#   give birth, but every female gives birth after 5 or 6 seasons.

# We include negative density dependent survival so if AFC population > 150, S is reduced so the population does not explode

# Maximum age = 87 seasons (equates to 30 years).

# Some SA4 (age 13 seasons = 4 yr) switch AFCs for natal dispersal. We code a 0.13 probability for SA4s to transition to a new AFC 
#    based on Bond et al. 2021 (JAE). We will code zero movements between MRC and N TNP SC2&3 (1,2,3,4,5,9) with barrier perturbations

# Include effects on S of covariates of predation pressure and rain.

# clear 
rm(list = ls())

# The following 3 lines are used if the code is run on a cluster, and can be uncommented for this purpose

# args=commandArgs(TRUE) 
# foldNumber <- as.integer(args[1]) 
# totalFold <- as.integer(args[2]) 


##------Core Model with reproduction, Aging, and Mortality loops---------

# Import a data frame of the individual giraffe characteristics in the initial population
InitPopGiraffe_DF=read.csv("InitPopGiraffe.csv")

# Initial population of 1000 by age class
# 601 AFs 16+
# 61 SA4  13:15
# 66 SA3  10:12
# 72 SA2  7:9
# 84 SA1  4:6
# 116 C   1:3

# Initial population of each AFC and SC
# AFC1 = W TNP; SC2       142
# AFC2 = SW MR; SC3       133
# AFC3 = N MR-MGCA; SC3   131
# AFC4 = SE MR; SC3       131
# AFC5 = NNW TNP; SC2     119
# AFC6 = SW TNP;  SC1      93
# AFC7 = Central TNP; SC1  88
# AFC8 = SSE TNP-LGCA; SC1 88
# AFC9 = NNE TNP-LGCA; SC2 75

inds <- InitPopGiraffe_DF

inds$AFSC <- as.character(inds$AFSC)
inds$age <- as.numeric(inds$age)
colnames(inds)[3] <- "group"

#########################################################################################

### SURVIVAL FUNCTIONS ###

# Calf: apply to individuals aged 1:3 seasons

survNeo <- function(age, age2, density, beta3, rain, com.effjuv, pred.eff){  
  
  beta0 = 1  # intercept, changes with density dependence: intercept-1 is the reduction in S with density dependence
  beta1 <- 0.3  # effect of age [Lee et al. 2017]
  beta2 <- -0.005  # effect of age2 [Lee et al. 2017]
  # beta3 is season of birth [Lee et al. 2017]
  # com.effjuv is community effect [Bond et al. 2021 JWM]
  # pred.eff is higher predation in short and long and lower predation in dry in TNP-no effect in MRC [Lee et al. 2016 E&E]
  
  eq.density=150
  if(density<eq.density){
    mu.surv <- exp(beta0 + beta1*age + beta2*age2 + beta3 - 0.061*rain) # rain [Bond et al. 2023]
  }else{
    
    mu.surv <- exp(beta0 -1 + beta1*age + beta2*age2 + beta3 - 0.061*rain) # -1 reduces S at higher density [Bonenfant et al. 2009]
  }
  
  prob.surv=mu.surv/(1+mu.surv)
  
  prob.surv=prob.surv + com.effjuv + pred.eff
  prob.surv[prob.surv>1]=1
  prob.surv[prob.surv<0]=0
  
  return(prob.surv)
  
}   

# Juvenile/Subadult: apply to individuals aged 4:12 seasons

survJuv <- function(age, age2, density, beta3, rain, com.effjuv){  
  
  beta0 = 1  # intercept, changes with density dependence: intercept-0.25 is the reduction in S with density dependence
  beta1 <- 0.3  # effect of age
  beta2 <- -0.005  # effect of age2
  # beta3 is season of birth
  # com.effjuv is community effect [Bond et al. 2021 JWM]
  
  eq.density=150
  if(density<eq.density){
    mu.surv <- exp(beta0 + beta1*age + beta2*age2 + beta3 - 0.061*rain) # rain [Bond et al. 2023]
  }else{
    
    mu.surv <- exp(beta0 -0.7 + beta1*age + beta2*age2 + beta3 - 0.061*rain)  # -0.25 reduces S at higher density
  }
  
  prob.surv=mu.surv/(1+mu.surv)
  
  prob.surv=prob.surv + com.effjuv
  prob.surv[prob.surv>1]=1
  prob.surv[prob.surv<0]=0
  
  return(prob.surv)
  
}    

#########################################################################################

#########################################
### BASELINE IBM SENSITIVITY ANALYSIS ###
#########################################

### INITIALIZE ###

sens.df=NULL

parameters=1:100 # number of simulations 

# sub=as.numeric(cut(parameters,totalFold)) # un-comment if running as clusters

# par.sub=parameters[sub==foldNumber]

par.sub=parameters

groups=c(1:9) # AF social community IDs  

seasons <- 150 # how many seasons of simulation

## COVARIATES

curr.seas.char=rep(c("short","long","dry"),seasons)

rain.obs=read.csv("rainfall_now.csv")  # load observed rainfall anomalies

# Vectors of community modifications to juvenile and adult survival

com.effjuv.vec=c(0.046,0.036,0.046,0.036,0.006,-0.074,-0.034,-0.044,-0.015) # adjust juv S so 6&8 do not crash and 9&7 are more stable
com.effaf.vec=c(0.005,0,0,0.0075,0.005,0.0075,0.01,0.005,0.01) # adjust ad S so high S pops do not explode

# Vector of community and seasonal specific modifications to juvenile survival due to predation

pred.eff.vec=c(-0.1,0,0,0,-0.1,-0.1,-0.1,-0.1,-0.1)  # predation higher in TNP in long and short rains [Lee et al. 2016 E&E]
pred.eff.vecV2=c(0.03,0,0,0,0.03,0.03,0.03,0.03,0.03)  # predation lower in TNP in dry 

### BEGIN #####

scen=c("control","survNeo","survJuv","survAd","Rep","Disp")

for(sc in 1:length(scen)){
  
  AB.DATA=NULL  # create database to hold abundances
  for(pu in 1:length(par.sub)){ # first loop: SIMULATIONS 
    
    print(paste("Simulation", par.sub[pu]))
    
    SIM.DATA=NULL  # Empty files to hold results for each simulation run in
    
    DATA.G=inds
    
    DATA.G=DATA.G[,c("ID","age","repCount","season","SOB","group")] 
    
    
    for(i in 1:(seasons-1)){ # second loop: seasons
      
      
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
        
        DATA$beta3[DATA$age%in%c(1:12)&DATA$SOB%in%"short"] <- 0.13
        DATA$beta3[DATA$age%in%c(1:12)&DATA$SOB%in%"dry"] <- 0.26
        DATA$beta3[DATA$age%in%c(1:12)&is.na(DATA$beta3)] <- 0
        
        # first season survival
        
        if(nrow(DATA[DATA$age%in%c(1:3),])>0){
          
          if(sc==2){
            mean.prob.neo=survNeo(DATA$age[DATA$age%in%c(1:3)],DATA$age[DATA$age%in%c(1:3)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(1:3)],rain,com.effjuv.vec[gg], ifelse(any(i==c(3,6,9,12,15,18,21,24,27,30,
                                                                                                                                                                                                       33,36,39,42,45,48,51,54,57,
                                                                                                                                                                                                       60,63,66,69,72,75,78,81,84,
                                                                                                                                                                                                       87,90,93,96,99,102,105,108,
                                                                                                                                                                                                       111,114,117,120,123,126,129,
                                                                                                                                                                                                       132,135,138,141,144,147)),pred.eff.vecV2[gg],pred.eff.vec[gg]))
            pert.prob.neo=mean.prob.neo-0.01
            pert.prob.neo[pert.prob.neo<0]=0
            
            DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
            
          }else{
            DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=survNeo(DATA$age[DATA$age%in%c(1:3)],DATA$age[DATA$age%in%c(1:3)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(1:3)],rain,com.effjuv.vec[gg], ifelse(any(i==c(3,6,9,12,15,18,21,24,27,30,
                                                                                                                                                                                                                                                                         33,36,39,42,45,48,51,54,57,
                                                                                                                                                                                                                                                                         60,63,66,69,72,75,78,81,84,
                                                                                                                                                                                                                                                                         87,90,93,96,99,102,105,108,
                                                                                                                                                                                                                                                                         111,114,117,120,123,126,129,
                                                                                                                                                                                                                                                                         132,135,138,141,144,147)),pred.eff.vecV2[gg],pred.eff.vec[gg])),size=1)
          }
          
          
        } 
        
        # juvenile survival
        
        if(nrow(DATA[DATA$age%in%c(4:12),])>0){
          
          if(sc==3){
            mean.prob.juv= survJuv(DATA$age[DATA$age%in%c(4:12)],DATA$age[DATA$age%in%c(4:12)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(4:12)],rain,com.effjuv.vec[gg])
            
            pert.prob.juv=mean.prob.juv-0.01
            pert.prob.juv[pert.prob.juv<0]=0
            
            DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
            
          }else{
            DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=survJuv(DATA$age[DATA$age%in%c(4:12)],DATA$age[DATA$age%in%c(4:12)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(4:12)],rain,com.effjuv.vec[gg]),size=1)
          }
          
          
          
        }
        
        # adult survival
        
        if(nrow(DATA[DATA$age%in%c(13:87),])>0){
          
          if(sc==4){
            
            mean.prob.ad= survAd(curr.seas.char[i],rain,nrow(DATA[DATA$age%in%c(1:87),]),com.effaf.vec[gg])
            
            pert.prob.ad=mean.prob.ad-0.01
            pert.prob.ad[pert.prob.ad<0]=0
            
            DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)
            
            
          }else{
            
            DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=survAd(curr.seas.char[i],rain,nrow(DATA[DATA$age%in%c(1:87),]),com.effaf.vec[gg]),size=1)
            
          }
          
          
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
        
        ## Update RepCount: repCount goes down by one each time step until 0, then goes to 5 and repeat
        
        if(sc==5){
          
          DATA$repCount[DATA$repCount%in%0] <-7 
        }else{
          
          DATA$repCount[DATA$repCount%in%0] <-6
        }
        
        DATA$repCount <- DATA$repCount-1
        
        ### DISPERSAL AFTER SURVIVAL ###
        
        #how many:
        
        DATA$disp=NA
        
        # Because dispersal is at the end of a time step, after the individual aged,
        # we set age of dispersal to 14 (they start as 14 seasons old in a new group)
        
        if(length(DATA$disp[DATA$age%in%14])>0){
          
          if(sc==6){
            DATA$disp[DATA$age%in%14]=rbinom(length(DATA$disp[DATA$age%in%14]),prob=0.12,size=1)
          }else{
            DATA$disp[DATA$age%in%14]=rbinom(length(DATA$disp[DATA$age%in%14]),prob=0.13,size=1)
          }
          
          
          if (any(DATA$disp[DATA$age%in%14]%in%1)) DATA$ID[DATA$disp%in%1]<-paste(DATA$ID[DATA$disp%in%1],"Disperser")
          
          DISP=rbind(DISP,DATA[DATA$disp%in%1,c(1:6)]) # save the ones that dispersed in dispersal data
          
          DATA=DATA[!DATA$disp%in%1,] # keep the ones that didn't disperse 
          
        }
        
        
        DATA=DATA[,c("ID","age","repCount","season","SOB","group")]
        
        DISP=DISP[,c("ID","age","repCount","season","SOB","group")]
        
        
        ## Add calves to the data
        
        if(exists("CALF.DATA")){
          
          DATA=rbind(DATA,data.frame(ID=paste(c(1:length(CALF.DATA$ID)),g.name,i,par.sub[pu]), # ID keeps name of natal group, season, and simulation
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
    
    SIM.DATA$run=par.sub[pu]
    
    DATA.sub=SIM.DATA
    DATA.sub$class=NA
    DATA.sub$class[DATA.sub$age%in%c(1:3)]="C"
    DATA.sub$class[DATA.sub$age%in%c(4:6)]="SA1"
    DATA.sub$class[DATA.sub$age%in%c(7:9)]="SA2"
    DATA.sub$class[DATA.sub$age%in%c(10:12)]="SA3"
    DATA.sub$class[DATA.sub$age%in%c(13:15)]="SA4"
    DATA.sub$class[DATA.sub$age>15]="A"
    
    data.agg=aggregate(ID~group+class+TIME.SIM+run, data=DATA.sub, function(x) length(unique(x)))
    
    AB.DATA=rbind(AB.DATA,data.agg)
    
  }
  
  AB.DATA$scen=scen[sc]
  sens.df=rbind(sens.df,AB.DATA[AB.DATA$TIME.SIM>120,]) # I only save the last 30 seasons for sensitivity comparisons
}


# The output of this simulation is saved as "sens.giraffe.base.csv" and is available here: 10.6084/m9.figshare.23587563
# and can be directly accessed on FigShare for plotting