
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

#########################################################################################

#####################
### SCENARIOS IBM ###
#####################

### INITIALIZE ###


pu.df=NULL

pu.sim=1:100

parameters=1:50 # number of simulations 

sub=as.numeric(cut(pu.sim,totalFold))

par.sub=pu.sim[sub==foldNumber]

groups=c(1:9) # AF social community IDs  

seasons <- 100 # how many seasons of simulation

## COVARIATES

curr.seas.char=rep(c("short","long","dry"),seasons)

rain.obs=read.csv("rainfall_now.csv")  # load observed rainfall anomalies 2000-2020

rain.10=read.csv("rainfall_10.csv")  # load observed rainfall anomalies
rain.25=read.csv("rainfall_25.csv")  # load observed rainfall anomalies

# Vectors of community modifications to juvenile and adult survival


### Scenario vectors

# EXPANDING HUMAN IMPACT FROM TOWNS  [covariate: Bond et al. 2021c]
# Vector for 2_town 
town.eff.vec=c(-0.00035,-0.00043,-0.00035,-0.00043,-0.00043,0,0,0,-0.00035)

# BLOCKING DISPERSAL scenario 3


# ANTI-POACHING EFECTS  [covariate: Lee et al. 2016a, Game Controlled Area diff in S]
# Vector for 4_decr_antipoach_TNP 
more.poachTP.eff.vec=c(-0.1,0,0,0,0,-0.1,-0.1,-0.1,-0.1)

# Vector for 5_decr_antipoach_MR
more.poachMR.eff.vec=c(0,-0.1,-0.1,-0.1,0,0,0,0,0)

# Vector for 6_decr_antipoach_all
more.poachAll.eff.vec=c(-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1)

# Vector for 7_incr_antipoach_TNP
less.poachTP.eff.vec=c(0.01,0,0,0,0,0.01,0.01,0.01,0.01)

# Vector for 8_incr_antipoach_MR
less.poachMR.eff.vec=c(0,0.01,0.01,0.01,0,0,0,0,0)

# Vector for 9_incr_antipoach_all
less.poachAll.eff.vec=c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)


# LOSING MIGRATORY ALTERNATIVE PREY  [Lee et al. 2016b]
# Vectors for 10_lose_migr = remove pred.eff.vec and .vecV2 and add: 
alt.prey.loss.eff.vec=c(-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03)
pred.eff.vecV3=c(0,0,0,0,0,0,0,0,0) 

# LOSING PREDATORS
# Vectors for 11_lose_predator = remove pred.eff.vec and .vec2 and add:
lion.loss.eff.vec=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)

# RESTORING MIGRATORY PREY AND PREDATORS, STRONGER PREDATOR SWAMPING EFFECT IN TNP IN DRY
pred.eff.vecV4=c(0.04,0,0,0,0.04,0.04,0.04,0.04,0.04)

# CLIMATE EFFECTS scenarios 12 and 13, 22-32


### BEGIN #####


scen.name=c("1_control","2_town","3_disp","4_decr_antiTNP","5_decr_antiMR","6_decr_antiALL",
            "7_incr_antiTNP", "8_incr_antiMR", "9_incr_antiALL", "10_lose_migr", "11_lose_predator",
            "12_rain10","13_rain25", "14_town_disp", "15_town_decr_antiTNP","16_town_decr_antiMR",
            "17_town_decr_antiALL","18_town_disp_decr_antiALL", "19_town_incr_antiTNP","20_town_incr_antiMR",
            "21_town_incr_antiALL", "22_rain10_predmigr","23_rain25_predmigr","24_rain10_incr_antiTNP",
            "25_rain10_incr_antiMR","26_rain10_incr_antiALL","27_rain25_incr_antiTNP","28_rain25_incr_antiMR",
            "29_rain25_incr_antiALL","30_rain25_incr_antiTNP_predmigr","31_rain25_incr_antiMR_predmigr",
            "32_rain25_incr_antiALL_predmigr")

scen=c(1,2,3,6,9,10,11,13,14,17,18,21,23,29,32)
for(ex in 1:length(par.sub)){
  
  print(paste("Parameter draw", par.sub[ex]))
  
  scen.df=NULL
  
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
                   rtruncnorm(1, a=(-0.99), b=0, mean = -0.074, sd = 0.002),rtruncnorm(1, a=(-0.99), b=0, mean = -0.034, sd = 0.002),rtruncnorm(1, a=(-0.99), b=0, mean = -0.044, sd = 0.002),rtruncnorm(1, a=(-0.99), b=0, mean = -0.015, sd = 0.002)) # [Bond et al. 2021b] 
  com.effaf.vec=c(rtruncnorm(1, a=0.001, b=0.99, mean = 0.005, sd = 0.0002),0,0,rtruncnorm(1, a=0.001, b=0.99, mean = 0.0075, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.005, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.0075, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.01, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.005, sd = 0.0002),rtruncnorm(1, a=0.001, b=0.99, mean = 0.01, sd = 0.0002)) # adjust ad S so high S pops do not explode
  
  # Vector of community and seasonal specific modifications to juvenile survival due to predation
  
  pred.eff.pu=rtruncnorm(1, a=(-0.99), b=0, mean = -0.1, sd = 0.06)
  pred.eff.vec=c(pred.eff.pu,0,0,0,pred.eff.pu,pred.eff.pu,pred.eff.pu,pred.eff.pu,pred.eff.pu)  # predation higher in TNP in long and short rains [Lee et al. 2016b]
  
  pred.eff.puV2=rtruncnorm(1, a=0.001, b=0.99, mean = 0.03, sd = 0.007)
  pred.eff.vecV2=c(pred.eff.puV2,0,0,0,pred.eff.puV2,pred.eff.puV2,pred.eff.puV2,pred.eff.puV2,pred.eff.puV2)  # predation lower in TNP in dry 
  
  for(sc in scen){
    
    AB.DATA=NULL  # create database to hold abundances
    for(pu in 1:length(parameters)){ # first loop: SIMULATIONS 
      
      print(paste("Simulation", parameters[pu]))
      
      SIM.DATA=NULL  # Empty files to hold results for each simulation run in
      
      DATA.G=inds
      
      DATA.G=DATA.G[,c("ID","age","repCount","season","SOB","group")] 
      
      
      for(i in 1:(seasons-1)){ # second loop: SEASONS
        
        DISP=NULL
        DATA.GROUP.TEMP=NULL
        
        if(curr.seas.char[i]=="short"){
          
          if(sc==12|sc==22|sc==24|sc==25|sc==26){  #10% more rain
            
            rain = sample(rain.10$rainfall[rain.10$season=="short"],1)
            
          }else if(sc==13|sc==23|sc==27|sc==28|sc==29
                   |sc==30|sc==31|sc==32){ #25% more rain
            
            rain = sample(rain.25$rainfall[rain.25$season=="short"],1)
            
          }else{
            
            rain = sample(rain.obs$rainfall[rain.obs$season=="short"],1)
          }
          
        }else if(curr.seas.char[i]=="long"){
          
          if(sc==12|sc==22|sc==24|sc==25|sc==26){  #10% more rain
            
            rain = sample(rain.10$rainfall[rain.10$season=="long"],1)
            
          }else if(sc==13|sc==23|sc==27|sc==28|sc==29|sc==30|sc==31|sc==32){ #25% more rain
            
            rain = sample(rain.25$rainfall[rain.25$season=="long"],1)
            
          }else{
            
            rain = sample(rain.obs$rainfall[rain.obs$season=="long"],1)
            
          }
          
        }else if(curr.seas.char[i]=="dry"){
          
          if(sc==12|sc==22|sc==24|sc==25|sc==26){  #10% more rain
            
            rain = sample(rain.10$rainfall[rain.10$season=="long"],1)
            
          }else if(sc==13|sc==23|sc==27|sc==28|sc==29|sc==30|sc==31|sc==32){ #25% more rain
            
            rain = sample(rain.25$rainfall[rain.25$season=="long"],1)
            
          }else {
            
            rain = sample(rain.obs$rainfall[rain.obs$season=="dry"],1)
            
          }
        }
        
        for(gg in 1:length(groups)){ # third loop: GROUP, which is social community
          
          g.name=groups[gg]
          
          ##### DEMOGRAPHIC RATES #####
          
          if(i==1){DATA<-DATA.G[DATA.G$group%in%g.name,] }else{DATA<-DATA.GROUP[DATA.GROUP$group%in%g.name,] }
          
          if(is.null(DATA)) next
          if(nrow(DATA)<1) next  
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
            
            # baseline survival
            
            mean.prob.neo=survNeo(DATA$age[DATA$age%in%c(1:3)],DATA$age[DATA$age%in%c(1:3)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(1:3)],rain,com.effjuv.vec[gg], 
                                  ifelse(any(i==c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147)),pred.eff.vecV2[gg],pred.eff.vec[gg]))
            
            # for alternative prey loss and predator loss scenarios
            
            mean.prob.neo.V2=survNeo(DATA$age[DATA$age%in%c(1:3)],DATA$age[DATA$age%in%c(1:3)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(1:3)],rain,com.effjuv.vec[gg],pred.eff.vecV3[gg])
            
            # for protecting alternative prey and predators
            
            mean.prob.neo.V3=survNeo(DATA$age[DATA$age%in%c(1:3)],DATA$age[DATA$age%in%c(1:3)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(1:3)],rain,com.effjuv.vec[gg], 
                                     ifelse(any(i==c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147)),pred.eff.vecV4[gg],pred.eff.vec[gg]))
            
            if(sc==2|sc==14){  # town effect 2_town
              
              pert.prob.neo=mean.prob.neo + town.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==4){ # more poaching in TNP - 4_decr_antipoach_TNP
              
              pert.prob.neo=mean.prob.neo + more.poachTP.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==5){ # more poaching in MR - 5_decr_antipoach_MR
              
              pert.prob.neo=mean.prob.neo + more.poachMR.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==6){ # more poaching in TNP and MR - 6_decr_antipoach_all
              
              pert.prob.neo=mean.prob.neo + more.poachAll.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==7|sc==24|sc==27){ # less poaching in TNP - 7_incr_antipoach_TNP
              
              pert.prob.neo=mean.prob.neo + less.poachTP.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==8|sc==25|sc==28){ # less poaching in MR - 8_incr_antipoach_MR
              
              pert.prob.neo=mean.prob.neo + less.poachMR.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==9|sc==26|sc==29){ # less poaching in TNP and MR - 9_incr_antipoach_all
              
              pert.prob.neo=mean.prob.neo + less.poachAll.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==10){ # lose migratory species - 10_lose_migr
              
              pert.prob.neo=mean.prob.neo.V2 + alt.prey.loss.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==11){ # lose predators - 11_lose_predator
              
              pert.prob.neo=mean.prob.neo.V2 + lion.loss.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==15){ # expand towns and more poaching in TNP
              
              pert.prob.neo=mean.prob.neo + town.eff.vec[gg] + more.poachTP.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==16){ # expand towns and more poaching in MRC
              
              pert.prob.neo=mean.prob.neo + town.eff.vec[gg] + more.poachMR.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==17|sc==18){ # expand towns and more poaching in TNP and MRC
              
              pert.prob.neo=mean.prob.neo + town.eff.vec[gg] + more.poachAll.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==19){ # expand towns and less poaching in TNP
              
              pert.prob.neo=mean.prob.neo + town.eff.vec[gg] + less.poachTP.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==20){ # expand towns and less poaching in MR
              
              pert.prob.neo=mean.prob.neo + town.eff.vec[gg] + less.poachMR.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==21){ # expand towns and less poaching in TNP and MR
              
              pert.prob.neo=mean.prob.neo + town.eff.vec[gg] + less.poachAll.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
              # The next scenarios pertain only to calves
              
            }else if(sc==22|23){ # Heavier or much heavier rainfall but protect migration and predators
              
              pert.prob.neo=mean.prob.neo.V3
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==30){ # Much heavier rainfall but protect migration and predators and improve antipoaching in TNP
              
              pert.prob.neo=mean.prob.neo.V3 + less.poachTP.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==31){ # Much heavier rainfall but protect migration and predators and improve antipoaching in MR
              
              pert.prob.neo=mean.prob.neo.V3 + less.poachMR.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
            }else if(sc==32){ # Much heavier rainfall but protect migration and predators and improve antipoaching in both TNP and MR
              
              pert.prob.neo=mean.prob.neo.V3 + less.poachAll.eff.vec[gg]
              pert.prob.neo[pert.prob.neo<0]=0
              pert.prob.neo[pert.prob.neo>1]=1
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=pert.prob.neo,size=1)
              
              
            }else{ # otherwise use baseline survival
              
              DATA$surv[DATA$age%in%c(1:3)]=rbinom(length(DATA$surv[DATA$age%in%c(1:3)]),prob=mean.prob.neo,size=1)  
            }
          }
          
          # juvenile survival
          
          if(nrow(DATA[DATA$age%in%c(4:12),])>0){
            
            # baseline survival
            
            mean.prob.juv=survJuv(DATA$age[DATA$age%in%c(4:12)],DATA$age[DATA$age%in%c(4:12)]^2,nrow(DATA[DATA$age%in%c(1:87),]),DATA$beta3[DATA$age%in%c(4:12)],rain,com.effjuv.vec[gg])
            
            if(sc==2|sc==14){  # town effect 2_town
              
              pert.prob.juv=mean.prob.juv + town.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==4){ # more poaching in TNP - 4_decr_antipoach_TNP
              
              pert.prob.juv=mean.prob.juv + more.poachTP.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==5){ # more poaching in MR - 5_decr_antipoach_MR
              
              pert.prob.juv=mean.prob.juv + more.poachMR.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==6){ # more poaching in TNP and MR - 6_decr_antipoach_all
              
              pert.prob.juv=mean.prob.juv + more.poachAll.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==7|sc==24|sc==27){ # less poaching in TNP - 7_incr_antipoach_TNP
              
              pert.prob.juv=mean.prob.juv + less.poachTP.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==8|sc==25|sc==28){ # less poaching in TNP - 4_incr_antipoach_MR
              
              pert.prob.juv=mean.prob.juv + less.poachTP.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==9|sc==26|sc==29){ # less poaching in TNP and MR - 9_incr_antipoach_all
              
              pert.prob.juv=mean.prob.juv + less.poachAll.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==10){ # lose migratory species - 10_lose_migr
              
              pert.prob.juv=mean.prob.juv
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==11){ # lose predators - 11_lose_predator
              
              pert.prob.juv=mean.prob.juv
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==15){ # expand towns and more poaching in TNP
              
              pert.prob.juv=mean.prob.juv + town.eff.vec[gg] + more.poachTP.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==16){ # expand towns and more poaching in MRC
              
              pert.prob.juv=mean.prob.juv + town.eff.vec[gg] + more.poachMR.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==17|sc==18){ # expand towns and more poaching in TNP and MRC
              
              pert.prob.juv=mean.prob.juv + town.eff.vec[gg] + more.poachAll.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==19){ # expand towns and less poaching in TNP
              
              pert.prob.juv=mean.prob.juv + town.eff.vec[gg] + less.poachTP.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==20){ # expand towns and less poaching in MR
              
              pert.prob.juv=mean.prob.juv + town.eff.vec[gg] + less.poachMR.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
            }else if(sc==21){ # expand towns and less poaching in TNP and MR
              
              pert.prob.juv=mean.prob.juv + town.eff.vec[gg] + less.poachAll.eff.vec[gg]
              pert.prob.juv[pert.prob.juv<0]=0
              pert.prob.juv[pert.prob.juv>1]=1
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=pert.prob.juv,size=1)
              
              # The next scenarios pertain only to calves
              
            }else{ # otherwise use baseline survival
              
              DATA$surv[DATA$age%in%c(4:12)]=rbinom(length(DATA$surv[DATA$age%in%c(4:12)]),prob=mean.prob.juv,size=1)
            }
            
          }
          
          # adult survival
          
          if(nrow(DATA[DATA$age%in%c(13:87),])>0){
            
            # baseline survival
            
            mean.prob.ad=survAd(curr.seas.char[i],rain,nrow(DATA[DATA$age%in%c(1:87),]),com.effaf.vec[gg])
            
            if(sc==2|sc==14){  # town effect 2_town
              
              pert.prob.ad=mean.prob.ad + town.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)  
              
            }else if(sc==4){ # more poaching in TNP - 4_decr_antipoach_TNP
              
              pert.prob.ad=mean.prob.ad + more.poachTP.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)  
              
            }else if(sc==5){ # more poaching in MR - 5_decr_antipoach_MR
              
              pert.prob.ad=mean.prob.ad + more.poachMR.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)  
              
            }else if(sc==6){ #  more poaching in TNP and MR - 6_decr_antipoach_all
              
              pert.prob.ad=mean.prob.ad + more.poachAll.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)  
              
            }else if(sc==7|sc==24|sc==27){ # less poaching in TNP - 7_incr_antipoach_TNP
              
              pert.prob.ad=mean.prob.ad + less.poachTP.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)  
              
            }else if(sc==8|sc==25|sc==28){ # less poaching in TNP - 8_incr_antipoach_MR
              
              pert.prob.ad=mean.prob.ad + less.poachMR.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)  
              
            }else if(sc==9|sc==26|sc==29){ # less poaching in TNP and MR - 9_incr_antipoach_all
              
              pert.prob.ad=mean.prob.ad + less.poachAll.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)  
              
            }else if(sc==15){ # expand towns and more poaching in TNP
              
              pert.prob.ad=mean.prob.ad + town.eff.vec[gg] + more.poachTP.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)
              
            }else if(sc==16){ # expand towns and more poaching in MRC
              
              pert.prob.ad=mean.prob.ad + town.eff.vec[gg] + more.poachMR.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)
              
            }else if(sc==17|sc==18){ # expand towns and more poaching in TNP and MRC
              
              pert.prob.ad=mean.prob.ad + town.eff.vec[gg] + more.poachAll.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)
              
            }else if(sc==19){ # expand towns and less poaching in TNP
              
              pert.prob.ad=mean.prob.ad + town.eff.vec[gg] + less.poachTP.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)
              
            }else if(sc==20){ # expand towns and less poaching in MR
              
              pert.prob.ad=mean.prob.ad + town.eff.vec[gg] + less.poachMR.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)
              
            }else if(sc==21){ # expand towns and less poaching in TNP and MR
              
              pert.prob.ad=mean.prob.ad + town.eff.vec[gg] + less.poachAll.eff.vec[gg]
              pert.prob.ad[pert.prob.ad<0]=0
              pert.prob.ad[pert.prob.ad>1]=1
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=pert.prob.ad,size=1)
              
              # The next scenarios pertain only to calves
              
            }else{ # otherwise use baseline survival
              
              DATA$surv[DATA$age%in%c(13:87)]=rbinom(length(DATA$surv[DATA$age%in%c(13:87)]),prob=mean.prob.ad,size=1)
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
          
          if(nrow(DATA)<1) next  
          
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
          
          if(is.null(DATA)) next
          if(nrow(DATA)<1) next  
          
          DATA.GROUP.TEMP = rbind(DATA.GROUP.TEMP,DATA)
          
        }
        
        DATA.GROUP=DATA.GROUP.TEMP # Assign new dataframe to loop through the groups
        
        # Now deal with dispersers. 
        # We randomly assign them to a different group name, done inside the group loop 
        if(!is.null(DISP)){
          
          if(nrow(DISP)>0){
            
            DISP$group.old=DISP$group
            if(sc==3|sc==14|sc==18){
              
              for(gg in 1:length(groups)){
                
                if(gg%in%c(2,3,4)){
                  groups.sub=c(2,3,4)
                  DISP$group[DISP$group.old%in%gg]=sample(groups.sub[groups.sub!=gg],length(DISP$group[DISP$group.old%in%gg]),replace = T)
                }else if(gg%in%c(1,5,6,7,8,9)){
                  
                  groups.sub=c(1,5,6,7,8,9)
                  DISP$group[DISP$group.old%in%gg]=sample(groups.sub[groups.sub!=gg],length(DISP$group[DISP$group.old%in%gg]),replace = T)
                  
                }
                
              }
            }else{
              
              for(gg in 1:length(groups)){
                
                DISP$group[DISP$group.old%in%gg]=sample(groups[groups!=gg],length(DISP$group[DISP$group.old%in%gg]),replace = T)
              }
            }
            
            
            DATA.GROUP=rbind(DATA.GROUP,DISP[,-which(colnames(DISP)%in%"group.old")])
          }
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
      
      
      if(nrow(DATA.sub[DATA.sub$TIME.SIM>78,])>0){
        
        data.agg=aggregate(ID~group+class+TIME.SIM+run, data=DATA.sub[DATA.sub$TIME.SIM>78,], function(x) length(unique(x)))
        
        # IBM.DATA=rbind(IBM.DATA,SIM.DATA)
        AB.DATA=rbind(AB.DATA,data.agg)
      }
    
      
    }
    
    AB.DATA$scen=scen.name[sc]
    scen.df=rbind(scen.df,AB.DATA[AB.DATA$TIME.SIM>1,])
    
  }
  
  scen.df$par.samp=par.sub[ex]
  pu.df=rbind(pu.df,scen.df)
}


#total abundance at time 99

tot.end.scen.pu=aggregate(ID~run+scen+par.samp,sum,data=pu.df[pu.df$TIME.SIM==99,])

#total abundance

tot=aggregate(ID~run+TIME.SIM+group+scen+par.samp,sum,data=pu.df)

# mean size 
mean.size.scen.pu=aggregate(ID~run+group+scen+par.samp,mean,data=tot)

#CV size
CV.size.scen.pu=aggregate(ID~run+group+scen,CV,data=tot)

# extinction
ext.df.scen.pu=aggregate(TIME.SIM~run+group+scen+par.samp,max,data=tot)

ext.df.scen.pu$ext=0
ext.df.scen.pu$ext[ext.df.scen.pu$TIME.SIM<99]=1

#Abundance through time

ab.time.scen.pu=aggregate(ID~run+TIME.SIM+scen+par.samp,sum,data=pu.df)

# The summary outputs of this simulation are saved as:

# tot.end.scen.pu.csv 
# ab.time.scen.pu.csv
# mean.size.scen.pu.csv 
# CV.size.scen.pu.csv 
# ext.df.scen.pu.csv

# and are available here: 10.6084/m9.figshare.23587563
# and can be directly accessed on FigShare for plotting
