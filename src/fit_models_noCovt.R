## Date: 20,July 2023
# Author: Dr Aminath Shausan 

####################################
#this program fits a Bayesian GLM with random effects to AMR incidences
# assuming incidences follow either a Poisson, negative binomial, and zero-inflated versions of them
# assume Knorr-Held (2000) type random effects with: LCAR (for spatial), RW1 or RW2 (for temporal) random effect priors  
# and 4-types of space-time interactions with iCAR as spatial prior
# exclude fixed effects (i.e no covariates are considered)
# Perform 5-fold 6-month ahead time series type cross validations
## model fitting code is adapted from: (Goicoa et al, 2018)  
###############################

.libPaths("~/Documents/R/r-libraries")
#install.packages("usethis", lib="~/Documents/R/r-libraries")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
# #to incrase size and get rid of error ' Error: vector memory exhausted (limit reached?)
# library(usethis)
# usethis::edit_r_environ()
options(digits=10)

rm(list = ls())
set.seed(963258)


library(INLA)
library(sf)
library(spdep) #to compute neighbours of regions
library(dplyr)


####################################
## read data and structure it for cross validation
## yrs 2017 to 2022 data is used for model training and validation
## yr 2023 (Jan -June) is used for forecasting

data<- read.csv('./data/data.csv') 
str(data) #date is a character
data$date <- as.Date(data$date, "%d/%m/%Y") #"%Y-%m-%d"

data <- data[order(data$date,data$region),] 
rownames(data)<-1:nrow(data)

#get data for train and validation period (2017 to 2022)
data.req <- data %>%
  filter(year >= 2017 & year <=2022)  
unique(data.req$year) #upto 2023


S <- length(unique(data.req$region)) #22 regions
T <- length(unique(data.req$date)) # total time points 72

#add spatial, temporal, and space-time indexes
data.req <- data.req %>%
  mutate( 
        ID.area=rep(1:S,T),
        ID.month=rep(1:T,each=S),
        ID.area.month=seq(1,T*S)
  )


#fold = 5 ## change fold from 1 to 5
T.fold <- as.integer(T-6*4) #change this to T-6*4, T-6*3, T-6*2, T-6*1, T-6*0 (fold 1 - fold 5)
data.valid <- data.req%>%
  filter(ID.month <=T.fold)

h.valid <- as.integer(T.fold-5) # 
data.valid <- data.valid%>%  #assume last 6 obs missing 
  mutate(obs = ifelse(ID.month%in%c(h.valid:T.fold), NA, obs),
         #E = ifelse(ID.month%in%c(h.valid:T.fold), NA, E),
         popn = ifelse(ID.month%in%c(h.valid:T.fold), NA, popn)
       )
 


##################################################
#define spatial structure matrix of LCAR
################################################
pathShapefile = './data/shapefile/'
map <- read_sf(dsn = './data', layer ='hotspots_24regions')
#remove Christmas and Cocos islands
map <- map %>%
  filter(!(SA3_NAME21 %in%  c('Christmas Island','Cocos (Keeling) Islands'))) %>%
  arrange(SA3_NAME21) #22 by 5

unique(map$SA3_NAME21) #22
nb <- poly2nb(map) #list of 22 (nbhood)
head(nb)
nb2INLA("map.adj", nb) #list of 271
g <- inla.read.graph(filename = "map.adj") #

###########
#model: log(theta_{it}) = \mu + \alpha_i + \gamma_t + \delta_{it}
## all Q matrices are structured variance-covariance matrices
#########
#define the spatial structure matrix as LCAR
Q.alpha <- matrix(0, g$n, g$n)  
for (i in 1:g$n){
  Q.alpha[i,i]=g$nnbs[[i]] #diagonal element is number of neighbors
  Q.alpha[i,g$nbs[[i]]]=-1 #non-diagonal elements are =-1
}
Q.alpha <- as(Q.alpha,"Matrix")
Q.Leroux <- diag(S)-K.alpha  

#define the temporal structure matrix of RW1/RW2
D1 <- diff(diag(T.fold),differences=1) # 
Q.gammaRW1 <- as(t(D1)%*%D1,"Matrix") #  

D2 <- diff(diag(T.fold),differences=2)
Q.gammaRW2 <- as(t(D2)%*%D2,"Matrix") # 

##  Define hyperprior distributions 
##	- Unif(0,Inf) for sigma (std)	of random effects						                           
##	- Unif(0,1) for rho (spatial smoothing parameter)		

lunif = "expression:
    a = 1;
    b = 1;
    beta = exp(theta)/(1+exp(theta));
    logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
    log_jacobian = log(beta*(1-beta));
    return(logdens+log_jacobian)"

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

#####################################
## set to compute  posterior distributions ##
compute.patterns <- TRUE  ## 

if(compute.patterns){
  source("./src/posterior_linr_combs.R")
}else{
  all.lc <- NULL
}

##########################################################
## Define formulas and fit the models 
########################################

strategy <- "simplified.laplace"

## formula for Type I interaction and RW1/RW2 prior for time ##
f1.RW1 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE, #sum-to-zero constraint
                             hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                        f(ID.month, model="rw1", constr=TRUE, 
                          hyper=list(prec=list(prior=sdunif))) +
                        f(ID.area.month, model="iid", constr=TRUE, 
                          hyper=list(prec=list(prior=sdunif)))

A.constr1.RW2 <-  matrix(rep(1:T.fold,each=S),1,S*T.fold) # 
f1.RW2 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                           hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                        f(ID.month, model="rw2", constr=TRUE,
                          hyper=list(prec=list(prior=sdunif))) +
                        f(ID.area.month, model="iid", constr=TRUE,
                          hyper=list(prec=list(prior=sdunif)),
                          extraconstr=list(A=A.constr1.RW2,e =rep(0,nrow(A.constr1.RW2)))) 

########################################
## formula for Type II interaction and RW1/RW2 prior for time ##

R2.RW1 <- kronecker(Q.gammaRW1,diag(S))   #structured temporal, but unstructured spatial
r.def2.RW1 <- S #rank deficiency
A.constr2.RW1  <- kronecker(matrix(1,1,T.fold),diag(S)) # 50 by 1050
A.constr2.RW1  <- A.constr2.RW1 [-1,] #delete 1st row (need this, otherwise caused error)

f2.RW1 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                            hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                         f(ID.month, model="rw1", constr=TRUE,
                            hyper=list(prec=list(prior=sdunif))) +
                        f(ID.area.month, model="generic0", Cmatrix=R2.RW1 , rankdef=r.def2.RW1,
                          constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                           extraconstr=list(A=A.constr2.RW1, e=rep(0,nrow(A.constr2.RW1)))
                          #extraconstr=list(A=A.constr, e=rep(0,S)) #this caused error 
                          )

R2.RW2  <- kronecker(Q.gammaRW2,diag(S))
r.def2.RW2  <- 2*S
A.constr2.RW2  <- kronecker(matrix(1,1,T.fold),diag(S))
A.constr2.RW2  <- A.constr2.RW2 [-1,] #delete 1st row (need this, otherwise caused error)

f2.RW2 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                             hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                        f(ID.month, model="rw2", constr=TRUE,
                           hyper=list(prec=list(prior=sdunif)),
                           extraconstr=list(A=matrix(1:T.fold,1,T.fold),e=0)) +
                        f(ID.area.month, model="generic0", Cmatrix=R2.RW2 , rankdef=r.def2.RW2 ,
                          constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                          extraconstr=list(A=A.constr2.RW2, e=rep(0,nrow(A.constr2.RW2)))
                           #extraconstr=list(A=A.constr, e=rep(0,S))
                          )
########################################
## formula for Type III interaction and RW1/RW2 prior for time ##

R3.RW1 <- kronecker(diag(T.fold),Q.alpha) #unstructured temporal, but structured spatial
r.def3.RW1 <- T.fold
A.constr3.RW1 <- kronecker(diag(T.fold),matrix(1,1,S))
A.constr3.RW1 <- A.constr3.RW1[-1,]#delete 1st row (need this otherwise caused error)

f3.RW1 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                       hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                    f(ID.month, model="rw1", constr=TRUE,
                       hyper=list(prec=list(prior=sdunif))) +
                   f(ID.area.month, model="generic0", Cmatrix=R3.RW1, rankdef=r.def3.RW1,
                    constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                    extraconstr=list(A=A.constr3.RW1, e=rep(0,nrow(A.constr3.RW1)))
                    #extraconstr=list(A=A.constr, e=rep(0,T))
                    )

R3.RW2 <- kronecker(diag(T.fold),Q.alpha)
r.def3.RW2 <- T.fold
A.constr3.RW2 <- kronecker(diag(T.fold),matrix(1,1,S))
A.constr3.RW2 <- A.constr3.RW2[-1,]

f3.RW2 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                       hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                     f(ID.month, model="rw1", constr=TRUE,
                        hyper=list(prec=list(prior=sdunif))) +
                     f(ID.area.month, model="generic0", Cmatrix=R3.RW2, rankdef=r.def3.RW2,
                       constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                      extraconstr=list(A=A.constr3.RW2, e=rep(0,nrow(A.constr3.RW2)))
                      )


########################################
## formula for Type IV interaction and RW1/RW2 prior for time ##

R4.RW1 <- kronecker(Q.gammaRW1,Q.alpha) #both structured rempotal and spatial
r.def4.RW1 <- S+T.fold-1 #rank deficiency
A1.RW1 <- kronecker(matrix(1,1,T.fold),diag(S))
A2.RW1 <- kronecker(diag(T.fold),matrix(1,1,S))
  
A.constr4.RW1 <- rbind(A1.RW1[-1,],A2.RW1[-1,]) 

f4.RW1 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                  f(ID.month, model="rw1", constr=TRUE,
                     hyper=list(prec=list(prior=sdunif))) +
                f(ID.area.month, model="generic0", Cmatrix=R4.RW1, rankdef=r.def4.RW1,
                  constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                  extraconstr=list(A=A.constr4.RW1, e=rep(0,nrow(A.constr4.RW1)))  
                ) 

R4.RW2 <- kronecker(Q.gammaRW2,Q.alpha)
r.def4.RW2 <- 2*S+T.fold-2 #rank deficiency
A1.RW2 <- kronecker(matrix(1,1,T.fold),diag(S))
A2.RW2 <- kronecker(diag(T.fold),matrix(1,1,S))

A.constr4.RW2 <- rbind(A1.RW2[-1,],A2.RW2[-1,]) #delete 1st row in each matrix

f4.RW2 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                      hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                  f(ID.month, model="rw2", constr=TRUE,
                     hyper=list(prec=list(prior=sdunif))) +
                  f(ID.area.month, model="generic0", Cmatrix=R4.RW2, rankdef=r.def4.RW2,
                    constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                     extraconstr=list(A=A.constr4.RW2, e=rep(0,nrow(A.constr4.RW2)))
                  ) # 
########################################
## Fit Poisson models for Types I to IV interaction: 
#summary.fitted.values: gives the estimates for incidence proportions (theta_{it})
#to get expected observed values: compute y_{it} = E_{it}*theta_{it}
###########################################
#Type 1 interactions
poiss1.RW1.valid <- inla(f1.RW1, family="poisson", data=data.valid, E=popn, #E=E, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                       #control.inla=list(strategy=strategy),
                        control.inla=list(strategy=strategy, tolerance.step = 1e-8) ####### ADD
                      #  control.family = list(hyper = list(size = list(initial = -1)))### ADD
                       )

poiss1.RW2.valid <- inla(f1.RW2, family="poisson", data=data.valid, E=popn,
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                         control.inla=list(strategy=strategy)
                         # control.inla=list(strategy=strategy, tolerance.step = 1e-8), ####### ADD
                          #control.family = list(hyper = list(size = list(initial = -1)))### ADD
                         )
summary(poiss1.RW2.valid)

#Type 2 interactions

poiss2.RW1.valid <- inla(f2.RW1, family="poisson", data=data.valid, E=popn,
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                         control.inla=list(strategy=strategy)
                         #control.inla=list(strategy=strategy, tolerance.step = 1e-8)
                         )

summary(poiss2.RW1.valid)

poiss2.RW2.valid <- inla(f2.RW2, family="poisson", data=data.valid, E=popn, #E=E,
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                         control.inla=list(strategy=strategy))
summary(poiss2.RW2.valid)

#Type 3 interactions
poiss3.RW1.valid <- inla(f3.RW1, family="poisson", data=data.valid, E=popn,#E=E,
                           control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                           control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                          # control.inla=list(strategy=strategy)
                         control.inla=list(strategy=strategy, tolerance.step = 1e-8)
                         )

summary(poiss3.RW1.valid)

poiss3.RW2.valid <- inla(f3.RW2, family="poisson", data=data.valid, E=popn,#E=E,
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                         #control.inla=list(strategy= strategy) #
                           control.inla=list(strategy=strategy, tolerance.step = 1e-8) ####### ADD
                          # control.family = list(hyper = list(size = list(initial = -1)))### ADD
)
summary(poiss3.RW2.valid)


#Type 4 interactions
poiss4.RW1.valid <- inla(f4.RW1, family="poisson", data=data.valid, E=popn,#E=E,
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                         control.inla=list(strategy=strategy))

summary(poiss4.RW1.valid)

poiss4.RW2.valid <- inla(f4.RW2, family="poisson", data=data.valid, E=popn,#E=E,
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                         control.inla=list(strategy=strategy))
summary(poiss4.RW2.valid)

#save models (change fold number as required)
poiss.Models.fld1 <- list(
                     poiss1.RW1.valid=poiss1.RW1.valid,
                     poiss1.RW2.valid=poiss1.RW2.valid,
                     poiss2.RW1.valid=poiss2.RW1.valid,
                     poiss2.RW2.valid=poiss2.RW2.valid, 
                     poiss3.RW1.valid=poiss3.RW1.valid,
                     poiss3.RW2.valid=poiss3.RW2.valid, 
                     poiss4.RW1.valid=poiss4.RW1.valid,
                     poiss4.RW2.valid=poiss4.RW2.valid 
                    )
save(poiss.Models.fld1, file= "./results/poiss.Models.fld1.Rdata")

########################################
## negative binomials for Types I to IV interaction

###########################################
#type 1 interactions: nb

nb1.RW1.valid<-  inla(f1.RW1,
                      family = "nbinomial", data = data.valid, E=popn, 
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                     # control.inla=list(strategy=strategy)
                     control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                     control.family = list(hyper = list(size = list(initial = -1))) 
)
 

summary(nb1.RW1.valid)

nb1.RW2.valid<-  inla(f1.RW2,  
                   family = "nbinomial", data = data.valid, E=popn, 
                   control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                   #control.inla=list(strategy=strategy)
                   control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                   control.family = list(hyper = list(size = list(initial = -1))) 
)

summary(nb1.RW2.valid)

#Type 2 interactions: nb
nb2.RW1.valid<-  inla(f2.RW1,     
                      family = "nbinomial", data = data.valid, E=popn, 
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      #control.inla=list(strategy=strategy)
                      control.inla=list(strategy=strategy, tolerance.step = 1e-8)
                    #  control.family = list(hyper = list(size = list(initial = -1)))
)

summary(nb2.RW1.valid)

nb2.RW2.valid<-  inla(f2.RW2,
                      family = "nbinomial", data = data.valid, E=popn, 
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                      control.family = list(hyper = list(size = list(initial = -1))) 
)

#summary(nb2.RW2.valid)

#Type 3 interactions: nb
nb3.RW1.valid<-  inla(f3.RW1,
                      family = "nbinomial", data = data.valid, E=popn, 
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                      control.family = list(hyper = list(size = list(initial = -1))) 
)

#summary(nb3.RW1.valid)

nb3.RW2.valid<-  inla(f3.RW2,
                      family = "nbinomial", data = data.valid, E=popn, 
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy, tolerance.step = 1e-8)  
                      #control.family = list(hyper = list(size = list(initial = -1))) 
)

#summary(nb3.RW2.valid)

#Type 4 interactions: nb
nb4.RW1.valid<-  inla(f4.RW1,
                      family = "nbinomial", data = data.valid, E=popn,#
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      #control.inla=list(strategy=strategy)
                      control.inla=list(strategy=strategy, tolerance.step = 1e-8)  
                      #control.family = list(hyper = list(size = list(initial = -1))) 
)

#summary(nb4.RW1.valid)

nb4.RW2.valid<-  inla(f4.RW2,     
                      family = "nbinomial", data = data.valid, E=popn, 
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy)
                      #control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                      #control.family = list(hyper = list(size = list(initial = -1))) 
)

#summary(nb4.RW2.valid)

#save models 
nb.Models.fld1 <- list( 
                     nb1.RW1.valid=nb1.RW1.valid,
                     nb1.RW2.valid=nb1.RW2.valid,
                     nb2.RW1.valid=nb2.RW1.valid,
                     nb2.RW2.valid=nb2.RW2.valid,
                     nb3.RW1.valid=nb3.RW1.valid,
                     nb3.RW2.valid=nb3.RW2.valid,
                     nb4.RW1.valid=nb4.RW1.valid,
                     nb4.RW2.valid=nb4.RW2.valid
                  )

save(nb.Models.fld1, file= "./results/nb.Models.fld5.Rdata")

########################################
## Zero-inflated Poisson (ZIP) for Types I to IV interaction
 
###########################################
#type 1 interactions: zip

zip1.RW1.valid<-  inla(f1.RW1,                   
                      family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      #control.inla=list(strategy=strategy)
                       control.inla=list(strategy=strategy, tolerance.step = 1e-8)  
                      # control.family = list(hyper = list(size = list(initial = -1))) 
)
 
#summary(zip1.RW1.valid)

zip1.RW2.valid<-  inla(f1.RW2,
                       family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.inla=list(strategy=strategy)
                      #  control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                       # control.family = list(hyper = list(size = list(initial = -1))) 
)
 
#summary(zip1.RW2.valid)

#Type 2 interactions: nb
zip2.RW1.valid<-  inla(f2.RW1,
                      family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy)
                      # control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                      # control.family = list(hyper = list(size = list(initial = -1))) 
)

#summary(zip2.RW1.valid)

zip2.RW2.valid<-  inla(f2.RW2,
                      family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy)
)

summary(zip2.RW2.valid)

#Type 3 interactions: nb
zip3.RW1.valid<-  inla(f3.RW1,
                      family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      #control.inla=list(strategy=strategy)
                       control.inla=list(strategy=strategy, tolerance.step = 1e-8)  
                      # control.family = list(hyper = list(size = list(initial = -1))) 
)

summary(zip3.RW1.valid)

zip3.RW2.valid<-  inla(f3.RW2,
                      family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      #control.inla=list(strategy=strategy)
                       control.inla=list(strategy=strategy, tolerance.step = 1e-8)
                      # control.family = list(hyper = list(size = list(initial = -1)))
)

#summary(zip3.RW2.valid)

#Type 4 interactions: nb
zip4.RW1.valid<-  inla(f4.RW1,
                      family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy)
)

#summary(zip4.RW1.valid)

zip4.RW2.valid<-  inla(f4.RW2,    #this failed
                      family = "zeroinflatedpoisson0", data = data.valid, E=popn,  
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.inla=list(strategy=strategy)
                      # control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                      # control.family = list(hyper = list(size = list(initial = -1))) 
)

summary(zip4.RW2.valid)

#save models 
zip.Models.fld1 <- list( 
                  zip1.RW1.valid=zip1.RW1.valid,
                  zip1.RW2.valid=zip1.RW2.valid,
                  zip2.RW1.valid=zip2.RW1.valid,
                  zip2.RW2.valid=zip2.RW2.valid,
                  zip3.RW1.valid=zip3.RW1.valid,
                  zip3.RW2.valid=zip3.RW2.valid,
                  zip4.RW1.valid=zip4.RW1.valid,
                  zip4.RW2.valid=zip4.RW2.valid
)

save(zip.Models.fld1, file=  "./results/zip.Models.fld1.Rdata")


########################################
## Zero-inflated Negative Binomials (ZINB) for Types I to IV interaction
 
###########################################
#type 1 interactions: zip
 
 
zinb1.RW1.valid<-  inla(f1.RW1,      
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn,  
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.fixed = list(prec.intercept = 1),
                       control.inla=list(strategy=strategy, cmin=0,tolerance.step = 1e-8),
                      control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb1.RW1.valid)

zinb1.RW2.valid<-  inla(f1.RW2, 
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                      control.fixed = list(prec.intercept = 1),
                       control.inla=list(strategy=strategy,cmin=0, tolerance.step = 1e-8),
                       control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb1.RW2.valid)

#Type 2 interactions: nb
zinb2.RW1.valid<-  inla(f2.RW1,
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.fixed = list(prec.intercept = 1),
                       control.inla=list(strategy=strategy, cmin=0,tolerance.step = 1e-8),
                      control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb2.RW1.valid)

zinb2.RW2.valid<-  inla(f2.RW2,
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.fixed = list(prec.intercept = 1),
                        control.inla=list(strategy=strategy, cmin=0,tolerance.step = 1e-8),
                        control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb2.RW2.valid)

#Type 3 interactions: nb
zinb3.RW1.valid<-  inla(f3.RW1,
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.fixed = list(prec.intercept = 1),
                       control.inla=list(strategy=strategy, cmin=0,tolerance.step = 1e-8),
                       control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb3.RW1.valid)

zinb3.RW2.valid<-  inla(f3.RW2, 
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.fixed = list(prec.intercept = 1),
                       control.inla=list(strategy=strategy, cmin=0,tolerance.step = 1e-8),
                       control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb3.RW2.valid)

#Type 4 interactions: nb
zinb4.RW1.valid<-  inla(f4.RW1,
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.fixed = list(prec.intercept = 1),
                       control.inla=list(strategy=strategy, cmin=0,tolerance.step = 1e-8),
                       control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb4.RW1.valid)

zinb4.RW2.valid<-  inla(f4.RW2,    
                       family = "zeroinflatednbinomial0", data = data.valid, E=popn, 
                       control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),  
                       control.fixed = list(prec.intercept = 1),
                       control.inla=list(strategy=strategy, cmin=0,tolerance.step = 1e-8),
                       control.family = list(hyper = list(size = list(initial = -1)))
)

summary(zinb4.RW2.valid)

#save models 
zinb.Models.fld1 <- list(
                     zinb1.RW1.valid=zinb1.RW1.valid,
                    zinb1.RW2.valid=zinb1.RW2.valid,
                    zinb2.RW1.valid=zinb2.RW1.valid,
                  zinb2.RW2.valid=zinb2.RW2.valid,
                 zinb3.RW1.valid=zinb3.RW1.valid,
                  zinb3.RW2.valid=zinb3.RW2.valid,
                zinb4.RW1.valid=zinb4.RW1.valid,
                  zinb4.RW2.valid=zinb4.RW2.valid
)

save(zinb.Models.fld1, file= "./results/zinb.Models.fld1.Rdata")


#########################################
#prepare forecasting dataframe (assume forecast for next 6 months)
############################
h.forecast <- 6 #nor of steps
Tf.start <- as.integer(T+1) #start index for Temporal ID
Tf.end <- as.integer(T+h.forecast) #end index for Temporal ID
Tf.int.start <- as.integer(T*S+1) #start index for interaction ID
Tf.int.end <- as.integer(Tf.end*S)

##df for just forecasting period
df.sub <- data.frame(year = rep(2023, h.forecast*S), 
                     month  = rep(1:6, each=S), 
                     obs= rep(NA, h.forecast*S),
                     popn= rep(NA, h.forecast*S),
                     p = rep(NA, h.forecast*S),
                     ID.area=rep(1:S,h.forecast),
                     ID.month=rep(Tf.start:Tf.end,each=S),
                     ID.area.month=seq(Tf.int.start,Tf.int.end)
)
data.forecast <- rbind(data.req[,c("year","month","obs","popn", "p", "ID.area","ID.month","ID.area.month" )], df.sub)  
T.new <- length(unique(data.forecast$ID.month))  
 
Q.alpha <- matrix(0, g$n, g$n) 
for (i in 1:g$n){
  Q.alpha[i,i]=g$nnbs[[i]] #diagonal element is number of neighbors
  Q.alpha[i,g$nbs[[i]]]=-1 #non-diagonal elements are =-1
}
Q.alpha <- as(Q.alpha,"Matrix")
Q.Leroux <- diag(S)-Q.alpha #diagonal= -(nor.neighbors of each region -1); off diagonal =1 (if i and j are neighbors)

#define the temporal structure matrix of RW1/RW2
D1 <- diff(diag(T.new),differences=1) # 
Q.gammaRW1 <- as(t(D1)%*%D1,"Matrix") #  

R2.RW1 <- kronecker(Q.gammaRW1,diag(S))  
r.def2.RW1 <- S #rank deficiency
A.constr2.RW1  <- kronecker(matrix(1,1,T.new),diag(S)) # 50 by 1050
A.constr2.RW1  <- A.constr2.RW1 [-1,] #delete 1st row (need this, otherwise caused error)

#formular for RW1 prior
f2.RW1 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                f(ID.month, model="rw1", constr=TRUE,
                hyper=list(prec=list(prior=sdunif))) +
                  f(ID.area.month, model="generic0", Cmatrix=R2.RW1 , rankdef=r.def2.RW1,
                constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                extraconstr=list(A=A.constr2.RW1, e=rep(0,nrow(A.constr2.RW1)))
     )

#forecast using Type II RW1: Poisson model (best model)
poiss2.RW1.forecast <- inla(f2.RW1, family="poisson", data=data.forecast, E = popn,  ,
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor = TRUE),
                         control.inla=list(strategy=strategy)
                         #control.inla=list(strategy=strategy, tolerance.step = 1e-8)
                         # control.family = list(hyper = list(size = list(initial = -1)))
                         )
summary(poiss2.RW1.forecast)
#plot(poiss2.RW1.forecast$marginals.fitted.values[[1]])

#forecast using Type II RW1: NB model 
nb2.RW1.forecast<-  inla(f2.RW1,     
                      family = "nbinomial", data = data.forecast, E = popn, 
                      control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor = TRUE),  
                      #control.inla=list(strategy=strategy)
                      control.inla=list(strategy=strategy, tolerance.step = 1e-8),  
                      control.family = list(hyper = list(size = list(initial = -1))) 
)

summary(nb2.RW1.forecast)
#plot(nb2.RW1.forecast$marginals.fitted.values[[1]])

#save models 
forecast.Models <- list(poiss2.RW1.forecast=poiss2.RW1.forecast,
                        nb2.RW1.forecast=nb2.RW1.forecast
)

save(forecast.Models, file=  "./results/forecast.Models.all.Rdata")

#######################################
# fit best model( Poiss: Type 2: RW1) with data upto end of 2022 (72 months)
##use this to get posterior estimates 
################################
#Type 2 interactions
 
Q.alpha <- matrix(0, g$n, g$n) #
for (i in 1:g$n){
  Q.alpha[i,i]=g$nnbs[[i]]  
  Q.alpha[i,g$nbs[[i]]]=-1  
}
Q.alpha <- as(Q.alpha,"Matrix")
Q.Leroux <- diag(S)-Q.alpha # 

#define the temporal structure matrix of RW1 
D1 <- diff(diag(T),differences=1) # 
Q.gammaRW1 <- as(t(D1)%*%D1,"Matrix") #  

R2.RW1 <- kronecker(Q.gammaRW1,diag(S))  
r.def2.RW1 <- S #rank deficiency
A.constr2.RW1  <- kronecker(matrix(1,1,T),diag(S)) #  
A.constr2.RW1  <- A.constr2.RW1 [-1,] #delete 1st row (need this, otherwise caused error)

f2.RW1 <- obs ~ f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
                  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
                f(ID.month, model="rw1", constr=TRUE,
                  hyper=list(prec=list(prior=sdunif))) +
                f(ID.area.month, model="generic0", Cmatrix=R2.RW1 , rankdef=r.def2.RW1,
                constr=TRUE, hyper=list(prec=list(prior=sdunif)),
                extraconstr=list(A=A.constr2.RW1, e=rep(0,nrow(A.constr2.RW1)))
  )


poiss2.RW1.best<- inla(f2.RW1, family="poisson", data=data.req, E=popn, 
                         control.predictor=list(compute=TRUE, cdf=c(log(1)), link =1),
                         control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                         lincomb=all.lc,  
                         control.inla=list(strategy=strategy)
                       )

summary(poiss2.RW1.best)

save(poiss2.RW1.best, file= "./results/poiss2.RW1.best2.Rdata")

#################################### modeling rates ###################################
 
