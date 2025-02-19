## Date: 20,July 2023
# Author: Dr Aminath Shausan 

####################################
#this program evaluates model fits generated from the file  'fit_models_noCovt.R'
###############################
.libPaths("~/Documents/R/r-libraries")
#install.packages("scales", lib="~/Documents/R/r-libraries")
options(digits=8)
rm(list = ls())
set.seed(963258)

#library(MASS)
library(INLA)
#library(inlatools)
library(sf)
library(spdep) #to compute neighbours of regions
#library("sf")
library(dplyr)
library(plotly)
library(tmap)
library(RColorBrewer)
library(ggplot2)

########################################################
#read data 
######################################
data<- read.csv('./data/data.csv') 
str(data) #date is a character
data$date <- as.Date(data$date, "%d/%m/%Y") #"%Y-%m-%d"

data <- data[order(data$date,data$region),] 
rownames(data)<-1:nrow(data)

#get data for the observed period (2017 to 2022)
data.req <- data %>%
  filter(year >= 2017 & year <=2022)  

S <- length(unique(data.req$region)) #22 regions
T <- length(unique(data.req$date))  # 72. months
data.req <- data.req %>%    #insert index cols
  mutate(ID.area=rep(1:S,T),
         ID.month=rep(1:T,each=S),
         ID.area.month=seq(1,T*S)
  )

########################################################
#load shape file
######################################
pathShapefile = './data/shapefile/'
map <- read_sf(dsn = './data', layer ='hotspots_24regions')
#remove Christmas and Cocos islands
map <- map %>%
  filter(!(SA3_NAME21 %in%  c('Christmas Island','Cocos (Keeling) Islands'))) %>%
  arrange(SA3_NAME21) #22 by 5

unique(map$SA3_NAME21) #22


#########################################
#plot fig2: map of avg observed AMR counts (or population)
####################################

overall.rate <- 10000*sum(data.req$obs, na.rm=T)/sum(data.req$popn, na.rm=T) #5.049
summary(data.req$rate )
summary(data.req$p )
df.spatial.avg <- data.req %>%
  group_by(region)%>%
  summarise(obs = sum(obs, na.rm = T),
            popn = sum(popn, na.rm = T),
            )%>%
  as.data.frame()
 
df.spatial.avg <- df.spatial.avg %>%
  mutate( avgObs = obs/72,
          avgPopn = (popn/72)/10000,
          rate = 10000*(df.spatial.avg$obs/df.spatial.avg$popn),
         # 
        ) %>%
  mutate(qavgObs=cut(avgObs, breaks = c(0, 5,10,20,30,40,50,Inf)),
         qavgPopn=cut(avgPopn, breaks = c(0, 1,5,8,10,12,16,Inf)),
         # qavgPopn=cut(avgPopn, breaks = c(0, 1,5,8,10,15,20,Inf)),
        #qpop=cut(popdiv, breaks = c(-Inf, 50,100,300,500,700,900,Inf))
         )

summary(df.spatial.avg$avgObs) 
summary(round(df.spatial.avg$avgPopn,0)) 
summary(round(df.spatial.avg$rate,0))

#add df.spatial.avg to map
map.rates <- map%>%
  mutate(avgObs = df.spatial.avg$avgObs,
         avgPopn = df.spatial.avg$avgPopn,
         qavgObs =df.spatial.avg$qavgObs,
         qavgPopn = df.spatial.avg$qavgPopn,
         rate = df.spatial.avg$rate
       )

map.rates$rate <- round(map.rates$rate, 0)


colr.pal <- brewer.pal(8,"RdYlGn")[8:1]
#values <- c(-Inf,5,10,15,20,25,Inf) #total obs
#values <- c(0,5,10,30,50,70,90,Inf) #total population
values.avgObs <-c(0, 5,10,20,30,40,50,Inf) 
values.avgPopn <- c(0, 1,5,8,10,12,16,Inf)
values.rate <- c(0, 5,9,15, 20,30,40,50,Inf) #c(0, 5,9,10,20,30,40,50,Inf)

# plot avg rate (spatial dimension)
fig.map.rate <- tm_shape(map.rates) +
  tm_polygons("rate", palette = colr.pal, breaks=values.rate, title="Average MRSA rate per region per 10,000 population",
              interval.closure="left", legend.reverse=F, legend.is.portrait = FALSE, border.col = "black",
              border.alpha = 0.5) +
  #tm_borders("black", alpha=.5) +
  tm_layout(legend.outside=T, legend.outside.position="bottom", legend.frame=F, legend.outside.size=0.5)+
  tm_layout(frame = FALSE, inner.margins = c(0.02, 0.08, 0.02, 0.02))
fig.map.rate

map.north <- map.rates %>%
  filter(SA3_NAME21%in%c("Darwin_Region", "Palmerston", "Litchfield"))
fig.north <- tm_shape(map.north) +
  tm_polygons("rate", palette = colr.pal, breaks=values.rate, legend.show = F, title=" ",
              interval.closure="left", legend.reverse=T, border.col = "black",
              border.alpha = 0.5) +
  tm_layout(frame = FALSE, outer.margins=c(0.02,0.01,0.02,0.01))

fig.north
rm(map.rates, map.north, df.rates)

# plot avg population (spatial dimension)
fig.map.avgPopn <- tm_shape(map.rates) +
  tm_polygons("avgPopn", palette = colr.pal, breaks=values.avgPopn, title="Average population per region (x 10,000)",
              legend.show=T, interval.closure="left", legend.reverse=F, legend.is.portrait = F) +
  tm_layout(legend.outside=T, legend.outside.position="bottom", legend.frame=F, legend.outside.size=0.5)+
  tm_layout(frame = FALSE, inner.margins = c(0.02, 0.08, 0.02, 0.02))
fig.map.avgPopn

fig.avgPopn.north <- tm_shape(map.north) +
  tm_polygons("avgPopn", palette = colr.pal, breaks=values.avgPopn, title=" ",
              legend.show=T, interval.closure="left", legend.reverse=F, legend.is.portrait = T) +
  tm_layout(legend.outside=T, legend.outside.position="bottom", legend.frame=F, legend.outside.size=0.5)+
  tm_layout(frame = FALSE, inner.margins = c(0.02, 0.08, 0.02, 0.02))
fig.avgPopn.north


########## regions with > 9  cases (overall rate) #######
map.rates[map.rates$rate > 9 ,] 

rm(df.spatial.avg) #, 
rm(map.rates)
rm(fig.map.avgObs)
rm(fig.map.avgPopn)
rm(fig.map.rate)
rm(fig_temp_obs_pop)
#########################################
#plot fig3: avg MRSA per time point
####################################

df.temp <- data.req %>%
 group_by(date) %>%
  # group_by(ID.month)%>%
  summarise(totalCases = sum(obs, na.rm=T),
            totalPopn = sum(popn, na.rm=T),
            avgObs =  sum(obs, na.rm=T)/22,
            avgPopn =  (sum(popn, na.rm=T)/22)/1000
            )
df.temp$avgRate = 10000*(df.temp$totalCases/df.temp$totalPopn)

summary(df.temp$avgRate)


x <- seq(1,72, 12)

#plot temporal trend of MARSA incident rate per 10,000 ppl
min <- min(df.temp$date) 
max <- max(df.temp$date) 
date.breaks = filter(df.temp, row_number() %% 12 == 1)  
fig_temp_obs_rate <- ggplot(df.temp, aes(x=as.Date(date, "%d/%m/%Y"))) + # aes(x=ID.month)
  geom_line(aes(y=avgRate))+
  geom_point(aes(y=avgRate))+
  scale_x_date(breaks = date.breaks$date, date_labels = "%m-%Y", limits = c(min, max), expand = c(0.01, 0.01))+
  labs(x="Time", y = "Average MRSA rate per month  \n per 10,000 population")
fig_temp_obs_rate


rm(df.temp)
rm(date.breaks)
rm(fig_temp_obs_rate)
##################################

##################################
# model choice and performance: compute avg from 5 folds 
## Manually Load model fits for each fold in the 'results' folder
#DIC, WAIC, MPL, are model choice: lower DIC and WAIC and larger MPL(ie lower log-likelihood = lower mlik values), the model is better
## The larger the MPL, the better the prediction
#################################
choice <- function(x){
   data.frame(
      DIC =x$dic$dic,       ## Deviance Information Criterion
      WAIC=x$waic$waic,   ## Watanabe-Akaike information criterion
      MPL = x$mlik[2]      #marginal predictive likelihood
  )
 
}

col.names.choice <- c("DIC1", "WAIC1", "MPL1", "DIC2", "WAIC2", "MPL2",
               "DIC3", "WAIC3", "MPL3", "DIC4", "WAIC4", "MPL4",
               "DIC5", "WAIC5", "MPL5"
)

poiss.choice <- cbind(do.call(rbind,lapply(poiss.Models.fld1, choice)), 
                      do.call(rbind,lapply(poiss.Models.fld2, choice)),
                      do.call(rbind,lapply(poiss.Models.fld3, choice)), 
                      do.call(rbind,lapply(poiss.Models.fld4, choice)), 
                      do.call(rbind,lapply(poiss.Models.fld5, choice))#[2:9,] #row1 corresponds to fit to all data
                      )

colnames(poiss.choice) <- col.names.choice
nb.choice <- cbind(do.call(rbind,lapply(nb.Models.fld1, choice)), 
                      do.call(rbind,lapply(nb.Models.fld2, choice)),
                      do.call(rbind,lapply(nb.Models.fld3, choice)), 
                      do.call(rbind,lapply(nb.Models.fld4, choice)), 
                      do.call(rbind,lapply(nb.Models.fld5, choice))#[2:9,] #row1 corresponds to fit to all data
)
colnames(nb.choice) <- col.names.choice


zip.choice <- cbind(do.call(rbind,lapply(zip.Models.fld1, choice)), 
                   do.call(rbind,lapply(zip.Models.fld2, choice)),
                   do.call(rbind,lapply(zip.Models.fld3, choice)), 
                   do.call(rbind,lapply(zip.Models.fld4, choice)), 
                   do.call(rbind,lapply(zip.Models.fld5, choice))#[2:9,] #row1 corresponds to fit to all data
)
colnames(zip.choice) <- col.names.choice
zinb.choice <- cbind(do.call(rbind,lapply(zinb.Models.fld1, choice)), 
                    do.call(rbind,lapply(zinb.Models.fld2, choice)),
                    do.call(rbind,lapply(zinb.Models.fld3, choice)), 
                    do.call(rbind,lapply(zinb.Models.fld4, choice)), 
                    do.call(rbind,lapply(zinb.Models.fld5, choice))#[2:9,] #row1 corresponds to fit to all data
)
colnames(zinb.choice) <- col.names.choice


df.choice <- rbind(poiss.choice, nb.choice,  zip.choice, zinb.choice) #combine
df.choice <- df.choice %>%                                         #compute avgs 
      mutate(avgDIC = round(rowMeans(subset(df.choice, select = c(DIC1, DIC2, DIC3, DIC4, DIC5))),4),
             avgWAIC= round(rowMeans(subset(df.choice, select = c(WAIC1, WAIC2, WAIC3, WAIC4, WAIC5))),4),
             avgMPL = round(rowMeans(subset(df.choice, select = c(MPL1, MPL2, MPL3, MPL4, MPL5))),4)
      )
write.csv(df.choice, file= "./results/df.choice.csv")   

########################################################
## check model performance using model choice matrics for each fold
#############################################

idx.fold <- c(as.integer(T-6*4), as.integer(T-6*3), as.integer(T-6*2), as.integer(T-6*1), as.integer(T-6*0)) #index for fold dfs

perform <- function(x,i){
  #i=1
    print(idx.fold[i])
    data.valid <- data.req%>%
      filter(ID.month <=idx.fold[i])
      data.valid$RR <- x$summary.fitted.values[, "0.5quant"] #estimated rate
      data.valid$obs.est <-  data.valid$RR*data.valid$E #expected observed value
  
      df_perform <- data.valid %>%
        filter(ID.month >=idx.fold[i]-5)%>%
        na.omit(df_perform, cols = "obs") 
 
    data.frame(RMSE = sqrt(mean((df_perform$obs - df_perform$obs.est)^2 , na.rm=TRUE)),
             MAE = mean(abs(df_perform$obs - df_perform$obs.est), na.rm=TRUE),
             r= cor(df_perform$obs, df_perform$obs.est,  method = 'pearson')
  ) 
}


col.names.perfm <- c("RMSE1", "MAE1", "r1", "RMSE2", "MAE2", "r2","RMSE3", "MAE3", "r3",
                      "RMSE4", "MAE4", "r4","RMSE5", "MAE5", "r5"
)
poiss.perfm  <- cbind(do.call(rbind,lapply(poiss.Models.fld1,i=1, perform)),
                do.call(rbind,lapply(poiss.Models.fld2,i=2, perform)), 
                do.call(rbind,lapply(poiss.Models.fld3,i=3, perform)), 
                do.call(rbind,lapply(poiss.Models.fld4,i=4, perform)) ,
                do.call(rbind,lapply(poiss.Models.fld5 ,i=5, perform))#[2:9,]
)
colnames(poiss.perfm) <- col.names.perfm

nb.perfm  <- cbind(do.call(rbind,lapply(nb.Models.fld1,i=1, perform)),
                      do.call(rbind,lapply(nb.Models.fld2,i=2, perform)), 
                      do.call(rbind,lapply(nb.Models.fld3,i=3, perform)), 
                      do.call(rbind,lapply(nb.Models.fld4,i=4, perform)) ,
                      do.call(rbind,lapply(nb.Models.fld5 ,i=5, perform))#[2:9,]
)
colnames(nb.perfm) <- col.names.perfm

zip.perfm  <- cbind(do.call(rbind,lapply(zip.Models.fld1,i=1, perform)),
                   do.call(rbind,lapply(zip.Models.fld2,i=2, perform)), 
                   do.call(rbind,lapply(zip.Models.fld3,i=3, perform)), 
                   do.call(rbind,lapply(zip.Models.fld4,i=4, perform)) ,
                   do.call(rbind,lapply(zip.Models.fld5 ,i=5, perform))#[2:9,]
)
colnames(zip.perfm) <- col.names.perfm
zinb.perfm  <- cbind(do.call(rbind,lapply(zinb.Models.fld1,i=1, perform)),
                    do.call(rbind,lapply(zinb.Models.fld2,i=2, perform)), 
                    do.call(rbind,lapply(zinb.Models.fld3,i=3, perform)), 
                    do.call(rbind,lapply(zinb.Models.fld4,i=4, perform)) ,
                    do.call(rbind,lapply(zinb.Models.fld5 ,i=5, perform))#[2:9,]
)
colnames(zinb.perfm) <- col.names.perfm



df.perfm <- rbind(poiss.perfm, nb.perfm,  zip.perfm, zinb.perfm)
df.perfm <- df.perfm %>%                                         #compute avgs 
  mutate(avgRMSE = round(rowMeans(subset(df.perfm, select = c(RMSE1, RMSE2, RMSE3, RMSE4, RMSE5))),4),
         avgMAE= round(rowMeans(subset(df.perfm, select = c(MAE1, MAE2, MAE3, MAE4, MAE5))),4),
         avgr = round(rowMeans(subset(df.perfm, select = c(r1, r2, r3, r4, r5))),4)
  )
 
write.csv(df.perfm, file= "./results/df.perfm.csv")  #save df




df.ferfm <- rbind(bin.perfm, nb.perfm, zibin.perfm)
df.ferfm <- df.ferfm %>%                                         #compute avgs 
  mutate(avgRMSE = round(rowMeans(subset(df.ferfm, select = c(RMSE1, RMSE2, RMSE3, RMSE4, RMSE5))),4),
         avgMAE= round(rowMeans(subset(df.ferfm, select = c(MAE1, MAE2, MAE3, MAE4, MAE5))),4),
         avgr = round(rowMeans(subset(df.ferfm, select = c(r1, r2, r3, r4, r5))),4)
  )

write.csv(df.ferfm, file= "./results/df.perfm.csv")  #save 


##################################
# dispersion checks: Poisson and NB models for Type2:RW1 
##Note:  If dispersion index for Poisson >1, implies overdispersion; 
# if dispersion parameter for NB is above 0, implies significant overdispersion
# for NB models: overdispersion is related to variance being higher than 'normal'.  
#             high value of 'size' gives lower overdispersion, so its like size = 1/overdispersion
#################################
#compute deviance index for poission models 

poiss.fld1 <- poiss.Models.fld1$poiss2.RW1.valid 
poiss.fld2 <- poiss.Models.fld2$poiss2.RW1.valid 
poiss.fld3 <- poiss.Models.fld3$poiss2.RW1.valid 
poiss.fld4 <- poiss.Models.fld4$poiss2.RW1.valid 
poiss.fld5 <- poiss.Models.fld5$poiss2.RW1.valid 

idx.fold <- c(as.integer(T-6*4), as.integer(T-6*3), as.integer(T-6*2), as.integer(T-6*1), as.integer(T-6*0)) #index for fold dfs
get.dev.index <- function(x,i){
  
       data.valid <- data.req%>%
        filter(ID.month <=idx.fold[i])
      
       y.hat  <- x$summary.fitted.values[, "0.5quant"] *data.valid$popn #estimated rate
       D <- sum(data.valid$obs*log(data.valid$obs/y.hat) -(data.valid$obs - y.hat),na.rm = T) # 
       df <- length(y.hat) - x$dic$p.eff #degreesof freedom
    dev.fold <- D/(df-1)  
    
    return(dev.fold)
}
dev.indx <- c(get.dev.index(poiss.fld1,1), get.dev.index(poiss.fld2,2),
              get.dev.index(poiss.fld3,3), get.dev.index(poiss.fld4,4), get.dev.index(poiss.fld5,5))
avg.dev.inx = mean(dev.indx) 

# get summary of dispersion parameter for NB models

nb.fld1 <- nb.Models.fld1$nb2.RW1.valid
nb.fld2 <- nb.Models.fld2$nb2.RW1.valid 
nb.fld3 <- nb.Models.fld3$nb2.RW1.valid 
nb.fld4 <- nb.Models.fld4$nb2.RW1.valid 
nb.fld5 <- nb.Models.fld5$nb2.RW1.valid 

size <-nb.fld1$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)` #marginal distribution of size param 

get.nb.overdispersion <- function(x){
         size =   x$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)` #get distribution of 1/overdispersion
         overdisp.dis <- inla.tmarginal(function(s) 1/s, size)    ##transform 1/overdispersion 
         return(overdisp.dis)
}

nb.overdisp <-  rbind(data.frame(inla.zmarginal(get.nb.overdispersion(nb.fld1))), 
                      data.frame(inla.zmarginal(get.nb.overdispersion(nb.fld2))),
                      data.frame(inla.zmarginal(get.nb.overdispersion(nb.fld3))),
                      data.frame(inla.zmarginal(get.nb.overdispersion(nb.fld4))),
                      data.frame(inla.zmarginal(get.nb.overdispersion(nb.fld5)))
) 

nb.overdisp <- cbind(fold = 1:5, nb.overdisp)

df.overdispersion <- cbind(nb.overdisp, poiss.dev = dev.indx)
round(colMeans(df.overdispersion[, c("mean", "quant0.025", "quant0.975", "poiss.dev")]),3)
write.csv(df.overdispersion, file=  "./results/df.overdispersion.csv", row.names = F)


##################################
# Posterior estimates of parameters from best model: Poisson  for Type2:RW1 
#################################

Model <- poiss2.RW1.best
summary(Model)
#Fixed effects:Intercept
Model$summary.fixed #gives summary of fixed effects 
Model$summary.hyperpar #gives summary of hyper parameters;(rho = Beta for ID.area)
  

marg.stdev.IDarea <- inla.tmarginal(function(x) x^(-1/2), #
                             Model$marginals.hyperpar$"Precision for ID.area") #

marg.stdev.IDmonth<- inla.tmarginal(function(x) x^(-1/2),
                                    Model$marginals.hyperpar$`Precision for ID.month`, 
                                    n = 2030L,
                                    h.diff = .Machine[["double.eps"]]^(1/2)) 

marg.stdev.IDarea.month <- inla.tmarginal(function(x) x^(-1/2),
                                    Model$marginals.hyperpar$`Precision for ID.area.month`) 

### get variance of the random effects:
marg.variance.area <- inla.tmarginal(function(x) 1/x,
                                Model$marginals.hyperpar$"Precision for ID.area") #
inla.emarginal(function(x) x, marg.variance.area) #
#expected value of variance
m <- inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area")) #
mm <- inla.emarginal(function(x) x^2, marg.variance.area)
sqrt(mm - m^2) #std of variance hyper para  

inla.zmarginal(marg.variance.area) #

marg.variance.month <- inla.tmarginal(function(x) 1/x,
                                     Model$marginals.hyperpar$"Precision for ID.month") #
marg.variance.area.month <- inla.tmarginal(function(x) 1/x,
                                      Model$marginals.hyperpar$"Precision for ID.area.month") #


par = c("mu","sigma2_S","rho_S","sigma2_T","sigma2_ST")

mean.model <- c(summary(Model)$fixed[1],
                inla.emarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area"),
                Model$summary.hyperpar[2,1],  
                inla.emarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.month"),
                inla.emarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area.month"))

sd.model <- c(summary(Model)$fixed[2],
              sqrt(inla.emarginal(function(x) 1/(x^2), Model$marginals.hyperpar$"Precision for ID.area")-mean.model[2]^2),
              Model$summary.hyperpar[2,2],
              sqrt(inla.emarginal(function(x) 1/(x^2), Model$marginals.hyperpar$"Precision for ID.month")-mean.model[4]^2),
              sqrt(inla.emarginal(function(x) 1/(x^2), Model$marginals.hyperpar$"Precision for ID.area.month")-mean.model[5]^2))

q1.model <- c(summary(Model)$fixed[3],
              inla.qmarginal(0.025, inla.tmarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area")),
                Model$summary.hyperpar[2,3],  
              inla.qmarginal(0.025, inla.tmarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.month")),
              inla.qmarginal(0.025, inla.tmarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area.month")))

q3.model <- c(summary(Model)$fixed[5],
              inla.qmarginal(0.975, inla.tmarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area")),
              Model$summary.hyperpar[2,5],  
              inla.qmarginal(0.975, inla.tmarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.month")),
              inla.qmarginal(0.975, inla.tmarginal(function(x) 1/x, Model$marginals.hyperpar$"Precision for ID.area.month")))

df.post.best <- data.frame( par,3, mean.model, sd.model, q1.model, q3.model) #summary of posteriors
print(df.post.best)

round(df.post.best[,3:6],3)
write.csv(df.post.best, file= "./results/df.post.best.csv", row.names = F)


##################################################
###compute Fractions of variance explained by random effects
##################################################
T <- length(unique(Model$.args$data$ID.month))
mu <- Model$summary.lincomb.derived$mean[1]
risks.S <- matrix(Model$summary.lincomb.derived$mean[2:(S+1)],S,T) # spatial random effect
risks.T <- t(matrix(Model$summary.lincomb.derived$mean[(S+2):(S+T+1)],T,S)) # temporal random effect
risks.ST <- matrix(Model$summary.lincomb.derived$'mean'[(S+T+2):(S+T+S*T+1)], S, T) # space-time random effect

varS <-var(as.vector(risks.S))
varT <- var(as.vector(risks.T))
varST <- var(as.vector(risks.ST))
round(100*c(varS,varT,varST)/(varS+varT+varST),4) #Fraction of variance for each random effect

check <- as.vector(risks.S) #

#####################################
# Fig 4: plot observed vs predicted rate
######################################
plot(x= Model$summary.fitted.values$mean, y=data.req$p , ylab="Observed", xlab="Estimated")
abline(a=0, b=1)

data.req$est.p <- Model$summary.fitted.values$mean
ggplot(data.req, aes(x=p, y= est.p)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept=0, slope=1) +
  labs(x='Actual rate', y='Predicted rate')

#####################################
# Moran's I test for best model
######################################
#function to compute Moran's I 
nb <- poly2nb(map)
colW <- nb2listw(nb, style="B", zero.policy=NULL)
z <-data.req$p - data.req$est.p
 
df.resid <- data.frame(region = data.req$region, monthNumber = data.req$ID.month,
                         residual = z)
df.resid <-  df.resid %>%
    mutate(residualFilled = if_else(is.na(residual), 0, residual))   #fill nan values with 0
  #df_resid <- na.omit(df_resid) 
  
df.morans <- df.resid %>% #compute moran's I for each time point
    group_by(monthNumber) %>%
    summarise(mstat = round(moran.test(x = residualFilled,listw = colW,
                                       na.action = na.pass, zero.policy=TRUE)$estimate[[1]],4),
              pvalue = round(moran.test(x = residualFilled,listw = colW,
                                        na.action = na.pass, zero.policy=TRUE)$p.value,4)
    )%>%
    data.frame()

#get number of time periods with significant spatial correlation
df.signf <- df.morans %>%
    filter(pvalue < 0.05) # 


##################################
# fitted  incidence rates from: best model: Type2:RW1 Poisson
## to get forecast models load the estimated model in: /results/forecast.Models.all.Rdata"
#################################
names(forecast.Models) #"poiss2.RW1.forecast" "nb2.RW1.forecast"   

df.frcst<- data

S <- length(unique(df.frcst$region)) #22 regions
T.forecast <- length(unique(df.frcst$date))  # 78. months

df.frcst <- df.frcst %>%    #insert index cols
  mutate(ID.area=rep(1:S,T.forecast),
         ID.month=rep(1:T.forecast,each=S),
         ID.area.month=seq(1,T.forecast*S)
  )%>%  
  filter(ID.month <=78) ##filter month nors <= last forecasted timepoint

#get estimated incidence rates from forecasted models 
frcst.model <- forecast.Models$poiss2.RW1.forecast   #  
summary(frcst.model) 
length(frcst.model$summary.fitted.values[, "0.5quant"]) #  

df.frcst$est.p <- frcst.model$summary.fitted.values[, "0.5quant"] #estimated rate ratio
df.frcst$est.rate <- 10000*(df.frcst$est.p)
summary(df.frcst$est.rate)
summary(df.frcst$rate)

##################################
# Fig 7: spatial distribution of estimated incidence rates from: best model 
#################################
color.pal <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(0,5,9,15,20,30,40,50,Inf)      

tmap_mode("plot")
#yr <- 2017
plt.est <- function(yr){
    df.est <- df.frcst %>%
        filter(year == yr)
    print(unique(df.est$year))
    df.est.sub <-  reshape(df.est[, c('region', 'month', 'est.rate') ],  
                 timevar = "month",
                 idvar = "region",
                 direction = "wide"
      )
    colnames(df.est.sub) <- c("region", seq(1,12))   #use seq(1,12) for yrs 2017 - 2022 and seq(1,6) for 2023

    map.est <- merge(map, df.est.sub,  by.x = "SA3_NAME21", by.y = "region")
    ###use seq(1,12) for yrs 2017 - 2022 and seq(1,6) for 2023
    fig.est <- tm_shape(map.est) +
      tm_polygons(col=paste(col = seq(1,12,by=1)), palette=color.pal,  breaks=values,  
            legend.show=T, 
            title=   paste0("Estimated rate per 10,000 population"),  #, in yr
            legend.width = 0.5,  
              legend.reverse=F,  legend.is.portrait = FALSE, interval.closure="left",
            border.col = "black",border.alpha = 0.5) + #style="fixed",
    tm_layout(legend.outside=T, legend.outside.position="top", legend.frame=F, legend.outside.size=0.1, #bottom
            panel.labels= month.name, #seq(1,12,by=1)
            outer.margins=c(0.02, 0.08, 0.02, 0.02)) + #0.02,0.01,0.02,0.01))0.0,0.0,0.0,0.0
     tm_facets(nrow=4, ncol=3)
      # tm_facets(nrow=2, ncol=3) #use this for 2023
}

fig.est <- plt.est(yr = 2017)  #change the year for each panel
print(fig.est)
#print(fig.est)

################################
#### get regions with rate > 9 cases per 10,000 population
################################
df.yr <- df.frcst %>%
  filter(year == 2023) %>% #do this for each yr 2017 to 2022
  filter(round(est.rate,0) > 9)
print(unique(df.yr$year))

get.high.inf.regs <- function(i){
  #for(i in 1:12){
  df.reg.high <- df.yr %>%
    filter(month == i)
  print(df.reg.high$region)
} 

#x <- c(1:12) ## for 2017 - 2022
x <- c(1:6) ## for 2023
high.reg <- sapply(x, get.high.inf.regs )
print(sort(unlist(Reduce(intersect, high.reg)))) #common values in all months
unique(sort(unlist(high.reg))) # unique values in all list


##################################
## Maps with posterior exceedence probabilities of rate ratio   ##
# P(theta > c) = 1- P(theta <= c), c= 0.0009 (for rate ratio)
#################################
summary(data.req$p)
#mean = 0.0009 *10000 = 9 average rate per 10000 popln
summary(df.frcst$est.p)

exc.prob <- sapply(frcst.model$marginals.fitted.values,
       FUN = function(marg){1-inla.pmarginal(q = 0.0009, marginal = marg)}) #length 1716 q = 1.2, q = 0.001
df.frcst$exc.prob <- exc.prob
round(summary(df.frcst$exc.prob),5)
hist(exc.prob)

color.pal <- brewer.pal(5,"RdYlGn")[5:1]
values <- c(0, 0.2, 0.5,0.8,1)

plt.exc.prob <- function(yr){
 
  df.prob <- df.frcst %>%
    filter(year == yr)
  print(unique(df.prob$year))
  df.prob.sub <-  reshape(df.prob[, c('region', 'month', 'exc.prob')],
                         timevar = "month",
                         idvar = "region",
                         direction = "wide"
  )
  
  colnames(df.prob.sub) <- c("region", seq(1,6))   #seq(1,12) for all yrs , seq(1,2) for 2023
  map.prob <- merge(map, df.prob.sub,  by.x = "SA3_NAME21", by.y = "region")
 
  fig.prob <- tm_shape(map.prob) +
    tm_polygons(col=paste(col = seq(1,6,by=1)), palette=color.pal,  breaks=values,  legend.show=T, #col = seq(1,12,by=1))
                legend.width = 0.5, title =   paste0("Posterior exceedance probability in ", yr),
               # title = expression(paste0("Exceedence probability risk ratio in ",  yr)),
                legend.reverse=F,  legend.is.portrait = FALSE, interval.closure="left",
               border.col = "black",border.alpha = 0.5) + #style="fixed",
    tm_layout(legend.outside=T, legend.outside.position="top", legend.frame=F, legend.outside.size=0.1,
              panel.labels= month.name, 
              outer.margins=c(0.02, 0.08, 0.02, 0.02)) + #0.02,0.01,0.02,0.01)) 0.0, 0.0, 0.0, 0.0
     tm_facets(nrow=2, ncol=3)
     #tm_facets(nrow=4, ncol=3)
}
fig.exc.prob.2023 <- plt.exc.prob(yr = 2023)
print(fig.exc.prob.2023)


#####check which regions are >0.8
df.prob <- df.frcst %>%
  filter(year == 2023) %>%
  filter(round(exc.prob,5) > 0.8)

##################################
#plot temporal trends of region specific rate
#################################

df.frcst$q1 <- frcst.model$summary.fitted.values[,"0.025quant"]
df.frcst$q2 <- frcst.model$summary.fitted.values[,"0.5quant"]
df.frcst$q3 <- frcst.model$summary.fitted.values[,"0.975quant"]


min <- min(df.frcst$date) ### 
max <- max(df.frcst$date) 
date.breaks = filter(df.frcst, row_number() %% 22 == 1) #use date column in this df to set breaks 
date.breaks = filter(date.breaks, row_number() %% 12 == 1)

fig.est.trend <- ggplot(df.frcst, aes(x= date) ) +  
  geom_line(aes(y = p, color = "Observed rate"))+ 
  geom_ribbon(aes(ymin=q1, ymax=q3, fill='95% Credible Interval'), alpha=0.5) +
  geom_line(aes(y = q2, color = "Predicted rate")) + 
  scale_color_manual("", values = c("red","black"))+
  scale_fill_manual("",values="grey")+
  scale_x_date(breaks = date.breaks$date, date_labels = "%m-%Y", limits = c(min, max))+ 
  facet_wrap(. ~ region, nrow = 7, ncol = 5, scales = "free")+
  theme(legend.position="bottom", legend.title=element_blank(), strip.text = element_text(size = 8))+ #,
  theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed",
         linewidth  = 0.1),panel.background = element_blank() )+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  labs(x="Time", y = " ")

fig.est.trend
 
##### trend for regions with high rates with high population or low data 
lst_rgn_high <- c("Barkly",  "Cairns_Region",
             "Charters Towers - Ayr - Ingham", 
             "Darwin_Region", "East Arnhem", "East Pilbara",
             "Far North",  "Katherine",
             "Kimberley",  
             "Outback - North",  
            "Townsville", "West Pilbara")
df_high_reg <-  df.frcst %>%
  filter(region %in% lst_rgn_high) 
print(unique(df_high_reg$region))
  
fig.est.trend_high <- ggplot(df_high_reg, aes(x= date) ) + 
   geom_line(aes(y = p, color = "Observed rate"))+
  geom_ribbon(aes(ymin=q1, ymax=q3, fill='95% Credible Interval'), alpha=0.5) +
  geom_line(aes(y = q2, color = "Predicted rate")) + 
  scale_color_manual("", values = c("red","black"))+
  scale_fill_manual("",values="grey")+
  scale_x_date(breaks = date.breaks$date, date_labels = "%m-%Y", limits = c(min, max))+ 
  facet_wrap(. ~ region, nrow = 5, ncol =3, scales = "free")+
  theme(legend.position="bottom", legend.title=element_blank(), strip.text = element_text(size = 8))+ #,
  theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed",
                                        linewidth  = 0.1),panel.background = element_blank() )+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  labs(x="Time", y = " ")

fig.est.trend_high

#### trend for some regions with low rates or low population
lst_rgn_low <- c("Alice Springs", "Bowen Basin - North", 
                  "Daly - Tiwi - West Arnhem",
                 "Innisfail - Cassowary Coast", 
                  "Litchfield", "Mackay",
                  "Palmerston", "Port Douglas - Daintree",
                 "Tablelands (East) - Kuranda",    "Whitsunday")
df_low_reg <-  df.frcst %>%
  filter(region %in% lst_rgn_low) 
print(unique(df_low_reg$region))

fig.est.trend_low <- ggplot(df_low_reg, aes(x= date) ) + 
  geom_line(aes(y = p, color = "Actual rate"))+ #round(10000*p,0)
  geom_ribbon(aes(ymin=q1, ymax=q3, fill='95% Credible Interval'), alpha=0.5) +
  geom_line(aes(y = q2, color = "Predicted rate")) + 
  scale_color_manual("", values = c("red","black"))+
  scale_fill_manual("",values="grey")+
  scale_x_date(breaks = date.breaks$date, date_labels = "%m-%Y", limits = c(min, max))+ 
  facet_wrap(. ~ region, nrow = 5, ncol =3, scales = "free")+
  theme(legend.position="bottom", legend.title=element_blank(), strip.text = element_text(size = 8))+ #,
  theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed",
                                        linewidth  = 0.1),panel.background = element_blank() )+  
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  labs(x="Time", y = " ")

fig.est.trend_low
####################################

 