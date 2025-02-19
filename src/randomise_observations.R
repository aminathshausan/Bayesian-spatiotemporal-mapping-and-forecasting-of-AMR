## Date: 19 Feb 2025
# Author: Dr Aminath Shausan 

####################################
#this program adds a random noise to observed AMR counts 
## this is done as due to ethical considerations, actual data can't be shared publicly 
###############################
.libPaths("~/Documents/R/r-libraries")
#install.packages("scales", lib="~/Documents/R/r-libraries")
options(digits=8)
rm(list = ls())
set.seed(963258)

#library(MASS)
#library(INLA)
#library(inlatools)
#library(sf)
#library(spdep) #to compute neighbours of regions
#library("sf")
library(dplyr)
library(plotly)
library(tmap)
library(RColorBrewer)
library(ggplot2)


########################################################
# add noise to raw data 
######################################
data<- read.csv('./data/data.csv') 
str(data) #date is a character
data$date <- as.Date(data$date, "%d/%m/%Y") #"%Y-%m-%d"

data <- data[order(data$date,data$region),] 
rownames(data)<-1:nrow(data)

#add noise to observed data
n <- length(data$obs)
data <- data %>%
        mutate(obs.noise = round(obs + runif(n, -1, -1) ))

#replace negative observation by 0
data$obs.noise <- ifelse(data$obs.noise < 0, 0, data$obs.noise)
 
#select rows with negatobe values repla
data[which(data$obs.noise < 0,arr.ind = TRUE),]

write.csv(data, file= "./data/data.csv")    


########################################################
# add noise to fitted values 
######################################
