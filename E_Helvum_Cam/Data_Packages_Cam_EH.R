#Packages
library("lubridate") 
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(latex2exp)
library(zoo)
library(minpack.lm)
library(rstan)
library(reshape2)
library(coda)
library(MCMCglmm)
library(bayestestR)
library(pracma)
library(reshape)
library(functional)
library(matrixcalc)
library(lhs)
set.seed(420)

### DATA
Cam_data <- read.csv("E_helvum_data.csv", header=TRUE)
Date <- as.Date(Cam_data$Dat_EC)
J_Date <- yday(as.Date(Cam_data$Dat_EC))
Week <- as.numeric(strftime(Date, format = "%V"))
Cam <- cbind(Cam_data, J_Date, Week)
#Cam <- subset(Cam, Cam$Age==c("A"))

##EBOV Seropositives by week
EBVagg_cases <- Cam %>% group_by(Week) %>% dplyr::summarize(EBVPos = sum(Res1GP.ZEBVkiss)) %>% as.data.frame()
Juv <- Cam %>% group_by(Week) %>% dplyr::summarize(Juv = sum(Juvinile)) %>% as.data.frame()
week_sample_count <- data.frame(table(Cam$Week))
EBVagg_cases <- cbind(EBVagg_cases, Juv$Juv, week_sample_count$Freq)
colnames(EBVagg_cases) <- c("Week", "EBVPos","Juv", "Total_Sampled")
EBVagg_cases$EBVPos_prop <- EBVagg_cases$EBVPos / EBVagg_cases$Total_Sampled
EBVagg_cases$Juv_prop <- EBVagg_cases$Juv/ EBVagg_cases$Total_Sampled

#0,1,2,6,7,8,11,12,15,16,17,20,21,22,23,26,27,30,31,32,33,34,35,36,39,51 missing
#missing 7 weeks starting week 30
#25 weeks total sampled

# plot(EBVagg_cases$Week, EBVagg_cases$EBVPos_prop)
# plot(EBVagg_cases$Week, EBVagg_cases$Juv_prop)

###Lactating females Week
lac <- data.frame(Cam$J_Date, Cam$Week, Cam$Lactation)
lac <- na.omit(lac)
colnames(lac) <- c("Day", "Week","Lac")
agg_ges <- lac %>% group_by(Week) %>% dplyr::summarize(ges_num = sum( Lac)) %>% as.data.frame()
week_sample_count <- data.frame(table(lac$Week))
agg_ges <- cbind(agg_ges, week_sample_count$Freq)
colnames(agg_ges) <- c("week", "ges", "Total_Sampled")
agg_ges$ges_Prop <- agg_ges$ges / agg_ges$Total_Sampled
#Freq is the total number of bats caught in each month, verified in the excel sheet
#plot(agg_ges$week,agg_ges$ges_Prop)

