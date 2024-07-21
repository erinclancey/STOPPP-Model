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
Eid_data <- read.csv("E_helvum_data.csv", header=TRUE)
Date <- as.Date(Eid_data$Dat_EC)
J_Date <- yday(as.Date(Eid_data$Dat_EC))
Week <- as.numeric(strftime(Date, format = "%V"))
Eid <- cbind(Eid_data, J_Date, Week)

##EBOV Seropositives by week
EBOVagg_cases <- Eid %>% group_by(Week) %>% dplyr::summarize(EBOVPos = sum(Res1GP.ZEBVkiss)) %>% as.data.frame()
Juv <- Eid %>% group_by(Week) %>% dplyr::summarize(Juv = sum(Juvinile)) %>% as.data.frame()
week_sample_count <- data.frame(table(Eid$Week))
EBOVagg_cases <- cbind(EBOVagg_cases, Juv$Juv, week_sample_count$Freq)
colnames(EBOVagg_cases) <- c("Week", "EBOVPos","Juv", "Total_Sampled")
EBOVagg_cases$EBOVPos_prop <- EBOVagg_cases$EBOVPos / EBOVagg_cases$Total_Sampled
EBOVagg_cases$Juv_prop <- EBOVagg_cases$Juv/ EBOVagg_cases$Total_Sampled


###Lactating females Week
lac <- data.frame(Eid$J_Date, Eid$Week, Eid$Lactation)
lac <- na.omit(lac)
colnames(lac) <- c("Day", "Week","Lac")
agg_lac <- lac %>% group_by(Week) %>% dplyr::summarize(lac_num = sum(Lac)) %>% as.data.frame()
week_sample_count <- data.frame(table(lac$Week))
agg_lac <- cbind(agg_lac, week_sample_count$Freq)
colnames(agg_lac) <- c("week", "lac", "Total_Sampled")
agg_lac$lac_Prop <- agg_lac$lac / agg_lac$Total_Sampled


