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
library(xtable) 
set.seed(420)

### DATA
Hyp_data <- read.csv("H_monstrosus_data.csv", header=TRUE)
Date <- format(as.Date(Hyp_data$Date_EC,format = "%m/%d/%Y"), "%Y-%m-%d")
J_Date <- yday(as.Date(Date))
Week <- as.numeric(strftime(Date, format = "%V"))
Hyp <- cbind(Hyp_data, J_Date, Week)


### EBOV Seropositives by week
EBVagg_cases <- Hyp %>% group_by(Week) %>% dplyr::summarize(EBVPos = sum(Res1GP.ZEBVkiss)) %>% as.data.frame()
week_sample_count <- data.frame(table(Hyp$Week))
EBVagg_cases <- cbind(EBVagg_cases, week_sample_count$Freq)
colnames(EBVagg_cases) <- c("Week", "EBVPos", "Total_Sampled")
EBVagg_cases$EBVPos_prop <- EBVagg_cases$EBVPos / EBVagg_cases$Total_Sampled


###Lactating females Week
Hyp_noNA <- na.omit(Hyp)
agg_preg <- Hyp_noNA %>% group_by(Week) %>% dplyr::summarize(Pregnant_num = sum(Lactation)) %>% as.data.frame()
week_sample_count <- data.frame(table(Hyp_noNA$Week))
agg_preg <- cbind(agg_preg, week_sample_count$Freq)
colnames(agg_preg) <- c("Week", "lac", "Total_Sampled")
agg_preg$lac_Prop <- agg_preg$lac / agg_preg$Total_Sampled



