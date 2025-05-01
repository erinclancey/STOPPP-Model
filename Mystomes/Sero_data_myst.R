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

Myst <- read.csv("MN serology data CMR Morogoro_EC.csv", header=TRUE)
Myst <- na.omit(Myst)

Myst$Antibody_status[Myst$Antibody_status == "neg_uncertain"] <- "neg"
Myst$Antibody_status[Myst$Antibody_status == "pos_uncertain"] <- "pos"
Myst$serpos[Myst$Antibody_status == "neg"] <- 0
Myst$serpos[Myst$Antibody_status == "pos"] <- 1


### DATA
Myst$Date_new <- format(as.Date(Myst$Date, format = "%d/%m/%Y"), "%Y-%m-%d")
Myst$J_date <- yday(as.Date(Myst$Date_new))
Myst$Week <- as.numeric(strftime(Myst$Date_new, format = "%V"))
Myst$Week_year <- strftime(Myst$Date_new, format = "%V-%Y")
Myst$Month <- as.numeric( format(as.Date(Myst$Date, format = "%d/%m/%Y"), "%m"))
Date <- as.numeric(format(as.Date(Myst$Date, format = "%d/%m/%Y"), "%Y%m%d"))
Myst$month_year <- strftime(Myst$Date_new, format = "%m-%Y")
Myst$Year <- strftime(Myst$Date_new, format = "%Y")



##EBV Seropositives by week
MORV_cases <- Myst %>% group_by(Week, Year) %>% summarize(MORVPos = sum(serpos)) %>% as.data.frame()
week_sample_count <- data.frame(table(Myst$Week, Myst$Year))
colnames(week_sample_count) <- c("Week", "Year", "Total_Sampled")
MORV_cases <- merge(MORV_cases, week_sample_count, by=c("Week","Year"))
MORV_cases$MORVPos_prop <- MORV_cases$MORVPos / MORV_cases$Total_Sampled

ggplot(MORV_cases, aes(x = Week, y = MORVPos_prop)) +
  geom_point(color = "darkorchid4") +
  facet_wrap( ~ Year ) + ylab("Proportion Seropositive")

##EBV Seropositives by week
MORV_cases <- Myst %>% group_by(Week) %>% summarize(MORVPos = sum(serpos)) %>% as.data.frame()
week_sample_count <- data.frame(table(Myst$Week))
colnames(week_sample_count) <- c("Week", "Total_Sampled")
MORV_cases <- merge(MORV_cases, week_sample_count, by=c("Week"))
MORV_cases$MORVPos_prop <- MORV_cases$MORVPos / MORV_cases$Total_Sampled

ggplot(MORV_cases, aes(x = Week, y = MORVPos_prop)) +
  geom_point(color = "darkorchid4")

library(npreg)
MORV_cases <- na.omit(MORV_cases)
t=MORV_cases$Week
x=MORV_cases$MORVPos
n=MORV_cases$Total_Sampled

seroprev <- x/n

m_spline_R <- ss(t,x/n, method = "REML", nknots=25)
p_spline_R <- predict(m_spline_R, t)$y
plot(t,x/n)
lines(t, p_spline_R, col = "blue", lwd = 2)
boot <- boot(m_spline_R, boot.dist=TRUE, R=1000)
plot(boot)
t <- seq(0,52, 0.26)
dist <- cbind(t,data.frame(t(boot$boot.dist)))
dist %>%
  ggplot(aes(x=t,y=X1))+geom_line()

long.dist <- melt(dist, id = c("t"))
long.dist %>%
  ggplot(aes(x=t,y=value,group=variable, color=ID))+
  geom_line(aes(group=variable))+
  scale_colour_gradient(low = "black", high = "#E69F00")+
  guides(color="none") + 
  theme_minimal()+theme_minimal()+ylab("Proportion Infected")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))
