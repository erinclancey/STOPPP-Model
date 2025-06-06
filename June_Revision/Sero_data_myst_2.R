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
library(dplyr)
setwd("~/STOPPP Model/STOPPP_Model/June_Revision")
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
Myst$month_date <- as.Date(paste0("15-", Myst$month_year), format = "%d-%m-%Y")



MORV_cases <- Myst %>% group_by(month_date) %>% summarize(MORVPos = sum(serpos)) %>% as.data.frame()
week_sample_count <- data.frame(table(Myst$month_date))
colnames(week_sample_count) <- c("month_date", "Total_Sampled")
MORV_cases <- merge(MORV_cases, week_sample_count, by=c("month_date"))
MORV_cases$MORVPos_prop <- MORV_cases$MORVPos / MORV_cases$Total_Sampled

ggplot(MORV_cases, aes(x = as.Date(month_date), y = MORVPos_prop)) +
  geom_point(color = "darkorchid4") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  geom_line()
  


df <- MORV_cases
df <- na.omit(df)
gamma=1/0.7
omega_a=0
#################################
day <- seq(1,365,1)
b_rate <- vector()
b_rate[1:110]=0
b_rate[111:255]=0.0012 * (day[111:255]-110)
b_rate[256:329]=0.0012 * (day[256:329]-255)
b_rate[330:365]=0
b_df <- data.frame(day, b_rate)

# Perform kernel smoothing
kspline <- ksmooth(MORV_cases$month_date, MORV_cases$MORVPos_prop, kernel = "normal", bandwidth = 50, n.points = 300)
R_t_k <- diff(kspline$y)
day <- yday(as.Date(kspline$x[-1]))

x_spline <- kspline$x[-1]
y_spline <- kspline$y[-1]
new_df <- data.frame(x_spline, y_spline,R_t_k, day)

kspline_df <- merge(b_df, new_df, by = "day")
kspline_df <- arrange(kspline_df,x_spline)


kspline_df$kspline_I <- (R_t_k + y_spline *(omega_a + kspline_df$b_rate ))/ gamma


# Plot with the smoothing line
ggplot(df, aes(x = as.Date(month_date), y = MORVPos_prop)) +
  geom_point(color = "darkorchid4") +
  geom_line(data = kspline_df, aes(x = x_spline, y = y_spline), color = "blue", linewidth = 1) +  
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Date", y = "Prop. Seropositive")

library(ggplot2)
library(scales)
library(latex2exp)  # For LaTeX-style annotations

A <- ggplot(df, aes(x = as.Date(month_date), y = MORVPos_prop)) +
  geom_point(aes(color = "Daily Data"), size = 2, alpha = 0.5, shape = 16) +  # Map color inside aes()
  geom_line(data = kspline_df, aes(x = x_spline, y = y_spline, color = "Predicted"), size = 1) +  
  #geom_line(data = kspline_df, aes(x = x_spline, y = kspline_I, color = "#CC79A7"), size = 1) +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(title = expression(paste(italic("M. natalensis "), "Seroprevalence")), x = "", y = "Prop. Seropositive") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.85,1),
    plot.title = element_text(hjust = 0)
  ) +
  scale_color_manual(
    values = c("Daily Data" = "#0072B2", "Predicted" = "#0072B2"),
    labels = c(TeX("Data"), TeX("Predicted $(\\hat{r}_a)$"))
  ) +
  guides(color = guide_legend(title = ""))
A
B <- ggplot(df, aes(x = as.Date(month_date), y = MORVPos_prop)) +
  geom_line(data = kspline_df, aes(x = x_spline, y = kspline_I, color = "Predicted"), size = 1) +  
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  labs(title = expression(paste(italic("M. natalensis "), "Predicted Prevalence")), x = "", y = "Prop. Infected") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.85,1),
    plot.title = element_text(hjust = 0)
  ) +
  scale_color_manual(
    values = c("Predicted" = "#CC79A7"),
    labels = c(TeX("Predicted $(\\hat{\\iota})$"))
  ) +
  guides(color = guide_legend(title = ""))

plot <- plot_grid(A,B,ncol = 1, rel_widths = c(1,1))
plot




