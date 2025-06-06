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
setwd("~/STOPPP Model/Mystomes_2")
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
View(Myst)


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

# Perform kernel smoothing
kspline <- ksmooth(MORV_cases$month_date, MORV_cases$MORVPos_prop, kernel = "normal", bandwidth = 75, n.points = 5000)

# Convert the output into a data frame
kspline_df <- data.frame(
  month_date = as.Date(kspline$x, origin = "1970-01-01"),
  MORVPos_prop = kspline$y
)

# Plot with the smoothing line
ggplot(df, aes(x = as.Date(month_date), y = MORVPos_prop)) +
  geom_point(color = "darkorchid4") +
  geom_line(data = kspline_df, aes(x = month_date, y = MORVPos_prop), color = "blue", linewidth = 1) +  
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Date", y = "Prop. Seropositive")

library(ggplot2)
library(scales)
library(latex2exp)  # For LaTeX-style annotations

ggplot(df, aes(x = as.Date(month_date), y = MORVPos_prop)) +
  geom_point(aes(color = "Daily Data"), size = 2, alpha = 0.5, shape = 16) +  # Map color inside aes()
  geom_line(data = kspline_df, aes(x = month_date, y = MORVPos_prop, color = "Predicted"), size = 1) +  
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = expression(paste(italic("M. natalensis "), "Seroprevalence")), x = "", y = "Prop. Seropositive") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.9, 0.9),
    plot.title = element_text(hjust = 0)
  ) +
  scale_color_manual(
    values = c("Daily Data" = "#0072B2", "Predicted" = "#0072B2"),
    labels = c(TeX("Data"), TeX("Predicted $(\\hat{r}_a)$"))
  ) +
  guides(color = guide_legend(title = ""))

