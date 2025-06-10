n_trials=20
############Daily Samples
sample_daily <- list()
for(i in 1:length(data.tau.list) ){
  t <- data.tau.list[[i]]$t
  I_sample <- vector()
  R_a_sample <- vector()
  for(k in 1:length(t)) {
    I_sample[k] <- with(data.tau.list[[i]], rbinom(1, n_trials, I[k]/N[k]))
    R_a_sample[k] <- with(data.tau.list[[i]], rbinom(1, n_trials, R_a[k]/N[k]))
  }
  df <- data.frame(t, I_sample, R_a_sample)
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("daily",length(t))
  sample_daily[[i]] <- df
}
##########Weekly Samples
sample_weekly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  week <- seq(1,length(df$t),7)
  df <- df[week,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("weekly",length(df$t))
  sample_weekly[[i]] <- df

}
##########Bi-Weekly Samples
sample_biweekly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  bi_week <- seq(1,length(df$t),14)
  df <- df[bi_week,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("biweekly",length(df$t))
  sample_biweekly[[i]] <- df
}
new_sample <- c(sample_daily,sample_weekly,sample_biweekly)    

omega_a=pars$omega_a
gamma=pars$gamma

set=2
df <- new_sample[[set]]
kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 45, n.points=500)
R_t_k <- diff(kspline$y)
#kspline_I <- (R_t_k + kspline$y[-1]*omega_a)/gamma

b_rates = with(pars, b_func(t=kspline$x, g=g, s=s, psi=psi, f=f)) 
kspline_I <- (R_t_k +
                 kspline$y[-1] *
                 (omega_a+b_rates[-1]))/gamma


t_spline <- kspline$x[-1]
spline_df <- data.frame(t_spline, kspline_I)

# par(mfrow = c(1,2), mar = c(5,5,1,1))
# plot(df$t,df$R_a_sample/n_trials,xaxt = "n",yaxt = "n",ylim=c(0,0.7),
#      xlab = "Time (Days)",
#      ylab = "Prop. Seropos.",
#      pch = 20,
#      col= "grey", cex.lab=1.5, cex.axis=1.5)
# axis(1, at = seq(0, max(df$t), by = 365))
# axis(2, at = seq(0, 0.8, by = 0.1))
# #title("(a) Low Amplitude Daily Sample", adj=0, cex.main=1.75)
# lines(df$t,kspline$y, lwd=2,col="#0072B2")
# lines(data.tau.list[[1]]$t[-1], data.tau.list[[1]]$R_a[-1]/data.tau.list[[1]]$N[-1],lwd=2,col="black")
# legend("topright", c(expression(r[a]),expression(hat(r)[a])),
#        col = c("black","#0072B2"), lty=1, cex=1.5)
# 
# plot(t_spline,kspline_I, ylim = c(0, 0.3),xaxt = "n",yaxt = "n",
#      xlab = "Time (Days)",
#      ylab = "Prop. Infected",
#      type="l",lwd=2,
#      col= "#CC79A7",cex.lab=1.5, cex.axis=1.5)
# axis(1, at = seq(0, max(df$t), by = 365))
# axis(2, at = seq(0, 0.3, by = 0.05))
# lines(data.tau.list[[1]]$t[-1], data.tau.list[[1]]$I[-1]/data.tau.list[[1]]$N[-1], col="black",lwd=2)
# legend("topright", c(expression(iota),expression(hat(iota))),
#        col = c("black","#CC79A7"), lty=1:1, cex=1.5)


library(ggplot2)

# Convert data to a tibble for ggplot compatibility
df$t <- as.numeric(df$t)  # Ensure time is numeric

# Create kspline plot data
kspline_df <- data.frame(t = kspline$x, kspline_y = kspline$y)
data_tau_df <- data.frame(
  t = data.tau.list[[1]]$t[-1],
  R_a = data.tau.list[[1]]$R_a[-1] / data.tau.list[[1]]$N[-1]
)

# First plot: Proportion seropositive

#####################

  p1 <- ggplot() +
    geom_point(data = df, aes(x = t/7, y = R_a_sample/n_trials, color = "Data"), alpha = 0.25, shape = 16) +
    geom_line(data = kspline_df, aes(x = t/7, y = kspline_y, color = "Predicted"), size = 1.5) +
    geom_line(data = data_tau_df, aes(x = t/7, y = R_a, color = "True"), size = 1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, 156, by = 52), labels = function(x) x) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
    labs(x = " ", y = "Prop. Seropos.") +
    theme_classic(base_size = 15) +
    theme(axis.text = element_text(size = 14)) +
    theme(legend.position = "right") +
    scale_color_manual(values = c("Data" = "#0072B2", "Predicted" = "#0072B2", "True" = "black"),
                       labels = c(TeX("Data"), TeX("Predicted $(\\hat{r}_a)$"), TeX("True $(r_a)$"))) +
    guides(color = guide_legend(title = ""))
  ########################################



# Second plot: Proportion infected
spline_df <- data.frame(t = t_spline, kspline_I = kspline_I)
data_tau_I_df <- data.frame(
  t = data.tau.list[[1]]$t[-1],
  I = data.tau.list[[1]]$I[-1] / data.tau.list[[1]]$N[-1]
)


library(pracma)  # Ensure the package is loaded

# Find peaks first
peaks <- findpeaks(kspline_I, minpeakheight = 0.1, minpeakdistance = 20, npeaks = 6, nups = 0, ndowns = 0)

# Create a subset based on detected peak times
peaks_df <- spline_df[spline_df$t %in% spline_df$t[peaks[,2]], ]

#####################
# Convert time variables to weeks
spline_df$t_weeks <- spline_df$t / 7
data_tau_I_df$t_weeks <- data_tau_I_df$t / 7

# Ensure peaks_df references the correct transformed time column
peaks_df$t_weeks <- peaks_df$t / 7

p2 <- ggplot() +
  geom_line(data = spline_df, aes(x = t_weeks, y = kspline_I, color = "Predicted"), size = 1.5) +
  geom_line(data = data_tau_I_df, aes(x = t_weeks, y = I, color = "True"), size = 1, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 156, by = 52), labels = function(x) x) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(-0.085, 0.5)) +
  labs(x = "Time (Weeks)", y = "Prop. Infected") +
  theme_classic(base_size = 15) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.position = "right") +
  scale_color_manual(values = c("Predicted" = "#CC79A7", "True" = "black"),
                     labels = c(TeX("Predicted $(\\hat{\\iota})$"), TeX("True $(\\iota)$"))) +
  guides(color = guide_legend(title = "")) +
  geom_vline(xintercept = peaks_df$t_weeks[1], linetype = 3, color = "#CC79A7", size = 1) +
  geom_vline(xintercept = peaks_df$t_weeks[2], linetype = 3, color = "#CC79A7", size = 1) +
  geom_vline(xintercept = peaks_df$t_weeks[3], linetype = 3, color = "#CC79A7", size = 1) +
  geom_vline(xintercept = peaks_df$t_weeks[4], linetype = 3, color = "#CC79A7", size = 1) +
  geom_vline(xintercept = peaks_df$t_weeks[5], linetype = 3, color = "#CC79A7", size = 1)
########################################

# Display side-by-side
library(patchwork)
p1 / p2

