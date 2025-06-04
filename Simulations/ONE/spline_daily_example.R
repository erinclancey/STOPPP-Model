# set=1
# df <- even.sample.list[[set]]
# kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 50, n.points=length(df$t) )
# R_t_k <- diff(kspline$y)
# kspline_I <- ((R_t_k + kspline$y[-1]*pars[even.sample.list[[set]]$pars_id[1],]$omega_a)/pars[even.sample.list[[set]]$pars_id[1],]$gamma)
# t_spline <- df$t[-1]
# spline_df <- data.frame(t_spline, kspline_I)
# par(mfrow = c(3,2), mar = c(5,5,5,1))
# plot(df$t,df$R_a_sample/n_trials, ylim = c(0, 0.7), 
#      xlab = "Time (Days)",
#      ylab = "Prop. Seropos.",
#      pch = 20,
#      col= "grey", cex.lab=1.75, cex.axis=1.75)
# title("(a) Low Amplitude Daily Sample", adj=0, cex.main=1.75)
# lines(df$t,kspline$y, lwd=2,col="black")
# lines(data.tau.list[[1]]$t[-1], data.tau.list[[set]]$R_a[-1]/data.tau.list[[set]]$N[-1],lwd=2,col="#0072B2")
# legend("topright", c(expression(r[a]),expression(hat(r)[a])),
#        col = c("#0072B2","black"), lty=1, cex=1.5)
# 
# plot(t_spline,kspline_I, ylim = c(0, 0.3), 
#      xlab = "Time (Days)",
#      ylab = "Prop. Infected",
#      type="l",lwd=2,
#      col= "black",cex.lab=1.75, cex.axis=1.75)
# lines(data.tau.list[[set]]$t[-1], data.tau.list[[set]]$I[-1]/data.tau.list[[set]]$N[-1], col="#E69F00",lwd=2)
# legend("topright", c(expression(iota),expression(hat(iota))),
#        col = c("#E69F00","black"), lty=1:1, cex=1.5)

###################################
library(ggplot2)

set=60
df <- even.sample.list[[set]]
kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 50, n.points=length(df$t) )
R_t_k <- diff(kspline$y)

kspline_I <- ((R_t_k +
                 kspline$y[-1]*
                 (pars[even.sample.list[[set]]$pars_id[1],]$omega_a))/
                pars[even.sample.list[[set]]$pars_id[1],]$gamma)

# kspline_I <- ((R_t_k +
#                  kspline$y[-1]*
#                  (pars[even.sample.list[[set]]$pars_id[1],]$omega_a+pars[even.sample.list[[set]]$pars_id[1],]$b_ave))/
#                 pars[even.sample.list[[set]]$pars_id[1],]$gamma)
# 
# kspline_I <- ((R_t_k +
#                  kspline$y[-1] *
#                  (pars[even.sample.list[[set]]$pars_id[1],]$omega_a+with(pars[set,], b_func(t=seq(-27,364,1), g=g, s=s, psi=psi, f=f)))  )/
#                 pars[even.sample.list[[set]]$pars_id[1],]$gamma)

t_spline <- df$t[-1]

# Create a new data frame for plotting
df_plot <- data.frame(
  Time = df$t,
  Seropositive = df$R_a_sample / n_trials
)

spline_df_plot <- data.frame(
  Time = t_spline,
  Prop_Infected = kspline_I
)

##########################################
A <- ggplot() +
  geom_point(data=df_plot, aes(x = Time, y = Seropositive, color = "Daily Data"), size = 2, alpha=0.3, shape=16)+ 
  geom_line(data = data.frame(Time = df$t, Value = kspline$y), 
            aes(x = Time, y = Value, color = "Predicted"), size = 1.5) +
  geom_line(data = data.frame(Time = data.tau.list[[1]]$t[-1], 
                              Value = data.tau.list[[set]]$R_a[-1] / data.tau.list[[set]]$N[-1]),
            aes(x = Time, y = Value, color = "True"), size = 1.2, alpha=1, linetype=1) +
  labs(title = "Seroprevalence", x = "", y = "Prop. Seropos.") +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0)) +
  scale_color_manual(values = c("Daily Data"="#0072B2", "Predicted" = "#0072B2", "True" = "black"),
                     labels = c(TeX("Data"),TeX("Predicted $(\\hat{r}_a)$"), TeX("True $(r_a)$"))) +
  guides(color = guide_legend(title = ""))


######################################
B <- ggplot() +
  geom_line(data = spline_df_plot, aes(x = Time, y = Prop_Infected, color = "Predicted"), size = 1.5) +
  geom_line(data = data.frame(Time = data.tau.list[[set]]$t[-1], 
                              Value = data.tau.list[[set]]$I[-1] / data.tau.list[[set]]$N[-1]),
            aes(x = Time, y = Value, color = "True"), size = 1.2, alpha=1) +
  labs(title = "Prevalence", x = "Time (Days)", y = "Prop. Infected") +
  theme_minimal() +
  theme(text = element_text(size = 18))+
  theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_color_manual(values = c("Predicted" = "#CC79A7", "True" = "black"),
                     labels = c(TeX("Predicted $(\\hat{\\iota})$"), TeX("True $(\\iota)$"))) +
  guides(color = guide_legend(title = ""))

plot <- plot_grid(A,B, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot

