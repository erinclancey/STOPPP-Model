
###################################
library(ggplot2)

set=60
df <- even.sample.list[[set]]
kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 50, n.points=length(df$t) )
R_t_k <- diff(kspline$y)

# kspline_I <- ((R_t_k +
#                  kspline$y[-1]*
#                  (pars[even.sample.list[[set]]$pars_id[1],]$omega_a))/
#                 pars[even.sample.list[[set]]$pars_id[1],]$gamma)

# kspline_I <- ((R_t_k +
#                  kspline$y[-1]*
#                  (pars[even.sample.list[[set]]$pars_id[1],]$omega_a+pars[even.sample.list[[set]]$pars_id[1],]$b_ave))/
#                 pars[even.sample.list[[set]]$pars_id[1],]$gamma)
# 
kspline_I <- ((R_t_k +
                 kspline$y[-1] *
                 (pars[even.sample.list[[set]]$pars_id[1],]$omega_a+with(pars[set,], b_func(t=seq(-27,364,1), g=g, s=s, psi=psi, f=f)))  )/
                pars[even.sample.list[[set]]$pars_id[1],]$gamma)

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

quantile <- subset(spline_df_plot, kspline_I>quantile(kspline_I, probs = 0.90))
min(quantile$Time)
max(quantile$Time)

##########################################
A <- ggplot() +
  geom_point(data=df_plot, aes(x = Time, y = Seropositive, color = "Daily Data"), size = 2, alpha=0.3, shape=16)+ 
  geom_line(data = data.frame(Time = data.tau.list[[1]]$t[-1], 
                              Value = data.tau.list[[set]]$R_a[-1] / data.tau.list[[set]]$N[-1]),
            aes(x = Time, y = Value, color = "True"), size = 1.2, alpha=1, linetype=1) +
  geom_line(data = data.frame(Time = df$t, Value = kspline$y), 
            aes(x = Time, y = Value, color = "Predicted"), size = 1.5) +

  labs(title = "Seroprevalence", x = "", y = "Prop. Seropos.") +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0)) +
  scale_color_manual(values = c("Daily Data"="#0072B2","True" = "black", "Predicted" = "#0072B2"),
                     labels = c(TeX("Data"),TeX("Predicted $(\\hat{r}_a)$"), TeX("True $(r_a)$"))) +
  guides(color = guide_legend(title = ""))


######################################
B <- ggplot() +
  geom_line(data = data.frame(Time = data.tau.list[[set]]$t[-1], 
                              Value = data.tau.list[[set]]$I[-1] / data.tau.list[[set]]$N[-1]),
            aes(x = Time, y = Value, color = "True"), size = 1.2, alpha=1) +
  geom_line(data = spline_df_plot, aes(x = Time, y = Prop_Infected, color = "Predicted"), size = 1.5) +

  labs(title = "Prevalence", x = "Time (Days)", y = "Prop. Infected") +
  theme_minimal() +
  theme(text = element_text(size = 18))+
  theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_color_manual(values = c("Predicted" = "#CC79A7", "True" = "black"),
                     labels = c(TeX("Predicted $(\\hat{\\iota})$"), TeX("True $(\\iota)$"))) +
  guides(color = guide_legend(title = ""))+
  geom_ribbon(data = subset(spline_df_plot, Time >= min(quantile$Time) & Time <= max(quantile$Time)), 
              aes(x = Time, ymin = 0, ymax = Prop_Infected), 
              fill = "#CC79A7", alpha = 0.5)+
  geom_vline(xintercept = pars[set,]$peak_I, linetype=3, size=1)

plot <- plot_grid(A,B, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot

