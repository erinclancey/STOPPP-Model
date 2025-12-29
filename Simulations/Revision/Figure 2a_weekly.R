# R code to produce Figure 2a
# Relies on Sim_epi_model.R and Sample_epi_model.R
###################################
set=465
df <- sample.list[[set]]
kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 70, n.points=500 )
R_t_k <- diff(kspline$y)
t <- kspline$x
t_spline <- kspline$x[-1]
kspline_I <- ((R_t_k +
                 kspline$y[-1] *
                 (pars[sample.list[[i]]$pars_id[1],]$omega_a+with(pars[sample.list[[i]]$pars_id[1],], b_func(t=t_spline, g=g, s=s, psi=psi, f=f)) )  )/
                pars[sample.list[[i]]$pars_id[1],]$gamma)
# Create a new data frame for plotting
df_plot <- data.frame(
  Time = df$t,
  Seropositive = df$R_a_sample / n_trials
)

spline_df_plot <- data.frame(
  Time =t_spline,
  Prop_Infected = kspline_I
)

quantile <- subset(spline_df_plot, kspline_I>quantile(kspline_I, probs = 0.90))
lower <- min(quantile$Time)/7
upper <- max(quantile$Time)/7

##########################################
A <- ggplot() +
  geom_point(data=df_plot, aes(x = Time/7, y = Seropositive, color = "Weekly Data"), size = 2, alpha=0.3, shape=16)+ 
  geom_line(data = data.frame(Time = data.tau.list.1_365[[set-300]]$t[-1]/7, 
                              Value = data.tau.list.1_365[[set-300]]$R_a[-1] / data.tau.list.1_365[[set-300]]$N[-1]),
            aes(x = Time, y = Value, color = "True"), size = 1.2, alpha=1, linetype=1) +
  geom_line(data = data.frame(Time = t/7, Value = kspline$y), 
            aes(x = Time, y = Value, color = "Predicted"), size = 1.5) +

  labs(title = "Seroprevalence", x = "", y = "Prop. Seropos.") +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(
    breaks = c(1, seq(6, 52*2, by = 6)) )+
  theme(plot.title = element_text(hjust = 0)) +
  scale_color_manual(values = c("Weekly Data"="#0072B2","True" = "black", "Predicted" = "#0072B2"),
                     labels = c(TeX("Data"),TeX("Predicted $(\\hat{r}_a)$"), TeX("True $(r_a)$"))) +
  guides(color = guide_legend(title = ""))


######################################
B <- ggplot() +
  geom_line(data = data.frame(Time = data.tau.list.1_365[[set-300]]$t[-1]/7, 
                              Value = data.tau.list.1_365[[set-300]]$I[-1] / data.tau.list.1_365[[set-300]]$N[-1]),
            aes(x = Time, y = Value, color = "True"), size = 1.2, alpha=1) +
  geom_line(data = spline_df_plot, aes(x = Time/7, y = Prop_Infected, color = "Predicted"), size = 1.5) +

  labs(title = "Prevalence", x = "Time (Weeks)", y = "Prop. Infected") +
  theme_minimal() +
  theme(text = element_text(size = 18))+
  theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(-0.02, 0.3)) +
  scale_x_continuous(
    breaks = c(1, seq(6, 52*2, by = 6)) )+
  scale_color_manual(values = c("Predicted" = "#CC79A7", "True" = "black"),
                     labels = c(TeX("Predicted $(\\hat{\\iota})$"), TeX("True $(\\iota)$"))) +
  guides(color = guide_legend(title = ""))+
  geom_ribbon(data = subset(spline_df_plot, Time/7 >= lower & Time/7 <= upper), 
              aes(x = Time/7, ymin = 0, ymax = Prop_Infected), 
              fill = "#CC79A7", alpha = 0.5)+
  geom_vline(xintercept = pars[set-300,]$peak_I_1/7, linetype=3, size=1)

plot <- plot_grid(A,B, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot

