set=60


Ra <- ggplot() +
  geom_line(data = pred_frame.list[[set]], aes(x = t, y = Ra_fit, group = ID, color = "Fitted"), 
            alpha = 1, linewidth = 0.75) +  # Keeps all blue lines
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = R_a/N, color="True"), 
            alpha=1,linewidth=0.75)+
  labs(title = "Seroprevalence", x = "Time (Days)", y = "Prop. Infected") +
  theme_minimal() +
  theme(text = element_text(size = 18))+ theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Fitted" = "#0072B2", "True" = "black"),
                     labels = c(TeX("Fitted $(\\hat{r}_a)$"), TeX("True $(r_a)$"))) +
  guides(color = guide_legend(title = ""))


I_fit <- ggplot() +
  geom_line(data = pred_frame.list[[set]], aes(x = t, y = I_pred_fit, group = ID, color = "Fitted"), 
            alpha = 1, linewidth = 0.75) +  # Keeps all blue lines
  geom_line(data = data.frame(Time = data.tau.list[[set]]$t[-1], 
                              Value = data.tau.list[[set]]$I[-1] / data.tau.list[[set]]$N[-1]),
            aes(x = Time, y = Value, color = "True"), 
            alpha = 1, linewidth = 0.75) +  # Keeps the black reference line
  labs(title = "Prevalence", x = "Time (Days)", y = "Prop. Infected") +
  theme_minimal() +
  theme(text = element_text(size = 18))+ theme(legend.position = c(0.9, 0.9))+
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_color_manual(values = c("Fitted" = "#CC79A7", "True" = "black"),
                     labels = c(TeX("Fitted $(\\hat{\\iota})$"), TeX("True $(\\iota)$"))) +
  guides(color = guide_legend(title = ""))

plot <- plot_grid(Ra,I_fit, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot
