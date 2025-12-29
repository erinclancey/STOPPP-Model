set=80

mk_legend <- ggplot(data = data.tau.list[[set]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
                      labels=c("s",expression(iota),expression(r[A]),expression(r[T])))+
  theme(plot.title = element_text(hjust = 0))
legend <- cowplot::get_legend(mk_legend)
grid.newpage()
grid.draw(legend) 

combo1 <- ggplot(data = data.tau.list[[set]], aes(x = t)) +
  theme_minimal() +
  ylab("Proportion") + xlab("") +
  
  # --- YEAR LABELS --- 
  annotate("text", x = (1 + 365) / 2, y = 0.95, label = "Year 1", size = 5, fontface = "bold", color = "grey30") + 
  annotate("text", x = (366 + max(data.tau.list[[set]]$t)) / 2, y = 0.95, label = "Year 2", size = 5, fontface = "bold", color = "grey30") +
  
  # --- BACKGROUND SHADING ---
  annotate("rect", xmin = 1, xmax = 365, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.3) +
  annotate("rect", xmin = 366, xmax = max(data.tau.list[[set]]$t),
           ymin = -Inf, ymax = Inf,
           fill = "lightyellow", alpha = 0.3) +
  # ---------------------------

geom_line(aes(y = S/N, color = "S"), linewidth = 0.75) +
  geom_line(aes(y = I/N, color = "I"), linewidth = 0.75) +
  geom_line(aes(y = R_a/N, color = "R_a"), linewidth = 0.75) +
  geom_line(aes(y = R_i/N, color = "R_i"), linewidth = 0.75) +
  
  scale_colour_manual(
    "",
    breaks = c("S","I","R_a","R_i"),
    values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
    labels = c("s", expression(iota), expression(r[A]), expression(r[T]))
  ) +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "none") +
  scale_x_continuous(
    breaks = c(1, seq(91.25, 365*2, by = 91.25))
  ) + ylim(0,1)

combo2 <- ggplot(data = data.tau.list[[100+set]], aes(x = t)) +
  theme_minimal() +
  ylab("Proportion") + xlab("") +
  
  annotate("rect", xmin = 1, xmax = 365, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.3) +
  annotate("rect", xmin = 366, xmax = max(data.tau.list[[100+set]]$t),
           ymin = -Inf, ymax = Inf,
           fill = "lightyellow", alpha = 0.3) +
  
  geom_line(aes(y = S/N, color = "S"), linewidth = 0.75) +
  geom_line(aes(y = I/N, color = "I"), linewidth = 0.75) +
  geom_line(aes(y = R_a/N, color = "R_a"), linewidth = 0.75) +
  geom_line(aes(y = R_i/N, color = "R_i"), linewidth = 0.75) +
  
  scale_colour_manual(
    "",
    breaks = c("S","I","R_a","R_i"),
    values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
    labels = c("s", expression(iota), expression(r[A]), expression(r[T]))
  ) +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "none") +
  scale_x_continuous(
    breaks = c(1, seq(91.25, 365*2, by = 91.25))
  )+ ylim(0,1)

combo3 <- ggplot(data = data.tau.list[[200+set]], aes(x = t)) +
  theme_minimal() +
  ylab("Proportion") + xlab("") +
  
  annotate("rect", xmin = 1, xmax = 365, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.3) +
  annotate("rect", xmin = 366, xmax = max(data.tau.list[[200+set]]$t),
           ymin = -Inf, ymax = Inf,
           fill = "lightyellow", alpha = 0.3) +
  
  geom_line(aes(y = S/N, color = "S"), linewidth = 0.75) +
  geom_line(aes(y = I/N, color = "I"), linewidth = 0.75) +
  geom_line(aes(y = R_a/N, color = "R_a"), linewidth = 0.75) +
  geom_line(aes(y = R_i/N, color = "R_i"), linewidth = 0.75) +
  
  scale_colour_manual(
    "",
    breaks = c("S","I","R_a","R_i"),
    values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
    labels = c("s", expression(iota), expression(r[A]), expression(r[T]))
  ) +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "none") +
  scale_x_continuous(
    breaks = c(1, seq(91.25, 365*2, by = 91.25))
  )+ ylim(0,1)

x.grob <- textGrob("Time (Days)", gp=gpar(fontsize=13))
plot <- plot_grid(combo1,combo2,combo3, ncol = 1, nrow = 3, rel_heights=c(1,1,1))
plot2 <- grid.arrange(arrangeGrob(plot, bottom = x.grob))
plot3 <- plot_grid(plot2, legend, ncol = 2, rel_widths = c(0.9, 0.1))
plot3

