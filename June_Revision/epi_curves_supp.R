low <- ggplot(data = data.tau.list[[set]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
                      labels=c("s",expression(iota),expression(r[A]),expression(r[T])))+
  theme(plot.title = element_text(hjust = 0))
legend <- cowplot::get_legend(low)
grid.newpage()
grid.draw(legend) 

low2 <- ggplot(data = data.tau.list[[set]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
                      labels=c("s",expression(iota),expression(r[A]),expression(r[T])))+
  theme(plot.title = element_text(hjust = 0))+theme(legend.position="none")


med <- ggplot(data = data.tau.list[[100+set]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
                      labels=c("s",expression(iota),expression(r[A]),expression(r[T])))+
  theme(plot.title = element_text(hjust = 0))+theme(legend.position="none")

high <- ggplot(data = data.tau.list[[200+set]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#CC79A7", "#0072B2", "#69b3a2"),
                      labels=c("s",expression(iota),expression(r[A]),expression(r[T])))+
  theme(plot.title = element_text(hjust = 0))+theme(legend.position="none")


x.grob <- textGrob("Time (Days)", gp=gpar(fontsize=13))
plot <- plot_grid(low2, med, high, ncol = 1, nrow = 3, rel_heights=c(1,1,1))
plot2 <- grid.arrange(arrangeGrob(plot, bottom = x.grob))
plot3 <- plot_grid(plot2, legend, ncol = 2, rel_widths = c(0.9, 0.1))
plot3

