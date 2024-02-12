set=905
p1_I <- pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=I_pred_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#E69F00")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Prop. Infected")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,365*5,100))+ylim(-0.15,0.3)+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = I/N), color="black", alpha=1,linewidth=0.75)+
  ggtitle("")
p1_Ra <- pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=Ra_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#0072B2")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Prop. Seropos.")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,365*5,100))+ylim(0,0.5)+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = R_a/N), color="black", alpha=1,linewidth=0.75)+
  ggtitle("(a) Low Amplitude 42-Day Week Cluster Sample")+theme(plot.title = element_text(size=15))

set=1006
p2_I <- pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=I_pred_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#E69F00")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Prop. Infected")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,365*5,100))+ylim(-0.15,0.3)+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = I/N), color="black", alpha=1,linewidth=0.75)+
  ggtitle("")
p2_Ra <- pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=Ra_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#0072B2")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Prop. Seropos.")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,365*5,100))+ylim(0,0.5)+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = R_a/N), color="black", alpha=1,linewidth=0.75)+
  ggtitle("(b) Medium Amplitude 42-Day Week Cluster Sample")+theme(plot.title = element_text(size=15))

set=1109
p3_I <- pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=I_pred_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#E69F00")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Prop. Infected")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,365*5,100))+ylim(-0.15,0.3)+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = I/N), color="black", alpha=1,linewidth=0.75)+
  ggtitle("")
p3_Ra <- pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=Ra_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#0072B2")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Prop. Seropos.")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,365*5,100))+ylim(0,0.52)+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = R_a/N), color="black", alpha=1,linewidth=0.75)+
  ggtitle("(c) High Amplitude 42-Day Week Cluster Sample")+theme(plot.title = element_text(size=15))

x.grob <- textGrob("Time (Days)", gp=gpar(fontsize=13))
plot <- plot_grid(p1_Ra,p1_I,p2_Ra,p2_I,p3_Ra,p3_I, ncol = 2, nrow = 3, rel_heights=c(1,1,1,1,1,1))
plot2 <- grid.arrange(arrangeGrob(plot, bottom = x.grob))
plot2
