
fit.posterior <- ra_posterior

omega_a=1/12.85714
gamma=1/1.428571
f=2/52
Week <- seq(0,51,1)
Week_new <- rep(NA, 52)
for(i in 1: length(Week)){
  if(Week[i]<16){Week_new[i]=Week[i]+52}else{
    Week_new[i]=Week[i]
  }
}

Week_new_sample <- rep(NA, length(agg_preg$Week))
for(i in 1: length(agg_preg$Week)){
  if(agg_preg$Week[i]<16){Week_new_sample[i]=agg_preg$Week[i]+52}else{
    Week_new_sample [i]=agg_preg$Week[i]
  }
}

week_cat_1 <- as.character(seq(16,51,2))
week_cat_2 <- as.character(seq(0,15,2))
week_cat <- c(week_cat_1,week_cat_2)


fit.posterior.long <- melt(fit.posterior , id=c("iter"))
modes <- fit.posterior.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode <- as.vector(modes$mode)
modes

# Re-sample the posterior
post_sample <- fit.posterior[sample(nrow(fit.posterior), 500, replace = TRUE), ] 
post_sample.long <- melt(post_sample, id=c("iter"))
modes <- post_sample.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode <- as.vector(modes$mode)
modes

cred_int=0.95
a_ci <- ci(post_sample$a_post, ci=cred_int, method = "HDI")
C1_ci <- ci(post_sample$C1_post, ci=cred_int, method = "HDI")
C2_ci <- ci(post_sample$C2_post, ci=cred_int, method = "HDI")
phi_ci <- ci(post_sample$phi_post, ci=cred_int, method = "HDI")


cred_post_sample <- subset(post_sample, a_ci$CI_low<a_post & a_post<a_ci$CI_high &
                             C1_ci$CI_low<C1_post & C1_post<C1_ci$CI_high &
                             C2_ci$CI_low<C2_post & C2_post<C2_ci$CI_high &
                             phi_ci$CI_low<phi_post & phi_post<phi_ci$CI_high)


#####
pred_list <- list()
  for(j in 1:length(cred_post_sample$iter)){
    I_pred_fit <- (ra_func_deriv(t=Week, C1=cred_post_sample$C1_post[j], a=cred_post_sample$a_post[j], 
                                 phi=cred_post_sample$phi_post[j], f=f)
                   + ra_func(t=Week, C1=cred_post_sample$C1_post[j],
                             a=cred_post_sample$a_post[j],
                             phi=cred_post_sample$phi_post[j], C2=cred_post_sample$C2_post[j], f=f)*omega_a)/gamma
    Ra_fit <- ra_func(t=Week, C1=cred_post_sample$C1_post[j], a=cred_post_sample$a_post[j], phi=cred_post_sample$phi_post[j], 
                      C2=cred_post_sample$C2_post[j], f=f)
    ID <- j
    df <- data.frame(Week, Week_new, ID, I_pred_fit, Ra_fit)
    pred_list[[j]] <- df
  }
pred_frame <- do.call(rbind, pred_list)



#effort <- agg_preg$Total_Sampled/650*10
p_I <- pred_frame %>%
  ggplot(aes(x=Week,y=I_pred_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#E69F00")+
  guides(color="none") + 
  theme_minimal()+theme_minimal()+ylab("Proportion Infected")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,51,2), limits = c(0,51))+scale_y_continuous(breaks=seq(-0.07,0.2,0.03), limits=c(-0.07,0.2))+
  #geom_vline(xintercept=agg_preg$Week,lwd=effort,colour="black", alpha=0.5)
  ggtitle(expression(paste("(b) ",italic("H. monstrosus "), "Predicted PCR Positive")))+theme(plot.title = element_text(size=15))

p_Ra <- pred_frame %>%
  ggplot(aes(x=Week,y=Ra_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#0072B2")+
  guides(color="none") + 
  theme_minimal()+theme_minimal()+ylab("Proportion Seropositive")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,51,2), limits = c(0,51))+scale_y_continuous(breaks=seq(0,0.5,0.04), limits=c(0,0.46))+
  ggtitle(expression(paste("(a) ", italic("H. monstrosus "), "Seropositive")))+theme(plot.title = element_text(size=15))
  
p_Ra <- p_Ra + geom_point(data = EBVagg_cases, aes(Week,EBVPos_prop), 
             fill = "black", color = "black", 
             size = 3, shape = 21)
  x.grob <- textGrob("Time (Weeks)", gp=gpar(fontsize=13))
  plot <- plot_grid(p_Ra, p_I, ncol = 1, nrow = 2, rel_heights=c(1,1))
  plot2 <- grid.arrange(arrangeGrob(plot, bottom = x.grob))


peak_list <- list()
for(j in 1:length(cred_post_sample$iter)){
  peak_list[[j]] <- sort(findpeaks(pred_list[[j]]$I_pred_fit, minpeakheight = 0.001,
                                   minpeakdistance = 20, npeaks = 2, nups = 0, ndowns = 0)[,2])
}
peak_vec_row <- data.frame(do.call(rbind, peak_list))
colnames(peak_vec_row) <- c("peak_1","peak_2")
peak_vec <- unlist(peak_list)
peaks <- peak_vec-1
peak_new <- rep(NA, length(peaks))
for(i in 1: length(peak_new)){
  if(peaks[i]<16){peak_new[i]=peaks[i]+52}else{
    peak_new[i]=peaks[i]
  }
}
label <- rep(NA, length(peaks))
for(i in 1: length(peak_new)){
  if(peak_new[i]<41){label[i]="peak 1"}else{
    label[i]="peak 2"
  }
}
peaks <- data.frame(peaks, peak_new, label)
colnames(peaks) <- c("peak_orig", "peak_new", "label")
peak_1 <- subset(peaks, label=="peak 1")
peak_2 <- subset(peaks, label=="peak 2")
peak1_ci <- ci(peak_1$peak_new, ci=cred_int, method = "HDI")
peak2_ci <- ci(peak_2$peak_new, ci=cred_int, method = "HDI")
mean(peak_1$peak_new)
mean(peak_2$peak_new)
peak1_mode <- peak_1 %>% summarise(mode = posterior.mode(mcmc(peak_new), adjust=1))
peak2_mode <- peak_2 %>% summarise(mode = posterior.mode(mcmc(peak_new), adjust=1))

peak1_P <- vector()
peak2_P <- vector()
peak1max <- vector()
peak2max <- vector()
min <- vector()
for(j in 1:length(cred_post_sample$iter)){
  min[j] <- min(pred_list[[j]]$Ra_fit)
  low_pred <- subset(pred_list[[j]], Week_new<40)
  high_pred <- subset(pred_list[[j]],  Week_new>=40)
  peak1_P[j] <- subset(low_pred , I_pred_fit==max(low_pred$I_pred_fit))$Week_new
  peak2_P[j] <- subset(high_pred, I_pred_fit==max(high_pred$I_pred_fit))$Week_new
  peak1max[j] <- max(low_pred$Ra_fit)
  peak2max[j] <- max(high_pred$Ra_fit)
}
ci1 <- ci(peak1_P, ci=cred_int, method = "HDI")
ci2 <- ci(peak2_P, ci=cred_int, method = "HDI")
mode1 <- Mode(peak1_P)
mode2 <- Mode(peak2_P)
peaks$peak_new_max <- c(peak1_P,peak2_P)
amp_1 <- peak1max - min
amp_2 <- peak2max - min
amp <- (amp_1+amp_2)/2
ci_amp <- ci(amp, method="HDI")
amp_bar <- mean(amp)

ci1 
ci2 
mode1 
mode2
amp_bar
ci_amp



p <- ggplot(peaks, aes(x=peak_new_max)) + theme_minimal()
p <- p + geom_density(adjust=0.5, color="black", fill="white", alpha=0)
dpb <- ggplot_build(p)
x1 <- min(which(dpb$data[[1]]$x >=ci1$CI_low))
x2 <- max(which(dpb$data[[1]]$x <=ci1$CI_high))
x3 <- min(which(dpb$data[[1]]$x >=ci2$CI_low))
x4 <- max(which(dpb$data[[1]]$x <=ci2$CI_high))
p <- p + geom_area(data=data.frame(x=dpb$data[[1]]$x[x1:x2],
                                   y=dpb$data[[1]]$y[x1:x2]),
                   aes(x=x, y=y), fill="#E7298A", alpha=0.5)
p <- p + geom_area(data=data.frame(x=dpb$data[[1]]$x[x3:x4],
                                   y=dpb$data[[1]]$y[x3:x4]),
                   aes(x=x, y=y), fill="#E7298A", alpha=0.5)
p <- p +  geom_histogram(aes(y=after_stat(density)), alpha=0.05, position="identity", binwidth=1, color="black") 
p <- p +  geom_vline(xintercept=mode1, color="black", size=0.75, linetype=2)+
  geom_vline(xintercept=mode2, color="black", size=0.75, linetype=2)
p <- p + scale_x_continuous(name="Time (Weeks)", breaks=seq(16,67,2),labels=week_cat, limits = c(16, 67))+ 
  scale_y_continuous(name="Predicted Peak Density", breaks = seq(0,0.12,0.01), limits=c(0,0.08))
p <- p + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14), 
               legend.position = "none")
p <- p + ggtitle(expression(paste(italic("H. monstrosus"))))+theme(plot.title = element_text(size=15))
plot(p)
