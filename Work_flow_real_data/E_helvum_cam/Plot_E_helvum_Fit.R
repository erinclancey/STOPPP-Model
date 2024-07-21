# R script to plot the fitted E. helvum data
##################################################
# name vectors and parameters
fit.posterior <- ra_posterior 
omega_a=1/12.85714
gamma=1/1.428571
f=1/52
Week <- seq(0,51,1)

# calculate modes from all posterior values
fit.posterior.long <- melt(fit.posterior , id=c("iter"))
modes <- fit.posterior.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode <- as.vector(modes$mode)

# Re-sample the posterior for plotting
post_sample <- fit.posterior[sample(nrow(fit.posterior), 500, replace = TRUE), ] 
post_sample.long <- melt(post_sample, id=c("iter"))
modes <- post_sample.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode <- as.vector(modes$mode)

# create credible intervals for plotting
cred_int=0.95
a_ci <- ci(post_sample$a_post, ci=cred_int, method = "HDI")
C1_ci <- ci(post_sample$C1_post, ci=cred_int, method = "HDI")
C2_ci <- ci(post_sample$C2_post, ci=cred_int, method = "HDI")
phi_ci <- ci(post_sample$phi_post, ci=cred_int, method = "HDI")
cred_post_sample <- subset(post_sample, a_ci$CI_low<a_post & a_post<a_ci$CI_high &
                             C1_ci$CI_low<C1_post & C1_post<C1_ci$CI_high &
                             C2_ci$CI_low<C2_post & C2_post<C2_ci$CI_high &
                             phi_ci$CI_low<phi_post & phi_post<phi_ci$CI_high)

# create fitted values for plotting
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
  df <- data.frame(Week,ID, I_pred_fit, Ra_fit)
  pred_list[[j]] <- df
}
pred_frame <- do.call(rbind, pred_list)

############################################################
# fit E. helvum data using the interpolation method
kspline <- ksmooth(EBOVagg_cases$Week, EBOVagg_cases$EBOVPos/EBOVagg_cases$Total_Sampled, "normal", 
                   bandwidth=20, n.points=100)
R_t_k <- diff(kspline$y)
kspline_I <- ((R_t_k + kspline$y[-1]*omega_a)/gamma)
R_spline <- kspline$y[-1]
t_spline <- kspline$x[-1]
spline_df <- data.frame(t_spline, kspline_I,R_spline)
# calculate peak week
floor(subset(spline_df, kspline_I==max(spline_df$kspline_I))$t_spline[1])

# fit E. helvum data using the model fitting method and calculate average amplitude
max <- vector()
min <- vector()
peak_list <- list()
for(j in 1:length(cred_post_sample$iter)){
  peak_list[[j]] <- subset(pred_list[[j]], I_pred_fit==max(pred_list[[j]]$I_pred_fit))$Week
  max[j] <- max(pred_list[[j]]$Ra_fit)
  min[j] <- min(pred_list[[j]]$Ra_fit)
}
peak_vec <- data.frame(do.call(rbind, peak_list))
colnames(peak_vec) <- c("peak")
peak_ci <- ci(peak_vec$peak, ci=cred_int, method = "HDI")
mean(peak_vec$peak)
peak_mode <- peak_vec %>% summarise(mode = posterior.mode(mcmc(peak), adjust=1))
# calculate peak mode and CI
peak_mode
peak_ci
amp <- max-min
ci_amp <- ci(amp, method="HDI")
amp_bar <- mean(amp)
#calculate avearge amplitude of the cycle
amp_bar

#########################################
# plot fitted curves for seroprevalenc and prevalence using both fitting methods

# plot seroprevalence
p_Ra <- pred_frame %>%
  ggplot(aes(x=Week,y=Ra_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#0072B2")+
  guides(color="none") + 
  theme_minimal()+theme_minimal()+ylab("Proportion Seropositive")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,51,3), limits = c(0,51))+scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
  ggtitle(expression(paste("(a) ",italic("E. helvum "), "Seropositive Model Fitting")))+theme(plot.title = element_text(size=15))
p_Ra <- p_Ra + geom_point(data = EBOVagg_cases, aes(Week,EBOVPos_prop), 
                          fill = "black", color = "black", 
                          size = 1, shape = 21)
p_Ra_spline <- ggplot()+
    geom_point(data = spline_df, aes(x=t_spline, y=R_spline), 
                          fill = "#0072B2", color = "black", 
                          size = 2, shape = 21)+ geom_line(data = spline_df, aes(x=t_spline, y=R_spline),color = "black")+
  theme_minimal()+theme_minimal()+ylab("Proportion Seropositive")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,51,3), limits = c(0,51))+scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
  ggtitle(expression(paste("(c) ",italic("E. helvum "), "Seropositive Interpolation")))+theme(plot.title = element_text(size=15))+
  geom_point(data = EBOVagg_cases, aes(Week,EBOVPos_prop), 
               fill = "black", color = "black", 
               size = 1, shape = 21)
#plot prevalence
p_I <- pred_frame %>%
  ggplot(aes(x=Week,y=I_pred_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#E69F00")+
  guides(color="none") + 
  theme_minimal()+theme_minimal()+ylab("Proportion Infected")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,51,3), limits = c(0,51))+scale_y_continuous(breaks=seq(-0.06,0.22,0.04), limits=c(-0.06,0.22))+
  ggtitle(expression(paste("(b) ",italic("E. helvum "), "PCR Positive Model Fitting")))+theme(plot.title = element_text(size=15))
x.grob <- textGrob("Time (Weeks)", gp=gpar(fontsize=13))
p_I_spline <- ggplot()+
  theme_minimal()+theme_minimal()+ylab("Proportion Infected")+xlab("")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  scale_x_continuous(breaks=seq(0,51,3), limits = c(0,51))+scale_y_continuous(breaks=seq(-0.06,0.22,0.04), limits=c(-0.06,0.22))+
  ggtitle(expression(paste("(d) ",italic("E. helvum "), "PCR Positive Interpolation")))+theme(plot.title = element_text(size=15))+
  geom_point(data = spline_df, aes(x=t_spline, y=kspline_I), fill = "#E69F00", color = "black", size = 2, shape = 21)

x.grob <- textGrob("Time (Weeks)", gp=gpar(fontsize=13))
plot <- plot_grid(p_Ra,p_Ra_spline, p_I, p_I_spline,ncol = 2, nrow = 2, rel_heights=c(1,1))
plot2 <- grid.arrange(arrangeGrob(plot, bottom = x.grob))
#8x6
# plot the distribution of predicted peaks
p <- ggplot(peak_vec, aes(x=peak)) + theme_minimal()
p <- p + geom_density(adjust=3, color="black", fill="white", alpha=0.5)
p <- p +  geom_histogram(aes(y=after_stat(density)), alpha=0.05, position="identity", binwidth=1, color="black") 
dpb <- ggplot_build(p)
x1 <- min(which(dpb$data[[1]]$x >=peak_ci$CI_low))
x2 <- max(which(dpb$data[[1]]$x <=peak_ci$CI_high))
p <- p + geom_area(data=data.frame(x=dpb$data[[1]]$x[x1:x2],
                                   y=dpb$data[[1]]$y[x1:x2]),
                   aes(x=x, y=y), fill="#66A61E", alpha=0.5)
p <- p +  geom_histogram(aes(y=after_stat(density)), alpha=0.05, position="identity", binwidth=1, color="black") 
p <- p +  geom_vline(xintercept=peak_mode$mode, color="black", linewidth=0.75, linetype=2)
p <- p + scale_x_continuous(name="Time (weeks)", breaks=seq(0,51,2), limits=c(0,51))+
  scale_y_continuous(name="Predicted Peak Density", breaks = seq(0,0.8,0.1), limits=c(0,0.8))
p <- p + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14), 
               legend.position = "none")
p <- p + ggtitle(expression(paste(italic("E. helvum"))))+theme(plot.title = element_text(size=15))
plot(p)



