cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p + scale_color_manual(values = cbp1)

set=905
pred_list=pred_frame.list[[set]]
pars_id=pred_list$pars_id[1]
#View(pred_list)
pred_list <- split(pred_list, pred_list$ID)
cred_int=0.95


## fit E. helvum data using the model fitting method and calculate average amplitude
peak_list <- list()
for(j in 1:500){
  peak_list[[j]] <- sort(findpeaks(pred_list[[j]]$I_pred_fit, minpeakheight = 0.001,
                                   minpeakdistance = 20, npeaks = 2, nups = 0, ndowns = 0)[,2])
}
peaks <- data.frame(do.call(rbind, peak_list))
colnames(peaks) <- c("peak_1","peak_2")

peak1_ci <- ci(peaks$peak_1, ci=cred_int, method = "HDI")
peak2_ci <- ci(peaks$peak_2, ci=cred_int, method = "HDI")
peak1_mode <- peaks %>% summarise(mode = posterior.mode(mcmc(peak_1), adjust=1))
peak2_mode <- peaks %>% summarise(mode = posterior.mode(mcmc(peak_2), adjust=1))
peak1_P <- vector()
peak2_P <- vector()
peak1max <- vector()
peak2max <- vector()
min <- vector()
for(j in 1:500){
  min[j] <- min(pred_list[[j]]$Ra_fit)
  low_pred <- subset(pred_list[[j]], t<150)
  high_pred <- subset(pred_list[[j]],  t>150)
  peak1_P[j] <- subset(low_pred , I_pred_fit==max(low_pred$I_pred_fit))$t
  peak2_P[j] <- subset(high_pred, I_pred_fit==max(high_pred$I_pred_fit))$t
  peak1max[j] <- max(low_pred$Ra_fit)
  peak2max[j] <- max(high_pred$Ra_fit)
}
ci1 <- ci(peak1_P, ci=cred_int, method = "HDI")
ci2 <- ci(peak2_P, ci=cred_int, method = "HDI")
mode1 <- Mode(peak1_P)
mode2 <- Mode(peak2_P)
peak1_mode <- summarise(mode = posterior.mode(mcmc(peak1_P), adjust=1))
peak2_mode <- peaks %>% summarise(mode = posterior.mode(mcmc(peak_2), adjust=1))
peaks_long <- data.frame(c(peak1_P,peak2_P))
colnames(peaks_long) <- c("peak_new_max")
amp_1 <- peak1max - min
amp_2 <- peak2max - min
amp <- (amp_1+amp_2)/2
ci_amp <- ci(amp, method="HDI")
amp_bar <- mean(amp)
# print values for modes, amplitude, and credible intervals
ci1 
ci2 
mode1 
mode2
amp_bar
ci_amp

mode1 <- posterior.mode(mcmc(peak1_P), adjust=2)
mode2 <- posterior.mode(mcmc(peak2_P), adjust=2)

# plot the distribution of predicted peaks
p <- ggplot(peaks_long , aes(x=peak_new_max)) + theme_minimal()
p <- p + geom_density(adjust=1, color="black", fill="white", alpha=0)
dpb <- ggplot_build(p)
x1 <- min(which(dpb$data[[1]]$x >=ci1$CI_low))
x2 <- max(which(dpb$data[[1]]$x <=ci1$CI_high))
x3 <- min(which(dpb$data[[1]]$x >=ci2$CI_low))
x4 <- max(which(dpb$data[[1]]$x <=ci2$CI_high))
p <- p + geom_area(data=data.frame(x=dpb$data[[1]]$x[x1:x2],
                                   y=dpb$data[[1]]$y[x1:x2]),
                   aes(x=x, y=y), fill="#CC79A7", alpha=0.5)
p <- p + geom_area(data=data.frame(x=dpb$data[[1]]$x[x3:x4],
                                   y=dpb$data[[1]]$y[x3:x4]),
                   aes(x=x, y=y), fill="#CC79A7", alpha=0.5)
p <- p +  geom_histogram(aes(y=after_stat(density)), alpha=0.05, position="identity", binwidth=5, color="black") 

p <- p +  geom_vline(xintercept=mode1, color="#CC79A7", size=0.75, linetype=1)+
  geom_vline(xintercept=mode2, color="#CC79A7", size=0.75, linetype=1)

p <- p +  geom_vline(xintercept=pars[pars_id,]$peak1_T, color="black", size=0.75, linetype=2)+
  geom_vline(xintercept=pars[pars_id,]$peak2_T, color="black", size=0.75, linetype=2)

p <- p + scale_x_continuous(name="Time (Days)", breaks=seq(-20,365,20), limits=c(-20,365))+
  scale_y_continuous(name="Predicted Peak Density", breaks = seq(0,0.12,0.01))

p <- p + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14), 
               legend.position = "none")
p <- p + ggtitle("Low Amplitude")+theme(plot.title = element_text(size=15))
plot(p)
