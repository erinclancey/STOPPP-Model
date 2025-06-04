############EVEN
peak_P <- rep(NA, length(cluster.sample.list))

peak_T <- rep(NA, length(cluster.sample.list))

pars_id <- rep(NA, length(cluster.sample.list))
peak_gamma <- rep(NA, length(cluster.sample.list))

type <- rep(NA, length(cluster.sample.list))

#bandwidth hs to be high (above 50) to capture the cycle. n.points is set to default max(100,length(x))
for(i in 1:length(cluster.sample.list)){
  df <- cluster.sample.list[[i]]
  kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth=75, n.points=500)
  
  # plot(df$t, df$R_a_sample/n_trials)
  # lines(kspline$x,kspline$y, col="red")
  
  R_t_k <- diff(kspline$y)
  rA_df <- data.frame(kspline$x, kspline$y)

  peak_gamma[i] <- floor(subset(rA_df, kspline.y==max(rA_df$kspline.y))$kspline.x[1]-
                           (1/pars[cluster.sample.list[[i]]$pars_id[1],]$gamma))
  
  kspline_I <- ((R_t_k + kspline$y[-1]*pars[cluster.sample.list[[i]]$pars_id[1],]$omega_a)/
                  pars[cluster.sample.list[[i]]$pars_id[1],]$gamma)
  t_spline <- kspline$x[-1]
  spline_df <- data.frame(t_spline, kspline_I)
  
  peak_P[i] <- floor(subset(spline_df, kspline_I==max(spline_df$kspline_I))$t_spline[1])
  
  pars_id[i] <- cluster.sample.list[[i]]$pars_id[1]
  type[i] <- cluster.sample.list[[i]]$type[1]
  peak_T[i] <- pars[cluster.sample.list[[i]]$pars_id[1],]$peak_I
}
#gamma=1/10 so subtract 10 to get the prediction
peaks_cluster_df <- data.frame(peak_P,peak_T,peak_gamma,pars_id,type)

bouts.list <- list()
for(j in 1:10){
  n_trials=20
  n_reps=10
  T_peak_sample <- vector()
  Pred_peak_sample <- vector()
  random_sample <- vector()
  gamma_sample <- vector()
  type <- vector()
  for(i in 1:length(cluster.sample.list)){
    data=data.tau.list[[cluster.sample.list[[i]]$pars_id[1]]]
    type[i]=peaks_cluster_df$type[i]
    Y_random=data[sample(nrow(data), 1, replace = FALSE), ] 
    Y_true <- subset(data, t==pars[cluster.sample.list[[i]]$pars_id[1],]$peak_I) 
    I_sample_T <- with(Y_true, rbinom(n_reps, n_trials, I/N))

    Y_pred <- subset(data, t==peaks_cluster_df$peak_P[i])

    I_sample_P <- with(Y_pred, rbinom(n_reps, n_trials, I/N))

    I_sample_R <- with(Y_random[1,], rbinom(n_reps, n_trials, I/N))

    Y_gamma <- subset(data, t==peaks_cluster_df$peak_gamma[i])

    I_sample_gamma <- with(Y_gamma, rbinom(n_reps, n_trials, I/N))
    
    I_sample_T[I_sample_T > 0] <- 1
    I_sample_P[I_sample_P > 0] <- 1
    I_sample_R[I_sample_R > 0] <- 1
    I_sample_gamma[I_sample_gamma > 0] <- 1
    T_peak_sample[i] = mean(I_sample_T)
    Pred_peak_sample[i] = mean(I_sample_P)
    random_sample[i] = mean(I_sample_R)
    gamma_sample[i] = mean(I_sample_gamma)
  }
  bouts.list[[j]] <- data.frame(type,T_peak_sample,Pred_peak_sample,gamma_sample,random_sample)
}
results <- do.call(rbind, bouts.list)
D.long <- melt(results, id=c("type"))
D.long$type <- factor(D.long$type, levels=c('42_random','3_day','42_week'))
row_names <- c('42_random'='Random Days','3_day'='3-Day Clusters','42_week'='Week Clusters')
ggplot(D.long, aes(x=variable, y=value, fill=variable))+ylab("Proportion Successful Sampling Bouts")+xlab("Sampling Time Point")+
  geom_boxplot(outlier.shape = NA)+scale_fill_brewer(palette="Paired",labels=c(TeX("True $(\\iota)$"), TeX("Predicted $(\\hat{\\iota})$"),TeX("Null ($1/\\gamma$)"),'Null (random)'))+
  theme(axis.text.x=element_blank())+ theme(legend.position="bottom",legend.text=element_text(size=16))+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),strip.text = element_text(
    size = 16), plot.title = element_text(hjust = 0.5, size=20)) +ggtitle("Kernal Smoothing")+
  facet_grid(rows = vars(type),labeller = labeller(.rows = row_names))

#10x12 

####Sim Stats Table
columns= c("Amplitude","Mean Amplitude","Sampling Design","Dist. to True") 
simple_stats = data.frame(matrix(nrow = 15, ncol = length(columns))) 
colnames(simple_stats) = columns
#############LOW
low <- subset(peaks_cluster_df, combo==1)
#Daily
daily <- subset(low, type=="daily")
simple_stats$Amplitude[c(1:5)]="Low"
simple_stats$`Mean Amplitude`[c(1:5)]=round(mean(subset(pars, combo==1)$amp_bar),2)
simple_stats$`Sampling Design`[c(1:5)]=c("Daily", "Weekly", "Bi-Weekly","Monthly","Bi-Monthly")
peak_true=c(daily$peak1_T, daily$peak2_T)
peak_pred=c(daily$peak1_P, daily$peak2_P)
simple_stats[1,4]= round(mean(abs(peak_true-peak_pred)),2)
#weekly
weekly <- subset(low, type=="weekly")
peak_true=c(weekly$peak1_T, weekly$peak2_T)
peak_pred=c(weekly$peak1_P, weekly$peak2_P)
simple_stats[2,4]= round(mean(abs(peak_true-peak_pred)),2)
#biweekly
biweekly <- subset(low, type=="biweekly")
peak_true=c(biweekly $peak1_T, biweekly $peak2_T)
peak_pred=c(biweekly $peak1_P, biweekly $peak2_P)
simple_stats[3,4]= round(mean(abs(peak_true-peak_pred)),2)
#monthly
monthly <- subset(low, type=="monthly")
peak_true=c(monthly$peak1_T, monthly$peak2_T)
peak_pred=c(monthly$peak1_P, monthly$peak2_P)
simple_stats[4,4]= round(mean(abs(peak_true-peak_pred)),2)
#bimonthly
bimonthly <- subset(low, type=="bimonthly")
peak_true=c(bimonthly$peak1_T, bimonthly$peak2_T)
peak_pred=c(bimonthly$peak1_P, bimonthly$peak2_P)
simple_stats[5,4]= round(mean(abs(peak_true-peak_pred)),2)
#############Medium
med <- subset(peaks_cluster_df, combo==2)
#Daily
daily <- subset(med, type=="daily")
simple_stats$Amplitude[c(6:10)]="Med"
simple_stats$`Mean Amplitude`[c(6:10)]=round(mean(subset(pars, combo==2)$amp_bar),2)
simple_stats$`Sampling Design`[c(6:10)]=c("Daily", "Weekly", "Bi-Weekly","Monthly","Bi-Monthly")
peak_true=c(daily$peak1_T, daily$peak2_T)
peak_pred=c(daily$peak1_P, daily$peak2_P)
simple_stats[6,4]= round(mean(abs(peak_true-peak_pred)),2)
#weekly
weekly <- subset(med, type=="weekly")
peak_true=c(weekly$peak1_T, weekly$peak2_T)
peak_pred=c(weekly$peak1_P, weekly$peak2_P)
simple_stats[7,4]= round(mean(abs(peak_true-peak_pred)),2)
#biweekly
biweekly <- subset(med, type=="biweekly")
peak_true=c(biweekly $peak1_T, biweekly $peak2_T)
peak_pred=c(biweekly $peak1_P, biweekly $peak2_P)
simple_stats[8,4]= round(mean(abs(peak_true-peak_pred)),2)
#monthly
monthly <- subset(med, type=="monthly")
peak_true=c(monthly$peak1_T, monthly$peak2_T)
peak_pred=c(monthly$peak1_P, monthly$peak2_P)
simple_stats[9,4]= round(mean(abs(peak_true-peak_pred)),2)
#bimonthly
bimonthly <- subset(med, type=="bimonthly")
peak_true=c(bimonthly$peak1_T, bimonthly$peak2_T)
peak_pred=c(bimonthly$peak1_P, bimonthly$peak2_P)
simple_stats[10,4]= round(mean(abs(peak_true-peak_pred)),2)
#############High
high <- subset(peaks_cluster_df, combo==3)
#Daily
daily <- subset(high, type=="daily")
simple_stats$Amplitude[c(11:15)]="High"
simple_stats$`Mean Amplitude`[c(11:15)]=round(mean(subset(pars, combo==3)$amp_bar),2)
simple_stats$`Sampling Design`[c(11:15)]=c("Daily", "Weekly", "Bi-Weekly","Monthly","Bi-Monthly")
peak_true=c(daily$peak1_T, daily$peak2_T)
peak_pred=c(daily$peak1_P, daily$peak2_P)
simple_stats[11,4]= round(mean(abs(peak_true-peak_pred)),2)
#weekly
weekly <- subset(high, type=="weekly")
peak_true=c(weekly$peak1_T, weekly$peak2_T)
peak_pred=c(weekly$peak1_P, weekly$peak2_P)
simple_stats[12,4]= round(mean(abs(peak_true-peak_pred)),2)
#biweekly
biweekly <- subset(high, type=="biweekly")
peak_true=c(biweekly $peak1_T, biweekly $peak2_T)
peak_pred=c(biweekly $peak1_P, biweekly $peak2_P)
simple_stats[13,4]= round(mean(abs(peak_true-peak_pred)),2)
#monthly
monthly <- subset(high, type=="monthly")
peak_true=c(monthly$peak1_T, monthly$peak2_T)
peak_pred=c(monthly$peak1_P, monthly$peak2_P)
simple_stats[14,4]= round(mean(abs(peak_true-peak_pred)),2)
#bimonthly
bimonthly <- subset(high, type=="bimonthly")
peak_true=c(bimonthly$peak1_T, bimonthly$peak2_T)
peak_pred=c(bimonthly$peak1_P, bimonthly$peak2_P)
simple_stats[15,4]= round(mean(abs(peak_true-peak_pred)),2)
simple_stats


print(xtable(simple_stats, type='latex'))




