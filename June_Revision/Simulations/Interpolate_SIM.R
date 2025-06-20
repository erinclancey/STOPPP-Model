############EVEN
peak_P <- rep(NA, length(sample.list))
peak_P_low <- rep(NA, length(sample.list))
peak_P_high <- rep(NA, length(sample.list))
peak_T_inCI <- rep(NA, length(sample.list))
peak_T <- rep(NA, length(sample.list))
pars_id <- rep(NA, length(sample.list))
peak_gamma <- rep(NA, length(sample.list))
type <- rep(NA, length(sample.list))

#bandwidth hs to be high (above 50) to capture the cycle. n.points is set to default max(100,length(x))
for(i in 1:length(sample.list)){
  df <- sample.list[[i]]
  kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth=50, n.points=500)
  
  R_t_k <- diff(kspline$y)
  rA_df <- data.frame(kspline$x, kspline$y)
  rA_df <- subset(rA_df, kspline.x>50)
  peak_gamma[i] <- floor(subset(rA_df, kspline.y==max(rA_df$kspline.y))$kspline.x[1]-
                           (1/pars[sample.list[[i]]$pars_id[1],]$gamma))
  t_spline <- kspline$x[-1]
  kspline_I <- ((R_t_k +
               kspline$y[-1] *
           (pars[sample.list[[i]]$pars_id[1],]$omega_a+with(pars[sample.list[[i]]$pars_id[1],], b_func(t=t_spline, g=g, s=s, psi=psi, f=f)) )  )/
                pars[sample.list[[i]]$pars_id[1],]$gamma)
  spline_df <- data.frame(t_spline, kspline_I)
  quantile <- subset(spline_df, kspline_I>quantile(kspline_I, probs = 0.90))
  peak_P_low[i] <- floor(min(quantile$t_spline))
  peak_P_high[i] <- floor(max(quantile$t_spline))
  
  peak_P[i] <- floor(subset(spline_df, kspline_I==max(spline_df$kspline_I))$t_spline[1])
  pars_id[i] <- sample.list[[i]]$pars_id[1]
  type[i] <- sample.list[[i]]$type[1]
  peak_T[i] <- pars[sample.list[[i]]$pars_id[1],]$peak_I
  peak_T_inCI[i] <- peak_P_low[i]<=peak_T[i] & peak_P_high[i]>=peak_T[i]
}
#gamma=1/10 so subtract 10 to get the prediction
peaks_df <- data.frame(peak_P,peak_P_low, peak_P_high, peak_T_inCI, peak_T,peak_gamma,pars_id,type)
peaks_df$CIsize <- peaks_df$peak_P_high-peaks_df$peak_P_low
peaks_df$dist <- abs(peaks_df$peak_P-peaks_df$peak_T)


bouts.list <- list()
for(j in 1:10){
  n_trials=20
  n_reps=10
  T_peak_sample <- vector()
  Pred_peak_sample <- vector()
  random_sample <- vector()
  gamma_sample <- vector()
  type <- vector()
  for(i in 1:length(sample.list)){
    data=data.tau.list[[sample.list[[i]]$pars_id[1]]]
    type[i]=peaks_df$type[i]
    Y_random=data[sample(nrow(data), 1, replace = FALSE), ] 
    Y_true <- subset(data, t==pars[sample.list[[i]]$pars_id[1],]$peak_I) 
    I_sample_T <- with(Y_true, rbinom(n_reps, n_trials, I/N))
    Y_pred <- subset(data, t==peaks_df$peak_P[i])
    I_sample_P <- with(Y_pred, rbinom(n_reps, n_trials, I/N))
    I_sample_R <- with(Y_random[1,], rbinom(n_reps, n_trials, I/N))
    Y_gamma <- subset(data, t==peaks_df$peak_gamma[i])
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
D.long$type <- factor(D.long$type, levels=c('daily','weekly','monthly', 'bimonthly'))
row_names <- c("daily"="Daily", "weekly"="Weekly","monthly"="Monthly", "bimonthly"="Bimonthly")
ggplot(D.long, aes(x=variable, y=value, fill=variable))+ylab("Proportion Successful Sampling Bouts")+xlab("Sampling Time Point")+
  geom_boxplot(outlier.shape = 20)+
  scale_fill_brewer(NULL,palette="Paired",
                    labels=c(TeX("True $(\\iota)$"), 
                             TeX("Predicted $(\\hat{\\iota})$"),
                             TeX("Null $(\\gamma^{-1})$   "),
                             TeX("Null $(random)$")))+
  theme(axis.text.x=element_blank())+ theme(legend.position="bottom",legend.text=element_text(size=16))+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),strip.text = element_text(
    size = 16),plot.title = element_text(hjust = 0, size=20)) +ggtitle("")+
  facet_grid(rows = vars(type),labeller = labeller(.rows = row_names))

 

####Sim Stats Table
columns= c("Sampling Design","Pred Dist. to True", "gamma Dist to True", "T in CI") 
simple_stats = data.frame(matrix(nrow = 4, ncol = length(columns))) 
colnames(simple_stats) = columns
#Daily
daily <- subset(peaks_df, type=="daily")
simple_stats$`Sampling Design`[c(1:4)]=c("Daily", "Weekly","Monthly","Bimonthly")
simple_stats[1,2]= with(daily, round(mean(abs(peak_T-peak_P)),2))
simple_stats[1,3]= with(daily, round(mean(abs(peak_T-peak_gamma)),2))
simple_stats[1,4]= round(length(daily$peak_T_inCI[daily$peak_T_inCI==TRUE]) / nrow(daily),2)

#weekly
weekly <- subset(peaks_df, type=="weekly")
simple_stats[2,2]= with(weekly, round(mean(abs(peak_T-peak_P)),2))
simple_stats[2,3]= with(weekly, round(mean(abs(peak_T-peak_gamma)),2))
simple_stats[2,4]= round(length(weekly$peak_T_inCI[weekly$peak_T_inCI==TRUE]) / nrow(weekly),2)

#monthly
monthly <- subset(peaks_df, type=="monthly")
simple_stats[3,2]= with(monthly, round(mean(abs(peak_T-peak_P)),2))
simple_stats[3,3]= with(monthly, round(mean(abs(peak_T-peak_gamma)),2))
simple_stats[3,4]= round(length(monthly$peak_T_inCI[monthly$peak_T_inCI==TRUE]) / nrow(monthly), 2)
#bimonthly
bimonthly <- subset(peaks_df, type=="bimonthly")
simple_stats[4,2]= with(bimonthly, round(mean(abs(peak_T-peak_P)),2))
simple_stats[4,3]= with(bimonthly, round(mean(abs(peak_T-peak_gamma)),2))
simple_stats[4,4]= round(length(bimonthly$peak_T_inCI[bimonthly$peak_T_inCI==TRUE]) / nrow(bimonthly),2)
simple_stats

print(xtable(simple_stats, type='latex'))




