############EVEN
peak1_P <- rep(NA, length(even.sample.list))
peak2_P <- rep(NA, length(even.sample.list))
peak1_T <- rep(NA, length(even.sample.list))
peak2_T <- rep(NA, length(even.sample.list))
pars_id <- rep(NA, length(even.sample.list))
peak1_gamma <- rep(NA, length(even.sample.list))
peak2_gamma <- rep(NA, length(even.sample.list))
type <- rep(NA, length(even.sample.list))
combo <- rep(NA, length(even.sample.list))
#bandwidth hs to be high (above 50) to capture the cycle. n.points is set to default max(100,length(x))
for(i in 1:length(even.sample.list)){
  df <- even.sample.list[[i]]
  kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 50)
  R_t_k <- diff(kspline$y)
  rA_df <- data.frame(kspline$x, kspline$y)
  low_rA <- subset(rA_df,  kspline.x>0 &  kspline.x<170)
  high_rA <- subset(rA_df,  kspline.x>170 &  kspline.x<360)
  peak1_gamma[i] <- floor(subset(low_rA, kspline.y==max(low_rA$kspline.y))$kspline.x[1]-(1/pars[even.sample.list[[i]]$pars_id[1],]$gamma))
  peak2_gamma[i] <- floor(subset(high_rA, kspline.y==max(high_rA$kspline.y))$kspline.x[1]-(1/pars[even.sample.list[[i]]$pars_id[1],]$gamma))
  kspline_I <- ((R_t_k + kspline$y[-1]*pars[even.sample.list[[i]]$pars_id[1],]$omega_a)/pars[even.sample.list[[i]]$pars_id[1],]$gamma)
  t_spline <- kspline$x[-1]
  spline_df <- data.frame(t_spline, kspline_I)
  low_pred <- subset(spline_df,  t_spline>0 &  t_spline<170)
  high_pred <- subset(spline_df,  t_spline>170 &  t_spline<360)
  peak1_P[i] <- floor(subset(low_pred, kspline_I==max(low_pred$kspline_I))$t_spline[1])
  peak2_P[i] <- floor(subset( high_pred, kspline_I==max(high_pred$kspline_I))$t_spline[1])
  pars_id[i] <- even.sample.list[[i]]$pars_id[1]
  type[i] <- even.sample.list[[i]]$type[1]
  combo[i] <-pars[even.sample.list[[i]]$pars_id[1],]$combo[1]
  peak1_T[i] <- pars[even.sample.list[[i]]$pars_id[1],]$peak1_T
  peak2_T[i] <- pars[even.sample.list[[i]]$pars_id[1],]$peak2_T
}
#gamma=1/10 so subtract 10 to get the prediction
peaks_even_df <- data.frame(peak1_P,peak2_P,peak1_T,peak2_T,peak1_gamma,peak2_gamma,pars_id,type,combo)

bouts.list <- list()
for(j in 1:10){
  n_trials=20
  n_reps=10
  T_peak_sample <- vector()
  Pred_peak_sample <- vector()
  random_sample <- vector()
  gamma_sample <- vector()
  combo <- vector()
  type <- vector()
  for(i in 1:length(even.sample.list)){
    data=data.tau.list[[even.sample.list[[i]]$pars_id[1]]]
    combo[i]=pars[even.sample.list[[i]]$pars_id[1],]$combo
    type[i]=peaks_even_df$type[i]
    Y_random=data[sample(nrow(data), 2, replace = FALSE), ] 
    Y_true1 <- subset(data, t==pars[even.sample.list[[i]]$pars_id[1],]$peak1_T) 
    Y_true2 <- subset(data, t==pars[even.sample.list[[i]]$pars_id[1],]$peak2_T) 
    I_sample_T1 <- with(Y_true1, rbinom(n_reps, n_trials, I/N))
    I_sample_T2 <- with(Y_true2, rbinom(n_reps, n_trials, I/N))
    Y_pred1 <- subset(data, t==peaks_even_df$peak1_P[i])
    Y_pred2 <- subset(data, t==peaks_even_df$peak2_P[i]) 
    I_sample_P1 <- with(Y_pred1, rbinom(n_reps, n_trials, I/N))
    I_sample_P2 <- with(Y_pred2, rbinom(n_reps, n_trials, I/N))
    I_sample_R1 <- with(Y_random[1,], rbinom(n_reps, n_trials, I/N))
    I_sample_R2 <- with(Y_random[2,], rbinom(n_reps, n_trials, I/N))
    Y_gamma1 <- subset(data, t==peaks_even_df$peak1_gamma[i])
    Y_gamma2 <- subset(data, t==peaks_even_df$peak2_gamma[i]) 
    I_sample_gamma1 <- with(Y_gamma1, rbinom(n_reps, n_trials, I/N))
    I_sample_gamma2 <- with(Y_gamma2, rbinom(n_reps, n_trials, I/N))
    I_sample_T <- c(I_sample_T1,I_sample_T2)
    I_sample_P <- c(I_sample_P1,I_sample_P2)
    I_sample_R <- c(I_sample_R1,I_sample_R2)
    I_sample_gamma <- c(I_sample_gamma1,I_sample_gamma2)
    I_sample_T[I_sample_T > 0] <- 1
    I_sample_P[I_sample_P > 0] <- 1
    I_sample_R[I_sample_R > 0] <- 1
    I_sample_gamma[I_sample_gamma > 0] <- 1
    T_peak_sample[i] = mean(I_sample_T)
    Pred_peak_sample[i] = mean(I_sample_P)
    random_sample[i] = mean(I_sample_R)
    gamma_sample[i] = mean(I_sample_gamma)
  }
  bouts.list[[j]] <- data.frame(combo,type,T_peak_sample,Pred_peak_sample,gamma_sample,random_sample)
}
results <- do.call(rbind, bouts.list)
D.long <- melt(results, id=c("combo", "type"))
D.long$type <- factor(D.long$type, levels=c('daily','weekly','biweekly','monthly','bimonthly'))
col_names <- c("1"="Low", "2"="Medium", "3"="High")
row_names <- c("daily"="Daily", "weekly"="Weekly", "biweekly"="Bi-Weekly","monthly"="Monthly","bimonthly"="Bi-Monthly")
ggplot(D.long, aes(x=variable, y=value, fill=variable))+ylab("Proportion Successful Sampling Bouts")+xlab("Sampling Time Point")+
  geom_boxplot(outlier.shape = NA)+scale_fill_brewer(NULL,palette="Dark2",labels=c(TeX("True $(\\iota)$"), TeX("Predicted $(\\hat{\\iota})$"),TeX("Null ($1/\\gamma$)"),'Null (random)'))+
  theme(axis.text.x=element_blank())+ theme(legend.position="bottom",legend.text=element_text(size=16))+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),strip.text = element_text(
    size = 16))+
  facet_grid(cols=vars(combo),rows = vars(type),labeller = labeller(.rows = row_names , .cols = col_names ))

#10x12 

####Sim Stats Table
columns= c("Amplitude","Mean Amplitude","Sampling Design","Dist. to True") 
simple_stats = data.frame(matrix(nrow = 15, ncol = length(columns))) 
colnames(simple_stats) = columns
#############LOW
low <- subset(peaks_even_df, combo==1)
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
med <- subset(peaks_even_df, combo==2)
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
high <- subset(peaks_even_df, combo==3)
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




