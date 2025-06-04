############EVEN
peak_P <- rep(NA, length(even.sample.list))

peak_T <- rep(NA, length(even.sample.list))

pars_id <- rep(NA, length(even.sample.list))
peak_gamma <- rep(NA, length(even.sample.list))

type <- rep(NA, length(even.sample.list))

#bandwidth hs to be high (above 50) to capture the cycle. n.points is set to default max(100,length(x))
for(i in 1:length(even.sample.list)){
  df <- even.sample.list[[i]]
  kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth=60, n.points=500)
  
  R_t_k <- diff(kspline$y)
  rA_df <- data.frame(kspline$x, kspline$y)

  peak_gamma[i] <- floor(subset(rA_df, kspline.y==max(rA_df$kspline.y))$kspline.x[1]-
                           (1/pars[even.sample.list[[i]]$pars_id[1],]$gamma))
  
  kspline_I <- ((R_t_k + kspline$y[-1]*pars[even.sample.list[[i]]$pars_id[1],]$omega_a)/
                  pars[even.sample.list[[i]]$pars_id[1],]$gamma)
  
  # kspline_I <- ((R_t_k +
  #                  kspline$y[-1] *
  #                  (pars[even.sample.list[[set]]$pars_id[1],]$omega_a+with(pars[set,], b_func(t=seq(-27,364,1), g=g, s=s, psi=psi, f=f)))  )/
  #                 pars[even.sample.list[[set]]$pars_id[1],]$gamma)
  
  # plot(df$t, df$R_a_sample/n_trials, ylim=c(-0.2,1))
  # lines(kspline$x,kspline$y, col="red")
  # lines(kspline$x[-1],kspline_I, col="blue")

  t_spline <- kspline$x[-1]
  spline_df <- data.frame(t_spline, kspline_I)
  
  peak_P[i] <- floor(subset(spline_df, kspline_I==max(spline_df$kspline_I))$t_spline[1])
  
  pars_id[i] <- even.sample.list[[i]]$pars_id[1]
  type[i] <- even.sample.list[[i]]$type[1]
  peak_T[i] <- pars[even.sample.list[[i]]$pars_id[1],]$peak_I
}
#gamma=1/10 so subtract 10 to get the prediction
peaks_even_df <- data.frame(peak_P,peak_T,peak_gamma,pars_id,type)

bouts.list <- list()
for(j in 1:10){
  n_trials=20
  n_reps=10
  T_peak_sample <- vector()
  Pred_peak_sample <- vector()
  random_sample <- vector()
  gamma_sample <- vector()
  type <- vector()
  for(i in 1:length(even.sample.list)){
    data=data.tau.list[[even.sample.list[[i]]$pars_id[1]]]
    type[i]=peaks_even_df$type[i]
    Y_random=data[sample(nrow(data), 1, replace = FALSE), ] 
    Y_true <- subset(data, t==pars[even.sample.list[[i]]$pars_id[1],]$peak_I) 
    I_sample_T <- with(Y_true, rbinom(n_reps, n_trials, I/N))
    Y_pred <- subset(data, t==peaks_even_df$peak_P[i])
    I_sample_P <- with(Y_pred, rbinom(n_reps, n_trials, I/N))
    I_sample_R <- with(Y_random[1,], rbinom(n_reps, n_trials, I/N))
    Y_gamma <- subset(data, t==peaks_even_df$peak_gamma[i])
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
  geom_boxplot(outlier.shape = NA)+scale_fill_brewer(NULL,palette="Paired",labels=c(TeX("True $(\\iota)$"), TeX("Predicted $(\\hat{\\iota})$"),TeX("Null ($1/\\gamma$)"),'Null (random)'))+
  theme(axis.text.x=element_blank())+ theme(legend.position="bottom",legend.text=element_text(size=16))+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),strip.text = element_text(
    size = 16),plot.title = element_text(hjust = 0, size=20)) +ggtitle("")+
  facet_grid(rows = vars(type),labeller = labeller(.rows = row_names))

#10x12 

####Sim Stats Table
columns= c("Sampling Design","Pred Dist. to True", "gamma Dist to True") 
simple_stats = data.frame(matrix(nrow = 4, ncol = length(columns))) 
colnames(simple_stats) = columns
#Daily
daily <- subset(peaks_even_df, type=="daily")
simple_stats$`Sampling Design`[c(1:4)]=c("Daily", "Weekly","Monthly","Bimonthly")
simple_stats[1,2]= with(daily, round(mean(abs(peak_T-peak_P)),2))
simple_stats[1,3]= with(daily, round(mean(abs(peak_T-peak_gamma)),2))

#weekly
weekly <- subset(peaks_even_df, type=="weekly")
simple_stats[2,2]= with(weekly, round(mean(abs(peak_T-peak_P)),2))
simple_stats[2,3]= with(weekly, round(mean(abs(peak_T-peak_gamma)),2))

#monthly
monthly <- subset(peaks_even_df, type=="monthly")
simple_stats[3,2]= with(monthly, round(mean(abs(peak_T-peak_P)),2))
simple_stats[3,3]= with(monthly, round(mean(abs(peak_T-peak_gamma)),2))
#bimonthly
bimonthly <- subset(peaks_even_df, type=="bimonthly")
simple_stats[4,2]= with(bimonthly, round(mean(abs(peak_T-peak_P)),2))
simple_stats[4,3]= with(bimonthly, round(mean(abs(peak_T-peak_gamma)),2))
simple_stats

print(xtable(simple_stats, type='latex'))




