########################Clustered Sampling
peak1_P <- rep(NA, length(cluster.sample.list))
peak2_P <- rep(NA, length(cluster.sample.list))
peak1_T <- rep(NA, length(cluster.sample.list))
peak2_T <- rep(NA, length(cluster.sample.list))
pars_id <- rep(NA, length(cluster.sample.list))
type <- rep(NA, length(cluster.sample.list))
combo <- rep(NA, length(cluster.sample.list))
for(i in 1:length(cluster.sample.list)){
  df <- cluster.sample.list[[i]]
  kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 90, n.points=length(df$t) )
  if(any(is.na(kspline$y))==TRUE){
    kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 110, n.points=length(df$t) )
  }else{kspline=kspline}
  R_t_k <- diff(kspline$y)
  kspline_I <- ((R_t_k + kspline$y[-1]*pars[cluster.sample.list[[i]]$pars_id[1],]$omega_a)/pars[cluster.sample.list[[i]]$pars_id[1],]$gamma)
  t_spline <- df$t[-1]
  spline_df <- data.frame(t_spline, kspline_I)
  low_pred <- subset(spline_df,  t_spline>0 &  t_spline<170)
  high_pred <- subset(spline_df,  t_spline>170 &  t_spline<360)
  if(length(low_pred$t_spline)!=0){peak1_P[i] <- subset(spline_df, kspline_I==max(low_pred$kspline_I))$t_spline}else{
    peak1_P[i] <- sample(seq(0,170,1),1)
  }
  if(length(high_pred$t_spline)!=0){peak2_P[i] <- subset(spline_df, kspline_I==max(high_pred$kspline_I))$t_spline}else{
    peak2_P[i] <- sample(seq(170,360,1),1)
  }
  pars_id[i] <- cluster.sample.list[[i]]$pars_id[1]
  type[i] <- cluster.sample.list[[i]]$type[1]
  combo[i] <-pars[cluster.sample.list[[i]]$pars_id[1],]$combo[1]
  peak1_T[i] <- pars[cluster.sample.list[[i]]$pars_id[1],]$peak1_T
  peak2_T[i] <- pars[cluster.sample.list[[i]]$pars_id[1],]$peak2_T
}
peaks_cluster_df <- data.frame(peak1_P,peak2_P,peak1_T,peak2_T,pars_id,type,combo)


bouts.list <- list()
for(j in 1:10){
  n_trials=20
  n_reps=10
  T_peak_sample <- vector()
  Pred_peak_sample <- vector()
  random_sample <- vector()
  combo <- vector()
  type <- vector()
  for(i in 1:length(cluster.sample.list)){
    data=data.tau.list[[cluster.sample.list[[i]]$pars_id[1]]]
    combo[i]=pars[cluster.sample.list[[i]]$pars_id[1],]$combo
    type[i]=peaks_cluster_df$type[i]
    Y_random=data[sample(nrow(data), 2, replace = FALSE), ] 
    Y_true1 <- subset(data, t==pars[cluster.sample.list[[i]]$pars_id[1],]$peak1_T) 
    Y_true2 <- subset(data, t==pars[cluster.sample.list[[i]]$pars_id[1],]$peak2_T) 
    I_sample_T1 <- with(Y_true1, rbinom(n_reps, n_trials, I/N))
    I_sample_T2 <- with(Y_true2, rbinom(n_reps, n_trials, I/N))
    Y_pred1 <- subset(data, t==peaks_cluster_df$peak1_P[i])
    Y_pred2 <- subset(data, t==peaks_cluster_df$peak2_P[i]) 
    I_sample_P1 <- with(Y_pred1, rbinom(n_reps, n_trials, I/N))
    I_sample_P2 <- with(Y_pred2, rbinom(n_reps, n_trials, I/N))
    I_sample_R1 <- with(Y_random[1,], rbinom(n_reps, n_trials, I/N))
    I_sample_R2 <- with(Y_random[2,], rbinom(n_reps, n_trials, I/N))
    I_sample_T <- c(I_sample_T1,I_sample_T2)
    I_sample_P <- c(I_sample_P1,I_sample_P2)
    I_sample_R <- c(I_sample_R1,I_sample_R2)
    I_sample_T[I_sample_T > 0] <- 1
    I_sample_P[I_sample_P > 0] <- 1
    I_sample_R[I_sample_R > 0] <- 1
    T_peak_sample[i] = mean(I_sample_T)
    Pred_peak_sample[i] = mean(I_sample_P)
    random_sample[i] = mean(I_sample_R)
  }
  bouts.list[[j]] <- data.frame(combo,type,T_peak_sample,Pred_peak_sample,random_sample)
}
results <- do.call(rbind, bouts.list)
D.long <- melt(results, id=c("combo", "type"))
D.long$type <- factor(D.long$type, levels=c('42_even','42_random','3_day','42_week'))
col_names <- c("1"="Low", "2"="Medium", "3"="High")
row_names <- c("42_even"="Even Days", "42_random"="Random Days", "3_day"="3-Day Clusters", "42_week"="Week Clusters")
ggplot(D.long, aes(x=variable, y=value, fill=variable))+ylab("Proportion Successful Sampling Bouts")+xlab("Sampling Time Point")+
  geom_boxplot(outlier.shape = NA)+scale_fill_brewer(NULL,palette="Dark2",labels=c('True', 'Predicted','Random'))+
  theme(axis.text.x=element_blank())+ theme(legend.position="bottom",legend.text=element_text(size=13))+
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),strip.text = element_text(
    size = 13))+
  ggtitle("(a) Interpolation")+theme(plot.title = element_text(size=18))+
  facet_grid(cols=vars(combo),rows = vars(type),labeller = labeller(.rows = row_names , .cols = col_names ))


####Sim Stats Table
columns= c("Amplitude","Mean Amplitude","Sampling Design","Dist. to True") 
cluster_stats = data.frame(matrix(nrow = 12, ncol = length(columns))) 
colnames(cluster_stats) = columns
#############LOW
low <- subset(peaks_cluster_df, combo==1)
#Even
even <- subset(low, type=="42_even")
cluster_stats$Amplitude[c(1:4)]="Low"
cluster_stats$`Mean Amplitude`[c(1:4)]=round(mean(subset(pars, combo==1)$amp_bar),2)
cluster_stats$`Sampling Design`[c(1:4)]=c("Even Days","Random Days","3-Day Clusters","Week Clusters")
peak_true=c(even$peak1_T, even$peak2_T)
peak_pred=c(even$peak1_P, even$peak2_P)
cluster_stats[1,4]= round(mean(abs(peak_true-peak_pred)),2)
#Random
random <- subset(low, type=="42_random")
peak_true=c(random$peak1_T, random$peak2_T)
peak_pred=c(random$peak1_P, random$peak2_P)
cluster_stats[2,4]= round(mean(abs(peak_true-peak_pred)),2)
#3-day clusters
three_day <- subset(low, type=="3_day")
peak_true=c(three_day $peak1_T, three_day $peak2_T)
peak_pred=c(three_day $peak1_P, three_day $peak2_P)
cluster_stats[3,4]= round(mean(abs(peak_true-peak_pred)),2)
#Week clusters
week <- subset(low, type=="42_week")
peak_true=c(week$peak1_T, week$peak2_T)
peak_pred=c(week$peak1_P, week$peak2_P)
cluster_stats[4,4]= round(mean(abs(peak_true-peak_pred)),2)
#############MED
med<- subset(peaks_cluster_df, combo==2)
#Even
even <- subset(med, type=="42_even")
cluster_stats$Amplitude[c(5:8)]="Med"
cluster_stats$`Mean Amplitude`[c(5:8)]=round(mean(subset(pars, combo==2)$amp_bar),2)
cluster_stats$`Sampling Design`[c(5:8)]=c("Even Days","Random Days","3-Day Clusters","Week Clusters")
peak_true=c(even$peak1_T, even$peak2_T)
peak_pred=c(even$peak1_P, even$peak2_P)
cluster_stats[5,4]= round(mean(abs(peak_true-peak_pred)),2)
#Random
random <- subset(med, type=="42_random")
peak_true=c(random$peak1_T, random$peak2_T)
peak_pred=c(random$peak1_P, random$peak2_P)
cluster_stats[6,4]= round(mean(abs(peak_true-peak_pred)),2)
#3-day clusters
three_day <- subset(med, type=="3_day")
peak_true=c(three_day $peak1_T, three_day $peak2_T)
peak_pred=c(three_day $peak1_P, three_day $peak2_P)
cluster_stats[7,4]= round(mean(abs(peak_true-peak_pred)),2)
#Week clusters
week <- subset(med, type=="42_week")
peak_true=c(week$peak1_T, week$peak2_T)
peak_pred=c(week$peak1_P, week$peak2_P)
cluster_stats[8,4]= round(mean(abs(peak_true-peak_pred)),2)
#############High
high <- subset(peaks_cluster_df, combo==3)
#Even
even <- subset(high, type=="42_even")
cluster_stats$Amplitude[c(9:12)]="High"
cluster_stats$`Mean Amplitude`[c(9:12)]=round(mean(subset(pars, combo==3)$amp_bar),2)
cluster_stats$`Sampling Design`[c(9:12)]=c("Even Days","Random Days","3-Day Clusters","Week Clusters")
peak_true=c(even$peak1_T, even$peak2_T)
peak_pred=c(even$peak1_P, even$peak2_P)
cluster_stats[9,4]= round(mean(abs(peak_true-peak_pred)),2)
#Random
random <- subset(high, type=="42_random")
peak_true=c(random$peak1_T, random$peak2_T)
peak_pred=c(random$peak1_P, random$peak2_P)
cluster_stats[10,4]= round(mean(abs(peak_true-peak_pred)),2)
#3-day clusters
three_day <- subset(high, type=="3_day")
peak_true=c(three_day $peak1_T, three_day $peak2_T)
peak_pred=c(three_day $peak1_P, three_day $peak2_P)
cluster_stats[11,4]= round(mean(abs(peak_true-peak_pred)),2)
#Week clusters
week <- subset(high, type=="42_week")
peak_true=c(week$peak1_T, week$peak2_T)
peak_pred=c(week$peak1_P, week$peak2_P)
cluster_stats[12,4]= round(mean(abs(peak_true-peak_pred)),2)

cluster_stats
print(xtable(cluster_stats, type='latex'))



