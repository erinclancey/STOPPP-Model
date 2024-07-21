############cluster
peak1_P <- rep(NA, length(cluster.sample.list))
peak2_P <- rep(NA, length(cluster.sample.list))
peak1_T <- rep(NA, length(cluster.sample.list))
peak2_T <- rep(NA, length(cluster.sample.list))
pars_id <- rep(NA, length(cluster.sample.list))
peak1_gamma <- rep(NA, length(cluster.sample.list))
peak2_gamma <- rep(NA, length(cluster.sample.list))
type <- rep(NA, length(cluster.sample.list))
combo <- rep(NA, length(cluster.sample.list))

for(i in 1:length(cluster.sample.list)){
  df <- cluster.sample.list[[i]]
  kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 50)
  if(any(is.na(kspline$y))==TRUE){
    kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 90)
  }else{kspline=kspline}
  R_t_k <- diff(kspline$y)
  rA_df <- data.frame(kspline$x, kspline$y)
  low_rA <- subset(rA_df,  kspline.x>0 &  kspline.x<170)
  high_rA <- subset(rA_df,  kspline.x>170 &  kspline.x<360)
  peak1_gamma[i] <- floor(subset(low_rA, kspline.y==max(low_rA$kspline.y))$kspline.x[1]-(1/pars[cluster.sample.list[[i]]$pars_id[1],]$gamma))
  peak2_gamma[i] <- floor(subset(high_rA, kspline.y==max(high_rA$kspline.y))$kspline.x[1]-(1/pars[cluster.sample.list[[i]]$pars_id[1],]$gamma))
  kspline_I <- ((R_t_k + kspline$y[-1]*pars[cluster.sample.list[[i]]$pars_id[1],]$omega_a)/pars[cluster.sample.list[[i]]$pars_id[1],]$gamma)
  t_spline <- kspline$x[-1]
  spline_df <- data.frame(t_spline, kspline_I)
  low_pred <- subset(spline_df,  t_spline>0 &  t_spline<170)
  high_pred <- subset(spline_df,  t_spline>170 &  t_spline<360)
  peak1_P[i] <- floor(subset(low_pred, kspline_I==max(low_pred$kspline_I))$t_spline[1])
  peak2_P[i] <- floor(subset( high_pred, kspline_I==max(high_pred$kspline_I))$t_spline[1])
  pars_id[i] <- cluster.sample.list[[i]]$pars_id[1]
  type[i] <- cluster.sample.list[[i]]$type[1]
  combo[i] <-pars[cluster.sample.list[[i]]$pars_id[1],]$combo[1]
  peak1_T[i] <- pars[cluster.sample.list[[i]]$pars_id[1],]$peak1_T
  peak2_T[i] <- pars[cluster.sample.list[[i]]$pars_id[1],]$peak2_T
}
#gamma=1/10 so subtract 10 to get the prediction
peaks_cluster_df <- data.frame(peak1_P,peak2_P,peak1_T,peak2_T,peak1_gamma,peak2_gamma,pars_id,type,combo)
#Warning occurs because not all weeks of the year are picked in week clusters. 
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
    Y_gamma1 <- subset(data, t==peaks_cluster_df$peak1_gamma[i])
    Y_gamma2 <- subset(data, t==peaks_cluster_df$peak2_gamma[i]) 
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
# 200 rows removed due to NA values. NA values can occur in cases when cluster sampling does not occur in the first half of the year randomly. 
D.long <- melt(results, id=c("combo", "type"))
D.long$type <- factor(D.long$type, levels=c('42_even','42_random','3_day','42_week'))
col_names <- c("1"="Low", "2"="Medium", "3"="High")
row_names <- c("42_even"="Even Days", "42_random"="Random Days", "3_day"="3-Day Clusters", "42_week"="Week Clusters")
ggplot(D.long, aes(x=variable, y=value, fill=variable))+ylab("Proportion Successful Sampling Bouts")+xlab("Sampling Time Point")+
  geom_boxplot(outlier.shape = NA)+scale_fill_brewer(NULL,palette="Dark2",labels=c(TeX("True $(\\iota)$"), TeX("Predicted $\\hat{\\iota}$"),TeX("Null ($1/\\gamma$)"),'Null (random)'))+
  theme(axis.text.x=element_blank())+ theme(legend.position="bottom",legend.text=element_text(size=14))+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),strip.text = element_text(
    size = 18))+ggtitle("(a) Interpolation")+theme(plot.title = element_text(size=30))+
  facet_grid(cols=vars(combo),rows = vars(type),labeller = labeller(.rows = row_names , .cols = col_names ))

#10x11


####Sim Stats Table
columns= c("Amplitude","Mean Amplitude","Sampling Design","Dist. to True") 
cluster_stats = data.frame(matrix(nrow = 12, ncol = length(columns))) 
colnames(cluster_stats) = columns
peaks_cluster_df <- na.omit(peaks_cluster_df)
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




