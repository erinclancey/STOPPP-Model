###########################STAN
post.list <- list()
for(i in 1:length(cluster.sample.list) ){   
  #Setup Stan Estimation
  t_stan=cluster.sample.list[[i]]$t
  n=rep(n_trials, length(t_stan))
  N=length(t_stan)
  ### Estimate the Seropositive Curve
  data_Ra=list(N=N,x=cluster.sample.list[[i]]$R_a_sample,n=n,t=t_stan)
  fit_sim_Ra = stan(model_code=ra_func.stan, data=data_Ra, iter=1000)
  theta_draws = rstan::extract(fit_sim_Ra)
  C1_post = theta_draws$C1
  a_post = theta_draws$a
  phi_post = theta_draws$phi
  C2_post = theta_draws$C2
  iter=seq(1,length(C1_post), 1)
  sim_ra_post = data.frame(iter, C1_post, a_post, phi_post, C2_post)
  type <- rep(cluster.sample.list[[i]]$type[1], 2000)
  pars_id <- rep(cluster.sample.list[[i]]$pars_id[1], 2000)
  sim_post <- cbind(sim_ra_post,type,pars_id)
  post.list[[i]] <- sim_post
  print(i)
}    

cred_sim_post <- list()
for(i in 1:length(post.list) ){
  C1_ci <- ci(post.list[[i]]$C1_post, ci=0.95, method = "HDI")
  a_ci <- ci(post.list[[i]]$a_post, ci=0.95, method = "HDI")
  phi_ci <- ci(post.list[[i]]$phi_post, ci=0.95, method = "HDI")
  C2_ci <- ci(post.list[[i]]$C2_post, ci=0.95, method = "HDI")
  cred_sim_post[[i]] <- subset(post.list[[i]], 
                               C1_ci$CI_low<C1_post & C1_post<C1_ci$CI_high &
                                 a_ci$CI_low<a_post & a_post<a_ci$CI_high &
                                 phi_ci$CI_low<phi_post & phi_post<phi_ci$CI_high &
                                 C2_ci$CI_low<C2_post & C2_post<C2_ci$CI_high )
  cred_sim_post[[i]] <- cred_sim_post[[i]][sample(nrow(cred_sim_post[[i]]), 500, replace = TRUE), ] 
}

pred_frame.list <- list()
for(i in 1:length(cred_sim_post) ){
  pred_list <- list()
  omega_a=pars[cluster.sample.list[[i]]$pars_id[1],]$omega_a
  gamma=pars[cluster.sample.list[[i]]$pars_id[1],]$gamma
  t <- data.tau.list[[1]]$t
  for(j in 1:length(cred_sim_post[[i]]$iter)){
    I_pred_fit <- (ra_func_deriv(t=t, C1=cred_sim_post[[i]]$C1_post[j], a=cred_sim_post[[i]]$a_post[j], 
                                 phi=cred_sim_post[[i]]$phi_post[j], f=f)
                   + ra_func(t=t, C1=cred_sim_post[[i]]$C1_post[j],
                             a=cred_sim_post[[i]]$a_post[j],
                             phi=cred_sim_post[[i]]$phi_post[j], C2=cred_sim_post[[i]]$C2_post[j], f=f)*(omega_a))/gamma
    Ra_fit <- ra_func(t=t, C1=cred_sim_post[[i]]$C1_post[j], a=cred_sim_post[[i]]$a_post[j], phi=cred_sim_post[[i]]$phi_post[j], 
                      C2=cred_sim_post[[i]]$C2_post[j], f=f)
    pars_id <- rep(cluster.sample.list[[i]]$pars_id[1], length(I_pred_fit))
    type <- rep(cluster.sample.list[[i]]$type[1], length(I_pred_fit))
    ID <- j
    pred_list[[j]] <- data.frame(ID,t,I_pred_fit,Ra_fit,pars_id,type)
  }
  pred_frame <- do.call(rbind, pred_list)
  pred_frame.list[[i]] <- pred_frame
  print(i)
}

# set=101
# #pairs(post.list[[set]][,6:9])
# hist(post.list[[set]]$a_post)
#   #ci(post.list[[set]]$a_post, ci=0.95, method = "HDI")
# pred_frame.list[[set]] %>%
#   ggplot(aes(x=t,y=I_pred_fit,group=ID, color=ID))+
#   geom_line(aes(group=ID))+
#   scale_colour_gradient(low = "black", high = "#E69F00")+
#   guides(color="none") +
#   theme_minimal()+theme_minimal()+ylab("Proportion Predicted I")+xlab("Day")+
#   theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))+
#   scale_x_continuous(breaks=seq(0,365*5,100))+scale_y_continuous(breaks=seq(-0.05,0.1,0.01))+
#   geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = I/N), color="black", alpha=1,linewidth=0.75)
# pred_frame.list[[set]] %>%
#   ggplot(aes(x=t,y=Ra_fit,group=ID, color=ID))+
#   geom_line(aes(group=ID))+
#   scale_colour_gradient(low = "black", high = "#0072B2")+
#   guides(color="none") +
#   theme_minimal()+theme_minimal()+ylab("Proportion Seropositive Fitted")+xlab("Day")+
#   theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))+
#   scale_x_continuous(breaks=seq(0,365*5,100))+scale_y_continuous(breaks=seq(0,0.5,0.02))+
#   geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = R_a/N), color="black", alpha=1,linewidth=0.75)

######
CI_list <- list()
for(i in 1:length(pred_frame.list)){
  peak1_P <- rep(NA, 500)
  peak2_P <- rep(NA, 500)
  peak1_gamma <- rep(NA, 500)
  peak2_gamma <- rep(NA, 500)
  for(j in 1:500){
    I_pred <- subset(pred_frame.list[[i]], ID==j)
    low_pred <- subset(I_pred,  t>0 &  t<170)
    high_pred <- subset(I_pred,  t>170 &  t<360)
    peak1_P[j] <- subset(I_pred, I_pred$I_pred_fit==max(low_pred$I_pred_fit))$t
    peak2_P[j] <- subset(I_pred, I_pred$I_pred_fit==max(high_pred$I_pred_fit))$t
    peak1_gamma[j] <- subset(I_pred, I_pred$Ra_fit==max(low_pred$Ra_fit))$t-10
    peak2_gamma[j] <- subset(I_pred, I_pred$Ra_fit==max(high_pred$Ra_fit))$t-10
  }
  peak1_ci <- ci(peak1_P, ci=0.95, method = "HDI")
  peak1_ci$Mode1 <- Mode(peak1_P)
  peak2_ci <- ci(peak2_P, ci=0.95, method = "HDI")
  peak2_ci$Mode2 <- Mode(peak2_P)
  peak1_ci$Mode1gamma <- Mode(peak1_gamma)
  peak2_ci$Mode2gamma <- Mode(peak2_gamma)
  ####
  CI_df <- cbind(peak1_ci, peak2_ci)
  CI_df <- CI_df[,-c(1,6)]
  colnames(CI_df) <- c("CI_low1", "CI_high1", "Mode1","Mode1gamma","CI_low2", "CI_high2", "Mode2","Mode2gamma")
  CI_df$peak1_T <- pars[cluster.sample.list[[i]]$pars_id[1],]$peak1_T
  CI_df$peak2_T <- pars[cluster.sample.list[[i]]$pars_id[1],]$peak2_T
  CI_df$pars_id <- pred_frame.list[[i]]$pars_id[1]
  CI_df$type <- pred_frame.list[[i]]$type[1]
  CI_df$combo <- pars[cluster.sample.list[[i]]$pars_id[1],]$combo
  CI_list[[i]] <- CI_df
  print(i)
}

#peaks are in row numbers not in t
peaks_frame <- do.call(rbind, CI_list)
peaks_frame$in_ci1 <- c(peaks_frame$CI_low1 <= peaks_frame$peak1_T & 
                          peaks_frame$CI_high1 >= peaks_frame$peak1_T )
peaks_frame$in_ci2 <- c(peaks_frame$CI_low2 <= peaks_frame$peak2_T & 
                          peaks_frame$CI_high2 >= peaks_frame$peak2_T )

bouts.list <- list()
for(j in 1:10){
  n_trials=20
  n_reps=10
  T_peak_sample <- vector()
  Pred_peak_sample <- vector()
  gamma_peak_sample <- vector()
  random_sample <- vector()
  combo <- vector()
  type <- vector()
  for(i in 1:length(cluster.sample.list)){
    data=data.tau.list[[cluster.sample.list[[i]]$pars_id[1]]]
    combo[i]=pars[cluster.sample.list[[i]]$pars_id[1],]$combo
    type[i]=peaks_frame$type[i]
    ############True 
    Y_true1 <- subset(data, t==peaks_frame$peak1_T[i]) 
    Y_true2 <- subset(data, t==peaks_frame$peak2_T[i]) 
    I_sample_T1 <- with(Y_true1, rbinom(n_reps, n_trials, I/N))
    I_sample_T2 <- with(Y_true2, rbinom(n_reps, n_trials, I/N))
    ########Pred
    Y_pred1 <- subset(data, t==peaks_frame$Mode1[i])
    Y_pred2 <- subset(data, t==peaks_frame$Mode2[i]) 
    I_sample_P1 <- with(Y_pred1, rbinom(n_reps, n_trials, I/N))
    I_sample_P2 <- with(Y_pred2, rbinom(n_reps, n_trials, I/N))
    #################gamma
    Y_gamma1 <- subset(data, t==peaks_frame$Mode1gamma[i])
    Y_gamma2 <- subset(data, t==peaks_frame$Mode2gamma[i]) 
    I_sample_gamma1 <- with(Y_gamma1, rbinom(n_reps, n_trials, I/N))
    I_sample_gamma2 <- with(Y_gamma2, rbinom(n_reps, n_trials, I/N))
    ################random
    Y_random=data[sample(nrow(data), 2, replace = FALSE), ]
    I_sample_R1 <- with(Y_random[1,], rbinom(n_reps, n_trials, I/N))
    I_sample_R2 <- with(Y_random[2,], rbinom(n_reps, n_trials, I/N))
    ##############
    I_sample_T <- c(I_sample_T1,I_sample_T2)
    I_sample_P <- c(I_sample_P1,I_sample_P2)
    I_sample_gamma <- c(I_sample_gamma1,I_sample_gamma2)
    I_sample_R <- c(I_sample_R1,I_sample_R2)
    I_sample_T[I_sample_T > 0] <- 1
    I_sample_P[I_sample_P > 0] <- 1
    I_sample_gamma[I_sample_gamma > 0] <- 1
    I_sample_R[I_sample_R > 0] <- 1
    T_peak_sample[i] = mean(I_sample_T)
    Pred_peak_sample[i] = mean(I_sample_P)
    gamma_peak_sample[i] = mean(I_sample_gamma)
    random_sample[i] = mean(I_sample_R)
  }
  bouts.list[[j]] <- data.frame(combo,type,T_peak_sample,Pred_peak_sample,gamma_peak_sample,random_sample)
}
results <- do.call(rbind, bouts.list)
D.long <- melt(results, id=c("combo", "type"))
D.long$type <- factor(D.long$type, levels=c('42_even','42_random','3_day','42_week'))
col_names <- c("1"="Low", "2"="Medium", "3"="High")
row_names <- c("42_even"="Even Days", "42_random"="Random Days", "3_day"="3-Day Clusters", "42_week"="Week Clusters")
ggplot(D.long, aes(x=variable, y=value, fill=variable))+ylab("Proportion Successful Sampling Bouts")+xlab("Sampling Time Point")+
  geom_boxplot(outlier.shape = NA)+scale_fill_brewer(NULL,palette="Dark2",labels=c(TeX("True $(\\iota)$"), TeX("Predicted $\\hat{\\iota}$"),TeX("Null ($1/\\gamma$)"),'Null (random)'))+
  theme(axis.text.x=element_blank())+ theme(legend.position="bottom",legend.text=element_text(size=14))+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),strip.text = element_text(
    size = 18))+ggtitle("(b) Model Fitting")+theme(plot.title = element_text(size=30))+
  facet_grid(cols=vars(combo),rows = vars(type),labeller = labeller(.rows = row_names , .cols = col_names ))

#10x11

####Sim Stats Table
columns= c("Amplitude","Mean Amplitude","Sampling Design","Dist. to True","95% CI Size","Prop. True in CI") 
lik_stats = data.frame(matrix(nrow = 12, ncol = length(columns))) 
colnames(lik_stats) = columns
#############LOW
low <- subset(peaks_frame, combo==1)
#Even
even <- subset(low, type=="42_even")
lik_stats$Amplitude[c(1:4)]="Low"
lik_stats$`Mean Amplitude`[c(1:4)]=round(mean(subset(pars, combo==1)$amp_bar),2)
lik_stats$`Sampling Design`[c(1:4)]=c("Even Days","Random Days","3-Day Clusters","Week Clusters")
peak_true=c(even$peak1_T, even$peak2_T)
peak_pred=c(even$Mode1, even$Mode2)
lik_stats[1,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[1,5]= round(mean(abs(c(even$CI_high1-even$CI_low1,even$CI_high2-even$CI_low2))),2)
lik_stats[1,6]= round((length(even$in_ci1[even$in_ci1==TRUE])+length(even$in_ci2[even$in_ci2==TRUE])) / (2*length(even$CI_low1)),2)
#Random
random <- subset(low, type=="42_random")
peak_true=c(random$peak1_T, random$peak2_T)
peak_pred=c(random$Mode1, random$Mode2)
lik_stats[2,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[2,5]= round(mean(abs(c(random$CI_high1-random$CI_low1,random$CI_high2-random$CI_low2))),2)
lik_stats[2,6]= round((length(random$in_ci1[random$in_ci1==TRUE])+length(random$in_ci2[random$in_ci2==TRUE])) / (2*length(random$CI_low1)),2)
#3-day clusters
three_day <- subset(low, type=="3_day")
peak_true=c(three_day $peak1_T, three_day $peak2_T)
peak_pred=c(three_day $Mode1, three_day $Mode2)
lik_stats[3,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[3,5]= round(mean(abs(c(three_day $CI_high1-three_day $CI_low1,three_day $CI_high2-three_day $CI_low2))),2)
lik_stats[3,6]= round((length(three_day $in_ci1[three_day $in_ci1==TRUE])+length(three_day $in_ci2[three_day $in_ci2==TRUE])) / (2*length(three_day $CI_low1)),2)
#Week clusters
week <- subset(low, type=="42_week")
peak_true=c(week$peak1_T, week$peak2_T)
peak_pred=c(week$Mode1, week$Mode2)
lik_stats[4,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[4,5]= round(mean(abs(c(week$CI_high1-week$CI_low1,week$CI_high2-week$CI_low2))),2)
lik_stats[4,6]= round((length(week$in_ci1[week$in_ci1==TRUE])+length(week$in_ci2[week$in_ci2==TRUE])) / (2*length(week$CI_low1)),2)
#############MED
med<- subset(peaks_frame, combo==2)
#Even
even <- subset(med, type=="42_even")
lik_stats$Amplitude[c(5:8)]="Med"
lik_stats$`Mean Amplitude`[c(5:8)]=round(mean(subset(pars, combo==2)$amp_bar),2)
lik_stats$`Sampling Design`[c(5:8)]=c("Even Days","Random Days","3-Day Clusters","Week Clusters")
peak_true=c(even$peak1_T, even$peak2_T)
peak_pred=c(even$Mode1, even$Mode2)
lik_stats[5,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[5,5]= round(mean(abs(c(even$CI_high1-even$CI_low1,even$CI_high2-even$CI_low2))),2)
lik_stats[5,6]= round((length(even$in_ci1[even$in_ci1==TRUE])+length(even$in_ci2[even$in_ci2==TRUE])) / (2*length(even$CI_low1)),2)
#Random
random <- subset(med, type=="42_random")
peak_true=c(random$peak1_T, random$peak2_T)
peak_pred=c(random$Mode1, random$Mode2)
lik_stats[6,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[6,5]= round(mean(abs(c(random$CI_high1-random$CI_low1,random$CI_high2-random$CI_low2))),2)
lik_stats[6,6]= round((length(random$in_ci1[random$in_ci1==TRUE])+length(random$in_ci2[random$in_ci2==TRUE])) / (2*length(random$CI_low1)),2)
#3-day clusters
three_day <- subset(med, type=="3_day")
peak_true=c(three_day $peak1_T, three_day $peak2_T)
peak_pred=c(three_day $Mode1, three_day $Mode2)
lik_stats[7,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[7,5]= round(mean(abs(c(three_day $CI_high1-three_day $CI_low1,three_day $CI_high2-three_day $CI_low2))),2)
lik_stats[7,6]= round((length(three_day $in_ci1[three_day $in_ci1==TRUE])+length(three_day $in_ci2[three_day $in_ci2==TRUE])) / (2*length(three_day $CI_low1)),2)
#Week clusters
week <- subset(med, type=="42_week")
peak_true=c(week$peak1_T, week$peak2_T)
peak_pred=c(week$Mode1, week$Mode2)
lik_stats[8,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[8,5]= round(mean(abs(c(week$CI_high1-week$CI_low1,week$CI_high2-week$CI_low2))),2)
lik_stats[8,6]= round((length(week$in_ci1[week$in_ci1==TRUE])+length(week$in_ci2[week$in_ci2==TRUE])) / (2*length(week$CI_low1)),2)

#############High
high <- subset(peaks_frame, combo==3)
#Even
even <- subset(high, type=="42_even")
lik_stats$Amplitude[c(9:12)]="High"
lik_stats$`Mean Amplitude`[c(9:12)]=round(mean(subset(pars, combo==3)$amp_bar),2)
lik_stats$`Sampling Design`[c(9:12)]=c("Even Days","Random Days","3-Day Clusters","Week Clusters")
peak_true=c(even$peak1_T, even$peak2_T)
peak_pred=c(even$Mode1, even$Mode2)
lik_stats[9,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[9,5]= round(mean(abs(c(even$CI_high1-even$CI_low1,even$CI_high2-even$CI_low2))),2)
lik_stats[9,6]= round((length(even$in_ci1[even$in_ci1==TRUE])+length(even$in_ci2[even$in_ci2==TRUE])) / (2*length(even$CI_low1)),2)
#Random
random <- subset(high, type=="42_random")
peak_true=c(random$peak1_T, random$peak2_T)
peak_pred=c(random$Mode1, random$Mode2)
lik_stats[10,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[10,5]= round(mean(abs(c(random$CI_high1-random$CI_low1,random$CI_high2-random$CI_low2))),2)
lik_stats[10,6]= round((length(random$in_ci1[random$in_ci1==TRUE])+length(random$in_ci2[random$in_ci2==TRUE])) / (2*length(random$CI_low1)),2)
#3-day clusters
three_day <- subset(high, type=="3_day")
peak_true=c(three_day $peak1_T, three_day $peak2_T)
peak_pred=c(three_day $Mode1, three_day $Mode2)
lik_stats[11,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[11,5]= round(mean(abs(c(three_day $CI_high1-three_day $CI_low1,three_day $CI_high2-three_day $CI_low2))),2)
lik_stats[11,6]= round((length(three_day $in_ci1[three_day $in_ci1==TRUE])+length(three_day $in_ci2[three_day $in_ci2==TRUE])) / (2*length(three_day $CI_low1)),2)
#Week clusters
week <- subset(high, type=="42_week")
peak_true=c(week$peak1_T, week$peak2_T)
peak_pred=c(week$Mode1, week$Mode2)
lik_stats[12,4]= round(mean(abs(peak_true-peak_pred)),2)
lik_stats[12,5]= round(mean(abs(c(week$CI_high1-week$CI_low1,week$CI_high2-week$CI_low2))),2)
lik_stats[12,6]= round((length(week$in_ci1[week$in_ci1==TRUE])+length(week$in_ci2[week$in_ci2==TRUE])) / (2*length(week$CI_low1)),2)
lik_stats


print(xtable(lik_stats, type='latex'))


