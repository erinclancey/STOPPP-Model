###########################STAN
post.list <- list()
for(i in 1:length(even.sample.list) ){   
  #Setup Stan Estimation
  t_stan=even.sample.list[[i]]$t
  n=rep(n_trials, length(t_stan))
  N=length(t_stan)
  ### Estimate the Seropositive Curve
  data_Ra=list(N=N,x=even.sample.list[[i]]$R_a_sample,n=n,t=t_stan)
  fit_sim_Ra = stan(model_code=ra_func.stan, data=data_Ra, iter=1000)
  theta_draws = rstan::extract(fit_sim_Ra)
  C1_post = theta_draws$C1
  a_post = theta_draws$a
  phi_post = theta_draws$phi
  C2_post = theta_draws$C2
  iter=seq(1,length(C1_post), 1)
  sim_ra_post = data.frame(iter, C1_post, a_post, phi_post, C2_post)
  type <- rep(even.sample.list[[i]]$type[1], 2000)
  pars_id <- rep(even.sample.list[[i]]$pars_id[1], 2000)
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
  omega_a=pars[even.sample.list[[i]]$pars_id[1],]$omega_a
  gamma=pars[even.sample.list[[i]]$pars_id[1],]$gamma
  t <- data.tau.list[[1]]$t
  for(j in 1:length(cred_sim_post[[i]]$iter)){
    I_pred_fit <- (ra_func_deriv(t=t, C1=cred_sim_post[[i]]$C1_post[j], a=cred_sim_post[[i]]$a_post[j], 
                                 phi=cred_sim_post[[i]]$phi_post[j], f=f)
                   + ra_func(t=t, C1=cred_sim_post[[i]]$C1_post[j],
                             a=cred_sim_post[[i]]$a_post[j],
                             phi=cred_sim_post[[i]]$phi_post[j], C2=cred_sim_post[[i]]$C2_post[j], f=f)*(omega_a))/gamma
    Ra_fit <- ra_func(t=t, C1=cred_sim_post[[i]]$C1_post[j], a=cred_sim_post[[i]]$a_post[j], phi=cred_sim_post[[i]]$phi_post[j], 
                      C2=cred_sim_post[[i]]$C2_post[j], f=f)
    pars_id <- rep(even.sample.list[[i]]$pars_id[1], length(I_pred_fit))
    type <- rep(even.sample.list[[i]]$type[1], length(I_pred_fit))
    ID <- j
    pred_list[[j]] <- data.frame(ID,t,I_pred_fit,Ra_fit,pars_id,type)
  }
  pred_frame <- do.call(rbind, pred_list)
  pred_frame.list[[i]] <- pred_frame
  print(i)
}

set=900
#pairs(post.list[[set]][,6:9])
hist(post.list[[set]]$phi_post)
  #ci(post.list[[set]]$a_post, ci=0.95, method = "HDI")
pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=I_pred_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#E69F00")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Proportion Predicted I")+xlab("Day")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))+
  scale_x_continuous(breaks=seq(0,365*5,100))+scale_y_continuous(breaks=seq(-0.05,0.1,0.01))+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = I/N), color="black", alpha=1,linewidth=0.75)
pred_frame.list[[set]] %>%
  ggplot(aes(x=t,y=Ra_fit,group=ID, color=ID))+
  geom_line(aes(group=ID))+
  scale_colour_gradient(low = "black", high = "#0072B2")+
  guides(color="none") +
  theme_minimal()+theme_minimal()+ylab("Proportion Seropositive Fitted")+xlab("Day")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))+
  scale_x_continuous(breaks=seq(0,365*5,100))+scale_y_continuous(breaks=seq(0,0.5,0.02))+
  geom_line(data=data.tau.list[[pred_frame.list[[set]]$pars_id[1]]], aes(x=t,y = R_a/N), color="black", alpha=1,linewidth=0.75)

######
CI_list <- list()
for(i in 1:length(pred_frame.list)){
  peak1_P <- rep(NA, 500)
  peak1_gamma <- rep(NA, 500)
  for(j in 1:500){
    I_pred <- subset(pred_frame.list[[i]], ID==j)
    peak1_P[j] <- subset(I_pred, I_pred$I_pred_fit==max(I_pred$I_pred_fit))$t
    peak1_gamma[j] <- subset(I_pred, I_pred$Ra_fit==max(I_pred$Ra_fit))$t-10
  }
  peak1_ci <- ci(peak1_P, ci=0.95, method = "HDI")
  peak1_ci$Mode1 <- Mode(peak1_P)
  peak1_ci$Mode1gamma <- Mode(peak1_gamma)
  ####
  CI_df <- peak1_ci
  CI_df <- CI_df[,-c(1,6)]
  colnames(CI_df) <- c("CI_low1", "CI_high1", "Mode1","Mode1gamma")
  CI_df$peak1_T <- pars[even.sample.list[[i]]$pars_id[1],]$peak_I
  CI_df$pars_id <- pred_frame.list[[i]]$pars_id[1]
  CI_df$type <- pred_frame.list[[i]]$type[1]
  CI_list[[i]] <- CI_df
  print(i)
}

#peaks are in row numbers not in t
peaks_frame <- do.call(rbind, CI_list)
peaks_frame$in_ci1 <- c(peaks_frame$CI_low1 <= peaks_frame$peak1_T & 
                          peaks_frame$CI_high1 >= peaks_frame$peak1_T )
View(peaks_frame)
#######################
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
    type[i]=peaks_frame$type[i]
    Y_random=data[sample(nrow(data), 1, replace = FALSE), ] 
    Y_true <- subset(data, t==pars[even.sample.list[[i]]$pars_id[1],]$peak_I) 
    I_sample_T <- with(Y_true, rbinom(n_reps, n_trials, I/N))
    Y_pred <- subset(data, t==peaks_frame$Mode1[i])
    I_sample_P <- with(Y_pred, rbinom(n_reps, n_trials, I/N))
    I_sample_R <- with(Y_random[1,], rbinom(n_reps, n_trials, I/N))
    Y_gamma <- subset(data, t==peaks_frame$Mode1gamma[i])
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
    size = 16),plot.title = element_text(hjust = 0, size=20)) +ggtitle("(b) Model Fitting")+
  facet_grid(rows = vars(type),labeller = labeller(.rows = row_names))


####Lik Stats Table
columns= c("Sampling Design","Pred Dist. to True", "gamma Dist to True", "CI size", "Prop. in True") 
lik_stats = data.frame(matrix(nrow = 4, ncol = length(columns))) 
colnames(lik_stats) = columns
#Daily
daily <- subset(peaks_frame, type=="daily")
lik_stats$`Sampling Design`[c(1:4)]=c("Daily", "Weekly","Monthly","Bimonthly")
lik_stats[1,2]= with(daily, round(mean(abs(peak1_T-Mode1)),2))
lik_stats[1,3]= with(daily, round(mean(abs(peak1_T-Mode1gamma)),2))
lik_stats[1,4]= round(mean(abs(c(daily$CI_high1-daily$CI_low1))),2)
lik_stats[1,5]= length(daily$in_ci1[daily$in_ci1==TRUE]) / length(daily$CI_low1)
#Random
#weekly
weekly <- subset(peaks_frame, type=="weekly")
lik_stats[2,2]= with(weekly, round(mean(abs(peak1_T-Mode1)),2))
lik_stats[2,3]= with(weekly, round(mean(abs(peak1_T-Mode1gamma)),2))
lik_stats[2,4]= round(mean(abs(c(weekly$CI_high1-weekly$CI_low1))),2)
lik_stats[2,5]= length(weekly$in_ci1[weekly$in_ci1==TRUE]) / length(weekly$CI_low1)

#monthly
monthly <- subset(peaks_frame, type=="monthly")
lik_stats[3,2]= with(monthly, round(mean(abs(peak1_T-Mode1)),2))
lik_stats[3,3]= with(monthly, round(mean(abs(peak1_T-Mode1gamma)),2))
lik_stats[3,4]= round(mean(abs(c(monthly$CI_high1-monthly$CI_low1))),2)
lik_stats[3,5]= length(monthly$in_ci1[monthly$in_ci1==TRUE]) / length(monthly$CI_low1)

#bimonthly
bimonthly <- subset(peaks_frame, type=="bimonthly")
lik_stats[4,2]= with(bimonthly, round(mean(abs(peak1_T-Mode1)),2))
lik_stats[4,3]= with(bimonthly, round(mean(abs(peak1_T-Mode1gamma)),2))
lik_stats[4,4]= round(mean(abs(c(bimonthly$CI_high1-bimonthly$CI_low1))),2)
lik_stats[4,5]= length(bimonthly$in_ci1[bimonthly$in_ci1==TRUE]) / length(bimonthly$CI_low1)
lik_stats

print(xtable(lik_stats, type='latex'))

