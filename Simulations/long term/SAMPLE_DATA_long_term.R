n_trials=20
############Daily Samples
sample_daily <- list()
for(i in 1:length(data.tau.list) ){
  t <- data.tau.list[[i]]$t
  I_sample <- vector()
  R_a_sample <- vector()
  for(k in 1:length(t)) {
    I_sample[k] <- with(data.tau.list[[i]], rbinom(1, n_trials, I[k]/N[k]))
    R_a_sample[k] <- with(data.tau.list[[i]], rbinom(1, n_trials, R_a[k]/N[k]))
  }
  df <- data.frame(t, I_sample, R_a_sample)
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("daily",length(t))
  sample_daily[[i]] <- df
}
##########Weekly Samples
sample_weekly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  week <- seq(1,length(df$t),7)
  df <- df[week,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("weekly",length(df$t))
  sample_weekly[[i]] <- df

}
##########Bi-Weekly Samples
sample_biweekly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  bi_week <- seq(1,length(df$t),14)
  df <- df[bi_week,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("biweekly",length(df$t))
  sample_biweekly[[i]] <- df
}
new_sample <- c(sample_daily,sample_weekly,sample_biweekly)    

omega_a=pars$omega_a
gamma=pars$gamma

set=1
df <- new_sample[[set]]
kspline <- ksmooth(df$t, df$R_a_sample/n_trials, "normal", bandwidth = 40, n.points=length(df$t) )
R_t_k <- diff(kspline$y)
kspline_I <- (R_t_k + kspline$y[-1]*omega_a)/gamma
t_spline <- df$t[-1]
spline_df <- data.frame(t_spline, kspline_I)
par(mfrow = c(1,2), mar = c(5,5,1,1))
plot(df$t,df$R_a_sample/n_trials,xaxt = "n",yaxt = "n",ylim=c(0,0.7),
     xlab = "Time (Days)",
     ylab = "Prop. Seropos.",
     pch = 20,
     col= "grey", cex.lab=1.5, cex.axis=1.5)
axis(1, at = seq(0, max(df$t), by = 365))
axis(2, at = seq(0, 0.8, by = 0.1))
#title("(a) Low Amplitude Daily Sample", adj=0, cex.main=1.75)
lines(df$t,kspline$y, lwd=2,col="black")
lines(data.tau.list[[1]]$t[-1], data.tau.list[[1]]$R_a[-1]/data.tau.list[[1]]$N[-1],lwd=2,col="#0072B2")
legend("topright", c(expression(r[a]),expression(hat(r)[a])),
       col = c("#0072B2","black"), lty=1, cex=1.5)

plot(t_spline,kspline_I, ylim = c(0, 0.3),xaxt = "n",yaxt = "n",
     xlab = "Time (Days)",
     ylab = "Prop. Infected",
     type="l",lwd=2,
     col= "black",cex.lab=1.5, cex.axis=1.5)
axis(1, at = seq(0, max(df$t), by = 365))
axis(2, at = seq(0, 0.3, by = 0.05))
lines(data.tau.list[[1]]$t[-1], data.tau.list[[1]]$I[-1]/data.tau.list[[1]]$N[-1], col="#E69F00",lwd=2)
legend("topright", c(expression(iota),expression(hat(iota))),
       col = c("#E69F00","black"), lty=1:1, cex=1.5)

