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
##########Even Samples
sample_weekly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  week <- seq(1,length(df$t),7)
  df <- df[week,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("weekly",length(df$t))
  sample_weekly[[i]] <- df
}
sample_biweekly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  bi_week <- seq(1,length(df$t),14)
  df <- df[bi_week,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("biweekly",length(df$t))
  sample_biweekly[[i]] <- df
}
sample_monthly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  month <- seq(1,length(df$t),30)
  df <- df[month,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("monthly",length(df$t))
  sample_monthly[[i]] <- df
}
sample_bimonthly <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  bi_month <- seq(1,length(df$t),60)
  df <- df[bi_month,]
  df$pars_id <- rep(i, length(df$t))
  df$type <- rep("bimonthly",length(df$t))
  sample_bimonthly[[i]] <- df
}
even.sample.list <- c(sample_daily,sample_weekly,sample_biweekly,sample_monthly,sample_bimonthly)    
#pars[sample.list[[180]]$pars_id[1],]$gamma
######################42 Day Samples
sample_42_even <- list()
sample_days=42
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  day <- seq(1,length(df$t),9)
  remove <- sample(seq(1,length(day),1),2,replace=FALSE)
  day <- day[-remove]
  df <- df[day,]
  df$pars_id <- rep(i, sample_days)
  df$type <- rep("42_even",sample_days)
  sample_42_even[[i]] <- df
}
sample_42_random <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  df <- df[sample(nrow(df), sample_days), ]
  df <- df[order(df$t),]
  df$pars_id <- rep(i, sample_days)
  df$type <- rep("42_random",sample_days)
  sample_42_random[[i]] <- df
}
sample_42_3_day <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  index <- seq(1:130)
  three <- sample(index,14,replace=FALSE)
  day <- three*3
  days <- c(sapply(day,function(x){seq(x,x+2)}))
  df <- df[days,]
  df <- df[order(df$t),]
  df$pars_id <- rep(i, sample_days)
  df$type <- rep("3_day", sample_days)
  sample_42_3_day[[i]] <- df
}
sample_42_week <- list()
for(i in 1:length(sample_daily) ){
  df <- sample_daily[[i]]
  index <- seq(1:55) #weeks in timeframe
  week <- sample(index,6,replace=FALSE)
  day <- week*7
  days <- c(sapply(day,function(x){seq(x,x+6)}))
  df_w <- df[days,]
  df <- df_w
  df <- df[order(df$t),]
  df$pars_id <- rep(i, sample_days)
  df$type <- rep("42_week", sample_days)
  sample_42_week[[i]] <- df
}
cluster.sample.list <- c(sample_42_even,sample_42_random,sample_42_3_day,sample_42_week)    
#pars[sample.list[[180]]$pars_id[1],]$gamma

