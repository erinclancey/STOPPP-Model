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
sample.list <- c(sample_daily,sample_weekly,sample_monthly,sample_bimonthly)    

