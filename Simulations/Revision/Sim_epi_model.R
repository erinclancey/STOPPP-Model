# R code to simulate Model 1
start_time <- Sys.time()
source("Functions and Packages_SIMS.R")

f=1/365

b_func_2 <- function(t) exp(-4*cos(pi*f*t-2.5)^2)
int=integrate(b_func_2, 1, 365)$value

pars1 <- expand.grid(rep=1:100,combo=2, N0=10000, beta=0.06e-03,omega_a=1/90, omega_i=1/730, gamma=1/10, mu=1/1095,
                     m=2/365, tau=1, g=1/730*365/int, s=4, psi=2.5, k=(1/730-1/1095)/10000)

pars2 <- expand.grid(rep=1:100,combo=2, N0=10000, beta=0.06e-03,omega_a=1/90, omega_i=1/365, gamma=1/10, mu=1/730,
                     m=2/365, tau=1, g=1/365*365/int, s=4, psi=2.5,k=(1/365-1/730)/10000)

pars3 <- expand.grid(rep=1:100,combo=3, N0=10000, beta=0.06e-03,omega_a=1/365, omega_i=1/365, gamma=1/10, mu=1/730,
                     m=2/365, tau=1, g=1/365*365/int, s=4, psi=2.5,k=(1/365-1/730)/10000)

pars <- rbind(pars1, pars2, pars3)

pars$R0 <- rep(NA, length(pars$rep))
pars$b_ave <- rep(NA, length(pars$rep))
for(i in 1:length(pars$rep)){
  pars$R0[i] <- calc_R0(with(pars[i,], c(beta=beta, mu=mu, k=k, gamma=gamma,g=g, s=s, psi=psi, f=f)))
  pars$b_ave[i] <- mean(with(pars[i,], b_func(t=seq(1,365,1), g=g, s=s, psi=psi, f=f)))
}

S_start = floor(with(pars,(b_ave + gamma)/beta))
I_start = floor(with(pars, ((b_ave*k - b_ave*beta + k*gamma + beta*mu)*(b_ave + omega_a)*(b_ave + omega_i)*(gamma*mu + mu^2 + gamma*omega_a+ mu*omega_a + gamma*omega_i + mu*omega_i + omega_a*omega_i))/
                       (-(b_ave^2*k*beta*gamma*mu) - b_ave*k*beta*gamma^2*mu - b_ave^2*k*beta*mu^2 - b_ave*k*beta*gamma*mu^2- b_ave^2*k*beta*gamma*omega_a - b_ave*k*beta*gamma^2*omega_a- b_ave^2*k*beta*mu*omega_a - 2*b_ave*k*beta*gamma*mu*omega_a - 
                          k*beta*gamma^2*mu*omega_a - b_ave*k*beta*mu^2*omega_a - k*beta*gamma*mu^2*omega_a - b_ave*k*beta*gamma*omega_a^2 - k*beta*gamma^2*omega_a^2 - b_ave*k*beta*mu*omega_a^2 - k*beta*gamma*mu*omega_a^2 - 
                          b_ave^2*k*beta*gamma*omega_i - b_ave*k*beta*gamma^2*omega_i - b_ave^2*k*beta*mu*omega_i - 2*b_ave*k*beta*gamma*mu*omega_i - k*beta*gamma^2*mu*omega_i - b_ave*k*beta*mu^2*omega_i - k*beta*gamma*mu^2*omega_i - b_ave^2*k*beta*omega_a*omega_i - 
                          3*b_ave*k*beta*gamma*omega_a*omega_i - 2*k*beta*gamma^2*omega_a*omega_i - 2*b_ave*k*beta*mu*omega_a*omega_i - 3*k*beta*gamma*mu*omega_a*omega_i - k*beta*mu^2*omega_a*omega_i - b_ave*k*beta*omega_a^2*omega_i - 2*k*beta*gamma*omega_a^2*omega_i - k*beta*mu*omega_a^2*omega_i - 
                          b_ave*k*beta*gamma*omega_i^2 - k*beta*gamma^2*omega_i^2 - b_ave*k*beta*mu*omega_i^2 - k*beta*gamma*mu*omega_i^2 - b_ave*k*beta*omega_a*omega_i^2 - 2*k*beta*gamma*omega_a*omega_i^2 - k*beta*mu*omega_a*omega_i^2 - 
                          k*beta*omega_a^2*omega_i^2)))
R_a_start = floor(with(pars, -((gamma*(b_ave*k - b_ave*beta + k*gamma + beta*mu)*(b_ave + omega_i))/(k*beta*(b_ave^2 + b_ave*gamma + b_ave*omega_a + gamma*omega_a + b_ave*omega_i + gamma*omega_i + omega_a*omega_i)))))
R_i_start = floor(with(pars,-((gamma*(b_ave*k - b_ave*beta + k*gamma + beta*mu)*omega_a)/(k*beta*(b_ave^2 + b_ave*gamma + b_ave*omega_a + gamma*omega_a + b_ave*omega_i + gamma*omega_i + omega_a*omega_i)))))
b_start = (S_start+I_start+R_a_start+R_i_start)*with(pars, b_func(t=0, g=g, s=s, psi=psi, f=f))

### Run Tau Leaping
data.tau.list.full <- list()
data.tau.list <- list()
for(i in 1:nrow(pars)){
  xstart <- with(pars, c(t=0,S=S_start[i],I=I_start[i],R_a=R_a_start[i],R_i=R_i_start[i],births=b_start[i]))
  data.tau <- Bat.sim.tau(xstart,with(pars, c(beta=beta[i], mu=mu[i], k=k[i], omega_a=omega_a[i], omega_i=omega_i[i], gamma=gamma[i], m=m[i], tau=tau[i],g=g[i], s=s[i], psi=psi[i], f=f)))
  data.tau$N <- with(data.tau , (S+I+R_a+R_i))
  data.tau.list.full[[i]] <-data.tau 
  data.tau.list[[i]] <- subset(data.tau, data.tau$t>365 & data.tau$t<365*3+1)
  print(i)
}

for(i in 1:length(data.tau.list)){
  data.tau.list[[i]]$t<- seq(1,365*2,1)
}

#calculate prevalence
pars$prev <- rep(NA, length(pars$rep))
pars$seroprev <- rep(NA, length(pars$rep))
for(i in 1:length(data.tau.list)){
  pars$prev[i] <- mean(data.tau.list[[i]]$I/data.tau.list[[i]]$N)
  pars$seroprev[i] <- mean(data.tau.list[[i]]$R_a/data.tau.list[[i]]$N)
}

#remove any datasets where no disease outbreak occurs
num_prev_zero <- as.numeric(rownames(pars[pars$prev==0,]))
if(length(num_prev_zero)!=0){data.tau.list <- data.tau.list[-num_prev_zero] }
pars <- pars[-pars$prev!=0,]

data.tau.list.1_365     <- lapply(data.tau.list, subset, t <= 365)
data.tau.list.366_730   <- lapply(data.tau.list, subset, t >= 366)
for(i in 1:length(data.tau.list.366_730)){
  data.tau.list.366_730[[i]]$t<- seq(1,365,1)
}


# calculate the timeing of the peaks
pars$peak_I_1 <- rep(NA, length(pars$rep))
pars$peak_Ra_1 <- rep(NA, length(pars$rep))
pars$peak_I_2 <- rep(NA, length(pars$rep))
pars$peak_Ra_2 <- rep(NA, length(pars$rep))
for(i in 1:length(data.tau.list)){
  pars$peak_Ra_1[i] <- data.tau.list.1_365[[i]][which.max(data.tau.list.1_365[[i]]$R_a/data.tau.list.1_365[[i]]$N),1] 
  pars$peak_I_1[i] <- data.tau.list.1_365[[i]][which.max(data.tau.list.1_365[[i]]$I/data.tau.list.1_365[[i]]$N),1]
  pars$peak_Ra_2[i] <- data.tau.list.366_730 [[i]][which.max(data.tau.list.366_730 [[i]]$R_a/data.tau.list.366_730 [[i]]$N),1] 
  pars$peak_I_2[i] <- data.tau.list.366_730 [[i]][which.max(data.tau.list.366_730 [[i]]$I/data.tau.list.366_730 [[i]]$N),1]
  
}

end_time <- Sys.time()
end_time - start_time 



