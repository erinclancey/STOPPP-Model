start_time <- Sys.time()
setwd("~/STOPPP Model/STOPPP_Model/Simulations/One_pulse")
source("Functions_Sims and Packages.R")


set.seed(520)
f=1/365

b_func_2 <- function(t) exp(-4*cos(pi*f*t-3)^2)
int=integrate(b_func_2, 0, 365)$value

pars1 <- expand.grid(rep=1:1,combo=2, N0=100000, beta=0.06e-04,omega_a=1/90, omega_i=1/730, gamma=1/10, mu=1/1095,
                     m=2/365, tau=1, g=1/730*365/int, s=4, psi=3, k=(1/730-1/1095)/100000)

pars2 <- expand.grid(rep=1:1,combo=2, N0=100000, beta=0.06e-04,omega_a=1/90, omega_i=1/365, gamma=1/10, mu=1/730,
                     m=2/365, tau=1, g=1/365*365/int, s=4, psi=3,k=(1/365-1/730)/100000)

pars3 <- expand.grid(rep=1:1,combo=3, N0=100000, beta=0.06e-04,omega_a=1/365, omega_i=1/365, gamma=1/10, mu=1/730,
                     m=2/365, tau=1, g=1/365*365/int, s=4, psi=3,k=(1/365-1/730)/100000)


pars <- rbind(pars1, pars2, pars3)

pars$R0 <- rep(NA, length(pars$rep))
pars$b_ave <- rep(NA, length(pars$rep))
for(i in 1:length(pars$rep)){
  pars$R0[i] <- calc_R0(with(pars[i,], c(beta=beta, mu=mu, k=k, gamma=gamma,g=g, s=s, psi=psi, f=f)))
  pars$b_ave[i] <- mean(with(pars[i,], b_func(t=seq(0,364,1), g=g, s=s, psi=psi, f=f)))
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
  data.tau.list[[i]] <- subset(data.tau, data.tau$t>=365*7 & data.tau$t<365*10)
  print(i)
}

for(i in 1:length(data.tau.list)){
  data.tau.list[[i]]$t<- seq(1,365*3,1)
}

##Data from the population has to stay in counts

ggplot(data = data.tau.list[[3]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("Time (Days)")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#E69F00", "#0072B2", "#69b3a2"),
                      labels=c("S","I","Ra","Ri"))+
  #scale_x_continuous(breaks=seq(0,52*5,10))+scale_y_continuous(breaks=seq(0,0.7,0.05))+
  theme(plot.title = element_text(hjust = 1), legend.position = c(0.98, 0.90))

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
length(data.tau.list)
length(pars$g)

# calculate the timeing of the peaks
pars$peak_I <- rep(NA, length(pars$rep))
pars$peak_Ra <- rep(NA, length(pars$rep))
for(i in 1:length(data.tau.list)){
  pars$peak_Ra[i] <-data.tau.list[[i]][which.max(data.tau.list[[i]]$R_a/data.tau.list[[i]]$N),1] 
  pars$peak_I[i] <-data.tau.list[[i]][which.max(data.tau.list[[i]]$I/data.tau.list[[i]]$N),1]
}

end_time <- Sys.time()
end_time - start_time 
   