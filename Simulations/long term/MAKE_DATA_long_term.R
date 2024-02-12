
setwd("~/Ebola/Bat Model/long term")
source("Functions_Packages_long_term.R")


f=1/365

# b_func_2 <- function(t) exp(-1*cos(pi*f*t-2)^2)
# int=integrate(b_func_2, 0, 365)$value
# pars <- expand.grid(rep=1:1,combo=3, N0=1000, beta=0.06e-2,omega_a=1/90, omega_i=1/300, gamma=1/10, mu=1/730,
#                      m=1/100, tau=1, g=1/365*365/int, s=1, psi=2,k=(1/365-1/730)/1000)

b_func_2 <- function(t) exp(-2*cos(pi*f*t-1)^2)
int=integrate(b_func_2, 0, 365)$value
pars <- expand.grid(rep=1:1,combo=2, N0=1000, beta=0.08e-02,omega_a=1/70, omega_i=1/300, gamma=1/7, mu=1/1095,
                     m=2/100, tau=1, g=1/730*365/int, s=2, psi=1, k=(1/730-1/1095)/1000)


pars$R0 <- rep(NA, length(pars$rep))
pars$b_ave <- rep(NA, length(pars$rep))
for(i in 1:length(pars$rep)){
  pars$R0[i] <- calc_R0(with(pars[i,], c(beta=beta, mu=mu, k=k, gamma=gamma,g=g, s=s, psi=psi, f=f)))
  pars$b_ave[i] <- mean(with(pars[i,], b_func(t=seq(0,364,1), g=g, s=s, psi=psi, f=f)))
}

pars

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
  data.tau.list[[i]] <- data.tau
  print(i)
}


##Data from the population has to stay in counts

ggplot(data = data.tau.list[[1]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("Time (Days)")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#E69F00", "#0072B2", "#69b3a2"),
                      labels=c("S","I","Ra","Ri"))+
  scale_x_continuous(breaks=seq(0,1095,365))+scale_y_continuous(breaks=seq(0,0.7,0.1))+
  theme(plot.title = element_text(hjust = 1), legend.position = c(0.98, 0.90))
