start_time <- Sys.time()
setwd("~/Ebola/Bat Model/STOPPP_Model/Simulations")
source("Functions_Sims and Packages.R")


set.seed(420)
f=2/365

b_func_2 <- function(t) exp(-6*cos(pi*f*t-1.25)^2)
int=integrate(b_func_2, 0, 365)$value
pars1 <- expand.grid(rep=1:100,combo=1, N0=100000, beta=0.06e-04,omega_a=1/90, omega_i=1/1095, gamma=1/10, mu=1/1825,
                     m=2/365, tau=1, g=1/1095*365/int, s=6, psi=1.25, k=(1/1095-1/1825)/100000)
b_func_2 <- function(t) exp(-4*cos(pi*f*t-1.25)^2)
int=integrate(b_func_2, 0, 365)$value
pars2 <- expand.grid(rep=1:100,combo=2, N0=100000, beta=0.06e-04,omega_a=1/90, omega_i=1/730, gamma=1/10, mu=1/1095,
                     m=2/365, tau=1, g=1/730*365/int, s=4, psi=1.25, k=(1/730-1/1095)/100000)

b_func_2 <- function(t) exp(-4*cos(pi*f*t-2)^2)
int=integrate(b_func_2, 0, 365)$value
pars3 <- expand.grid(rep=1:100,combo=3, N0=100000, beta=0.06e-04,omega_a=1/90, omega_i=1/365, gamma=1/10, mu=1/730,
                     m=2/365, tau=1, g=1/365*365/int, s=4, psi=2,k=(1/365-1/730)/100000)


pars <- rbind(pars1,pars2,pars3)

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
  data.tau.list[[i]] <- subset(data.tau, data.tau$t>365*9-29 & data.tau$t<365*10)
  print(i)
}

for(i in 1:length(data.tau.list)){
  data.tau.list[[i]]$t<- c(seq(-28,-1),seq(0,364,1))
}

##Data from the population has to stay in counts

ggplot(data = data.tau.list[[1]], aes(x = t))+ theme_minimal()+ ylab("Proportion")+xlab("Time (Days)")+
  geom_line(aes(y = S/N, color="S"), alpha=1, linewidth=0.75) +
  geom_line(aes(y = I/N, color="I"), alpha=1,linewidth=0.75) +
  geom_line(aes(y = R_a/N, color="R_a"), alpha=1, linewidth=0.75)+
  geom_line(aes(y = R_i/N, color="R_i"), alpha=1, linewidth=0.75)+
  scale_colour_manual("", breaks = c("S","I","R_a","R_i"),values = c("black", "#E69F00", "#0072B2", "#69b3a2"),
                      labels=c("S","I","Ra","Ri"))+
  #scale_x_continuous(breaks=seq(0,52*5,10))+scale_y_continuous(breaks=seq(0,0.7,0.05))+
  theme(plot.title = element_text(hjust = 1), legend.position = c(0.98, 0.90))
prev <- vector()
seroprev <- vector()
for(i in 1:length(data.tau.list)){
  prev[i] <- mean(data.tau.list[[i]]$I/data.tau.list[[i]]$N)
  seroprev[i] <- mean(data.tau.list[[i]]$R_a/data.tau.list[[i]]$N)
}
pars$prev <- prev
pars$seroprev <- seroprev

num_prev_zero <- as.numeric(rownames(pars[pars$prev==0,]))
if(length(num_prev_zero)!=0){data.tau.list <- data.tau.list[-num_prev_zero] }
pars <- pars[-pars$prev!=0,]
length(data.tau.list)
length(pars$g)

pars$amp_bar <- rep(NA,length(pars$g))
for(i in 1:length(data.tau.list)){
  Ra_prop <- data.tau.list[[i]]$R_a/data.tau.list[[i]]$N
  peak_matrix <- findpeaks(Ra_prop, minpeakheight = 0.0001,
                           minpeakdistance = 20, npeaks = 2, nups = 0, ndowns = 0)
  if(nrow(peak_matrix)<2 || is.null(peak_matrix)==TRUE){pars$amp_bar[i]=NA}else{
    peak_matrix <- peak_matrix[order(peak_matrix[,2],decreasing=FALSE),]
    if(abs(peak_matrix[1,2]-peak_matrix[1,4]) <= abs(peak_matrix[2,2]-peak_matrix[2,3])){
      trough <- peak_matrix[2,3]+1}else{trough <- peak_matrix[1,4]+1}
    #Returns a matrix where each row represents one peak found. 
    #The first column gives the height, the second the position/index where the maximum is reached, the third and forth the indices of where the peak begins and ends — in the sense of where the pattern starts and ends.
    amp_1 <- peak_matrix[1,1] - Ra_prop[trough]
    amp_2 <- peak_matrix[2,1] - Ra_prop[trough]
    pars$amp_bar[i] <- (amp_1+amp_2)/2
    print(i)}
}
pars_na <- as.numeric(rownames(pars[is.na(pars$amp_bar)==TRUE,]))
pars <- na.omit((pars))
amp_bar_neg <- as.numeric(rownames(pars[pars$amp_bar<0,]))
remove <- c(amp_bar_neg, pars_na)
if(length(remove)!=0){data.tau.list <- data.tau.list[-remove] }
pars <- pars[-pars$amp_bar<0,]
length(data.tau.list)
length(pars$g)

row.names(pars) <- 1:nrow(pars)
#View(pars)

# low <- subset(pars, combo==1)
# mean(low$seroprev)
# mean(low$prev)
# mean(low$amp_bar)
# med<- subset(pars, combo==2)
# mean(med$seroprev)
# mean(med$prev)
# hi <- subset(pars, combo==3)
# mean(hi$seroprev)
# mean(hi$prev)

for(i in 1:length(data.tau.list)){
low_true <- subset(data.tau.list[[i]], t>0 & t<170)
high_true <- subset(data.tau.list[[i]], t>190 & t<340)
pars$peak1_T[i] <- subset(low_true, I/N==max(low_true$I/low_true$N))$t
pars$peak2_T[i]  <- subset(high_true, I/N==max(high_true$I/high_true$N))$t
}

end_time <- Sys.time()
end_time - start_time 