setwd("~/Ebola/Bat Model/STOPPP_Model/H_monstrosus_Cam")
source("Functions_H.monstrosus_Cam.R")
source("Data_Packages_Cam_HM.R")

t=agg_preg$Week
x=agg_preg$lac
n=agg_preg$Total_Sampled
N=length(t)


data=list(N=N,x=x,n=n,t=t)
fit_preg = stan(model_code=b_func.stan.week, data=data, iter=10000)

print(fit_preg, probs=c(0.1, 0.9))
theta_draws = rstan::extract(fit_preg)
g_post=theta_draws$g
s_post=theta_draws$s
psi_post = theta_draws$psi
iter=seq(1,length(psi_post), 1)
b.func_posterior = data.frame(iter, g_post, s_post, psi_post)
plotpostre = ggplot(b.func_posterior, aes(x=psi_post)) +
  geom_histogram(bins=20, color="gray")
plotpostre

#Calculate posterior modes
b.func_posterior.long <- melt(b.func_posterior, id=c("iter"))
modes.preg <- b.func_posterior.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode.preg <- as.vector(modes.preg$mode)
modes.preg

ci(g_post, ci=0.95, method = "HDI")
ci(s_post, ci=0.95, method = "HDI")
ci(psi_post, ci=0.95, method = "HDI")


ggplot(agg_preg, aes(x=Week, y=lac_Prop))+geom_point()+
  scale_x_continuous(breaks=seq(0,51,2), limits=c(0,51)) +xlab("Time (Weeks)")+
  scale_y_continuous(breaks=seq(0,0.7, 0.1),limits=c(0,0.7))+
  theme_minimal()+ylab("Proportion Lactating")+
  ggtitle(expression(paste(italic("H. monstrosus "), "Semi-annual Birth Pulse")))+theme(plot.title = element_text(size=15))+
  stat_function(fun=function(t=Week, g=mode.preg[1],
                             s=mode.preg[2],
                             psi=mode.preg[3]) g*exp(-s*cos(2*pi/52*t-psi)^2), color="black")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))


#weaning is 97.458 days or 13.9 weeks
#I was using 2 weeks for the preg to birth interval
# max longevity is 4462.125 days or 12.225 years

#Death rates by day
#effective lifespan needs to be updated in simulation
N=5000
mu=1/4462
k=(1/2231-1/4462)/N
b=N*k+mu
b
1/2231
#convert b (logevity to weeks to calcualte g_new below)


g_new <- rep(NA, length(b.func_posterior$iter))
# 15 week offset to include birth and weaning
psi_new <- b.func_posterior$psi_post-(15*pi/52)

for (i in 1:length(g_new)){
  b_func_2 <- function(t) exp(-b.func_posterior$s_post[i]*cos(2*pi/52*t-psi_new[i])^2)
  int=integrate(b_func_2, 0, 52)$value
  g_new[i]=0.003137606*52/int
}

b.func_posterior <- cbind(b.func_posterior, g_new, psi_new)

## The period has to be divided by 52 and not 51 otherwise the last day of the year will be identical to the first. We want the cycle to restart on day 52 not end on day 52. 
#Have to integrate over the full cycle which restarts on week 52.
brate <- b_func(t=seq(0,51,1), g=b.func_posterior$g_new[1], s= b.func_posterior$s_post[1] , psi=b.func_posterior$psi_new[1],f=1/52)
#0.003137606, 318.7143 weeks
mean(brate)

