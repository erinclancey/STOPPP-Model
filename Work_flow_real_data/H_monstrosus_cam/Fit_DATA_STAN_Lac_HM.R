setwd("~/STOPPP Model/STOPPP_Model/H_monstrosus_Cam")
source("Functions_H.monstrosus_Cam.R")
source("Data_Packages_Cam_HM.R")

t=agg_lac$Week
x=agg_lac$lac
n=agg_lac$Total_Sampled
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

#Calculate posterior modes
b.func_posterior.long <- melt(b.func_posterior, id=c("iter"))
modes.preg <- b.func_posterior.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode.preg <- as.vector(modes.preg$mode)

ggplot(agg_lac, aes(x=Week, y=lac_Prop))+geom_point()+
  scale_x_continuous(breaks=seq(0,51,2), limits=c(0,51)) +xlab("Time (Weeks)")+
  scale_y_continuous(breaks=seq(0,0.7, 0.1),limits=c(0,0.7))+
  theme_minimal()+ylab("Proportion Lactating")+
  ggtitle(expression(paste(italic("H. monstrosus "), "Semi-annual Birth Pulse")))+theme(plot.title = element_text(size=15))+
  stat_function(fun=function(t=Week, g=mode.preg[1],
                             s=mode.preg[2],
                             psi=mode.preg[3]) g*exp(-s*cos(2*pi/52*t-psi)^2), color="black")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))