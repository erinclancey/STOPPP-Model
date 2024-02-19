setwd("~/Ebola/Bat Model/STOPPP_Model/E_Helvum_Cam")
source("Data_Packages_Cam_EH.R")
source("Functions_E_helvum.R")

# Fit E_helvum lactation data with STAN
t=agg_lac$week
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

ggplot(agg_lac, aes(x=week, y=lac_Prop))+geom_point(color="black",fill="black", size = 2, shape = 21)+
  scale_x_continuous(breaks=seq(0,51,2), limits=c(0,51))+
  scale_y_continuous(breaks=seq(0,0.8, 0.1))+
  ggtitle(expression(paste(italic("E. helvum "), "Annual Birth Pulse")))+theme(plot.title = element_text(size=15))+
  theme_minimal()+ylab("Proportion Lactating")+xlab("Time (Weeks)")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  stat_function(fun=function(t=Week, g=mode.preg[1],
                             s=mode.preg[2],
                             psi=mode.preg[3]) g*exp(-s*cos(1*pi/52*t-(psi-pi/52))^2), color="black")


