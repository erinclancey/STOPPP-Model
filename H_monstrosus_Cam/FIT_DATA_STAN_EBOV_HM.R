setwd("~/Ebola/Bat Model/STOPPP_Model/H_monstrosus_Cam")
source("Functions_H.monstrosus_Cam.R")
source("Data_Packages_Cam_HM.R")

t=EBVagg_cases$Week
x=EBVagg_cases$EBVPos
n=EBVagg_cases$Total_Sampled
N=length(t)
data=list(N=N,x=x,n=n,t=t)

#Now that I am fitting the new version of the Peel Function I am not getting any amplitude

fit_EBV = stan(model_code=ra_func.stan.week, data=data, iter=10000)

print(fit_EBV, probs=c(0.1, 0.9))
theta_draws = rstan::extract(fit_EBV)
a_post = theta_draws$a
C1_post = theta_draws$C1
C2_post = theta_draws$C2
phi_post = theta_draws$phi
iter=seq(1,length(a_post), 1)
ra_posterior = data.frame(iter, a_post, C1_post, C2_post, phi_post)
plotpostre = ggplot(ra_posterior, aes(x=phi_post)) +
  geom_histogram(bins=20, color="gray")
plotpostre



#Calculate posterior modes
ra_posterior.long <- melt(ra_posterior, id=c("iter"))
modes <- ra_posterior.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode <- as.vector(modes$mode)
modes

ci(a_post, ci=0.95, method = "HDI")
ci(C1_post, ci=0.95, method = "HDI")
ci(C2_post, ci=0.95, method = "HDI")
ci(phi_post, ci=0.95, method = "HDI")

ggplot(EBVagg_cases, aes(x=Week, y=EBVPos_prop))+geom_point()+scale_x_continuous(breaks=seq(0,51,2))+
  theme_minimal()+ylab("Proportion Seropositive")+
  stat_function(fun=function(t=Week, a=mode[1],
                             C1=mode[2],
                             C2=mode[3],
                             phi=mode[4]) -C1*exp(-a*cos(2*pi/52*t-phi)^2)+C2, color="blue")


# omega_a=1/12.85714
# gamma=1/1.428571
# f=2/52
# kspline <- ksmooth(EBVagg_cases$Week, EBVagg_cases$EBVPos_prop, "box", 
#                    bandwidth = 10, n.points=length(EBVagg_cases$Week) )
# plot(EBVagg_cases$Week,EBVagg_cases$EBVPos_prop)
# lines(EBVagg_cases$Week,kspline$y)
# R_t_k <- diff(kspline$y)
# I_new_k <- (R_t_k + kspline$y[-1]*omega_a)/gamma
# t_spline <- EBVagg_cases$Week[-1]
# lines(t_spline,I_new_k, col="red")
# 
# Ispline <- ksmooth(t_spline, I_new_k , "box", 
#                    bandwidth = 20, n.points=length(t_spline) )
# lines(t_spline, Ispline$y, col="blue")
# plot(t_spline, I_new_k , col="blue")
