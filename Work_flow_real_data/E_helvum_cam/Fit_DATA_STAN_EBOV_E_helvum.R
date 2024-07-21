setwd("~/STOPPP-Model/Work_flow_real_data/E_helvum_cam")
source("Data_Packages_Cam_EH.R")
source("Functions_E_helvum.R")

# Fit E_helvum serology data with STAN
t=EBOVagg_cases$Week
x=EBOVagg_cases$EBOVPos
n=EBOVagg_cases$Total_Sampled
N=length(t)
data=list(N=N,x=x,n=n,t=t)

fit_EBOV = stan(model_code=ra_func.stan.week, data=data, iter=10000)

print(fit_EBOV, probs=c(0.1, 0.9))
theta_draws = rstan::extract(fit_EBOV)
a_post = theta_draws$a
C1_post = theta_draws$C1
C2_post = theta_draws$C2
phi_post = theta_draws$phi
iter=seq(1,length(a_post), 1)
ra_posterior = data.frame(iter, a_post, C1_post, C2_post, phi_post)
