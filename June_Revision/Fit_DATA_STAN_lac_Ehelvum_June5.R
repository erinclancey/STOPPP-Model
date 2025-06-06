setwd("~/STOPPP Model/STOPPP_Model/E_Helvum_Cam")
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

#Calculate posterior modes and plot the fitted function
b.func_posterior.long <- melt(b.func_posterior, id=c("iter"))
modes.preg <- b.func_posterior.long %>% group_by(variable) %>% summarise(mode = posterior.mode(mcmc(value), adjust=1))
mode.preg <- as.vector(modes.preg$mode)

ggplot(agg_lac, aes(x=week, y=lac_Prop))+geom_point(color="black",fill="black", size = 2, shape = 21)+
  scale_x_continuous(breaks=seq(0,51,2), limits=c(0,51))+
  scale_y_continuous(breaks=seq(0,0.8, 0.1))+
  ggtitle(expression(paste(italic("E. helvum "), "Annual Birth Pulse")))
+theme(plot.title = element_text(size=15))+
  theme_minimal()+ylab("Proportion Lactating")+xlab("Time (Weeks)")+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  stat_function(fun=function(t=Week, g=mode.preg[1],
                             s=mode.preg[2],
                             psi=mode.preg[3]) g*exp(-s*cos(1*pi/52*t-(psi-pi/52))^2), color="black")


b_func_helvum <- function(t) exp(-mode.preg[2]*cos(pi*1/52*seq(0,51,1)-mode.preg[3])^2)
int=integrate(b_func_2, 0, 52)$value

omega_a=1/12.85714
gamma=1/1.428571
f=1/52
g=1/(7969/7)*52 / int

EBOVagg_cases$b <- b_func(t=EBOVagg_cases$Week, g=g, s=mode.preg[2], psi=mode.preg[3], f=1/52)


df <- EBOVagg_cases
kspline <- ksmooth(df$Week, df$EBOVPos_prop, "normal", bandwidth = 17, n.points=500)
R_t_k <- diff(kspline$y)
kspline_I <- (R_t_k +
                 kspline$y[-1] *
                 (omega_a + b_func(t=kspline$x[-1], g=g, s=mode.preg[2], psi=mode.preg[3], f=1/52) ) )/ gamma

t_spline <- kspline$x[-1]

# Create a new data frame for plotting
df_plot <- data.frame(
  Time = df$Week,
  Seropositive = df$EBOVPos_prop
)

spline_df_plot <- data.frame(
  Time = t_spline,
  Prop_Sero=kspline$y[-1],
  Prop_Infected = kspline_I
)

quantile <- subset(spline_df_plot, kspline_I>quantile(kspline_I, probs = 0.90))
min <- floor(min(quantile$Time))
max <- ceiling(max(quantile$Time))
peak <- subset(spline_df_plot, Prop_Infected==max(Prop_Infected))
peak_Ra <- subset(spline_df_plot, Prop_Sero==max(Prop_Sero))
##########################################
A <- ggplot() +
  geom_point(data=df_plot, aes(x = Time, y = Seropositive, color = "Predicted"), size = 2, alpha=0.5, shape=16) +
  geom_line(data = spline_df_plot, 
            aes(x = Time, y = Prop_Sero, color = "Weekly Data"), size = 1.5) +
  
  labs(title = expression(paste(italic("E. helvum "), "Seroprevalence")), x = "", y = "Prop. Seropositive") +
  theme_minimal() +
  theme(text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  theme(legend.position = "right")+
  scale_y_continuous(limits = c(0, 1)) + scale_x_continuous(breaks=seq(0,51,2), limits=c(0,51))+
  theme(plot.title = element_text(hjust = 0)) +
 scale_color_manual(values = c("Weekly Data"="#0072B2","Predicted" = "#0072B2"),
                     labels = c(TeX("Data"),TeX("Predicted $(\\hat{r}_a)$")) ) +
  guides(color = guide_legend(title = ""))

A
######################################
B <- ggplot() +
  geom_line(data = spline_df_plot, aes(x = Time, y = Prop_Infected, color = "Predicted"), size = 1.5)+
  labs(title = expression(paste(italic("E. helvum "), "Predicted Prevalence")), x = "Week of Year", y = "Prop. Infected") +
  theme_minimal() +
  theme(text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13))+
  theme(legend.position = "right")+
  scale_y_continuous(limits = c(0, 0.25)) +  scale_x_continuous(breaks=seq(0,51,2), limits=c(0,51))+
  scale_color_manual(values = c("Predicted" = "#CC79A7"),
                     labels = c(TeX("Predicted $(\\hat{\\iota})$"))) +
  guides(color = guide_legend(title = ""))+
  geom_ribbon(data = subset(spline_df_plot, Time >= min & Time <= max), 
              aes(x = Time, ymin = 0, ymax = Prop_Infected), 
              fill = "#CC79A7", alpha = 0.5)+
  geom_vline(xintercept = floor(peak$Time), linetype=3, size=1)

plot <- plot_grid(A,B, ncol = 1, nrow = 2, rel_heights=c(1,1))
plot

# E. helvum is a relatively long-lived species, with a mean life expectancy estimates across the colonies ranging from 2.3 to 6.8 years, as can be seen after converting the model results to annual survival probabilities
                #Can survival analyses detect hunting pressure in a highly connected species? Lessons from straw-coloured fruit bats
