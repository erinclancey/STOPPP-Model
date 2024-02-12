#Packages
library("lubridate") 
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(latex2exp)
library(zoo)
library(minpack.lm)
library(rstan)
library(reshape2)
library(coda)
library(MCMCglmm)
library(bayestestR)
library(pracma)
library(reshape)
library(functional)
library(matrixcalc)
library(lhs)
library(xtable) 
set.seed(420)

## R Mathematical Functions
b_func <- function(t,g,s,psi,f) g*exp(-s*cos(pi*f*t-psi)^2)
ra_func <- function(t,C1,a,phi,f, C2) -C1*exp(-a*cos(pi*f*t-phi)^2)+C2
ra_func_deriv <- function(t,C1,a,phi,f, C2) -4*C1*a*f*pi*exp(-a*cos(pi*f*t-phi)^2)*cos(pi*f*t-phi)*sin(pi*f*t-phi)

##Function to Calculate R0
calc_R0 <- function(params){
  beta <- params["beta"]
  mu <- params["mu"]
  k <- params["k"]
  omega_a <- params["omega_a"]
  omega_i <- params["omega_i"]
  gamma <- params["gamma"]
  m <- params["m"]
  tau <- params["tau"]
  g <- params["g"]
  s <- params["s"]
  psi <- params["psi"]
  f<- params["f"]
  t <- seq(0,364,1)
  b <- b_func(t,g,s,psi,f)
  R_0 <- beta/(gamma+b) * ((b-mu)/k)
  mean(R_0)
}

## Function for the Tau Leaping Algorithm and simulation
Bat.onestep.tau <- function (x, params) {
  t <- x[1]
  S <- x[2]
  I <- x[3]
  R_a <- x[4]
  R_i <- x[5]
  births <- x[6]
  N <- S+I+R_a+R_i
  beta <- params["beta"]
  mu <- params["mu"]
  k <- params["k"]
  omega_a <- params["omega_a"]
  omega_i <- params["omega_i"]
  gamma <- params["gamma"]
  m <- params["m"]
  tau <- params["tau"]
  g <- params["g"]
  s <- params["s"]
  psi <- params["psi"]
  f <- params["f"]
  b <- b_func(t,g,s,psi,f)
  
  events <- c(
    birth = rpois(1,b*N),
    death_S = rpois(1,S*(mu + k*N)),
    death_I = rpois(1,I*(mu + k*N)),
    death_R_a = rpois(1,R_a*(mu + k*N)),
    death_R_i = rpois(1,R_i*(mu + k*N)),
    infection = rpois(1,beta*S*I),
    migrant = rpois(1, m),
    recovery_a = rpois(1,gamma*I),
    recovery_i = rpois(1,omega_a*R_a),
    wane_i = rpois(1,omega_i*R_i)
  )
  
  transitions <- data.frame( 
    S_birth = c(1,0,0,0)*events[1],
    S_death = c(-1,0,0,0)*events[2],
    I_death = c(0,-1,0,0)*events[3],
    R_a_death = c(0,0,-1,0)*events[4],
    R_i_death = c(0,0,0,-1)*events[5],
    infect = c(-1,1,0,0)*events[6],
    migrate = c(0,1,0,0)*events[7],
    recov_a = c(0,-1,1,0)*events[8],
    recov_i = c(0,0,-1,1)*events[9],
    wane_i = c(1,0,0,-1)*events[10]
  )
  x[6]=events[1]
  c(x[1:5] + c(tau, rowSums(transitions)), x[6])
}

Bat.sim.tau <- function (x, params, maxstep = 1000000) {
  simdata <- array(dim=c(maxstep+1,6))
  colnames(simdata) <- names(x)
  simdata[1,] <- x
  for(i in 1:100000000){
    simdata[i,] <- x <- Bat.onestep.tau(x,params)
    if(x["S"]<0){x["S"]=0}
    if(x["I"]<0){x["I"]=0}
    if(x["R_a"]<0){x["R_a"]=0}
    if(x["R_i"]<0){x["R_i"]=0}
    if( (i>maxstep) | (x["t"] > 365*10)){break}
  }
  as.data.frame(simdata[1:i,])
}

## Stan functions

ra_func.stan =
  "data {
  int<lower=0> N;   // number of data items
  int<lower=0> n[N];    // number of trials in i
  int<lower=0> x[N];   // success on trial n
  vector[N] t;               // time
}
parameters {
  real<lower=0, upper=1> C1;    // Peak 
  real<lower=0, upper=1> a ;   // Gaussian flatness
  real<lower=0, upper=pi()> phi;    // phase shift 
  real<lower=C1, upper=1+C1/exp(a)> C2;    // verticle shift 
}
model {
  x ~ binomial(n,  C2 - C1*exp (-a* pow (cos (-phi + pi()*t*2/365), 2) ) ); // likelihood
}"





