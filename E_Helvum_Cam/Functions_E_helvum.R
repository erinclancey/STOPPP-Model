library (rstan)

## R Mathematical Functions

b_func <- function(t,g,s,psi,f) g*exp(-s*cos(pi*f*t-psi)^2)
ra_func <- function(t,C1,a,phi,f, C2) -C1*exp(-a*cos(pi*f*t-phi)^2)+C2
ra_func_deriv <- function(t,C1,a,phi,f, C2) -2*C1*a*f*pi*exp(-a*cos(pi*f*t-phi)^2)*cos(pi*f*t-phi)*sin(pi*f*t-phi)


## Stan functions

b_func.stan =
  "data {
  int<lower=0> N;   // number of data items
  int<lower=0> n[N];    // number of trials in i
  int<lower=0> x[N];   // success on trial n
  vector[N] t;               // time
}

parameters {
  real<lower=0, upper=1> g;    // g Amplitude 
  real<lower=0, upper=positive_infinity()> s;   // Gaussian flatness
  real<lower=0, upper=pi()> psi;    // phase shift 
}

model {
  x ~ binomial(n,  g*exp (-s* pow (cos (-psi + pi()*t*1/365), 2) ) ); // likelihood
}"

b_func.stan.week =
  "data {
  int<lower=0> N;   // number of data items
  int<lower=0> n[N];    // number of trials in i
  int<lower=0> x[N];   // success on trial n
  vector[N] t;               // time
}

parameters {
  real<lower=0, upper=1> g;    // g Amplitude 
  real<lower=0, upper=positive_infinity()> s;   // Gaussian flatness
  real<lower=-0.5*pi(), upper=0.5*pi()> psi;    // phase shift 
}

model {
  x ~ binomial(n,  g*exp (-s* pow (cos (-psi + pi()*t*1/52), 2) ) ); // likelihood
}"

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
  x ~ binomial(n,  C2 - C1*exp (-a* pow (cos (-phi + pi()*t*1/365), 2) ) ); // likelihood
}"

ra_func.stan.week =
  "data {
  int<lower=0> N;   // number of data items
  int<lower=0> n[N];    // number of trials in i
  int<lower=0> x[N];   // success on trial n
  vector[N] t;               // time
}

parameters {
  real<lower=0, upper=1> C1;    // Peak 
  real<lower=0, upper=1> a ;   // Gaussian flatness
  real<lower=-0.5*pi(), upper=pi()*0.5> phi;    // phase shift 
  real<lower=C1, upper=1+C1/exp(a)> C2;    // verticle shift 
}

model {
  x ~ binomial(n,  C2 - C1*exp (-a* pow (cos (-phi + pi()*t*1/52), 2) ) ); // likelihood
}"
