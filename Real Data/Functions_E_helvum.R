
library (rstan)

## R Mathematical Functions to fit E. helvum data

b_func <- function(t,g,s,psi,f) g*exp(-s*cos(pi*f*t-psi)^2)
## Stan functions

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
