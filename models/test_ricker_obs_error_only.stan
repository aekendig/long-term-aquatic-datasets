data {
  int<lower=1> n;
  vector[n] y;
  real mu0;
}

parameters {
  // growth rate
  real b0;
  // -growth rate/K
  //real<upper=0> b1;
  real<lower=0> K;
  // observation error
  real<lower=0> sigma_irreg;
  // process error
  //real<lower=0> sigma_level;
}

transformed parameters {
  // log abundance
  vector[n] mu;
  // back-transformed abundance
  vector[n] yhat;
  
  // population growth
  mu[1] = mu0;
  for(t in 2:n){
    //mu[t] = mu[t-1] + b0 + b1 * exp(mu[t-1]);
    mu[t] = mu[t-1] + b0 - (b0/K) * exp(mu[t-1]);
  }
  // back-transformed abundance
  yhat = exp(mu);
}

model {
  // prior distributions
  sigma_irreg ~ normal(0, 1);
  b0 ~ normal(0, 10);
  K ~ normal(20, 10);
  //sigma_level ~ normal(0, 1);
  
  // process model
  //for(t in 2:n){
    //mu[t] ~ normal(mu[t-1] + b0 + b1 * exp(mu[t-1]), sigma_level);
  //}
  
  // observation model
  for(t in 1:n){
    y[t] ~ normal(yhat, sigma_irreg);
  }
  
}
