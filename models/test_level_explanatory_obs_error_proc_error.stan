data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
}

parameters {
  // level
  vector[n] mu;
  // explanatory
  real beta;
  // observation error
  real<lower=0> sigma_irreg;
  // process error
  real<lower=0> sigma_level;
}

transformed parameters {
  vector[n] yhat;
  yhat = mu + beta * x;
}

model {
  for(t in 2:n)
    mu[t] ~ normal(mu[t-1], sigma_level);
    
  y ~ normal(yhat, sigma_irreg);
}
