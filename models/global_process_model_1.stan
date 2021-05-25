// Global process model without herbicide

data {
  int<lower=0> TT; // length of time series
  vector[TT] y; // observations
  vector[TT] x; // independent variable
}

parameters {
  real beta0; // growth rate
  real beta1; // -growth rate/K
  real<lower=0> sigma; // error
}

model {
  y ~ normal(beta0 + beta1 * x, sigma);
}
