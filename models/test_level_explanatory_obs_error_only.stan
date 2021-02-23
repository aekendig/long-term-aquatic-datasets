data {
  int<lower=1> n;
  vector[n] y;
  vector[n] x;
}

parameters {
  // level
  real mu;
  // explanatory
  real beta;
  // observation error
  real<lower=0> sigma_irreg;
}

transformed parameters {
  vector[n] yhat;
  yhat = mu + beta * x;
}

model {
  y ~ normal(yhat, sigma_irreg);
}
