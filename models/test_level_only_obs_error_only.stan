data {
  int<lower=1> n;
  vector[n] y;
}
parameters {
  // level
  real mu;
  // observation error
  real<lower=0> sigma;
}
model {
  y ~ normal(mu, sigma);
}
