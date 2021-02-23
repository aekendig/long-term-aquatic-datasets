data {
  int<lower=1> n;
  vector[n] y;
}

parameters {
  // level
  vector[n] mu;
  // slope
  vector[n-1] v;
  // process error slope
  real<lower=0> sigma_drift;
  // observation error
  real<lower=0> sigma_irreg;
  // process error level
  real<lower=0> sigma_level;
}

transformed parameters {
  vector[n] yhat;
  yhat = mu;
}

model {
  // slope model
  v[1] ~ normal(0, sigma_drift);
  for(t in 2:n-1)
    v[t] ~ normal(v[t-1], sigma_drift);
    
  // level model
  mu[1] ~ normal(y[1], sigma_level);
  for(t in 2:n)
    mu[t] ~ normal(mu[t-1] + v[t-1], sigma_level);

  // data model
  y ~ normal(yhat, sigma_irreg);

  // priors for errors
  sigma_drift ~ student_t(4, 0, 1);
  sigma_irreg ~ student_t(4, 0, 1);
  sigma_level ~ student_t(4, 0, 1);
}

