data {
  int<lower=1> n; // number of years
  vector[n] y; // observations
  vector[n] J; // months past regular survey time
  int<lower=0> x0; // initial population
}

parameters {
  vector[n] x; // estimated true values
  real br; // intrinsic growth rate
  real<lower=0> sigma_proc; // process error
  real<lower=0> sigma_obs; // observation error
  real<lower=0> sigma_j; // error due to irrugular suvey
}

model {
  // priors
  sigma_proc ~ student_t(3, 0, 1);
  sigma_obs ~ student_t(3, 0, 1);
  sigma_j ~ student_t(3, 0, 1);
  br ~ normal(0, 1);
  
  // data model
  x[1] ~ normal(x0, sigma_proc);
  
  for(i in 2:n){
    x[i] ~ normal(x[i-1] + br, sigma_proc);
  }
  
  // process model
  for(i in 1:n){
    y[i] ~ normal(x[i], sigma_obs + sigma_j * J);
  }
}
