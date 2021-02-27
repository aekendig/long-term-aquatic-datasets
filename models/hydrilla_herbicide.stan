data {
  int<lower=1> n; // number of observations
  //int<lower=1> t; // number of years
  vector[n] y; // observations
  //int<lower=1> intval[n]; // time interval of observations
  //vector[n] j; // proportion of time interval
  //vector[t] h; // herbicide used in year t
  vector[n] h; // herbicide used in year t
  int<lower=0> z0; // initial population
}

parameters {
  //vector[t] z; // estimated true values
  vector[n] z; // estimated true values
  real bH; // coefficient for herbicide effect
  real<lower=0> sigma_obs; // observation error
  real<lower=0> sigma_proc; // process error
}

model {
  // priors
  sigma_obs ~ student_t(3, 0, 1);
  sigma_proc ~ student_t(3, 0, 1);
  bH ~ normal(0, 1);
  
  // data model
  z[1] ~ normal(z0, sigma_proc);
  
  //for(i in 2:t){
  for(i in 2:n){
    z[i] ~ normal(z[i-1] + bH * h[i-1], sigma_proc);
  }
  
  // process model
  for(i in 1:n){
    //y[i] ~ normal(j[i] * z[intval[i]-1] + (1 - j[i]) * z[intval[i]], sigma_obs);
    y[i] ~ normal(z[i], sigma_obs);
  }
}
