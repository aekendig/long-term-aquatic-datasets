data {
  int<lower=1> n; // number of years
  vector[n] y; // observations
  vector[n] sigmaJ; // error due to months past regular survey time
  vector[n] H; // herbicides applied
  real x0; // initial population estimate
  real ymis; // average value for missing values
  real br_mean; // mean prior for growth rate
  real br_sd; // sd prior for growth rate
  real obs_mean; // mean prior for growth rate
  real obs_sd; // sd prior for growth rate
  real proc_mean; // mean prior for growth rate
  real proc_sd; // sd prior for growth rate
}

parameters {
  vector[n] x; // estimated true values
  vector[n] y_surv; // y if it were measured during the same month every year
  real br; // intrinsic growth rate
  real bH; // herbicide effect
  real<lower=0> sigma_proc; // process error
  real<lower=0> sigma_obs; // observation error
}

model {
  // priors
  sigma_proc ~ student_t(3, proc_mean, proc_sd);
  sigma_obs ~ student_t(3, obs_mean, obs_sd);
  br ~ normal(br_mean, br_sd);
  bH ~ normal(0, 1);
  
  // data model
  x[1] ~ normal(x0, sigma_proc);
  
  for(i in 2:n){
    x[i] ~ normal(x[i-1] + br - bH * H[i], sigma_proc);
  }
  
  // process model
  for(i in 1:n){
    if(y[i] != -99){
      y[i] ~ normal(y_surv[i], sigmaJ[i]);
      y_surv[i] ~ normal(x[i], sigma_obs);
    } else {
      y_surv[i] ~ normal(ymis, 1);
    }
  }
}
