data {
  int<lower=1> n; // number of years
  int<lower=1> l; // number of lakes
  
  matrix[n,l] y; // log-transformed area covered
  matrix[n,l] sigmaJ; // error due to months past regular survey time
  matrix[n,l] H; // herbicides applied
  
  real x0[l]; // initial population estimate
  real ymis[l]; // average value for missing values
  
  real sr_mean; // mean prior for growth rate
  real sr_sd; // sd prior for growth rate
  real obs_mean; // mean prior for observational error
  real obs_sd; // sd prior for observational error
  real proc_mean; // mean prior for process error
  real proc_sd; // sd prior for process error
}

parameters {
  matrix[n,l] x; // estimated true values
  matrix[n,l] y_surv; // y if it were measured during the same month every year
  real bH; // herbicide effect
  real rmean; // overall growth rate
  real<lower=0> sigma_proc; // process error
  real<lower=0> sigma_r; // growth rate error
  real<lower=0> sigma_obs; // observation error
  vector[l] br; // intrinsic growth rate (random effects) vector
}

model {
  // priors
  bH ~ normal(0, 1);
  rmean ~ normal(0, 1);
  sigma_proc ~ student_t(3, proc_mean, proc_sd);
  sigma_r ~ student_t(3, sr_mean, sr_sd);
  sigma_obs ~ student_t(3, obs_mean, obs_sd);
  
  // data model
  br ~ normal(rmean, sigma_r);
  
  for(i in 1:l){
    x[1,i] ~ normal(x0[i], sigma_proc);
    
    for(t in 2:n){
      x[t,i] ~ normal(x[t-1,i] + br[i] - bH * H[t,i], sigma_proc);
    }
  }
  
  // process model
  for(i in 1:l){
      for(t in 1:n){
        if(y[t,i] != -99){
          y[t,i] ~ normal(y_surv[t,i], sigmaJ[t,i]);
          y_surv[t,i] ~ normal(x[t,i], sigma_obs);
          } else {
            y_surv[t,i] ~ normal(ymis[i], 1);
          }
      }
    }

