data {
  int<lower=1> n; // number of years
  vector[n] y; // observations
  vector[n] sigmaJ; // error due to months past regular survey time
  real x0; // initial population
}

parameters {
  vector[n] x; // estimated true values
  vector[n] y_surv; // y if it were measured during the same month every year
  real br; // intrinsic growth rate
  real<lower=0> sigma_proc; // process error
  real<lower=0> sigma_obs; // observation error
}

model {
  // priors
  sigma_proc ~ student_t(3, 0, 1);
  sigma_obs ~ student_t(3, 0, 1);
  br ~ normal(0, 1);
  
  // data model
  x[1] ~ normal(x0, sigma_proc);
  
  for(i in 2:n){
    x[i] ~ normal(x[i-1] + br, sigma_proc);
  }
  
  // process model
  for(i in 1:n){
    if(y[i] != -99){
      y[i] ~ normal(y_surv[i], sigmaJ[i]);
      y_surv[i] ~ normal(x[i], sigma_obs);
    } else {
      y_surv[i] ~ normal(0, 1);
    }
  }
}
