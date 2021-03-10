data {
  int<lower=1> n; // number of years
  int<lower=1> l; // number of lakes
  
  matrix[n,l] y; // log-transformed area covered
  matrix[n,l] sigmaJ; // error due to months past regular survey time
  
  real x0[l]; // initial population estimate
  real ymis[l]; // average value for missing values
}

parameters {
  matrix[n,l] x; // estimated true values
  matrix[n,l] y_surv; // y if it were measured during the same month every year
  real<lower=0> sigma_proc; // process error
  real<lower=0> sigma_obs; // observation error
  real br_hat; // overall intrinsic growth rate
  real<lower=0> sigma_r; // growth rate error
  real eta[l]; //re-parameterization of growth rate
}

transformed parameters {
  real br[l];
  for (i in 1:l){
    br[i] = br_hat + sigma_r * eta[i];
  }
}

model {
  // priors
  sigma_proc ~ student_t(3, 0, 1);
  sigma_obs ~ student_t(3, 0, 1);
  br_hat ~ normal(0, 1);
  sigma_r ~ student_t(3, 0, 1);
  eta ~ normal(0, 1);
  
  // data model
  for(i in 1:l){
    x[1,i] ~ normal(x0[i], sigma_proc);
    
    for(t in 2:n){
      x[t,i] ~ normal(x[t-1,i] + br[i], sigma_proc);
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
}
