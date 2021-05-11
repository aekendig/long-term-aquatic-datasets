// Multivariate autoregressive model
// https://nwfsc-timeseries.github.io/atsa-labs/sec-marss-fitting-with-stan.html
// modified for unique growth rates by group
// modified standard error priors

data {
  int<lower=0> TT; // length of time series
  int<lower=0> N; // number of groups; rows of y
  int<lower=0> n_pos; // number of non-NA values in y
  int<lower=0> col_indx_pos[n_pos]; // col index of non-NA vals
  int<lower=0> row_indx_pos[n_pos]; // row index of non-NA vals
  vector[n_pos] y; // non-NA observations
  vector[N] y0; // first observation for each group
}

parameters {
  vector[N] x0; // initial states
  vector[N] u; // population growth rates
  real<lower=0> sd_q; // process error
  real<lower=0> sd_r[N]; // observation error - one for each group
  vector[N] pro_dev_tilde[TT]; // refed as pro_dev_tilde[TT,N]; latent process error
}
transformed parameters {
  vector[N] x[TT]; // refed as x[TT,N]; array of TT objects, each a vector of length N
  vector[N] pro_dev[TT]; // refed as pro_dev[TT,N]; process error for each time/group (w_t)
  for(i in 1:N){ // cycle through groups
    pro_dev[1,i] = sd_q*pro_dev_tilde[1,i];
    x[1,i] = x0[i] + u[i] + pro_dev[1,i]; // x at t=1 = x at t=0 + growth + process error
    for(t in 2:TT) { // cycle through time points
      pro_dev[t,i] = sd_q*pro_dev_tilde[t,i];
      x[t,i] = x[t-1,i] + u[i] + pro_dev[t,i]; // x at t = x at t-1 + growth + process error
    }
  }
}
model {
  sd_q ~ student_t(3, 0, 1); // general prior for process error
  for(i in 1:N){ // cycle through group
    x0[i] ~ normal(y0[i], 2); // prior for initial depends on first observation
    sd_r[i] ~ student_t(3, 0, 1); // general prior for observation error
    for(t in 1:TT){ // cycle through time
      pro_dev_tilde[t,i] ~ normal(0, 1); // latent process error
    }
  }
  u ~ normal(0, 2); // general prior for growth rate
  for(i in 1:n_pos){ // cycle through observations
    y[i] ~ normal(x[col_indx_pos[i], row_indx_pos[i]], sd_r[row_indx_pos[i]]); // map y values onto x matrix (flipped axes), x value is mean, group-level observation error is standard error
  }
}
//generated quantities {
//  vector[n_pos] log_lik; // log-likelihood value for each observation
//  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[col_indx_pos[n], row_indx_pos[n]], sd_r[row_indx_pos[n]]);
//}
