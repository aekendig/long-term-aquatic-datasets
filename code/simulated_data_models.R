#### info ####

# goal: simulate data to test models


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(rstan)

# avoid recompiling unchanged Stan program
rstan_options(auto_write = TRUE)

# do in parallel
options(mc.cores = parallel::detectCores())

# model settings
l <- 10 # number of lakes
tt <- 18  # number of years
seed <- 220 # random number seed


#### global process model ####

# n_t+1 = n_t x exp[u x (1 - n_t/k)]
# log(n_t+1/n_t) = beta0 + beta1 x n_t + e_t
# beta0 = growth rate (u)
# beta1 = growth rate/carrying capacity (-u/k)
# non-hierarchical

# extend time series because there is only one population
tt2 <- tt * 4

# process error variance
df_t <- 8 # degrees of freedom
x <- seq(-2, 2, length.out = 30)
p <- dt(x, df = df_t)
plot(x, p, type = "l")
x2 <- rt(n = 100, df = df_t)
hist(x2)

sd_norm <- 0.05
p2 <- dnorm(x, sd = sd_norm)
plot(x, p2, type = "l")
x3 <- rnorm(n = 100, sd = sd_norm)
hist(x3)

set.seed(seed)
# q <- abs(rt(n = 1, df = df_t)) # default mean = 0, variance = df/(df-2) for df > 2
# variances too large with t-distribution
q <- abs(rnorm(n = 1, sd = sd_norm))

# growth rate (should be between ~ -1 and 1)
u <- 0.1

# carrying capacity
k <- 0.8

# initial population size
x0 <- 0.05

# process error by year
set.seed(seed)
w <- rnorm(n = tt2, mean = 0, sd = q)

# log population size
n1 <- rep(NA, tt2)
n2 <- rep(NA, tt2)
n1[1] <- u - (u/k) * x0 + log(x0)
n2[1] <- u - (u/k) * x0 + log(x0) + w[1]
for(i in 2:tt2){
  n1[i] <- u - (u/k) * exp(n1[i-1]) + n1[i-1]
  n2[i] <- u - (u/k) * exp(n2[i-1]) + n2[i-1] + w[i]
}

# dataframe 
df <- tibble(time = 1:tt2, logN = n2) %>%
  mutate(N = exp(logN),
         logNDiff = lead(logN) - logN,
         error = w,
         N_raw = exp(n1))

# visualize
ggplot(df, aes(time, N)) +
  geom_line() +
  geom_line(aes(y = N_raw), color = "red") +
  geom_line(aes(y = error), color = "blue")

ggplot(df, aes(N, logNDiff)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")

summary(lm(logNDiff ~ N, data = df))
# should have intercept = u, slope = -u/k

# model file
mod_global1 <- 'models/global_process_model_1.stan'

# remove NA's from data (last year)
dat_global1 <- filter(df, !is.na(logNDiff))

# data
input_global1 <- list(TT = max(dat_global1$time),
                      y = dat_global1$logNDiff,
                      x = dat_global1$N)

# fit model
fit_global1 <- stan(file = mod_global1, 
                    data = input_global1, 
                    iter = 4000, chains = 3, thin = 1, warmup = 1000)

# summary
print(fit_global1)


#### global process model + herbicides ####

# n_t+1 = n_t x exp[u x (1 - [n_t + h x h_t]/k)]
# log(n_t+1/n_t) = beta0 + beta1 x n_t + beta2 x h_t + e_t
# beta0 = growth rate (u)
# beta1 = growth rate/carrying capacity (-u/k)
# beta2 = herbicide effect (-uh/k)

# make longer to see if estimates improve
tt3 <- tt2 * 4

# herbicide effect
h <- 1

# herbicide time series
ht <- runif(n = tt3)

# make errors smaller to see if estimates improve (no)
set.seed(seed)
q2 <- abs(rnorm(n = 1, sd = 0.01))

# process error by year
set.seed(seed)
w <- rnorm(n = tt3, mean = 0, sd = q2)

# log population size
n1 <- rep(NA, tt3)
n2 <- rep(NA, tt3)
n1[1] <- u - (u/k) * x0 - (u*h/k) * ht[1] + log(x0)
n2[1] <- u - (u/k) * x0 - (u*h/k) * ht[1] + log(x0) + w[1]
for(i in 2:tt3){
  n1[i] <- u - (u/k) * exp(n1[i-1]) - (u*h/k) * ht[i] + n1[i-1]
  n2[i] <- u - (u/k) * exp(n2[i-1]) - (u*h/k) * ht[i] + n2[i-1] + w[i]
}

# dataframe 
df <- tibble(time = 1:tt3, logN = n2) %>%
  mutate(N = exp(logN),
         logNDiff = lead(logN) - logN,
         error = w,
         N_raw = exp(n1),
         herbicide = ht)

# visualize
ggplot(df, aes(time, N)) +
  geom_line() +
  geom_line(aes(y = N_raw), color = "red") +
  geom_line(aes(y = error), color = "blue")

ggplot(df, aes(N, logNDiff)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(df, aes(herbicide, logNDiff + u*N/k)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")
# not sure how to visualize herbicide effect

herb_mod1 <- lm(logNDiff ~ N + herbicide, data = df)
summary(herb_mod1)
u_est <- coef(herb_mod1)[1]
k_est <- -u_est/coef(herb_mod1)[2]
h_est <- coef(herb_mod1)[3] * -k_est / u_est
# still off by a lot, even with smaller errors and longer time series

# okay to use proportions?
plot(herb_mod1)