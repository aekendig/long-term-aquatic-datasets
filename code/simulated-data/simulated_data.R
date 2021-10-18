#### info ####

# goal: simulate data to test pre-herbicide model


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# model settings
l <- 10 # number of lakes
tt <- 18  # number of years
seed <- 220 # random number seed


#### global process model ####

# n_t+1 = n_t x exp[u x (1 - n_t/k)]
# log(n_t+1/n_t) = beta0 + beta1 x n_t + e_t
# beta0 = growth rate (u)
# beta1 = growth rate/carrying capacity (-u/k)

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
# tails are too fat on t-distribution, large variances
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
# correct estimates when w is removed

# save
write_csv(df, "intermediate-data/simulated_data_pre_herbicide_global_model.csv")


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
# model may not be identifiable

# okay to use proportions?
plot(herb_mod1)


#### process model ####

# process error sd
set.seed(seed)
q <- abs(rt(n = 1, df = 3)) # default mean = 0, variance = 3

# growth rate by lake
set.seed(seed)
u <- rnorm(n = l, mean = 0, sd = 2)

# process error by year and lake
set.seed(seed)
w <- matrix(rnorm(l*tt, mean = 0, sd = q),
            nrow = tt,
            ncol = l)

# initial population size
set.seed(seed)
x0 <- rnorm(n = l, mean = log(50), sd = 10)

# process values
x <- matrix(NA, nrow = tt, ncol = l)
x[1, ] <- x0 + u + w[1, ]
for(i in 2:tt){
  x[i, ] <- x[(i - 1), ] + u + w[(i), ]
}


#### observation model ####

# observation error sd
set.seed(seed)
r <- abs(rt(n = l, df = 3))

# observation error
v <- matrix(NA, nrow = tt, ncol = l)
set.seed(seed)
for(i in 1:l){
  v[ , i] <- rnorm(n = tt, mean = 0, sd = r[i])
}

# observation values
y <- x + v


#### format data ####

# extract indices
y_indx <- which(!is.na(y), arr.ind = T)

# generate NA's 
na_val <- 10
set.seed(seed)
na_indx <- y_indx[sample(l*tt, size = na_val), ]
y[na_indx] <- NA

# transpose
y2 <- t(y)

# make dataframe
y3 <- data.frame(y2)

# round all values to 3 decimal points
y4 <- round(y3, digits = 3)

# save
write_csv(y4, "intermediate-data/simulated_data_pre_herbicide_model.csv")
