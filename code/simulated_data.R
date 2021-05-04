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
seed <- 20 # random number seed


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

# save
write_csv(y3, "intermediate-data/simulated_data_pre_herbicide_model.csv")
