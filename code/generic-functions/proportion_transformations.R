# functions to transform data to account for 0's and 1's
# Douma and Weedon 2019
# for cases where values cannot exactly equal 0 or 1
transform01 <- function(x) {
  n <- sum(!is.na(x))
  (x * (n - 1) + 0.5) / n
}

backtransform01 <- function(x) {
  n <- sum(!is.na(x))
  (x * n - 0.5) / (n - 1)
}  

# for car function when adjust > 0
# https://stackoverflow.com/questions/23845283/logit-transformation-backwards
inv_logit_adjust <- function(f,a) {
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

# standardize proportion adjustment
prop_adjust <- 0.001

# inverse logit
logit2prob <- function(x) {
  exp(x)/(1+exp(x))
}