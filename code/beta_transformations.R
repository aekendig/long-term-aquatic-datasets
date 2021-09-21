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