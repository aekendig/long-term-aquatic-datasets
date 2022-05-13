# function to mean from cut_ functions
cut_mean <- function(x) {
  
  out <- str_extract_all(x, "-?[0-9.]+") %>%
    unlist() %>%
    as.numeric() %>%
    mean()
  
}