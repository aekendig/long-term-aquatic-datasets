# function to find longest time interval with the most lakes
time_int_fun2 <- function(year1){
  
  dat <- nat_fwc %>% 
    filter(TaxonName == nat_fwc$TaxonName[1] & GSYear < 2020) %>%
    full_join(inv_fwc2 %>%
                select(PermanentID, GSYear) %>%
                unique()) # all possible surveys
  
  dat2 <- dat %>%
    filter(GSYear >= year1 & is.na(Detected)) %>% # select missing years
    group_by(PermanentID) %>%
    summarise(year2 = min(GSYear)) %>% # identify first year missing data
    ungroup() %>%
    full_join(dat %>%
                select(PermanentID) %>%
                unique()) %>% # add all lakes (in case some had no missing data)
    mutate(year2 = replace_na(year2, max(dat$GSYear) + 1),   # assign last year (add one to count that year)
           years = year2 - year1) %>%
    group_by(years) %>%
    count() %>% # summarise number of lakes per timespan
    ungroup() %>%
    expand_grid(tibble(years_out = 0:(max(dat$GSYear) + 1 - year1))) %>% # expand for every years_out value
    filter(years >= years_out) %>%
    group_by(years_out) %>%
    summarise(lakes = sum(n))
  
  return(dat2)
}