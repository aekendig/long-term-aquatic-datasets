# function to find longest time interval of a species being present with the most lakes
time_int_fun <- function(year1, taxon, dat_in){
  
  # permanent ID's in dataset for each year
  perm_yr <- dat_in %>%
    select(PermanentID) %>%
    unique() %>%
    expand_grid(GSYear = min(dat_in$GSYear):max(dat_in$GSYear))
  
  # filter for taxon
  # add NA's for missing years
  dat <- dat_in %>%
    filter(CommonName == taxon) %>%
    full_join(perm_yr)
  
  dat2 <- dat %>%
    filter(GSYear >= year1 & is.na(SpeciesAcres)) %>% # select missing years
    group_by(PermanentID) %>%
    summarise(year2 = min(GSYear)) %>% # identify first year missing data
    ungroup() %>%
    full_join(dat %>%
                select(PermanentID) %>%
                unique()) %>% # add all lakes (in case some had no missing data)
    mutate(year2 = replace_na(year2, max(dat$GSYear) + 1),   # assign last year (add one to count that year)
           years = year2 - year1) %>%
    expand_grid(tibble(years_out = 0:(max(dat$GSYear) + 1 - year1))) %>% # expand for every years_out value
    filter(years >= years_out) %>%
    group_by(years_out) %>%
    mutate(waterbodies = n_distinct(PermanentID)) %>%
    ungroup() %>%
    mutate(data_points = years_out * waterbodies) 
  
  return(dat2)
}
