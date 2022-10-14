# function to find longest time interval of monitoring with the most lakes
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

# function to find longest time interval of monitoring in each lake
time_cont_fun <- function(yearT, permID, dat_in){
  
  # select waterbody
  dat <- dat_in %>%
    filter(PermanentID == permID) 
  
  # add NA's for missing years
  dat2 <- dat %>%
    full_join(tibble(GSYear = min(dat$GSYear):max(dat$GSYear)))

  # identify years missing data closes to focal year
  dat3 <- dat2 %>%
    filter(is.na(SpeciesAcres)) %>% # select missing years
    transmute(Year = GSYear,
              YearsToT = yearT - Year, # years to focal year
              YearDir = if_else(YearsToT > 0, "before", "after")) %>% # before or after focal year
    group_by(YearDir) %>%
    mutate(YearsMin = min(abs(YearsToT))) %>% # years with missing data closes to focal year
    ungroup() %>%
    filter(abs(YearsToT) == YearsMin) # select for closest years

  # add min/max year if none are missing data
  if(!("before" %in% dat3$YearDir)){

    dat3 <- dat3 %>%
      add_row(Year = min(dat$GSYear) - 1, YearDir = "before")

  }

  if(!("after" %in% dat3$YearDir)){

    dat3 <- dat3 %>%
      add_row(Year = max(dat$GSYear) + 1, YearDir = "after")

  }
  
  return(dat3)
}

# function to find longest time interval of monitoring with the most lakes
time_int_qual_fun <- function(year1, taxon, quarter, dat_in){
  
  # permanent ID's in dataset for each year
  perm_yr <- dat_in %>%
    select(PermanentID) %>%
    unique() %>%
    expand_grid(Year = min(dat_in$Year):max(dat_in$Year))
  
  # filter for taxon & quarter
  # add NA's for missing years
  # count quality values per year and waterbody
  dat <- dat_in %>%
    filter(CommonName == taxon) %>%
    full_join(perm_yr) %>%
    group_by(PermanentID, Year, PercCovered) %>%
    summarize(QualityValues = sum(!is.na(Quarter))) %>%
    ungroup()
  
  dat2 <- dat %>%
    filter(Year >= year1 & (is.na(PercCovered) | QualityValues < 4)) %>% # select missing years
    group_by(PermanentID) %>%
    summarise(year2 = min(Year)) %>% # identify first year missing data
    ungroup() %>%
    full_join(dat %>%
                select(PermanentID) %>%
                unique()) %>% # add all lakes (in case some had no missing data)
    mutate(year2 = replace_na(year2, max(dat$Year) + 1),   # assign last year (add one to count that year)
           years = year2 - year1) %>%
    expand_grid(tibble(years_out = 0:(max(dat$Year) + 1 - year1))) %>% # expand for every years_out value
    filter(years >= years_out) %>%
    group_by(years_out) %>%
    mutate(waterbodies = n_distinct(PermanentID)) %>%
    ungroup() %>%
    mutate(data_points = years_out * waterbodies) 
  
  return(dat2)
}

