# function to remove duplicates from abundance plant data
rem_abun_dups <- function(dat) {
  
  if (nrow(dat) == 1){
    dat2 <- dat
  } else if (var(dat$AreaCovered_ha) > 0){
    dat2 <- dat %>%
      filter(AreaCovered_ha == max(AreaCovered_ha)) # choose survey with maximum area
  } else if (var(dat$AreaCovered_ha) == 0){
    dat2 <- dat %>%
      mutate(DaysDiff = abs(lakeO_area_days - Days)) %>%
      filter(DaysDiff == min(DaysDiff)) %>% # choose survey closes to max. abundance time
      select(-DaysDiff)
  }
  
  return(dat2)
}