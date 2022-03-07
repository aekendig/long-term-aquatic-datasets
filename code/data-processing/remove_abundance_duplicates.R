# function to remove duplicates from abundance plant data
rem_abun_dups <- function(dat) {
  
  if (sum(!is.na(dat$AreaCovered_ha)) == 0){ # all values are NA
    
    dat2 <- dat[1,]
    
  } else {
    
    dat1 <- dat %>% filter(!is.na(AreaCovered_ha)) # remove NA values
    
    if (nrow(dat1) == 1){  # only one value
    
    dat2 <- dat1
    
    } else  if (var(dat1$AreaCovered_ha) > 0){ # values differ from one another
      
    dat2 <- dat1 %>%
      filter(AreaCovered_ha == max(AreaCovered_ha)) # choose survey with maximum area
    
    } else { # values are the same
      
    dat2 <- dat1 %>%
      mutate(DaysDiff = abs(lakeO_area_days - Days)) %>%
      filter(DaysDiff == min(DaysDiff)) %>% # choose survey closes to max. abundance time
      select(-DaysDiff)
    
    }
  }
  
  return(dat2)
}