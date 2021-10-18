quality_abundance_dataset <- function(qual_dat, abu_dat){
  
  # combine treatment and invasion datasets
  abu_qual <- abu_dat %>%
    left_join(qual_dat)# only include data for lakes and years when surveys occurred
  
  return(abu_qual)
  
}

