quality_dataset <- function(qual_dat, min_n){
  
  # summarize quality data by growing season year
  qual_dat2 <- qual_dat %>%
    mutate(GSYear = case_when(Month >= 4 ~ Year,
                              Month < 4 ~ Year - 1)) %>%
    group_by(PermanentID, GSYear) %>%
    summarise(across(.cols = c(TP_ug_L, TN_ug_L, CHL_ug_L, Secchi, Color_Pt_Co_Units, Cond_uS_cm), 
                     ~ mean(.x, na.rm = TRUE)),
              n_qual = n()) %>%
    ungroup() %>%
    filter(n_qual >= min_n) # at least min_n measurements per year
  
  return(qual_dat2)

}