growth_model <- function(dat, treat_col, mod_name){
  
  # assign treatment column
  treat_col <- enquo(treat_col)
  dat <- dat %>%
    mutate(Trt = !!treat_col)
  
  # fit model
  abu_mod <- glmmTMB(LogPropCovered ~ Trt * SurveyorExperienceB + PrevPropCoveredAdjCS + (1|PermanentID) + (1|GSYear), family = gaussian(), data = dat)  
    
  # save model
  save(abu_mod, file = paste0("output/", mod_name, "_model.rda"))
  
  # return model object
  return(abu_mod)
}