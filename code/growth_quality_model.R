qual_model <- function(dat, lag_interval, mod_name){
  
  # assign treatment column
  dat <- dat %>%
    select(starts_with(lag_interval), c("LogPropCovered", "PrevPropCoveredAdjCS", "PermanentID", "GSYear"))
  
  # # fit model
  # abu_mod <- glmmTMB(LogPropCovered ~ Trt * SurveyorExperienceB + PrevPropCoveredAdjCS + (1|PermanentID) + (1|GSYear), family = gaussian(), data = dat)  
  # 
  # # save model
  # save(abu_mod, file = paste0("output/", mod_name, "_model.rda"))
  # 
  # return model object
  return(dat)
}