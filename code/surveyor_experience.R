surv_exp <- function(dat){
  
  # surveys per surveyor
  surveyor <- dat %>%
    select(Surveyor, AreaOfInterest, AreaOfInterestID, PermanentID, SurveyDate) %>%
    unique() %>%
    mutate(Survey = 1) %>%
    arrange(SurveyDate, AreaOfInterestID) %>%
    group_by(Surveyor) %>%
    mutate(SurveyorExperience = cumsum(Survey)) %>%
    ungroup() %>%
    select(-Survey) %>%
    filter(!is.na(Surveyor)) %>%
    mutate(SurveyorExperienceCS = (SurveyorExperience - mean(SurveyorExperience)) / sd(SurveyorExperience),
           SurveyorExperienceB = cut_number(SurveyorExperience, n = 3))
  
  return(surveyor)
}
