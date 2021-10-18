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
  
  # surveyor bins
  surveyor_bins <- surveyor %>%
    select(SurveyorExperienceB) %>%
    unique() %>%
    rowwise() %>%
    mutate(Min = str_remove_all(SurveyorExperienceB, "\\[|\\]|\\(|\\)") %>%
             str_split(",") %>%
             unlist() %>% first() %>%
             as.numeric(),
           Max = str_remove_all(SurveyorExperienceB, "\\[|\\]|\\(|\\)") %>%
             str_split(",") %>%
             unlist() %>% last() %>%
             as.numeric())
  
  # categorize surveyor bins
  levels(surveyor_bins$SurveyorExperienceB) <- c("low", "medium", "high")
  
  # make wide
  surveyor_bins %<>%
    pivot_wider(names_from = SurveyorExperienceB,
                values_from = c(Min, Max),
                names_glue = "{SurveyorExperienceB}_{.value}")
  
  return(list(surveyor, surveyor_bins))
}
