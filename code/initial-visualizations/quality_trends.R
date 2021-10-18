#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
qual <- read_csv("intermediate-data/LW_quality_formatted.csv",
                 col_types = list(PermanentID = col_character(),
                                  JoinNotes = col_character()))


#### edit data ####

# summarize across stations
qual2 <- qual %>%
  group_by(County_LW, Lake, PermanentID, ShapeArea, Date, Month, Day, Year) %>%
  summarise(across(.cols = c(TP_ug_L, TN_ug_L, CHL_ug_L, Secchi, Color_Pt_Co_Units, Cond_uS_cm), 
                   ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(Lake_Year = paste0(County_LW, Lake, Year))


#### within-year trends ####

ggplot(qual2, aes(x = Month, y = TP_ug_L)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "line", fun = mean)
# peaks in April
# lowest in January

ggplot(qual2, aes(x = Month, y = TN_ug_L)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "line", fun = mean)
# peaks in April
# lowest in July

ggplot(qual2, aes(x = Month, y = CHL_ug_L)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "line", fun = mean)
# peaks in September
# lowest in January

ggplot(qual2, aes(x = Month, y = Secchi)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "line", fun = mean)
# peaks in January
# lowest in October

ggplot(qual2, aes(x = Month, y = Color_Pt_Co_Units)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "line", fun = mean)
# peaks in September
# lowest in April/May

ggplot(qual2, aes(x = Month, y = Cond_uS_cm)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "line", fun = mean)
# peaks in December
# lowest in August


#### sample sizes ####

ggplot(qual2, aes(x = Month)) +
  geom_bar()
# highest in March
# lowest in December
# relatively consistent throughout the year

ggplot(qual2, aes(x = Month)) +
  geom_bar() +
  facet_wrap(~ Year)
# low in December may be because of 2019
# not all samples processed?