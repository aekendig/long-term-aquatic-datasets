lw_qual <- read_csv("intermediate-data/LW_quality_formatted.csv")
wa_qual <- read_csv("intermediate-data/water_atlas_quality_formatted.csv")

# water quality
qual <- lw_qual %>%
  full_join(wa_qual) %>%
  filter(QualityMetric == "Secchi_ft") %>%
  group_by(PermanentID, GSYear) %>%
  summarize(QualityValue = mean(QualityValue),
            MonthsSampled = mean(MonthsSampled))