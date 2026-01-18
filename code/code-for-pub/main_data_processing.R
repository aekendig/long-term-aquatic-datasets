#### description ####

# processes raw data into dataset for target and methods analyses
# make Figures 2 and 4, Table S1


#### set-up ####

# load packages
library(tidyverse)
library(ggtext)
library(patchwork)

# import data
plants <- read_csv("intermediate-data/FWC_plant_formatted.csv")
mgmt <- read_csv("intermediate-data/FWC_management_formatted.csv")
loc_dat <- read_csv("intermediate-data/FWC_lakes_formatted.csv")

# figure settings
source("code/code-for-pub/figure_settings.R")


#### format data ####

# filter plant dataset for lakes
# remove problematic survey
# select plant survey data within management timeframe
# select years where all species were surveyd
# select waterbodies with at least three years of data
plants2 <- plants %>%
  filter(WaterbodyType == "Lake" & Outlier == 0) %>%
  cross_join(mgmt %>%
               summarize(MinTreatmentYear = min(TreatmentYear),
                         MaxTreatmentYear = max(TreatmentYear))) %>%
  filter(SurveyYear >= MinTreatmentYear & SurveyYear <= MaxTreatmentYear) %>%
  select(-c(MinTreatmentYear, MaxTreatmentYear)) %>%
  group_by(AreaOfInterestID) %>%
  mutate(SurveyYears = n_distinct(SurveyYear)) %>%
  ungroup() %>%
  filter(SurveyYears >= 3)

# waterbodies surveyed
# range of survey years
waterbodies <- plants2 %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County, 
           FType, WaterbodyAcres) %>%
  summarize(MinSurveyYear = min(SurveyYear),
            MaxSurveyYear = max(SurveyYear),
            .groups = "drop")

# select waterbodies from management data with surveys
# select years during survey span
mgmt2 <- mgmt %>%
  inner_join(waterbodies) %>%
  filter(TreatmentYear >= MinSurveyYear & TreatmentYear <= MaxSurveyYear) %>%
  mutate(TaxonName = replace_na(TaxonName, "unknown"),
         Species = replace_na(Species, "unknown"),
         Target = case_when(Species %in% c("Eichhornia crassipes", 
                                           "Pistia stratiotes", 
                                           "Floating Plants (Eichhornia and Pistia)") ~ 
                              "Float",
                            Species == "Hydrilla verticillata" ~ "Hydr",
                            TRUE ~ "Other")) 

# management from new dataset (with details)
mgmt_new <- mgmt2 %>%
  filter(CtrlSet == "new") %>%
  select(-c(MinSurveyYear, MaxSurveyYear))

# get surveys within new data time frame
plants_new <- plants2 %>%
  cross_join(mgmt_new %>%
               summarize(MinTreatmentYear = min(TreatmentYear),
                         MaxTreatmentYear = max(TreatmentYear))) %>%
  filter(SurveyYear >= MinTreatmentYear & SurveyYear <= MaxTreatmentYear) %>%
  select(-c(MinTreatmentYear, MaxTreatmentYear)) %>%
  group_by(AreaOfInterestID) %>%
  mutate(SurveyYears = n_distinct(SurveyYear)) %>%
  ungroup() %>%
  filter(SurveyYears >= 3)

# waterbodies surveyed
# range of survey years
waterbodies_new <- plants_new %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County, 
           FType, WaterbodyAcres) %>%
  summarize(MinSurveyYear = min(SurveyYear),
            MaxSurveyYear = max(SurveyYear),
            .groups = "drop")

# select waterbodies from management data with surveys
# select years during survey span
# make methods categories
# get survey years and prop treated for summaries
mgmt_new2 <- mgmt_new %>%
  inner_join(waterbodies_new) %>%
  filter(TreatmentYear >= MinSurveyYear & TreatmentYear <= MaxSurveyYear) %>%
  mutate(Method = case_when(MethodHerbicide == "yes" & Contact == 1 ~ "Con", # contact
                            MethodHerbicide == "yes" & Contact == 0 ~ "Sys", # systemic
                            MethodHerbicide == "no" ~ "Non", # non-herbicide
                            TRUE ~ "Unk"), # unknown
         SurveyYears = MaxSurveyYear - MinSurveyYear + 1,
         PropTreated = if_else(TotalAcres <= WaterbodyAcres, # Some exceed waterbody acres. Not sure why.
                               TotalAcres / WaterbodyAcres,
                               1))


#### management target dataset ####

# richness
# choose the survey with the highest richness for a given waterbody and year
plant_target_sum <- plants2 %>%
  filter(IsDetected == "Yes" & !(TaxonName %in% c("Hydrilla verticillata",
                                                  "Eichhornia crassipes", 
                                                  "Pistia stratiotes"))) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, 
           SurveyYear, SurveyDate) %>%
  summarize(NativeRichness = sum(Origin == "Native"),
            .groups = "drop") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, 
           SurveyYear) %>%
  mutate(MaxNativeRichness = max(NativeRichness)) %>%
  ungroup() %>%
  filter(NativeRichness == MaxNativeRichness) %>% 
  select(-MaxNativeRichness)

# mean invasive cover (PAC)
# choose the survey with the highest abundance for a given waterbody and year
inv_target_sum <- plants2 %>%
  filter(TaxonName %in% c("Hydrilla verticillata", "Eichhornia crassipes", 
                          "Pistia stratiotes")) %>%
  mutate(Invasive = case_when(TaxonName == "Hydrilla verticillata" ~ "Hydr",
                              TaxonName == "Eichhornia crassipes" ~ "Hyac", 
                              TaxonName == "Pistia stratiotes" ~ "Lett")) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Invasive,
           WaterbodyAcres, SurveyYear) %>%
  summarize(SpeciesAcres = max(SpeciesAcres),
            .groups = "drop") %>%
  pivot_wider(names_from = Invasive,
              values_from = SpeciesAcres,
              names_glue = "{Invasive}Acres") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
  summarize(HydrCov = mean(100 * HydrAcres / WaterbodyAcres),
            FloatCov = mean(100 * (HyacAcres + LettAcres) / WaterbodyAcres),
            .groups = "drop") %>%
  mutate(HydrCovC = HydrCov - mean(HydrCov),
         FloatCovC = FloatCov - mean(FloatCov))

# summarize by waterbody and target
mgmt_target_sum <- mgmt2 %>%
  mutate(SurveyYears = MaxSurveyYear - MinSurveyYear + 1,
         PropTreated = if_else(TotalAcres <= WaterbodyAcres, # Some exceed waterbody acres. Not sure why.
                               TotalAcres / WaterbodyAcres,
                               1)) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Target) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  full_join(waterbodies %>% # add in waterbodies without management
              select(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
              expand_grid(Target = c("Float", "Hydr", "Other"))) %>%
  mutate(TrtFreq = replace_na(TrtFreq, 0),
         TrtArea = replace_na(TrtArea, 0)) %>%
  pivot_wider(names_from = Target,
              values_from = c(TrtFreq, TrtArea),
              names_glue = "{Target}{.value}") %>%
  mutate(HydrTrtFreqC = HydrTrtFreq - mean(HydrTrtFreq),
         FloatTrtFreqC = FloatTrtFreq - mean(FloatTrtFreq),
         OtherTrtFreqC = OtherTrtFreq - mean(OtherTrtFreq),
         HydrTrtAreaC = HydrTrtArea - mean(HydrTrtArea),
         FloatTrtAreaC = FloatTrtArea - mean(FloatTrtArea),
         OtherTrtAreaC = OtherTrtArea - mean(OtherTrtArea))

# latitude
lat_target_sum <- loc_dat %>% 
  filter(AreaOfInterestID %in% mgmt_target_sum$AreaOfInterestID) %>% 
  mutate(LatitudeC = Latitude - mean(Latitude)) 

# combine data
# create a time variable
target_dat <- plant_target_sum %>%
  full_join(inv_target_sum) %>%
  full_join(mgmt_target_sum) %>%
  left_join(lat_target_sum %>% 
              select(AreaOfInterestID, Latitude, LatitudeC)) %>% 
  mutate(Time = SurveyYear - min(SurveyYear),
         SurveyDay = yday(SurveyDate))

# for individual species
target_taxa_dat <- plants2 %>%
  filter(!(TaxonName %in% c("Hydrilla verticillata", "Eichhornia crassipes", 
                            "Pistia stratiotes"))) %>%
  select(AreaOfInterestID, SurveyYear, TaxonName, Origin, Habitat, 
         IsDetected) %>% 
  inner_join(target_dat) %>%
  mutate(Detected = if_else(IsDetected == "Yes", 1, 0))

# save
write_csv(target_dat, 
          "intermediate-data/FWC_plant_management_target_analysis_formatted.csv")
write_csv(target_taxa_dat, 
          "intermediate-data/FWC_plant_management_target_taxa_analysis_formatted.csv")


#### management target figure ####

# summarize management across groups
mgmt_target_year_sum <- mgmt2 %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, 
           TreatmentYear, Target) %>%
  count(TreatmentYear, Target) %>%
  mutate(Target = fct_recode(Target,
                             "*Hydrilla verticillata*" = "Hydr",
                             "floating plants" = "Float",
                             "other plants" = "Other"))

# visualize
mgmt_time_fig <- ggplot(mgmt_target_year_sum, aes(x = TreatmentYear, y = n, color = Target)) +
  geom_line() + 
  geom_point() +
  scale_color_manual(values = col_pal) +
  labs(x = "Year", y = "Number of waterbodies with management records") +
  def_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.15),
        legend.text = element_markdown(),
        axis.text.x = element_text(hjust = 0.7),
        axis.title.y = element_text(hjust = 0))

# expand on other species
mgmt_target_other_sum <- mgmt2 %>%
  filter(!(Species %in% c("Eichhornia crassipes", "Pistia stratiotes", 
                          "Floating Plants (Eichhornia and Pistia)",
                          "Hydrilla verticillata"))) %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, 
           TreatmentYear, TaxonName) %>%
  mutate(TaxonName = case_when(TaxonName %in% c("Mixed Grasses", "Trees (blocking navigation)",
                                                "Tussocks", "Other", "unknown") ~ "unknown taxa",
                               str_detect(TaxonName, "AHRES") ~ "unknown taxa",
                               TRUE ~ TaxonName) %>%
           str_replace_all(" spp.| spp| sp.| sp", " spp.") %>%
           str_replace("Limnobium spp.ngia", "Limnobium spongia") %>% # not sure how to exclude from above
           str_remove_all(", sub| \\(exotic\\)| \\(other natives\\)| \\(other\\)|, natives|, emersed|, sub/floating") %>%
           str_remove_all("\\/floating"),
         Genus = word(TaxonName, 1, 1)) %>%
  count(TaxonName) %>%
  mutate(Taxon = if_else(TaxonName != "unknown taxa",
                         paste0("*", TaxonName, "*"),
                         TaxonName),
         Taxon = fct_reorder(Taxon, n))

# visualize
mgmt_other_fig <- mgmt_target_other_sum %>%
  filter(n >= 20) %>%
  ggplot(aes(y = Taxon, x = n)) +
  geom_col(fill = col_pal[2]) +
  labs(y = "Other plants", x = "Number of recorded management events") +
  def_theme_paper +
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_text(hjust = 1))

# save figures
mgmt_target_fig <- mgmt_time_fig + mgmt_other_fig +
  plot_layout(nrow = 1, widths = c(1, 0.7)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(",
                  tag_suffix = ")",
                  title = "Figure 2") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0),
        plot.title = element_text(size = 12))

ggsave("output/figure_2.tiff", mgmt_target_fig, 
       device = "tiff", width = 18, height = 9.5,
       units = "cm", dpi = 600)


#### management methods dataset ####

# richness
plant_methods_sum <- plants_new %>%
  filter(IsDetected == "Yes" & !(TaxonName %in% c("Hydrilla verticillata", 
                                                  "Eichhornia crassipes", 
                                                  "Pistia stratiotes")) &
           AreaOfInterestID %in% mgmt_new2$AreaOfInterestID) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, SurveyYear, SurveyDate) %>%
  summarize(NativeRichness = sum(Origin == "Native"),
            .groups = "drop") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, SurveyYear) %>%
  mutate(MaxNativeRichness = max(NativeRichness)) %>%
  ungroup() %>%
  filter(NativeRichness == MaxNativeRichness) %>% 
  select(-MaxNativeRichness)

# invasive cover
inv_methods_sum <- plants_new %>%
  filter(TaxonName %in% c("Hydrilla verticillata", 
                          "Eichhornia crassipes", 
                          "Pistia stratiotes")&
           AreaOfInterestID %in% mgmt_new2$AreaOfInterestID) %>%
  mutate(Invasive = case_when(TaxonName == "Hydrilla verticillata" ~ "Hydr", 
                              TaxonName == "Eichhornia crassipes" ~ "Hyac", 
                              TaxonName == "Pistia stratiotes" ~ "Lett")) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Invasive,
           WaterbodyAcres, SurveyYear) %>%
  summarize(SpeciesAcres = max(SpeciesAcres),
            .groups = "drop") %>%
  pivot_wider(names_from = Invasive,
              values_from = SpeciesAcres,
              names_glue = "{Invasive}Acres") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
  summarize(HydrPAC = mean(100 * HydrAcres / WaterbodyAcres),
            FloatPAC = mean(100 * (HyacAcres + LettAcres) / WaterbodyAcres),
            .groups = "drop") %>%
  mutate(HydrPACc = HydrPAC - mean(HydrPAC),
         FloatPACc = FloatPAC - mean(FloatPAC))

# summed area treated and frequency for each method
mgmt_methods_sum <- mgmt_new2 %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Method) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  pivot_wider(names_from = Method,
              values_from = c(TrtArea, TrtFreq),
              names_glue = "{.value}{Method}") %>%
  mutate(across(.cols = starts_with("Trt"), 
                .fns = ~replace_na(.x, 0)),
         TrtAreaConC = TrtAreaCon - mean(TrtAreaCon),
         TrtAreaSysC = TrtAreaSys - mean(TrtAreaSys),
         TrtAreaNonC = TrtAreaNon - mean(TrtAreaNon),
         TrtFreqConC = TrtFreqCon - mean(TrtFreqCon),
         TrtFreqSysC = TrtFreqSys - mean(TrtFreqSys),
         TrtFreqNonC = TrtFreqNon - mean(TrtFreqNon)) 

# same as above, but combine systemic and contact into herbicide
mgmt_herb_sum <- mgmt_new2 %>%
  filter(MethodHerbicide == "yes") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
  summarize(TrtFreqHerb = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtAreaHerb = mean(100 * PropTreated),
            .groups = "drop") %>%
  full_join(mgmt_new2 %>%
              distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County)) %>%
  mutate(across(.cols = starts_with("Trt"), 
                .fns = ~replace_na(.x, 0)),
         TrtAreaHerbC = TrtAreaHerb - mean(TrtAreaHerb),
         TrtFreqHerbC = TrtFreqHerb - mean(TrtFreqHerb)) 

# average management month
mgmt_methods_month_sum <- mgmt_new2 %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
  summarize(TrtMonth = mean(TreatmentMonth),
            .groups = "drop") %>%
  mutate(TrtMonthStd = TrtMonth - mean(TrtMonth))

# latitude
lat_methods_sum <- loc_dat %>% 
  filter(AreaOfInterestID %in% mgmt_methods_sum$AreaOfInterestID) %>% 
  mutate(LatitudeC = Latitude - mean(Latitude)) 

# combine data
methods_dat <- plant_methods_sum %>%
  full_join(inv_methods_sum) %>%
  full_join(mgmt_methods_sum %>% 
               full_join(mgmt_herb_sum) %>%
               full_join(mgmt_methods_month_sum)) %>%
  left_join(lat_methods_sum %>% 
              select(AreaOfInterestID, Latitude, LatitudeC)) %>% 
  mutate(Time = SurveyYear - min(SurveyYear))

# taxa dat for specific methods
methods_taxa_dat <- plants_new %>%
  filter(!(TaxonName %in% c("Hydrilla verticillata", "Eichhornia crassipes", 
                            "Pistia stratiotes"))) %>%
  select(AreaOfInterestID, SurveyYear, TaxonName, Origin, Habitat, 
         IsDetected) %>% 
  mutate(Detected = if_else(IsDetected == "Yes", 1, 0)) %>%
  inner_join(methods_dat)

# save
write_csv(methods_dat, 
          "intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
write_csv(methods_taxa_dat, 
          "intermediate-data/FWC_plant_management_methods_taxa_analysis_formatted.csv")


#### management methods figure ####

# management methods proportion
mgmt_methods_prop <- mgmt_new2 %>%
  count(TreatmentYear, Method, Target) %>%
  group_by(TreatmentYear, Target) %>%
  mutate(Prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(Method = fct_recode(Method,
                             "contact herbicide" = "Con",
                             "systemic herbicide" = "Sys",
                             "non-herbicide" = "Non",
                             "unknown" = "Unk") %>%
           fct_relevel("contact herbicide", "systemic herbicide"),
         Target = fct_recode(Target,
                             "Floating plants" = "Float",
                             "*Hydrilla verticillata*" = "Hydr",
                             "Other plants" = "Other") %>%
           fct_relevel("*Hydrilla verticillata*"))

# management methods figure
mgmt_methods_fig <- mgmt_methods_prop %>%
  filter(Method != "unknown") %>%
  ggplot(aes(x = as.character(TreatmentYear),
             y = Prop, fill = Method)) +
  geom_col(color = "white", linewidth = 0.2) +
  facet_wrap(~ Target) +
  scale_fill_manual(values = col_pal) +
  labs(x = "Year", y = "Proportion of recorded management events",
       title = "Figure 4") +
  def_theme_paper +
  theme(strip.text = element_markdown(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.114, 0.7),
        legend.background = element_rect(fill = "white", color = NA),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.title = element_text(size = 12))

ggsave("output/figure_4.tiff", mgmt_methods_fig, device = "tiff",
       width = 18, height = 9.5, units = "cm", dpi = 600)


#### smaller management target dataset ####

mgmt_methods_sum <- mgmt_new2 %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Method) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  pivot_wider(names_from = Method,
              values_from = c(TrtArea, TrtFreq),
              names_glue = "{.value}{Method}") %>%
  mutate(across(.cols = starts_with("Trt"), 
                .fns = ~replace_na(.x, 0)),
         TrtAreaConC = TrtAreaCon - mean(TrtAreaCon),
         TrtAreaSysC = TrtAreaSys - mean(TrtAreaSys),
         TrtAreaNonC = TrtAreaNon - mean(TrtAreaNon),
         TrtFreqConC = TrtFreqCon - mean(TrtFreqCon),
         TrtFreqSysC = TrtFreqSys - mean(TrtFreqSys),
         TrtFreqNonC = TrtFreqNon - mean(TrtFreqNon)) 

# management target (for the new dataset)
mgmt_target_sum_new <- mgmt_new2 %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Target) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  full_join(mgmt_new2 %>% # add in all waterbodies, all targets
              distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
              expand_grid(Target = c("Float", "Hydr", "Other"))) %>%
  mutate(TrtFreq = replace_na(TrtFreq, 0),
         TrtArea = replace_na(TrtArea, 0)) %>%
  pivot_wider(names_from = Target,
              values_from = c(TrtFreq, TrtArea),
              names_glue = "{Target}{.value}") %>%
  mutate(HydrTrtFreqC = HydrTrtFreq - mean(HydrTrtFreq),
         FloatTrtFreqC = FloatTrtFreq - mean(FloatTrtFreq),
         OtherTrtFreqC = OtherTrtFreq - mean(OtherTrtFreq),
         HydrTrtAreaC = HydrTrtArea - mean(HydrTrtArea),
         FloatTrtAreaC = FloatTrtArea - mean(FloatTrtArea),
         OtherTrtAreaC = OtherTrtArea - mean(OtherTrtArea))

# latitude
lat_methods_sum_new <- loc_dat %>% 
  filter(AreaOfInterestID %in% mgmt_target_sum_new$AreaOfInterestID) %>% 
  mutate(LatitudeC = Latitude - mean(Latitude)) 

# combine (same plant data as method daaset)
target_dat_new <- plant_methods_sum %>%
  full_join(inv_methods_sum) %>%
  full_join(mgmt_target_sum_new) %>%
  left_join(lat_methods_sum_new %>% 
              select(AreaOfInterestID, Latitude, LatitudeC)) %>% 
  mutate(Time = SurveyYear - min(SurveyYear))

# save
write_csv(target_dat_new, 
          "intermediate-data/FWC_plant_management_new_target_analysis_formatted.csv")


#### values for text ####

# proportion events with fluridone
mgmt_new2 %>% 
  filter(str_detect(ControlMethod, "Fluridone") & 
           Species == "Hydrilla verticillata") %>% 
  nrow() /
  mgmt_new2 %>% 
  filter(Species == "Hydrilla verticillata") %>% 
  nrow()

# systemic herbicides
mgmt_new2 %>% 
  filter(Method == "Sys") %>% 
  count(ControlMethod) %>% 
  arrange(desc(n))

# contact herbicides
mgmt_new2 %>% 
  filter(Method == "Con") %>% 
  count(ControlMethod) %>% 
  arrange(desc(n))

# non-herbicides
mgmt_new2 %>% 
  filter(Method == "Non") %>% 
  count(ControlMethod) %>% 
  arrange(desc(n))

# proportion mechanical
mgmt_new2 %>% 
  filter(ControlMethod != "Grass Carp" & 
           Method == "Non") %>% 
  nrow() /
  mgmt_new2 %>% 
  filter(Method == "Non") %>% 
  nrow()

# unknown
mgmt_new2 %>% 
  filter(Method == "Unk") %>% 
  count(ControlMethod) %>% 
  arrange(desc(n))

# summarize taxa by habitat
plants2 %>%
  distinct(TaxonName, Origin, Habitat) %>%
  count(Origin, Habitat)

plants_new %>%
  distinct(TaxonName, Origin, Habitat) %>%
  count(Origin, Habitat)

# table of taxa in surveys
plants2 %>%
  distinct(TaxonName, Origin, Habitat) %>%
  full_join(plants_new %>%
              distinct(TaxonName, Origin, Habitat)) %>%
  mutate(Habitat = tolower(Habitat) %>%
           fct_recode("emergent" = "emersed") %>%
           fct_relevel("submersed"),
         Origin = tolower(Origin) %>%
           fct_recode("non-native" = "exotic") %>%
           fct_relevel("native")) %>%
  arrange(Origin, Habitat, TaxonName) %>%
  write_csv("output/surveyed_taxa_table.csv")
