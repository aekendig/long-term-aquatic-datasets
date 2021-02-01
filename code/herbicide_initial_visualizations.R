#### info ####

# goal: visualize management activity


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")


#### export data for QGIS ####

ctrl_gis <- ctrl %>%
  select(AreaOfInterestID, longitude, latitude, AreaOfInterest) %>%
  unique()

write_csv(ctrl_gis, "intermediate-data/FWC_Herbicide_Coordinates.csv")


#### edit data ####

# isolate herbicides
unique(ctrl$ControlMethod)
unique(ctrl$HerbicideType)
ctrl %>%
  group_by(HerbicideType) %>%
  summarise(Control = paste(unique(ControlMethod), collapse = ", ")) %>%
  data.frame()
# HerbicideType column has errors

# missing control method
mis_ctrl <- ctrl %>%
  filter(is.na(ControlMethod)) %>%
  select_if(~ !all(is.na(.)))

trtID <- ctrl %>%
  select(TreatmentID, HerbicideID:HerbicideType) %>%
  unique()

left_join(mis_ctrl, trtID) %>%
  data.frame()
# unable to determine control method using treatment ID

# non-herbicide methods
unique(ctrl$ControlMethod)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# unique AreaOfInterest
unique(ctrl$AreaOfInterest) %>%
  length() # 446

# unique AreaOfInterest ID
unique(ctrl$AreaOfInterestID) %>%
  length() # 454

# duplicate AreaOfInterest names
ctrl_AOI <- ctrl %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterest)) %>%
  filter(dup == T) %>%
  select(AreaOfInterest) %>%
  unique()

# counties per name
ctrl_county <- ctrl %>%
  select(AreaOfInterest, AreaOfInterestID, County) %>%
  unique() %>%
  group_by(AreaOfInterest) %>%
  summarise(IDs = length(AreaOfInterestID),
            Counties = length(County))

ctrl_county %>%
  filter(IDs != Counties)
# county is unique like AreaOfInterestID for each AreaOfInterest

ctrl_county %>%
  filter(Counties > 1)
# same lakes as ctrl_AOI

# same lake in multiple counties or different lakes?
ctrl %>%
  filter(AreaOfInterest %in% ctrl_AOI$AreaOfInterest) %>%
  select(AreaOfInterest, AreaOfInterestID, County) %>%
  unique() %>%
  arrange(AreaOfInterest)
# Alligator: different
# Butler: different
# Cypress: different
# Jackson: different
# Minnehaha: different
# Trout: different
# Withlacoochee River: different

# check for missing data
sum(is.na(ctrl$County)) # 0
sum(is.na(ctrl$AreaOfInterest)) # 0

# managed species
sort(unique(ctrl$species))

# indicate herbicide treatments
# make Area ID a factor
ctrl2 <- ctrl %>%
  mutate(HerbicideF = ifelse(is.na(ControlMethod) |
                               (ControlMethod %in% non_herb),
                             "no",
                             "yes"),
         ControlF = ifelse(!is.na(ControlMethod), "yes", "no"),
         AOIF = as.factor(ifelse(AreaOfInterest %in% ctrl_AOI$AreaOfInterest, 
                                 paste(AreaOfInterest, County, sep = "_"), 
                                 AreaOfInterest)),
         AOI_group = cut(as.numeric(AOIF), breaks = 8),
         AOIF = fct_rev(AOIF),
         date = as.Date(BeginDate, "%m/%d/%y"),
         managed_species = case_when(species == "Hydrilla verticillata" ~ "hydrilla",
                                     str_detect(species, "Eichhornia") == T | str_detect(species, "Pistia") == T ~ "floating",
                                     TRUE ~ "other"),
         hydrilla = ifelse(managed_species == "hydrilla", 1, 0),
         floating = ifelse(managed_species == "floating", 1, 0))

# groups for pdf
AOI_grp = sort(unique(ctrl2$AOI_group))

# years of dataset
ctrl_years = tibble(date = c(min(ctrl2$date), max(ctrl2$date))) %>%
  expand_grid(ctrl2 %>%
                select(AOIF, AOI_group) %>%
                unique())


#### visualizations ####

# management over time
pdf("output/management_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_sub <- ctrl2 %>%
    filter(AOI_group == AOI_grp[i] & ControlF == "yes")
  
  ctrl_years_sub <- ctrl_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_sub,
               aes(x = date,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          xlab("Date (points = management action)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# hydrilla management over time
pdf("output/hydrilla_management_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_sub <- ctrl2 %>%
    filter(AOI_group == AOI_grp[i] & ControlF == "yes" & hydrilla == 1)
  
  ctrl_years_sub <- ctrl_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_sub,
               aes(x = date,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          xlab("Date (points = hydrilla management action)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# hydrilla management over time
pdf("output/hydrilla_herbicide_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_sub <- ctrl2 %>%
    filter(AOI_group == AOI_grp[i] & HerbicideF == "yes" & hydrilla == 1)
  
  ctrl_years_sub <- ctrl_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_sub,
               aes(x = date,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          xlab("Date (points = hydrilla herbicide application)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# floating plant management over time
pdf("output/floating_management_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_sub <- ctrl2 %>%
    filter(AOI_group == AOI_grp[i] & ControlF == "yes" & floating == 1)
  
  ctrl_years_sub <- ctrl_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_sub,
               aes(x = date,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          xlab("Date (points = floating plant management action)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# hydrilla management over time
pdf("output/floating_herbicide_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_sub <- ctrl2 %>%
    filter(AOI_group == AOI_grp[i] & HerbicideF == "yes" & floating == 1)
  
  ctrl_years_sub <- ctrl_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_sub,
               aes(x = date,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          xlab("Date (points = hydrilla herbicide application)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# managed species
pdf("output/management_lake_species.pdf")
ctrl2 %>%
  filter(ControlF == "yes") %>%
  group_by(managed_species) %>%
  summarise(units = length(unique(AOIF))) %>%
  ggplot(aes(x = managed_species, y = units)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = units), vjust = -0.25, size = 3) +
  ylab("Water bodies") +
  xlab("Managed species") +
  theme_bw()

ctrl2 %>%
  filter(ControlF == "yes") %>%
  group_by(managed_species) %>%
  summarise(units = length(unique(paste(AOIF, date)))) %>%
  ggplot(aes(x = managed_species, y = units)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = units), vjust = -0.25, size = 3) +
  ylab("Management actions") +
  xlab("Managed species") +
  theme_bw()

dev.off()

# management actions per lake
pdf("output/management_action_histogram.pdf")
ctrl2 %>%
  filter(ControlF == "yes") %>%
  group_by(AOIF, managed_species) %>%
  summarise(actions = n()) %>%
  ggplot(aes(x = actions, fill = managed_species)) +
  geom_histogram(binwidth = 1, color = "black") +
  scale_x_log10() +
  ylab("Number of water bodies") +
  xlab("Number of management actions (2010-2020)") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))

ctrl2 %>%
  filter(HerbicideF == "yes") %>%
  group_by(AOIF, managed_species) %>%
  summarise(actions = n()) %>%
  ggplot(aes(x = actions, fill = managed_species)) +
  geom_histogram(binwidth = 1, color = "black") +
  scale_x_log10() +
  ylab("Number of water bodies") +
  xlab("Number of herbicide applications (2010-2020)") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
dev.off()

# number of lakes
pdf("output/herbicide_lake_counts.pdf")
ctrl2 %>%
  group_by(AOIF) %>%
  summarise(control = paste(unique(HerbicideF), collapse = "_" ) %>%
              recode("yes_no" = "both", 
                     "no_yes" = "both",
                     "yes" = "herbicide only",
                     "no" = "non-herbicide only")) %>%
  mutate(control = fct_relevel(control, "herbicide only",
                               "non-herbicide only")) %>%
  group_by(control) %>%
  summarise(units = length(unique(AOIF))) %>%
  ggplot(aes(x = control, y = units)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = units), vjust = -0.25, size = 3) +
  ylab("Water bodies") +
  xlab("Control method") +
  theme_bw()
dev.off()
