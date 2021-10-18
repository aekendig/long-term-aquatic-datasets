#### info ####

# goal: visualize management activity from 97/98 to 2000


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(fst)

# import data
ctrl_old <- read_fst("original-data/PrePMARS_IPMData.fst")


#### edit data ####

# species
unique(ctrl_old$Species)

# two AOI columns
ctrl_old %>%
  filter(as.character(AOI) != as.character(AreaOfInterest)) %>%
  select(AOI, AreaOfInterest) %>%
  unique()
# not identical, AreaOfInterest seems more accurate

ctrl_old %>%
  filter(is.na(AreaOfInterest))
# 4 
# also missing county
# only one record of interest (floating plants)
# all ID's are zero

ctrl_old %>%
  filter(is.na(AOI))
# 0

# isolate missing data
(ctrl_old_miss <- ctrl_old %>%
    filter(is.na(County) | is.na(AreaOfInterest)) %>%
    data.frame())

# check rest of data for AOI
ctrl_old %>%
  filter(AOI %in% ctrl_old_miss$AOI)

ctrl_old %>%
  filter(AreaOfInterest %in% ctrl_old_miss$AOI)
# no other information

# duplicate AreaOfInterest names
ctrl_old_AOI <- ctrl_old %>%
  mutate(AreaOfInterest = as.character(AreaOfInterest),
         AreaOfInterest = ifelse(is.na(AreaOfInterest), 
                                 as.character(AOI), 
                                 AreaOfInterest)) %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterest)) %>%
  filter(dup == T) %>%
  select(AreaOfInterest) %>%
  unique()

# same lake in multiple counties or different lakes?
ctrl_old %>%
  filter(AreaOfInterest %in% ctrl_old_AOI$AreaOfInterest) %>%
  select(AreaOfInterest, AreaOfInterestID, County) %>%
  unique() %>%
  arrange(AreaOfInterest)
# Administration: different (AOI = "General lakes")
# Alligator: different
# Butler: different
# Francis: different
# Jackson: different
# Minnehaha: different
# Silver: different
# Trout: different
# Withlacoochee River: different

# make Area ID a factor
ctrl_old2 <- ctrl_old %>%
  filter(!is.na(AreaOfInterest)) %>%
  mutate(AreaOfInterest = as.character(AreaOfInterest),
         AOIF = as.factor(ifelse(AreaOfInterest %in% ctrl_old_AOI$AreaOfInterest, 
                                 paste(AreaOfInterest, County, sep = "_"), 
                                 AreaOfInterest)),
         AOI_group = cut(as.numeric(AOIF), breaks = 8),
         AOIF = fct_rev(AOIF),
         managed_species = case_when(Species == "Hydrilla verticillata" ~ "hydrilla",
                                     str_detect(Species, "Eichhornia") == T | str_detect(Species, "Pistia") == T ~ "floating",
                                     TRUE ~ "other"),
         hydrilla = ifelse(managed_species == "hydrilla", 1, 0),
         floating = ifelse(managed_species == "floating", 1, 0))

# check species grouping
ctrl_old2 %>%
  group_by(managed_species) %>%
  summarise(species = paste(unique(Species), collapse = ", "))

# groups for pdf
AOI_grp = sort(unique(ctrl_old2$AOI_group))

# years of dataset
ctrl_old_years = tibble(Year = c(min(ctrl_old2$Year), max(ctrl_old2$Year))) %>%
  expand_grid(ctrl_old2 %>%
                select(AOIF, AOI_group) %>%
                unique())


#### visualizations ####

# management over time
pdf("output/old_management_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_old_sub <- ctrl_old2 %>%
    filter(AOI_group == AOI_grp[i])
  
  ctrl_old_years_sub <- ctrl_old_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_old_sub,
               aes(x = Year,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_old_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          scale_x_continuous(breaks = 1998:2010) +
          xlab("Year (points = management action)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# hydrilla management over time
pdf("output/old_hydrilla_management_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_old_sub <- ctrl_old2 %>%
    filter(AOI_group == AOI_grp[i] & hydrilla == 1)
  
  ctrl_old_years_sub <- ctrl_old_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_old_sub,
               aes(x = Year,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_old_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          scale_x_continuous(breaks = 1998:2010) +
          xlab("Year (points = hydrilla management action)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# floating management over time
pdf("output/old_floating_management_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_old_sub <- ctrl_old2 %>%
    filter(AOI_group == AOI_grp[i] & floating == 1)
  
  ctrl_old_years_sub <- ctrl_old_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_old_sub,
               aes(x = Year,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = ctrl_old_years_sub,
                    size = 0.3,
                    aes(color = AOIF)) +
          scale_x_continuous(breaks = 1998:2010) +
          xlab("Year (points = floating plant management action)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
}
dev.off()

# managed species
pdf("output/old_management_lake_species.pdf")
ctrl_old2 %>%
  group_by(managed_species) %>%
  summarise(units = length(unique(AOIF))) %>%
  ggplot(aes(x = managed_species, y = units)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = units), vjust = -0.25, size = 3) +
  ylab("Water bodies") +
  xlab("Managed species") +
  theme_bw()

ctrl_old2 %>%
  group_by(managed_species) %>%
  summarise(units = length(unique(paste(AOIF, Year)))) %>%
  ggplot(aes(x = managed_species, y = units)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = units), vjust = -0.25, size = 3) +
  ylab("Management actions") +
  xlab("Managed species") +
  theme_bw()

dev.off()

# management actions per lake
pdf("output/old_management_action_histogram.pdf")
ctrl_old2 %>%
  group_by(AOIF, managed_species) %>%
  summarise(actions = n()) %>%
  ggplot(aes(x = actions, fill = managed_species)) +
  geom_histogram(binwidth = 1, color = "black") +
  scale_x_log10() +
  ylab("Number of water bodies") +
  xlab("Number of management actions (1998-2010)") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
dev.off()


#### output data as csv ####
write_csv(ctrl_old, "original-data/PrePMARS_IPMData.csv")
