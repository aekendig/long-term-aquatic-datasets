#### info ####

# goal: visualize herbicide activity from 97/98 to 2000


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(fst)

# import data
mgmt <- read_fst("original-data/PrePMARS_IPMData.fst")


#### edit data ####

# species
unique(mgmt$Species)

# two AOI columns
mgmt %>%
  filter(as.character(AOI) != as.character(AreaOfInterest)) %>%
  select(AOI, AreaOfInterest) %>%
  unique()
# not identical, AreaOfInterest seems more accurate

mgmt %>%
  filter(is.na(AreaOfInterest))
# 4 

mgmt %>%
  filter(is.na(AOI))
# 0

# duplicate AreaOfInterest names
mgmt_AOI <- mgmt %>%
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
mgmt %>%
  filter(AreaOfInterest %in% mgmt_AOI$AreaOfInterest) %>%
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
mgmt2 <- mgmt %>%
  mutate(AreaOfInterest = as.character(AreaOfInterest),
         AOIF = as.factor(ifelse(AreaOfInterest %in% mgmt_AOI$AreaOfInterest, 
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
mgmt2 %>%
  group_by(managed_species) %>%
  summarise(species = paste(unique(Species), collapse = ", "))

# groups for pdf
AOI_grp = sort(unique(mgmt2$AOI_group))

# years of dataset
mgmt_years = tibble(Year = c(min(mgmt2$Year), max(mgmt2$Year))) %>%
  expand_grid(mgmt2 %>%
                select(AOIF, AOI_group) %>%
                unique())


#### visualizations ####

# management over time
pdf("output/old_management_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  mgmt_sub <- mgmt2 %>%
    filter(AOI_group == AOI_grp[i])
  
  mgmt_years_sub <- mgmt_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = mgmt_sub,
               aes(x = Year,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = mgmt_years_sub,
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
  
  mgmt_sub <- mgmt2 %>%
    filter(AOI_group == AOI_grp[i] & hydrilla == 1)
  
  mgmt_years_sub <- mgmt_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = mgmt_sub,
               aes(x = Year,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = mgmt_years_sub,
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
  
  mgmt_sub <- mgmt2 %>%
    filter(AOI_group == AOI_grp[i] & floating == 1)
  
  mgmt_years_sub <- mgmt_years %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = mgmt_sub,
               aes(x = Year,
                   y = AOIF,
                   color = AOIF)) +
          geom_point() +
          geom_line(data = mgmt_years_sub,
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
mgmt2 %>%
  group_by(managed_species) %>%
  summarise(units = length(unique(AOIF))) %>%
  ggplot(aes(x = managed_species, y = units)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = units), vjust = -0.25, size = 3) +
  ylab("Water bodies") +
  xlab("Managed species") +
  theme_bw()

mgmt2 %>%
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
mgmt2 %>%
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
write_csv(mgmt, "original-data/PrePMARS_IPMData.csv")
