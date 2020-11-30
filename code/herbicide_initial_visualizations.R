#### info ####

# goal: visualize herbicide activity


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")


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

# duplicate AreaOfInterest names
AOI <- ctrl %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterest)) %>%
  filter(dup == T)

# indicate herbicide treatments
# make Area ID a factor
ctrl2 <- ctrl %>%
  mutate(HerbicideF = ifelse(is.na(ControlMethod) |
                               (ControlMethod %in% non_herb),
                             "no",
                             "yes"),
         AOIF = as.factor(ifelse(AreaOfInterest %in% AOI$AreaOfInterest, paste(AreaOfInterest, AreaOfInterestID, sep = "_"), AreaOfInterest)),
         AOI_group = cut(as.numeric(AOIF), breaks = 8),
         AOIF = fct_rev(AOIF),
         date = as.Date(BeginDate, "%m/%d/%y"))

AOI_grp = sort(unique(ctrl2$AOI_group))


#### visualizations ####

# herbicide use over time
pdf("output/herbicide_lake_time_series.pdf")
for(i in 1:length(AOI_grp)){
  
  ctrl_sub <- ctrl2 %>%
    filter(AOI_group == AOI_grp[i])
  
  print(ggplot(data = ctrl_sub,
               aes(x = date,
                   y = AOIF,
                   color = AOIF)) +
          geom_line(size = 0.5) +
          geom_point(data = filter(ctrl_sub, HerbicideF == "yes")) +
          xlab("Date (points = herbicide application)") +
          ylab("Water body") +
          theme_bw() +
          theme(legend.position = "none"))
  
}
dev.off()
# relevant characteristics: time series length, application intensity

# time series length
pdf("output/herbicide_lake_time_histogram.pdf")
ctrl2 %>%
  group_by(AOIF) %>%
  summarise(days = max(date)-min(date)) %>%
  mutate(years = as.numeric(days)/365) %>%
  ggplot(aes(x = years)) +
  geom_histogram(binwidth = 1) +
  ylab("Number of water bodies") +
  xlab("Length of dataset (years)") +
  theme_bw()
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
