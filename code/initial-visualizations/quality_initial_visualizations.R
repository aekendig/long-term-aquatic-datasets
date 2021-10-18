#### info ####

# goal: visualize water quality data


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")


#### edit data ####

# check for missing location data
sum(is.na(qual$Lake))
sum(is.na(qual$County))

# add columns
qual2 <- qual %>%
  mutate(date = as.Date(paste(Month, Day, Year, sep = "-"), "%m-%d-%Y"),
         lake_county = paste(Lake, County, sep = "_") %>%
           as.factor(),
         secchi_combined = ifelse(is.na(SECCHI_ft), SECCHI_2, SECCHI_ft) %>%
           tolower(),
         secchi_bottom = ifelse(str_detect(secchi_combined, "bottom") == T, 1, 0),
         secchi_weeds = ifelse(str_detect(secchi_combined, "weeds") == T, 1, 0),
         secchi = case_when(secchi_combined == "." ~ NA_character_, 
                            secchi_combined == "weeds" ~ NA_character_,
                            secchi_combined == "weeds (surface)" ~ NA_character_,
                            secchi_combined == "bottom" ~ NA_character_,
                            TRUE ~ secchi_combined) %>%
           parse_number(),
         lake_group = cut(as.numeric(lake_county), breaks = 30),
         lake_county = fct_rev(lake_county))

# summarize by lake
qual3 <- qual2 %>%
  group_by(lake_county, lake_group, date) %>%
  summarize(TP = sum(!is.na(TP_ug_L)),
            TN = sum(!is.na(TN_ug_L)),
            chlorophyll = sum(!is.na(CHL_ug_L)),
            turbidity = sum(!is.na(secchi)),
            color = sum(!is.na(Color_Pt_Co_Units)),
            conductivity = sum(!is.na(Cond_uS_cm))) %>%
  pivot_longer(cols = -c(lake_county:date), 
               names_to = "measurements", 
               values_to = "replicates") %>%
  filter(replicates > 0)

lake_grp = sort(unique(qual3$lake_group))


#### Visualizations ####

# lake quality sampling over time
pdf("output/lake_quality_time_series.pdf")
for(i in 1:length(lake_grp)){
  
  qual_sub <- qual3 %>%
    filter(lake_group == lake_grp[i])
  
  print(ggplot(data = qual_sub %>%
                 select(lake_county:date) %>%
                 unique(),
               aes(x = date,
                   y = lake_county,
                   color = lake_county)) +
          geom_line(size = 0.5) +
          geom_point(data = qual_sub, aes(shape = measurements), alpha = 0.5) +
          scale_shape_manual(values = 0:6) +
          xlab("Date") +
          ylab("Lake (Lake_County)") +
          theme_bw() +
          theme(plot.margin = margin(t = 1, r = 1, b = 20, l = 1, unit = "pt"),
                legend.position = c(0.4, -0.08),
                legend.direction = "horizontal",
                legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt")) +
          guides(color = "none", shape = guide_legend(nrow = 1, override.aes = list(alpha = 1))))
  
}
dev.off()

# time series length
pdf("output/lake_quality_time_histogram.pdf")
qual3 %>%
  select(lake_county:date) %>%
  unique() %>%
  group_by(lake_county) %>%
  summarise(days = max(date)-min(date)) %>%
  mutate(years = as.numeric(days)/365) %>%
  ggplot(aes(x = years)) +
  geom_histogram(binwidth = 1) +
  ylab("Number of lakes") +
  xlab("Length of dataset (years)") +
  theme_bw()
dev.off()


# number of lakes
pdf("output/lake_quality_counts.pdf")
qual3 %>%
  group_by(lake_county) %>%
  summarise(measurements = length(unique(measurements))) %>%
  group_by(measurements) %>%
  summarise(units = length(unique(lake_county))) %>%
  ggplot(aes(x = measurements, y = units)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = units), vjust = -0.25, size = 3) +
  ylab("Lakes") +
  xlab("Number of quality variables") +
  theme_bw()
dev.off()


# Total P
# removes two values of zero
pdf("output/lake_phosphorus_histogram.pdf")
qual2 %>%
  group_by(lake_county, date, Month, Year) %>%
  summarise(TP = mean(TP_ug_L, na.rm = T)) %>%
  ungroup() %>%
  mutate(month = month.name[Month] %>%
           reorder(Month)) %>%
  filter(TP > 0) %>%
  ggplot(aes(x = TP, fill = month)) +
  geom_histogram() +
  ylab("Lake-year combinations") +
  xlab("Total P (ug/L)") +
  scale_x_log10() +
  theme_bw()
dev.off()

# Total N
# remove two values of zero (same as TP, error in data?)
pdf("output/lake_nitrogen_histogram.pdf")
qual2 %>%
  group_by(lake_county, date, Month, Year) %>%
  summarise(TN = mean(TN_ug_L, na.rm = T)) %>%
  ungroup() %>%
  mutate(month = month.name[Month] %>%
           reorder(Month)) %>%
  filter(TN > 0) %>%
  ggplot(aes(x = TN, fill = month)) +
  geom_histogram() +
  ylab("Lake-year combinations") +
  xlab("Total N (ug/L)") +
  scale_x_log10() +
  theme_bw()
dev.off()

# Chlorophyll
# add 0.1 to all values (minimum = 0.1)
pdf("output/lake_chlorophyll_histogram.pdf")
qual2 %>%
  group_by(lake_county, date, Month, Year) %>%
  summarise(Chl = mean(CHL_ug_L, na.rm = T)) %>%
  ungroup() %>%
  mutate(month = month.name[Month] %>%
           reorder(Month),
         Chl2 = Chl + 0.1) %>%
  filter(!is.na(Chl)) %>%
  ggplot(aes(x = Chl2, fill = month)) +
  geom_histogram() +
  ylab("Lake-year combinations") +
  xlab("Chlorophyll + 0.1 (ug/L)") +
  scale_x_log10() +
  theme_bw()
dev.off()

# Secchi
pdf("output/lake_secchi_histogram.pdf")
qual2 %>%
  group_by(lake_county, date, Month, Year) %>%
  summarise(secchi = mean(secchi, na.rm = T)) %>%
  ungroup() %>%
  mutate(month = month.name[Month] %>%
           reorder(Month)) %>%
  filter(!is.na(secchi)) %>%
  ggplot(aes(x = secchi, fill = month)) +
  geom_histogram() +
  ylab("Lake-year combinations") +
  xlab("Water clarity (ft)") +
  theme_bw()
dev.off()

# color
# add 1 to visualize on log scale
# minimum value other than 0 is 1
pdf("output/lake_color_histogram.pdf")
qual2 %>%
  group_by(lake_county, date, Month, Year) %>%
  summarise(color = mean(Color_Pt_Co_Units, na.rm = T)) %>%
  ungroup() %>%
  mutate(month = month.name[Month] %>%
           reorder(Month),
         color2 = color +1) %>%
  filter(!is.na(color)) %>%
  ggplot(aes(x = color2, fill = month)) +
  geom_histogram() +
  ylab("Lake-year combinations") +
  xlab("Water color + 1 (Pt-Co)") +
  scale_x_log10() +
  theme_bw()
dev.off()

# conductivity
pdf("output/lake_conductivity_histogram.pdf")
qual2 %>%
  group_by(lake_county, date, Month, Year) %>%
  summarise(conductivity = mean(Cond_uS_cm, na.rm = T)) %>%
  ungroup() %>%
  mutate(month = month.name[Month] %>%
           reorder(Month)) %>%
  filter(!is.na(conductivity)) %>%
  ggplot(aes(x = conductivity, fill = month)) +
  geom_histogram() +
  ylab("Lake-year combinations") +
  xlab("Conductivity (uS/cm)") +
  scale_x_log10() +
  theme_bw()
dev.off()