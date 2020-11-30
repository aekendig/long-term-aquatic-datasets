#### info ####

# goal: visualize herbicide activity


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")


#### edit data ####

# add columns
qual2 <- qual %>%
  mutate(date = as.Date(paste(Month, Day, Year, sep = "-"), "%m-%d-%Y"),
         lake_county = paste(Lake, County, sep = "_") %>%
           as.factor(),
         secchi_combined = ifelse(is.na(SECCHI_ft), SECCHI_2, SECCHI_ft),
         lake_group = cut(as.numeric(lake_county), breaks = 30),
         lake_county = fct_rev(lake_county))

# summarize by lake
qual3 <- qual2 %>%
  group_by(lake_county, lake_group, date) %>%
  summarize(TP = sum(!is.na(TP_ug_L)),
            TN = sum(!is.na(TN_ug_L)),
            chlorophyll = sum(!is.na(CHL_ug_L)),
            turbidity = sum(!is.na(secchi_combined)),
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