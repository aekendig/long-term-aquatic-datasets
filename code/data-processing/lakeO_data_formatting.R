#### info ####

# goal: Estimate change in plant cover over a year. Assume that most of the change in cover is due to plant growth, not reproduction/mortality. The growth rate should be in units of vegetation/(vegetation x days) so that it can be applied to lakes of any size. Used standardized areas to achieve this.


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(lme4)
library(scales)

# import data
lakeo <- read_csv("original-data/Lake_O_Helicopter_Data_1989_2020_by_Section.csv")
gis <- read_csv("intermediate-data/gis_fwc_lakewatch_fwri.csv",
                col_types = list(wkt_geom = col_character(),
                                 AOI = col_character(),
                                 Lake = col_character()))

# figure settings
source("code/settings/figure_settings.R")

#### edit data ####

# lake area
lake_area <- gis %>%
  filter(str_detect(AreaOfInterest, "Okeechobee") == T) %>%
  mutate(Area_ha = ShapeArea * 100) %>% # convert lake area from km-squared to hectares
  pull(Area_ha)

# make long
lakeo2 <- lakeo %>%
  mutate(Date = as.Date(Date, format = "%d-%b-%Y"),
         Year = year(Date),
         Month = month(Date),
         Day = day(Date),
         MonthDay = as.Date(paste("2020", Month, Day, sep = "-")), # 2020 is arbitrary
         Total = TorreyKreamer + RittaHalifax + EastWallCootBay + WestWallWhiddenObservation +
           FisheatingBay + HarneyIndianP + IndianPPierceKissimmeeCodys + KingsBar +
           KissimmeeTaylorEagleBay + TaylorLock7Chancey) %>%
  pivot_longer(cols = -c(Date, Year, Month, Day, MonthDay, Season, LakeRegulationSchedule, Agency, LakeLevel),
               names_to = "Section",
               values_to = "AreaCovered_acres") %>%
  mutate(AreaCovered_ha = AreaCovered_acres * 0.405,
         log_AreaCovered_ha = log(AreaCovered_ha + 0.001),
         PropCovered = AreaCovered_ha / lake_area)

# leap year?
lakeo2 %>%
  filter(Month == 2 & Day == 29)

# subset data
# days number
# standardized area
totDat <- lakeo2 %>%
  filter(Section == "Total" & !is.na(AreaCovered_ha)) %>%
  mutate(Yearf = case_when(Month >= 4 ~ as.character(Year),
                           Month < 4 ~ as.character(Year - 1)),
         MonthDay = case_when(Month >= 4 ~ as.Date(paste("2020", Month, Day, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                              Month < 4 ~ as.Date(paste("2021", Month, Day, sep = "-"))),
         Days = as.numeric(MonthDay - as.Date("2020-04-01")), # days away from 4/1 of the same yearf
         Days = Days / 364,
         AreaCovered_s = (AreaCovered_ha - mean(AreaCovered_ha)) / sd(AreaCovered_ha))

# data attributes
range(totDat$Year)


#### figure ####

# facet by section
ggplot(lakeo2, aes(MonthDay, AreaCovered_ha)) +
  geom_line(aes(color = as.factor(Year))) +
  geom_smooth(color = "black") +
  facet_wrap(~ Section, scales = "free") +
  theme(legend.position = "none")

# total only
lakeo2 %>%
  filter(Section == "Total") %>%
  ggplot(aes(MonthDay, AreaCovered_ha)) +
  geom_vline(xintercept = as.Date("2020-04-01")) +
  geom_line(aes(color = as.factor(Year))) +
  geom_smooth(color = "black") +
  theme(legend.position = "none")
# growing season seems to start around April/May

totDat %>%
  ggplot(aes(MonthDay, AreaCovered_s)) +
  geom_line(aes(color = Yearf)) +
  geom_smooth(color = "black") +
  theme(legend.position = "none")

totDat %>%
  ggplot(aes(MonthDay, PropCovered)) +
  geom_line(aes(color = Yearf)) +
  geom_smooth(color = "black") +
  theme(legend.position = "none")

totDat %>%
  ggplot(aes(Days)) +
  geom_histogram(binwidth = 0.1)

range(totDat$Days)


#### explore data ####

# total acreage zero for a streak within a year
lakeo2 %>%
  filter(Section == "Total" & AreaCovered_ha == 0)
# looked at original data and it does look like there was a crash

# range of dates
range(totDat$MonthDay) # 4/1 to 3/30


#### model fit ####

# models
lake0_area_mod <- lmer(AreaCovered_s ~ Days + I(Days^2) + (1 | Yearf), data = totDat)
summary(lake0_area_mod)

lake0_prop_mod <- lmer(PropCovered ~ Days + I(Days^2) + (1 | Yearf), data = totDat)
summary(lake0_prop_mod)

# model fit
fitDat <- totDat %>%
  select(Days, MonthDay) %>%
  unique() %>%
  mutate(AreaCovered_s = predict(lake0_area_mod, newdata = ., re.form = NA),
         PropCovered = predict(lake0_prop_mod, newdata = ., re.form = NA))


#### figure ####

pdf("output/lakeO_intraannual_veg_area_change_mod.pdf", width = 3.5, height = 3.5)
ggplot(totDat, aes(MonthDay, AreaCovered_s)) +
  geom_line(aes(color = Yearf)) +
  geom_line(data = fitDat, color = "black", size = 2) +
  labs(x = "Month", y = "Standardized area covered") +
  scale_x_date(labels = date_format("%b")) +
  def_theme_paper +
  theme(legend.position = "none")
dev.off()

pdf("output/lakeO_intraannual_veg_prop_change_mod.pdf", width = 3.5, height = 3.5)
ggplot(totDat, aes(MonthDay, PropCovered)) +
  geom_line(aes(color = Yearf)) +
  geom_line(data = fitDat, color = "black", size = 2) +
  labs(x = "Month", y = "Proportion of lake covered") +
  scale_x_date(labels = date_format("%b")) +
  def_theme_paper +
  theme(legend.position = "none")
dev.off()


#### days for other lakes ####

dayDat <- tibble(MonthDay = seq(as.Date("2020-04-01"), as.Date("2021-03-31"), by = 1)) %>%
  mutate(Days = as.numeric((MonthDay - min(MonthDay))),
         Days = Days / 364)


#### high numbers ####

totDat %>%
  filter(AreaCovered_s > 2) %>%
  select(Yearf) %>%
  unique()


#### output ####

# model
save(lake0_area_mod, file = "output/lakeO_intraannual_veg_area_change_mod.rda")
save(lake0_prop_mod, file = "output/lakeO_intraannual_veg_prop_change_mod.rda")
write_csv(totDat, "intermediate-data/lakeO_intraannual_veg_change_data.csv")
write_csv(dayDat, "intermediate-data/LakeO_day_data_for_model.csv")
