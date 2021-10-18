# extract information from Lake Okeechobee model

# load model
load("output/lakeO_intraannual_veg_area_change_mod.rda")
load("output/lakeO_intraannual_veg_prop_change_mod.rda")

# extract coefficients
lakeO_area_beta1 <- coef(summary(lake0_area_mod))[2, "Estimate"] # days
lakeO_area_beta2 <- coef(summary(lake0_area_mod))[3, "Estimate"] # days^2
lakeO_prop_beta1 <- coef(summary(lake0_prop_mod))[2, "Estimate"] # days
lakeO_prop_beta2 <- coef(summary(lake0_prop_mod))[3, "Estimate"] # days^2

# date of max abundance (-b/2a)
lakeO_area_days <- -lakeO_area_beta1 / (2 * lakeO_area_beta2)
lakeO_prop_days <- -lakeO_prop_beta1 / (2 * lakeO_prop_beta2)
# should be equal

# data for days
dayDat <- read_csv("intermediate-data/LakeO_day_data_for_model.csv")