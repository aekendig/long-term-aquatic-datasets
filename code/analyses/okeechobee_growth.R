# extract information from Lake Okeechobee model

# load model
load("output/lakeO_intraannual_veg_change_mod.rda")

# extract coefficients
lakeO_beta1 <- coef(summary(lake0_mod))[2, "Estimate"] # days
lakeO_beta2 <- coef(summary(lake0_mod))[3, "Estimate"] # days^2

# date of max abundance (-b/2a)
lakeO_days <- -lakeO_beta1 / (2 * lakeO_beta2)

# data for days
dayDat <- read_csv("intermediate-data/LakeO_day_data_for_model.csv")