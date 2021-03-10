#### info ####

# goal: discrete logistic model simulation


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# figure settings
def_theme <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.text = element_text(size = 10, color="black"),
        strip.background = element_blank())


#### create data ####
Kval = 30

rdat <- tibble(r = c(0.1, 0.25, 0.5, 1),
               rx = c(13, 13.5, 15, 20),
               ry = c(1, 2.5, 1, 1)) %>%
  mutate(rtext = paste("r = ", r, sep = ""))

dat <- tibble(r = rep(rdat$r, each = 30),
              time = rep(1:30, 4)) %>%
  left_join(rdat) %>%
  mutate(n = ifelse(time == 1, 1, NA_real_),
         rtext = ifelse(time != 15, NA_character_, rtext))

for(i in rdat$r){
  for(t in 2:30){
    dat <- dat %>%
      mutate(n = case_when(r == i & time == t ~ lag(n) * exp(r * (1 - lag(n) / Kval)),
                           TRUE ~ n))
  }
}


#### figure ####
pdf("output/discrete_logistic_simulation.pdf", width = 4, height = 4)
ggplot(dat, aes(x = time, y = n, color = as.factor(r))) +
  geom_line() +
  geom_text(aes(x = rx, y = n+ry, label = rtext)) +
  xlab("Time") +
  ylab("Population size") +
  scale_color_viridis_d() +
  def_theme +
  theme(legend.position = "none")
dev.off()