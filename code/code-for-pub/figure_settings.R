# load package
library(khroma)

# figure settings
def_theme_paper <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 9, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.background = element_blank(),
        strip.text = element_text(size = 9, color="black"),
        strip.background = element_blank(),
        plot.title = element_text(size = 9))

paper_text_size = 2.5

col_pal <- colour("muted")(7)[c(3, 5, 7)]
