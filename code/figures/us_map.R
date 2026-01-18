# Load required libraries
library(ggplot2)
library(maps)

# Get map data for the contiguous US states
us_states <- map_data("state")

# Create a variable to identify Florida
us_states$fill_color <- ifelse(us_states$region == "florida", "black", "white")

# Create the map
us_map <- ggplot(us_states, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = fill_color), color = "black", linewidth = 0.1) +
  scale_fill_identity() +
  coord_map("albers", lat0 = 30, lat1 = 40) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# Calculate aspect ratio of the map data
# The contiguous US is roughly 1.6:1 (width:height)
lon_range <- diff(range(us_states$long))
lat_range <- diff(range(us_states$lat))
aspect_ratio <- lon_range / lat_range

# figure dimensions
fig_height = 1
fig_width = fig_height * aspect_ratio

# save
ggsave("output/us_map.tiff", us_map, device = "tiff", dpi = 600, 
       height = fig_height, width = fig_width, units = "in")
