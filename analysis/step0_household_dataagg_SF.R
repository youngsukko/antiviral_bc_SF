# ==============================================================
# SF Household Location Map
# Background : maptiles — CartoDB.Positron
# Points     : one dot per household, semi-transparent
# Data       : RTI synthetic population, California 2019
#
# Required packages:
#   sf, maptiles, tidyterra, ggplot2, dplyr, arrow
# ==============================================================

library(dplyr)
library(arrow)
library(ggplot2)
library(sf)
library(maptiles)
library(tidyterra)

# --------------------------------------------------------------
# 0. Bounding box — San Francisco
# --------------------------------------------------------------
SF_LON_MIN <- -122.52
SF_LON_MAX <- -122.35
SF_LAT_MIN <-   37.70
SF_LAT_MAX <-   37.82

# --------------------------------------------------------------
# 1. Load and filter RTI synthetic population
# --------------------------------------------------------------
df_persons    <- read_parquet("data/RTI_household/persons/CA_2019_persons.parquet")
df_households <- read_parquet("data/RTI_household/households/CA_2019_households.parquet")

df_merged_sf <- df_persons %>%
  left_join(df_households, by = "hh_id") %>%
  filter(
    lon_4326 >= SF_LON_MIN, lon_4326 <= SF_LON_MAX,
    lat_4326 >= SF_LAT_MIN, lat_4326 <= SF_LAT_MAX
  )

cat(sprintf("Households plotted: %d\n", nrow(df_merged_sf)))

# --------------------------------------------------------------
# 2. Download background tiles (CartoDB.Positron — minimal style)
# --------------------------------------------------------------
pts_sf <- st_as_sf(df_merged_sf,
                   coords = c("lon_4326", "lat_4326"),
                   crs    = 4326)

tiles <- get_tiles(
  x        = pts_sf,
  provider = "CartoDB.Positron",
  zoom     = 13,
  crop     = TRUE
)

# --------------------------------------------------------------
# 3. Plot
# --------------------------------------------------------------
p <- ggplot() +
  
  # Basemap tiles
  geom_spatraster_rgb(data = tiles) +
  
  # One dot per household
  geom_point(
    data  = df_merged_sf,
    aes(x = lon_4326, y = lat_4326),
    color = "#2B3FBF",   # deep blue
    size  = 0.1,        # decrease if too dense; increase if too sparse
    alpha = 0.05         # semi-transparent: dense areas appear darker naturally
  ) +
  
  coord_sf(expand = FALSE) +
  
  labs(
    # title    = "San Francisco — Household Locations (2019)",
    # subtitle = "RTI synthetic population · each dot = one household",
    x = NULL,
    y = NULL
  ) +
  
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    # plot.title      = element_text(hjust = 0.5, face = "bold",
    #                                size  = 13,  margin = margin(b = 4)),
    # plot.subtitle   = element_text(hjust = 0.5, color = "grey45",
    #                                size  = 9,   margin = margin(b = 8)),
    plot.margin     = margin(12, 12, 12, 12)
  )

# --------------------------------------------------------------
# 4. Save
# --------------------------------------------------------------
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

ggsave("figures/SF_households_map.png",
       plot = p, width = 7, height = 7, dpi = 500, units = "in")

# ggsave("figures/sf_households_map.pdf",
#        plot = p, width = 7, height = 7, units = "in")

# ============================================================
# Save processed data
# ============================================================
write_parquet(df_merged_sf, "output/df_hh_sf.parquet")