library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
library(leaflet)
library(magrittr)
library(arrow)
library(sf)        # st_as_sf(), st_bbox()
library(maptiles)  # get_tiles(): downloads map tiles without API key
library(tidyterra) # geom_spatraster_rgb(): renders tile SpatRaster in ggplot2
# If not installed:
# install.packages(c("sf", "maptiles", "tidyterra"))

df_persons    <- read_parquet("data/RTI_household/persons/CA_2019_persons.parquet")
df_households <- read_parquet("data/RTI_household/households/CA_2019_households.parquet")

df_merged <- df_persons %>%
  left_join(df_households, by = "hh_id")

df_merged_sf <- df_merged %>%
  filter(
    lon_4326 >= -122.52 & lon_4326 <= -122.35,
    lat_4326 >= 37.70   & lat_4326 <= 37.82
  )

# ============================================================
# Interactive leaflet map (unchanged)
# ============================================================
m <- df_merged_sf %>%
  leaflet() %>%
  addTiles() %>%
  setView(lng = -122.44, lat = 37.76, zoom = 12) %>%
  addCircleMarkers(
    lng         = ~lon_4326,
    lat         = ~lat_4326,
    radius      = 0.03,
    color       = "blue",
    stroke      = FALSE,
    fillOpacity = 0.02,
    popup       = ~paste("Age:", agep, "<br> Household ID:", hh_id)
  )

m  # display in viewer

# ============================================================
# Static map — maptiles + ggplot2 (no API key, no PhantomJS)
# ============================================================

# Convert points to sf object (WGS84)
pts_sf <- st_as_sf(df_merged_sf,
                   coords = c("lon_4326", "lat_4326"),
                   crs    = 4326)

# Download background tiles for the bounding box
# Provider options (all free, no key):
#   "OpenStreetMap"          — standard OSM road map
#   "CartoDB.Positron"       — minimal light/grey style
#   "CartoDB.DarkMatter"     — dark minimal style
#   "Esri.WorldShadedRelief" — terrain/relief feel
tiles <- get_tiles(
  x        = pts_sf,
  provider = "CartoDB.Positron",   # <- change style here
  zoom     = 13,
  crop     = TRUE
)

# Build ggplot — square canvas enforced by coord_sf(expand = FALSE)
p_map <- ggplot() +
  geom_spatraster_rgb(data = tiles) +          # basemap tiles
  stat_density_2d(
    data          = df_merged_sf,
    aes(x = lon_4326, y = lat_4326, fill = after_stat(level)),
    geom          = "polygon",
    contour       = TRUE,
    bins          = 50,
    alpha         = 0.2
  ) +
  scale_fill_gradientn(
    colours = c("#FFFF00", "#FFA500", "#FF4500", "#CC0000"),  # yellow -> orange -> red
    name    = "Density"
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title    = "San Francisco — Household Locations (2019)",
    subtitle = "RTI synthetic population",
    x        = NULL,
    y        = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    plot.margin   = margin(10, 10, 10, 10),
    legend.position      = "right",
    legend.title         = element_text(size = 9),
    legend.text          = element_text(size = 8)
  )

# ============================================================
# Save as PNG and PDF — square (7 x 7 inches), no external tools needed
# ============================================================
if (!dir.exists("output"))  dir.create("output",  recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

ggsave("figures/sf_households_map.png",
       plot   = p_map,
       width  = 7,
       height = 7,
       dpi    = 300,
       units  = "in")

ggsave("figures/sf_households_map.pdf",
       plot   = p_map,
       width  = 7,
       height = 7,
       units  = "in")

cat("Saved: figures/sf_households_map.png\n")
cat("Saved: figures/sf_households_map.pdf\n")

# ============================================================
# Save processed data
# ============================================================
write_parquet(df_merged_sf, "output/df_hh_sf.parquet")