# ==============================================================================
# scp_visualize_network_sf_plot.R
#
# Purpose:
#   Load pre-built SF network (net_data_sf.rds) and generate all figures.
#   No network construction — visualization only.
#
# Input:
#   output/net_data_sf.rds  : saved from scp_build_network_sf.R
#
# Outputs (figures/ folder):
#   plot_sf_hh_degree.png        : household contact degree distribution
#   plot_sf_comm_degree.png      : community contact degree distribution
#   plot_sf_network_sample.png   : synthetic example network (100 HH)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(igraph)
library(ggforce)   # geom_mark_ellipse()

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------
output_dir <- "output"
figure_dir <- "figures"
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load network data
# ------------------------------------------------------------------------------
rds_path <- file.path(output_dir, "net_data_sf.rds")
cat(sprintf("Loading: %s\n", rds_path))
net_data <- readRDS(rds_path)

cat(sprintf("  Population      : %d\n", nrow(net_data$nodes)))
cat(sprintf("  Households      : %d\n", net_data$n_households))
cat(sprintf("  Household edges : %d\n", nrow(net_data$hh_edges)))
cat(sprintf("  Community edges : %d\n", nrow(net_data$comm_edges)))


# ==============================================================================
# Figure 1 — Household degree distribution
# ==============================================================================
hh_degree <- data.frame(
  person_id = c(net_data$hh_edges$from, net_data$hh_edges$to)
) %>%
  count(person_id, name = "degree") %>%
  right_join(data.frame(person_id = net_data$nodes$person_id), by = "person_id") %>%
  mutate(degree = ifelse(is.na(degree), 0L, degree))

p_hh_deg <- ggplot(hh_degree, aes(x = degree)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count)) * 100),
    binwidth = 1, fill = "steelblue", color = "white", alpha = 0.85
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))+
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    # title    = "Household Contact Degree Distribution",
    # subtitle = sprintf("mean = %.2f | N = %d people",
    #                    mean(hh_degree$degree), nrow(hh_degree)),
    x = "Number of household contacts",
    y = "Proportion of individuals (%)"
  ) +
  theme_minimal(base_size = 12, base_family = "sans")+
  theme(plot.title = element_text(face = "bold"))+
  coord_cartesian(xlim = c(0, 10))

path_hh_deg <- file.path(figure_dir, "contact_hh_degree.png")
ggsave(path_hh_deg, p_hh_deg, width = 5, height = 4, dpi = 500)
cat(sprintf("Saved: %s\n", path_hh_deg))


# ==============================================================================
# Figure 2 — Community degree distribution
# ==============================================================================
comm_degree <- data.frame(
  person_id = c(net_data$comm_edges$from, net_data$comm_edges$to)
) %>%
  count(person_id, name = "degree") %>%
  right_join(data.frame(person_id = net_data$nodes$person_id), by = "person_id") %>%
  mutate(degree = ifelse(is.na(degree), 0L, degree))

p_comm_deg <- ggplot(comm_degree, aes(x = degree)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count)) * 100),
    binwidth = 1, fill = "tomato", color = "white", alpha = 0.85
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))+
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    # title    = "Community Contact Degree Distribution",
    # subtitle = sprintf("mean = %.2f | N = %d people",
    #                    mean(comm_degree$degree), nrow(comm_degree)),
    x = "Number of community contacts",
    y = "Proportion of individuals (%)"
  ) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(plot.title = element_text(face = "bold"))+
  coord_cartesian(xlim = c(0, 30))

path_comm_deg <- file.path(figure_dir, "contact_nonhh_degree.png")
ggsave(path_comm_deg, p_comm_deg, width = 5, height = 4, dpi = 500)
cat(sprintf("Saved: %s\n", path_comm_deg))


# ==============================================================================
# Figure 3 — Synthetic example network (100 households, same contact model)
# Note that this is not relevant with the simulation at all! (just visualization)
# ==============================================================================
library(ggforce)
library(maptiles)
library(tidyterra)
library(sf)

set.seed(123)

# --- Bounding box ---
CENTER_LON <- -122.4064
CENTER_LAT <-   37.7952
DELTA      <-   0.0005

# --- Extract all members from 20 sampled households ---
# Step 1: find households in the area (use hh representative coordinate)
hh_in_area <- df_merged_sf %>%
  filter(
    lon_4326 >= CENTER_LON - DELTA, lon_4326 <= CENTER_LON + DELTA,
    lat_4326 >= CENTER_LAT - DELTA, lat_4326 <= CENTER_LAT + DELTA
  ) %>%
  distinct(hh_id) %>%
  slice_sample(n = 20)

# Step 2: pull ALL members of those households
# jitter radius in degrees (~2-3m) so members of same HH don't overlap
JITTER_DEG <- 0.00003

syn_nodes <- df_merged_sf %>%
  filter(hh_id %in% hh_in_area$hh_id) %>%
  group_by(hh_id) %>%
  mutate(
    # household center = mean coordinate of members
    hh_cx = mean(lon_4326),
    hh_cy = mean(lat_4326),
    n_members = n(),
    member_idx = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    # evenly space members around HH center in a small circle
    angle     = (2 * pi * (member_idx - 1)) / n_members,
    x = hh_cx + ifelse(n_members > 1, JITTER_DEG * cos(angle), 0),
    y = hh_cy + ifelse(n_members > 1, JITTER_DEG * sin(angle), 0),
    person_id = row_number(),
    household = as.integer(factor(hh_id)),
    age_group = pmin(floor(agep / 5) + 1L, 16L)
  ) %>%
  select(person_id, household, age_group, x, y)

N_vis <- nrow(syn_nodes)
cat(sprintf("Total people: %d | Households: %d\n", N_vis, n_distinct(syn_nodes$household)))

# --- Household edges ---
syn_hh_edges <- do.call(rbind, lapply(unique(syn_nodes$household), function(hh) {
  members <- syn_nodes$person_id[syn_nodes$household == hh]
  if (length(members) < 2) return(NULL)
  pairs <- combn(members, 2)
  data.frame(from = pairs[1, ], to = pairs[2, ], stringsAsFactors = FALSE)
}))

# --- Community edges (Prem Poisson sampling) ---
prem_syn <- matrix(0.3, nrow = 16, ncol = 16)
diag(prem_syn) <- 3.0
prem_syn[3:4, 3:4] <- prem_syn[3:4, 3:4] + 2.5
prem_syn[5:9, 5:9] <- prem_syn[5:9, 5:9] + 1.5
prem_syn[3:4, 5:9] <- prem_syn[3:4, 5:9] + 0.8
prem_syn[5:9, 3:4] <- prem_syn[5:9, 3:4] + 0.8

ag_members_vis <- lapply(1:16, function(ag)
  syn_nodes$person_id[syn_nodes$age_group == ag])

comm_from <- integer(0); comm_to <- integer(0)
for (i in seq_len(N_vis)) {
  ag_i <- syn_nodes$age_group[i]
  hh_i <- syn_nodes$household[i]
  for (ag_j in 1:16) {
    k <- rpois(1, prem_syn[ag_i, ag_j])
    if (k == 0) next
    candidates <- ag_members_vis[[ag_j]]
    candidates <- candidates[
      candidates != i &
        syn_nodes$household[candidates] != hh_i &
        candidates > i
    ]
    if (length(candidates) == 0) next
    sampled   <- sample(candidates, min(k, length(candidates)), replace = FALSE)
    comm_from <- c(comm_from, rep(i, length(sampled)))
    comm_to   <- c(comm_to, sampled)
  }
}
syn_comm_edges <- data.frame(from = comm_from, to = comm_to,
                             stringsAsFactors = FALSE)

# --- Edge coordinate tables ---
node_xy <- syn_nodes %>% select(person_id, x, y)

hh_edge_plot <- if (!is.null(syn_hh_edges) && nrow(syn_hh_edges) > 0) {
  syn_hh_edges %>%
    left_join(node_xy, by = c("from" = "person_id")) %>% rename(x_start = x, y_start = y) %>%
    left_join(node_xy, by = c("to"   = "person_id")) %>% rename(x_end   = x, y_end   = y)
} else NULL

comm_edge_plot <- if (nrow(syn_comm_edges) > 0) {
  syn_comm_edges %>%
    left_join(node_xy, by = c("from" = "person_id")) %>% rename(x_start = x, y_start = y) %>%
    left_join(node_xy, by = c("to"   = "person_id")) %>% rename(x_end   = x, y_end   = y)
} else NULL

# --- Background map tiles (zoom 18 = building level) ---
pts_sf <- st_as_sf(syn_nodes, coords = c("x", "y"), crs = 4326)

bbox_expanded <- st_bbox(
  c(xmin = min(syn_nodes$x) - JITTER_DEG * 4,
    xmax = max(syn_nodes$x) + JITTER_DEG * 8,
    ymin = min(syn_nodes$y) - JITTER_DEG * 4,
    ymax = max(syn_nodes$y) + JITTER_DEG * 8),
  crs = st_crs(4326)
) %>% st_as_sfc()

tiles <- get_tiles(
  x        = bbox_expanded,
  provider = "CartoDB.Positron",
  zoom     = 18,
  crop     = TRUE
)

# --- Age group palette ---
age_pal        <- colorRampPalette(c("#4575b4","#91bfdb","#fee090","#fc8d59","#d73027"))(16)
age_labels_vis <- paste0(seq(0, 75, 5), "-", seq(4, 79, 5), "y")

# --- Plot ---
p_net <- ggplot() +
  # Basemap
  geom_spatraster_rgb(data = tiles) +
  # Community edges
  {if (!is.null(comm_edge_plot))
    geom_segment(data = comm_edge_plot,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 color = "gray50", alpha = 0.35, linewidth = 0.35)} +
  # Household bubble
  geom_mark_ellipse(
    data        = syn_nodes,
    aes(x = x, y = y, group = factor(household)),
    fill        = "gray88", color = "gray55",
    alpha       = 0.4, expand = unit(3, "mm"),
    show.legend = FALSE
  ) +
  # Household edges
  {if (!is.null(hh_edge_plot))
    geom_segment(data = hh_edge_plot,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 color = "gray20", alpha = 0.8, linewidth = 0.7)} +
  # Nodes
  geom_point(
    data  = syn_nodes,
    aes(x = x, y = y, fill = factor(age_group)),
    shape = 21, size = 5.0, color = "white", stroke = 0.5
  ) +
  scale_fill_manual(values = age_pal, labels = age_labels_vis, name = "Age group") +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 3))) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 12, face = "bold")
  )

path_net <- file.path(figure_dir, "contact_network_sample.png")
ggsave(path_net, p_net, width = 7, height = 5, dpi = 500)
cat(sprintf("Saved: %s\n", path_net))


# ==============================================================================
# Figure 4 — Prem contact matrices (3-panel: school / work / other)
#             Home excluded (household structure already captured in network)
#             Independent color scale per panel; large text for readability
# ==============================================================================
library(readxl)
library(tidyr)
library(patchwork)   # arrange 3 independent ggplots in one row

prem_dir     <- "data/Prem_contact"
prem_country <- "United States of America"

file_map_3 <- c(
  School           = "MUestimates_school_2.xlsx",
  Work             = "MUestimates_work_2.xlsx",
  "Other locations" = "MUestimates_other_locations_2.xlsx"
)

# Age group labels: 0-4y, 5-9y, ..., 75+y
age_labels <- c(paste0(seq(0, 70, 5), "-", seq(4, 74, 5), "y"), "75+y")

# Helper: load one matrix -> long data frame
load_prem_long <- function(fpath, comp_label) {
  mat <- as.matrix(read_excel(fpath, sheet = prem_country, col_names = FALSE))
  storage.mode(mat) <- "numeric"
  as.data.frame(mat) %>%
    setNames(age_labels) %>%
    mutate(age_i = factor(age_labels, levels = age_labels)) %>%
    pivot_longer(-age_i, names_to = "age_j", values_to = "contact_rate") %>%
    mutate(age_j = factor(age_j, levels = age_labels),
           component = comp_label)
}

# Helper: build one heatmap panel with its own color scale
make_prem_panel <- function(df, title) {
  ggplot(df, aes(x = age_j, y = age_i)) +
    geom_tile(aes(fill = contact_rate), color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", contact_rate)),
              size = 2.6, color = "gray15") +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      name      = "Mean daily\ncontacts"
    ) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_y_discrete(limits = rev(age_labels)) +   # age 0-4y at top
    labs(title = title,
         x     = "Age of contact",
         y     = "Age of respondent") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.x  = element_text(size = 9),
      axis.text.y  = element_text(size = 9),
      axis.title   = element_text(size = 11),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      panel.grid   = element_blank()
    )
}

# Build 3 independent panels
panels <- mapply(
  function(label, fname) {
    make_prem_panel(
      load_prem_long(file.path(prem_dir, fname), label),
      title = label
    )
  },
  label = names(file_map_3),
  fname = unname(file_map_3),
  SIMPLIFY = FALSE
)

# Combine into one row with patchwork
p_prem <- wrap_plots(panels, nrow = 1) +
  plot_annotation(
    # title    = "Prem Contact Matrices — United States of America",
    # subtitle = "Mean number of daily contacts by age group of respondent (row) and contact (column)\nHome matrix excluded — household contacts captured via RTI household data",
    theme    = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10,
                                   margin = margin(b = 10))
    )
  )

path_prem <- file.path(figure_dir, "contact_prem_matrix.png")
ggsave(path_prem, p_prem, width = 19, height = 6, dpi = 500)
cat(sprintf("Saved: %s\n", path_prem))


# ==============================================================================
# Figure 5 — San Francisco age distribution (5-year bands, Prem-aligned)
# ==============================================================================

age_labels_5y <- c(paste0(seq(0, 70, 5), "-", seq(4, 74, 5)), "75+")

age_dist <- net_data$nodes %>%
  mutate(age_group_5y = factor(
    age_labels_5y[pmin(floor(agep / 5) + 1L, 16L)],
    levels = age_labels_5y
  )) %>%
  count(age_group_5y, name = "n") %>%
  mutate(pct = n / sum(n) * 100)

p_age <- ggplot(age_dist, aes(x = age_group_5y, y = pct)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            vjust = -0.4, size = 3.0, color = "gray30") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    # title    = "San Francisco Population Age Distribution (2019)",
    # subtitle = sprintf("RTI synthetic population · N = %s people",
    #                    format(sum(age_dist$n), big.mark = ",")),
    x = "Age group",
    y = "Share of population (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray45", size = 9,
                                 margin = margin(b = 8)),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

path_age <- file.path(figure_dir, "SF_age_dist.png")
ggsave(path_age, p_age, width = 7, height = 5, dpi = 500)
cat(sprintf("Saved: %s\n", path_age))