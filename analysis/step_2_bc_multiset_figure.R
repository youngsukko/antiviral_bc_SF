# ==============================================================================
# scp_plot_heatmap.R
#
# Figure 1: 4-panel heatmap (one panel per intervention level),
#           delay averaged across all four delay conditions.
#           X-axis : intervention start day (0, 7, 14, 21, 28)
#           Y-axis : antiviral efficacy multiplier (0.1 → 1.0, ascending)
#           Fill   : mean_infected (viridis)
#           Label  : P(epidemic) — red + bold if > 25%, white/black by luminance
#
# Figure 2: Grouped bar chart of P(epidemic) by delay × intervention,
#           averaged across e_max_mult and intervention_start.
#
# Typography: all text uses sans-serif at base_size = 9pt.
#             geom_text size = 3.2 (mm units, ≈ 9pt) to match theme base_size.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)

summary_path <- "output/scenario_summary.csv"
figure_dir   <- "figures"
if (!dir.exists(figure_dir)) dir.create(figure_dir)

BASE_SIZE  <- 9      # pt — applies to all theme text
GEOM_SIZE  <- 3.2    # mm — geom_text/annotate, equivalent to ~9pt


# ==============================================================================
# [Section 1] Load data
# ==============================================================================

df <- read.csv(summary_path, stringsAsFactors = FALSE) %>%
  mutate(
    intervention = recode(intervention, "Mild" = "Minimum"),
    intervention = factor(intervention,
                          levels = c("Minimum", "Moderate", "Intensive", "Maximum")),
    e_max_mult   = round(e_max_mult, 1),
    delay_label  = factor(paste0("0\u2013", delay_max, "d"),
                          levels = paste0("0\u2013", c(1,2,3,4), "d"))
  )


# ==============================================================================
# [Section 2] Figure 1 — Heatmap (delay-averaged)
# ==============================================================================

hm_df <- df %>%
  group_by(intervention, antiviral_start, e_max_mult) %>%
  summarise(
    mean_infected = mean(mean_infected),
    p_epidemic    = mean(p_epidemic),
    .groups       = "drop"
  ) %>%
  mutate(
    e_max_fac = factor(e_max_mult, levels = sort(unique(e_max_mult))),
    av_start  = factor(antiviral_start, levels = sort(unique(antiviral_start))),
    pep_label = paste0(round(p_epidemic * 100), "%")
  )

fill_range <- range(hm_df$mean_infected)

hm_df <- hm_df %>%
  mutate(
    fill_norm  = (mean_infected - fill_range[1]) / (fill_range[2] - fill_range[1]),
    text_color = case_when(
      p_epidemic > 0.25 ~ "#CC2222",
      fill_norm  < 0.55 ~ "white",
      TRUE              ~ "black"
    )
  )

p1 <- ggplot(hm_df, aes(x = av_start, y = e_max_fac)) +
  geom_tile(aes(fill = mean_infected), color = NA) +
  geom_text(aes(label    = pep_label,
                color    = text_color,
                fontface = ifelse(p_epidemic > 0.25, "bold", "plain")),
            size   = GEOM_SIZE,
            family = "sans") +
  scale_fill_viridis_c(
    name   = "Mean infections",
    option = "viridis",
    guide  = guide_colorbar(
      title.position = "right",
      title.hjust    = 0.5,
      barwidth       = unit(6, "pt"),
      barheight      = unit(80,   "pt"),
      ticks          = TRUE,
      ticks.colour   = "white"
    )
  ) +
  scale_color_identity() +
  facet_wrap(~ intervention, nrow = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Intervention start timing since index onset",
    y = "Antiviral efficacy multiplier"
  ) +
  theme_minimal(base_size = BASE_SIZE, base_family = "sans") +
  theme(
    panel.spacing     = unit(3, "pt"),
    panel.border      = element_rect(color = "grey75", fill = NA, linewidth = 0.3),
    panel.grid        = element_blank(),
    strip.text        = element_text(size = BASE_SIZE, face = "bold",
                                     family = "sans", margin = margin(b = 4)),
    axis.text.x       = element_text(size = BASE_SIZE, color = "grey40", family = "sans"),
    axis.text.y       = element_text(size = BASE_SIZE, color = "grey40", family = "sans"),
    axis.title.x      = element_text(size = BASE_SIZE, family = "sans", margin = margin(t = 5)),
    axis.title.y      = element_text(size = BASE_SIZE, family = "sans", margin = margin(r = 5)),
    axis.ticks        = element_blank(),
    legend.position   = "right",
    legend.title      = element_text(size = BASE_SIZE, family = "sans", angle = 270, hjust = 0.5),
    legend.text       = element_text(size = BASE_SIZE, family = "sans"),
    legend.margin     = margin(l = 6),
    plot.margin       = margin(6, 6, 0, 6),
    aspect.ratio = 2
  )


# ==============================================================================
# [Section 3] Figure 2 — Bar chart: P(epidemic) by delay x intervention
# ==============================================================================

bar_df <- df %>%
  group_by(intervention, delay_label) %>%
  summarise(
    p_epidemic = mean(p_epidemic),
    .groups    = "drop"
  )

delay_colors <- c("0\u20131d" = "#B5D4F4",
                  "0\u20132d" = "#378ADD",
                  "0\u20133d" = "#185FA5",
                  "0\u20134d" = "#0C447C")

p2 <- ggplot(bar_df, aes(x = intervention, y = p_epidemic, fill = delay_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.9) +
  geom_hline(yintercept = 0.25, linetype = "dashed",
             color = "#CC2222", linewidth = 0.5) +
  annotate("text", x = 0.55, y = 0.265, label = "25%",
           size = GEOM_SIZE, family = "sans", color = "#CC2222", hjust = 0) +
  scale_fill_manual(values = delay_colors, name = "Logistical delay") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = NULL, y = "P(epidemic)") +
  theme_minimal(base_size = BASE_SIZE, base_family = "sans") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    axis.text.x        = element_text(size = BASE_SIZE, color = "grey30", family = "sans"),
    axis.text.y        = element_text(size = BASE_SIZE, color = "grey40", family = "sans"),
    axis.title.y       = element_text(size = BASE_SIZE, family = "sans", margin = margin(r = 5)),
    axis.ticks         = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = BASE_SIZE, family = "sans"),
    legend.text        = element_text(size = BASE_SIZE, family = "sans"),
    legend.key.size    = unit(8, "pt"),
    plot.margin        = margin(0, 6, 6, 6)
  )


# ==============================================================================
# [Section 4] Combine and save
# ==============================================================================

combined <- p1 / p2 +
  plot_layout(heights = c(2.8, 1)) +
  plot_annotation(
    theme = theme(plot.margin = margin(8, 8, 8, 8))
  )

path_out <- file.path(figure_dir, "heatmap_scenario_grid.png")
ggsave(path_out, combined,
       width  = 10,
       height = 8,
       dpi    = 180)

cat(sprintf("Saved: %s\n", path_out))