# Stage 1: Descriptive Summary Analysis
# EW Guinean Lowland Forests, maxdist=2000m, maxtime=2d
# Reads interpolated TIFs from Apr25_Tri2nb_Run_nh_2000_2/YYYY/ folders
# Run this script from inside the 2000,2 directory
#setwd("/Users/sri/Downloads/Africa Forest Fires 2/2000,2")

#install.packages(c("terra", "sf", "tidyverse", "tidyterra", "scales", "patchwork", "viridis"))
#install.packages("dplyr")
#install.packages("tidyterra")
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(scales)
library(patchwork)
library(viridis)

# ---- Paths and setup --------------------------------------------------------

base_dir <- "./Apr25_Tri2nb_Run_nh_2000_2"
out_dir  <- "./2000,2_analysis_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

years <- 2013:2024

# Read the ecoregion boundary for map overlays
eco <- st_read("EW_Guinean_Lowland_Forests.shp", quiet = TRUE) %>%
  st_transform(st_crs("ESRI:102022"))

# Pixel area in hectares (375m resolution, so each pixel is 375x375 m)
pixel_ha <- (375 * 375) / 10000

# ---- Step 1: Stack all interpolated DOY and FRP rasters ----------------------

# We use the interpolated layers because they represent full burned area
# including the gaps filled inside fire polygons
cat("Loading raster stacks...\n")

doy_files <- file.path(base_dir, years, paste0("ugfinterpdoy_", years, ".tif"))
frp_files <- file.path(base_dir, years, paste0("ugfinterpfrp_", years, ".tif"))

doy_stack <- rast(doy_files)
frp_stack <- rast(frp_files)

names(doy_stack) <- as.character(years)
names(frp_stack) <- as.character(years)

# ---- Step 2: Build the annual summary table ----------------------------------

cat("Computing annual statistics...\n")

summary_rows <- lapply(seq_along(years), function(i) {
  yr       <- years[i]
  doy_lyr  <- doy_stack[[i]]
  frp_lyr  <- frp_stack[[i]]
  
  doy_vals <- values(doy_lyr, na.rm = TRUE)
  frp_vals <- values(frp_lyr, na.rm = TRUE)
  
  # Burned area is just the count of non-NA pixels times pixel size
  n_burned    <- length(doy_vals)
  burned_ha   <- n_burned * pixel_ha
  
  tibble(
    year          = yr,
    burned_ha     = round(burned_ha, 1),
    n_fire_pixels = n_burned,
    mean_doy      = round(mean(doy_vals),   1),
    median_doy    = round(median(doy_vals), 1),
    sd_doy        = round(sd(doy_vals),     1),
    earliest_doy  = round(min(doy_vals),    0),
    latest_doy    = round(max(doy_vals),    0),
    mean_frp_mw   = round(mean(frp_vals),   2),
    median_frp_mw = round(median(frp_vals), 2),
    total_frp_mw  = round(sum(frp_vals),    1),
    max_frp_mw    = round(max(frp_vals),    2)
  )
})

summary_df <- bind_rows(summary_rows)

# Save the summary table
write_csv(summary_df, file.path(out_dir, "annual_summary_2000_2.csv"))
cat("Summary table saved.\n")
print(summary_df)

# ---- Step 3: Map 1 - Mean burn DOY across all years -------------------------
# This shows WHERE in the ecoregion fires happen earliest vs latest on average

cat("Making mean DOY map...\n")

mean_doy <- mean(doy_stack, na.rm = TRUE)

p_mean_doy <- ggplot() +
  geom_spatraster(data = mean_doy) +
  geom_sf(data = eco, fill = NA, color = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    option    = "inferno",
    name      = "Mean burn\nday of year",
    na.value  = "transparent",
    direction = -1
  ) +
  labs(
    title    = "Mean Burn Day of Year (2013 to 2024)",
    subtitle = "EW Guinean Lowland Forests | maxdist=2000m",
    caption  = "Interpolated VIIRS SNPP active fire detections"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption  = element_text(size = 9,  hjust = 1,   color = "grey50"),
    legend.position      = "right",
    legend.title         = element_text(size = 10),
    plot.background      = element_rect(fill = "white", color = NA),
    plot.margin          = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "map_mean_burn_doy.png"),
       p_mean_doy, width = 10, height = 7, dpi = 200)

# ---- Step 4: Map 2 - Mean FRP across all years ------------------------------
# Brighter = more intense fires on average at that location

cat("Making mean FRP map...\n")

mean_frp <- mean(frp_stack, na.rm = TRUE)

p_mean_frp <- ggplot() +
  geom_spatraster(data = mean_frp) +
  geom_sf(data = eco, fill = NA, color = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    option   = "plasma",
    name     = "Mean FRP\n(MW)",
    na.value = "transparent"
  ) +
  labs(
    title    = "Mean Fire Radiative Power (2013 to 2024)",
    subtitle = "EW Guinean Lowland Forests | maxdist=2000m",
    caption  = "Interpolated VIIRS SNPP active fire detections"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption  = element_text(size = 9,  hjust = 1,   color = "grey50"),
    legend.position      = "right",
    legend.title         = element_text(size = 10),
    plot.background      = element_rect(fill = "white", color = NA),
    plot.margin          = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "map_mean_frp.png"),
       p_mean_frp, width = 10, height = 7, dpi = 200)

# ---- Step 5: Map 3 - Total burned area per pixel (fire frequency proxy) -----
# How many years did each pixel burn? Warm colors = burned more often

cat("Making fire frequency map...\n")

# Count how many years had a valid (non-NA) DOY value at each pixel
n_years_burned <- app(doy_stack, fun = function(x) sum(!is.na(x)))

p_freq <- ggplot() +
  geom_spatraster(data = n_years_burned) +
  geom_sf(data = eco, fill = NA, color = "white", linewidth = 0.6) +
  scale_fill_viridis_c(
    option   = "turbo",
    name     = "Years\nburned",
    na.value = "transparent",
    breaks   = 1:12,
    limits   = c(1, 12)
  ) +
  labs(
    title    = "Number of Years Burned per Pixel (2013 to 2024)",
    subtitle = "EW Guinean Lowland Forests | maxdist=2000m",
    caption  = "Based on interpolated VIIRS fire detections"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption  = element_text(size = 9,  hjust = 1,   color = "grey50"),
    legend.position      = "right",
    legend.title         = element_text(size = 10),
    plot.background      = element_rect(fill = "white", color = NA),
    plot.margin          = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "map_fire_frequency.png"),
       p_freq, width = 10, height = 7, dpi = 200)

# ---- Step 6: Map 4 - Small multiples of burned area by year -----------------
# One map per year so you can see spatial patterns shifting over time

cat("Making annual facet maps...\n")

# We will plot the DOY raster for each year as a faceted panel
# Stack already has names = year labels so tidyterra picks them up nicely
p_facet <- ggplot() +
  geom_spatraster(data = doy_stack) +
  geom_sf(data = eco, fill = NA, color = "white", linewidth = 0.4) +
  scale_fill_viridis_c(
    option   = "inferno",
    name     = "Burn DOY",
    na.value = "transparent",
    direction = -1
  ) +
  facet_wrap(~lyr, ncol = 4) +
  labs(
    title    = "Annual Burn Day of Year Maps (2013 to 2024)",
    subtitle = "EW Guinean Lowland Forests | maxdist=2000m",
    caption  = "Interpolated VIIRS SNPP active fire detections"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "grey40"),
    plot.caption     = element_text(size = 8,  hjust = 1,   color = "grey50"),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "map_annual_doy_facets.png"),
       p_facet, width = 16, height = 12, dpi = 200)

# ---- Step 7: Chart 1 - Burned area over time --------------------------------
# A clean line chart showing how total burned area changed year to year

cat("Making burned area trend chart...\n")

p_ba <- ggplot(summary_df, aes(x = year, y = burned_ha)) +
  geom_line(color = "#E84855", linewidth = 1.2) +
  geom_point(color = "#E84855", size = 3.5, shape = 21,
             fill = "white", stroke = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF6B35",
              fill = "#FF6B35", alpha = 0.15, linewidth = 0.8) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Total Burned Area per Year (2013 to 2024)",
    subtitle = "EW Guinean Lowland Forests | maxdist=2000m",
    x        = "Year",
    y        = "Burned area (ha)",
    caption  = "Interpolated VIIRS SNPP | shaded band is 95% confidence interval"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption    = element_text(size = 9,  hjust = 1,   color = "grey50"),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(10, 15, 10, 10)
  )

ggsave(file.path(out_dir, "chart_burned_area_trend.png"),
       p_ba, width = 10, height = 6, dpi = 200)

# ---- Step 8: Chart 2 - Mean burn DOY over time ------------------------------
# Are fires happening earlier or later in the season?

cat("Making burn DOY trend chart...\n")

p_doy <- ggplot(summary_df, aes(x = year, y = mean_doy)) +
  geom_ribbon(aes(ymin = mean_doy - sd_doy,
                  ymax = mean_doy + sd_doy),
              fill = "#4CC9F0", alpha = 0.25) +
  geom_line(color = "#4361EE", linewidth = 1.2) +
  geom_point(color = "#4361EE", size = 3.5, shape = 21,
             fill = "white", stroke = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "#7209B7",
              linewidth = 0.8, linetype = "dashed") +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(
    breaks = seq(0, 365, by = 30),
    labels = function(x) {
      months <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")
      m <- ceiling(x / 30.4)
      ifelse(m >= 1 & m <= 12, months[m], as.character(x))
    }
  ) +
  labs(
    title    = "Mean Burn Day of Year Over Time (2013 to 2024)",
    subtitle = "Shaded band shows plus/minus one standard deviation across pixels",
    x        = "Year",
    y        = "Mean burn day of year",
    caption  = "EW Guinean Lowland Forests | maxdist=2000m | dashed line is linear trend"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption     = element_text(size = 9,  hjust = 1,   color = "grey50"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(10, 15, 10, 10)
  )

ggsave(file.path(out_dir, "chart_mean_doy_trend.png"),
       p_doy, width = 10, height = 6, dpi = 200)

# ---- Step 9: Chart 3 - Mean FRP over time -----------------------------------
# Is fire intensity going up or down?

cat("Making FRP trend chart...\n")

p_frp <- ggplot(summary_df, aes(x = year, y = mean_frp_mw)) +
  geom_line(color = "#F72585", linewidth = 1.2) +
  geom_point(color = "#F72585", size = 3.5, shape = 21,
             fill = "white", stroke = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#FF9F1C",
              fill = "#FF9F1C", alpha = 0.15, linewidth = 0.8) +
  scale_x_continuous(breaks = years) +
  labs(
    title    = "Mean Fire Radiative Power Over Time (2013 to 2024)",
    subtitle = "EW Guinean Lowland Forests | maxdist=2000m",
    x        = "Year",
    y        = "Mean FRP (MW)",
    caption  = "Interpolated VIIRS SNPP | shaded band is 95% confidence interval"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption     = element_text(size = 9,  hjust = 1,   color = "grey50"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(10, 15, 10, 10)
  )

ggsave(file.path(out_dir, "chart_mean_frp_trend.png"),
       p_frp, width = 10, height = 6, dpi = 200)

# ---- Step 10: Combined summary panel ----------------------------------------
# One image with all three trend charts stacked so you can compare them easily

cat("Making combined summary panel...\n")

p_combined <- p_ba / p_doy / p_frp +
  plot_annotation(
    title   = "Stage 1 Summary: EW Guinean Lowland Forests (2013 to 2024)",
    caption = "maxdist=2000m, maxtime=2d | Interpolated VIIRS SNPP detections",
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 1, color = "grey50"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave(file.path(out_dir, "chart_combined_summary_panel.png"),
       p_combined, width = 11, height = 16, dpi = 200)

# ---- Done -------------------------------------------------------------------

cat("\nStage 1 complete! All outputs saved to:", out_dir, "\n\n")
cat("Files produced:\n")
cat("  annual_summary_2000_2.csv         (the full stats table)\n")
cat("  map_mean_burn_doy.png             (where fires burn earliest/latest)\n")
cat("  map_mean_frp.png                  (where fires are most intense)\n")
cat("  map_fire_frequency.png            (how many years each pixel burned)\n")
cat("  map_annual_doy_facets.png         (one map per year side by side)\n")
cat("  chart_burned_area_trend.png       (total burned area over time)\n")
cat("  chart_mean_doy_trend.png          (fire timing shift over time)\n")
cat("  chart_mean_frp_trend.png          (fire intensity over time)\n")
cat("  chart_combined_summary_panel.png  (all three charts in one image)\n")