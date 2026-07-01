library(terra)
library(tidyverse)
library(sf)
library(Kendall)

# ---- Step 1: Load reserves — top 6 largest by area ----
reserves <- st_read("reserves new/Reserves_Prof_Dissolved/Ghana_Reserves.shp") %>%
  st_transform(st_crs("ESRI:102022"))

# Calculate area and keep top 6
reserves$area_km2 <- as.numeric(expanse(vect(reserves), unit = "km"))
reserves <- reserves %>%
  arrange(desc(area_km2)) %>%
  slice(1:6)

cat("Top 6 largest reserves:\n")
print(reserves %>% st_drop_geometry() %>% select(Name, area_km2))

# ---- Step 2: Extract fire density per reserve per year ----
years <- 2013:2024
results <- list()

for(yr in years) {
  cat("Processing year:", yr, "\n")
  
  fire_file <- paste0("ugf_activedoy_", yr, ".tif")
  if(!file.exists(fire_file)) {
    cat("  MISSING:", fire_file, "\n")
    next
  }
  
  fire_r <- rast(fire_file)
  
  for(i in 1:nrow(reserves)) {
    res_name <- reserves$Name[i]
    res_vect <- vect(reserves[i, ])
    
    fire_crop    <- crop(fire_r, res_vect)
    fire_mask    <- mask(fire_crop, res_vect)
    fire_count   <- sum(!is.na(values(fire_mask)))
    area_km2     <- as.numeric(expanse(res_vect, unit = "km"))
    fire_density <- (fire_count / area_km2) * 100
    
    results[[length(results) + 1]] <- data.frame(
      year         = yr,
      reserve      = res_name,
      fire_count   = fire_count,
      area_km2     = round(area_km2, 1),
      fire_density = round(fire_density, 3)
    )
  }
}

fire_df <- bind_rows(results)
write.csv(fire_df, "reserve_fire_density_top6.csv", row.names = FALSE)
cat("Data saved to reserve_fire_density_top6.csv\n")
print(fire_df)

# ---- Step 3: Mann-Kendall trend test ----
pvals <- fire_df %>%
  group_by(reserve) %>%
  summarise(
    p_value = round(MannKendall(fire_density)$sl, 3),
    .groups = "drop"
  )

fire_df_plot <- left_join(fire_df, pvals, by = "reserve") %>%
  mutate(p_label = paste0("P = ", sprintf("%.3f", p_value)))

# ---- Step 4: Numeric summary table ----
summary_table <- fire_df_plot %>%
  group_by(reserve) %>%
  summarise(
    area_km2      = round(unique(area_km2), 1),
    mean_density  = round(mean(fire_density), 3),
    max_density   = round(max(fire_density), 3),
    min_density   = round(min(fire_density), 3),
    trend_p_value = unique(p_value),
    significant   = ifelse(unique(p_value) < 0.05, "YES", "no"),
    .groups = "drop"
  ) %>%
  arrange(desc(area_km2))

cat("\n=== Summary Table ===\n")
print(summary_table)
write.csv(summary_table, "reserve_trend_summary_top6.csv", row.names = FALSE)

# ---- Step 5: Plot — 2 rows x 3 columns ----
ggplot(fire_df_plot, aes(x = year, y = fire_density)) +
  geom_point(size = 2, color = "black") +
  geom_smooth(method = "loess", color = "steelblue",
              fill = "grey70", alpha = 0.4, linewidth = 1) +
  geom_text(data = fire_df_plot %>% distinct(reserve, p_label),
            aes(label = p_label),
            x = -Inf, y = Inf,
            hjust = -0.1, vjust = 1.4,
            size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ reserve, scales = "free_y", ncol = 3) +
  labs(title    = "Annual Fire Density — Top 6 Largest Forest Reserves, Ghana (2013–2024)",
       subtitle = "Blue line = LOESS trend | P-value from Mann-Kendall trend test",
       x        = NULL,
       y        = "Fires / 100 km²") +
  scale_x_continuous(breaks = c(2013, 2016, 2019, 2022, 2024)) +
  theme_bw() +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
        strip.text       = element_text(face = "bold", size = 9),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y      = element_text(size = 8),
        panel.grid.minor = element_blank())

ggsave("reserve_fire_trends_top6.png", width = 12, height = 7, 
       dpi = 300, limitsize = FALSE)
cat("Plot saved to reserve_fire_trends_top6.png\n")