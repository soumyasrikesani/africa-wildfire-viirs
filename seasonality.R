# =============================================================================
# Time Series EDA — VIIRS Active Fire Data (2012–2024)
# 5 plots displayed separately, all 13 years present in each
# Purpose: Establish temporal baselines that justify future validation processes
# Requires: bbox already loaded (Sections 1–2 of main EDA script)
# =============================================================================

library(data.table)
library(ggplot2)
library(scales)
library(lubridate)

# -----------------------------------------------------------------------------
# STEP 1 — Load all 13 years, bbox filter, extract monthly summaries
# Only reads what is needed: date, frp, confidence, daynight
# One year at a time — raw table freed after each year
# -----------------------------------------------------------------------------

years        <- 2013:2024
monthly_list <- list()
annual_list  <- list()

for (yr in years) {
  
  infile <- paste0("africa_fires_", yr, ".csv")
  
  if (!file.exists(infile)) {
    cat(sprintf("  [SKIP] %s not found\n", infile))
    next
  }
  
  cat(sprintf("Reading %d ...\n", yr))
  
  raw <- fread(infile,
               select = c("latitude", "longitude", "frp",
                          "confidence", "daynight", "acq_date"))
  
  # Bbox filter — numeric, before any sf object is created
  sub <- raw[longitude >= bbox["xmin"] & longitude <= bbox["xmax"] &
               latitude  >= bbox["ymin"] & latitude  <= bbox["ymax"]]
  
  # Keep all confidence levels consistent with main script
  sub <- sub[confidence %in% c("h", "n", "l")]
  
  # Parse dates
  sub[, acq_date := as.Date(acq_date)]
  sub[, year     := year(acq_date)]
  sub[, month    := lubridate::month(acq_date)]
  sub[, yearmon  := as.Date(paste0(yr, "-", month, "-01"))]
  
  # Monthly summary per year
  monthly_sum <- sub[, .(
    n_fires    = .N,
    frp_median = median(frp, na.rm = TRUE),
    frp_mean   = mean(frp,   na.rm = TRUE),
    frp_total  = sum(frp,    na.rm = TRUE),
    pct_high   = round(100 * sum(confidence == "h") / .N, 1)
  ), by = .(year, month, yearmon)][order(yearmon)]
  
  monthly_list[[as.character(yr)]] <- monthly_sum
  
  # Annual summary
  annual_list[[as.character(yr)]] <- data.table(
    year       = yr,
    n_total    = nrow(sub),
    frp_median = median(sub$frp, na.rm = TRUE),
    frp_mean   = mean(sub$frp,   na.rm = TRUE),
    frp_total  = sum(sub$frp,    na.rm = TRUE),
    pct_high   = round(100 * sum(sub$confidence == "h") / nrow(sub), 1)
  )
  
  rm(raw, sub); gc()
}

# Combine
monthly_all <- rbindlist(monthly_list)
annual_all  <- rbindlist(annual_list)

# Z-score and anomaly flag on annual counts — used across multiple plots
annual_all[, zscore  := round((n_total - mean(n_total)) / sd(n_total), 2)]
annual_all[, anomaly := fcase(
  zscore >  1,  "Above average",
  zscore < -1,  "Below average",
  default =     "Normal"
)]

# Ensure month is ordered factor for clean x-axis labels
month_labels <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")
monthly_all[, month_label := factor(month_labels[month], levels = month_labels)]

cat("\n--- Annual summary table ---\n")
print(annual_all)
cat("\n--- Monthly table (first 24 rows) ---\n")
print(head(monthly_all, 24))

# =============================================================================
# PLOT 1 — Monthly fire count time series (line, all 13 years as one timeline)
# Shows the raw temporal signal: seasonality + inter-annual variation together
# Validation purpose: establishes the expected seasonal cycle baseline
# =============================================================================

p1 <- ggplot(monthly_all, aes(x = yearmon, y = n_fires)) +
  
  # Shaded fire season bands — highlight peak months (Jul–Oct typical W Africa)
  # Adjust month range if your data shows a different peak
  geom_rect(
    data = data.frame(
      xmin = as.Date(paste0(years, "-07-01")),
      xmax = as.Date(paste0(years, "-10-31")),
      ymin = -Inf, ymax = Inf
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#fff3cd", alpha = 0.4
  ) +
  
  geom_line(color = "#3a7abf", linewidth = 0.7, alpha = 0.9) +
  geom_point(color = "#3a7abf", size = 0.8, alpha = 0.7) +
  
  # Annotate each year on the x-axis for easy reading
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand      = expansion(mult = 0.01)
  ) +
  scale_y_continuous(labels = comma) +
  
  labs(
    title    = "Plot 1 — Monthly Fire Detection Count (2012–2024)",
    subtitle = "Each point = one month  |  Yellow bands = typical fire season (Jul–Oct)\nValidation: baseline seasonal cycle for comparing model outputs",
    x        = "Date",
    y        = "Monthly detections",
    caption  = "VIIRS SNPP  |  Confidence: h + n + l"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p1)
cat("Plot 1 displayed: Monthly fire count time series\n")

# =============================================================================
# PLOT 2 — Heatmap: Year × Month fire count
# The single most informative temporal EDA for multi-year fire data
# Rows = years, columns = months, fill = detection count
# Validation purpose: shows whether the fire season timing is stable enough
# to use a fixed seasonal window in validation or anomaly detection
# =============================================================================

p2 <- ggplot(monthly_all,
             aes(x = month_label, y = factor(year), fill = n_fires)) +
  
  geom_tile(color = "white", linewidth = 0.5) +
  
  # Label each cell with count — helps spot data gaps (cells near 0)
  geom_text(aes(label = ifelse(n_fires > 0, comma(n_fires), "—")),
            size = 2.5, color = "grey20") +
  
  scale_fill_gradient2(
    low      = "#f7fbff",
    mid      = "#6baed6",
    high     = "#c0392b",
    midpoint = median(monthly_all$n_fires, na.rm = TRUE),
    labels   = comma,
    name     = "Detections"
  ) +
  
  scale_y_discrete(limits = rev) +   # most recent year at top
  
  labs(
    title    = "Plot 2 — Fire Detection Heatmap: Year × Month (2012–2024)",
    subtitle = "Red = peak activity  |  White = low/no activity\nValidation: confirms fire season window; '—' or near-zero cells = potential data gaps",
    x        = "Month",
    y        = "Year",
    caption  = "VIIRS SNPP  |  Confidence: h + n + l"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text     = element_text(size = 9),
    legend.position = "right",
    panel.grid    = element_blank()
  )

print(p2)
cat("Plot 2 displayed: Year x Month heatmap\n")

# =============================================================================
# PLOT 3 — Annual median FRP trend (line + points + smoothed trend)
# Shows fire INTENSITY over time, independent of detection count
# Validation purpose: FRP is the interpolated variable in the main pipeline —
# a stable FRP baseline justifies using multi-year mean as a reference surface
# for validating interpolated rasters
# =============================================================================

p3 <- ggplot(annual_all, aes(x = year, y = frp_median)) +
  
  # 13-year median reference line
  geom_hline(yintercept = median(annual_all$frp_median),
             linetype = "dashed", color = "grey50", linewidth = 0.7) +
  
  annotate("text",
           x     = min(annual_all$year),
           y     = median(annual_all$frp_median) * 1.015,
           label = paste0("13-yr median FRP = ",
                          round(median(annual_all$frp_median), 1), " MW"),
           hjust = 0, size = 3, color = "grey40") +
  
  # Smoothed trend
  geom_smooth(method  = "loess", span = 0.7,
              color   = "#c0392b", fill = "#f5b7b1",
              linewidth = 1, alpha = 0.25, se = TRUE) +
  
  # Raw annual values
  geom_line(color = "#3a7abf", linewidth = 0.8) +
  geom_point(aes(size = n_total, color = anomaly), alpha = 0.9) +
  
  scale_color_manual(
    values = c("Above average" = "#c0392b",
               "Normal"        = "#3a7abf",
               "Below average" = "#2d8e5e"),
    name   = "Count anomaly"
  ) +
  scale_size_continuous(range = c(2, 7), labels = comma,
                        name  = "Annual detections") +
  scale_x_continuous(breaks = years) +
  
  labs(
    title    = "Plot 3 — Annual Median FRP Trend (2012–2024)",
    subtitle = "Point size = total detections  |  Red LOESS = smoothed intensity trend\nValidation: stable FRP baseline justifies using multi-year mean as reference surface",
    x        = "Year",
    y        = "Median FRP (MW)",
    caption  = "VIIRS SNPP  |  Confidence: h + n + l"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

print(p3)
cat("Plot 3 displayed: Annual median FRP trend\n")

# =============================================================================
# PLOT 4 — Seasonal profile ribbons: all years overlaid on one calendar
# X-axis = month (1–12), Y-axis = detection count
# Each year is one line; the ribbon shows min–max range across all years
# Validation purpose: the ribbon width shows how much year-to-year variation
# exists within each month — wide ribbons = high inter-annual variability =
# validation metrics should account for seasonal timing uncertainty
# =============================================================================

# Min/max ribbon across all years per month
ribbon <- monthly_all[, .(
  ymin   = min(n_fires,    na.rm = TRUE),
  ymax   = max(n_fires,    na.rm = TRUE),
  ymean  = mean(n_fires,   na.rm = TRUE)
), by = month][order(month)]

ribbon[, month_label := factor(month_labels[month], levels = month_labels)]

p4 <- ggplot() +
  
  # Min–max ribbon across all years
  geom_ribbon(data = ribbon,
              aes(x = as.numeric(month_label),
                  ymin = ymin, ymax = ymax),
              fill = "#aed6f1", alpha = 0.4) +
  
  # Mean line
  geom_line(data = ribbon,
            aes(x = as.numeric(month_label), y = ymean),
            color = "#2980b9", linewidth = 1.2, linetype = "solid") +
  
  # Individual year lines — thin, coloured by year
  geom_line(data = monthly_all,
            aes(x     = as.numeric(month_label),
                y     = n_fires,
                color = factor(year),
                group = factor(year)),
            linewidth = 0.5, alpha = 0.65) +
  
  scale_color_manual(
    values = setNames(
      colorRampPalette(c("#1a5276", "#2e86c1", "#85c1e9",
                         "#f39c12", "#e74c3c", "#922b21"))(length(years)),
      as.character(years)
    ),
    name = "Year"
  ) +
  
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_y_continuous(labels = comma) +
  
  labs(
    title    = "Plot 4 — Seasonal Fire Profile: All Years Overlaid (2012–2024)",
    subtitle = "Thin lines = individual years  |  Blue line = 13-year monthly mean\nBlue ribbon = min–max range across all years per month\nValidation: wide ribbon months have high inter-annual variability",
    x        = "Month",
    y        = "Monthly detections",
    caption  = "VIIRS SNPP  |  Confidence: h + n + l"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

print(p4)
cat("Plot 4 displayed: Seasonal profile ribbons — all years overlaid\n")

# =============================================================================
# PLOT 5 — Cumulative fire detections per year (all years on one chart)
# X-axis = day of year (1–365), Y-axis = cumulative count
# Each year is one line
# Validation purpose: lines that diverge early indicate the fire season started
# earlier or later than usual — critical for validating whether burned area
# polygons from the pipeline correctly capture the full fire season
# A line that plateaus early = fire season ended early or data gap
# =============================================================================

# Rebuild from scratch — need day-of-year and cumulative count per year
# Load coordinates again — just date and year this time
cumul_list <- list()

for (yr in years) {
  
  infile <- paste0("africa_fires_", yr, ".csv")
  if (!file.exists(infile)) next
  
  raw <- fread(infile, select = c("latitude", "longitude",
                                  "acq_date", "confidence"))
  
  sub <- raw[longitude >= bbox["xmin"] & longitude <= bbox["xmax"] &
               latitude  >= bbox["ymin"] & latitude  <= bbox["ymax"]]
  sub <- sub[confidence %in% c("h", "n", "l")]
  
  sub[, acq_date := as.Date(acq_date)]
  sub[, doy      := yday(acq_date)]
  
  # Daily counts
  daily <- sub[, .N, by = doy][order(doy)]
  
  # Fill missing DOYs with 0 so cumsum is continuous
  all_doys   <- data.table(doy = 1:365)
  daily_full <- merge(all_doys, daily, by = "doy", all.x = TRUE)
  daily_full[is.na(N), N := 0]
  
  daily_full[, cumulative := cumsum(N)]
  daily_full[, year        := yr]
  
  cumul_list[[as.character(yr)]] <- daily_full
  rm(raw, sub, daily, daily_full, all_doys); gc()
}

cumul_all <- rbindlist(cumul_list)

p5 <- ggplot(cumul_all,
             aes(x     = doy,
                 y     = cumulative,
                 color = factor(year),
                 group = factor(year))) +
  
  geom_line(linewidth = 0.75, alpha = 0.85) +
  
  # Mark DOY 180 (late June) and DOY 300 (late October) as season markers
  geom_vline(xintercept = c(180, 300),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  
  annotate("text", x = 181, y = Inf,
           label = "DOY 180\n(late Jun)",
           hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
  
  annotate("text", x = 301, y = Inf,
           label = "DOY 300\n(late Oct)",
           hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
  
  scale_color_manual(
    values = setNames(
      colorRampPalette(c("#1a5276", "#2e86c1", "#85c1e9",
                         "#f39c12", "#e74c3c", "#922b21"))(length(years)),
      as.character(years)
    ),
    name = "Year"
  ) +
  
  scale_x_continuous(breaks = c(1, 60, 120, 180, 240, 300, 365),
                     labels = c("Jan 1","Mar 1","May 1",
                                "Jun 29","Aug 28","Oct 27","Dec 31")) +
  scale_y_continuous(labels = comma) +
  
  labs(
    title    = "Plot 5 — Cumulative Fire Detections by Day of Year (2012–2024)",
    subtitle = "Each line = one year  |  Steep slope = active fire season\nA line that plateaus early = short season or data gap\nValidation: divergence between lines reveals inter-annual timing differences",
    x        = "Day of Year",
    y        = "Cumulative detections",
    caption  = "VIIRS SNPP  |  Confidence: h + n + l"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

print(p5)
cat("Plot 5 displayed: Cumulative fire detections by DOY\n")

cat("\n=== ALL 5 TIME SERIES PLOTS DISPLAYED ===\n")
cat("Use the Plots pane arrows to navigate between plots.\n")
cat("\nPlot summary:\n")
cat("  Plot 1 — Monthly count time series     : raw temporal signal 2012-2024\n")
cat("  Plot 2 — Year x Month heatmap          : fire season stability\n")
cat("  Plot 3 — Annual median FRP trend       : fire intensity over time\n")
cat("  Plot 4 — Seasonal profiles overlaid    : inter-annual variability per month\n")
cat("  Plot 5 — Cumulative DOY curves         : fire season timing & length\n")