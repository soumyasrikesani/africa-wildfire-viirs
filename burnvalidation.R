##############################################################################
# Burn Probability Classification using Random Forest
# Version 3 - Visualizations display in RStudio Plots panel
#
# Inputs needed:
#   - ugf_interpdoy_<year>.tif     (interpolated DOY from Loop #2)
#   - ugf_interpfrp_<year>.tif     (interpolated FRP from Loop #2)
#   - ugf_fires_<year>.shp         (fire clusters from Loop #1)
#   - ghana_burnMask_2016.tif      (exported from QGIS)
#   - ghana_dNBR_2016.tif          (exported from QGIS)
#   - EW_Guinean_Lowland_Forests.shp
##############################################################################

library(terra)
library(sf)
library(tidyverse)
library(ranger)       # install.packages("ranger")
library(caret)        # install.packages("caret")
library(pROC)         # install.packages("pROC")
library(ggplot2)

# ─────────────────────────────────────────────────────────────
# SECTION 0: MEMORY MANAGEMENT
# ─────────────────────────────────────────────────────────────

setwd("/Users/sri/Downloads/Africa Forest Fires 2")

terraOptions(memfrac = 0.5,
             tempdir = "/Users/sri/Downloads/Africa Forest Fires 2/temp")

dir.create("./temp", showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────
# SECTION 1: CONFIGURATION
# ─────────────────────────────────────────────────────────────

geo_name    <- "ugf"
val_year    <- 2016
train_years <- 2013:2015

# ─────────────────────────────────────────────────────────────
# SECTION 2: LOAD REFERENCE DATA
# ─────────────────────────────────────────────────────────────

cat("Loading reference data...\n")

geo_area <- st_read("EW_Guinean_Lowland_Forests.shp") %>%
  st_transform(st_crs("ESRI:102022"))

burn_mask <- rast("ghana_burnMask_2016.tif")
dnbr_rast <- rast("ghana_dNBR_2016.tif")

burn_mask_aea <- project(burn_mask, "ESRI:102022", method = "near")
dnbr_aea      <- project(dnbr_rast, "ESRI:102022", method = "bilinear")

# ─────────────────────────────────────────────────────────────
# MEMORY FIX: Reduce raster resolution
# ─────────────────────────────────────────────────────────────

cat("Original raster size:", ncell(burn_mask_aea), "pixels\n")

burn_mask_aea <- aggregate(burn_mask_aea, fact = 10, fun = "modal")
dnbr_aea      <- aggregate(dnbr_aea,      fact = 10, fun = "mean")

cat("Reduced raster size:", ncell(burn_mask_aea), "pixels\n")

# ─────────────────────────────────────────────────────────────
# SECTION 3: BUILD FEATURE STACK FUNCTION
# Uses interpolated DOY and FRP from viirs_fire_interpolation.R
# ─────────────────────────────────────────────────────────────

build_feature_stack <- function(year, ref_raster) {
  
  doy_file <- paste0(geo_name, "_interpdoy_", 2016, ".tif")
  frp_file <- paste0(geo_name, "_interpfrp_", 2016, ".tif")
  
  if (!file.exists(doy_file) | !file.exists(frp_file)) {
    message("Missing interpolated raster files for year: ", year)
    return(NULL)
  }
  
  doy_rast <- rast(doy_file)
  frp_rast <- rast(frp_file)
  
  doy_rast <- resample(doy_rast, ref_raster, method = "bilinear")
  frp_rast <- resample(frp_rast, ref_raster, method = "bilinear")
  
  shp_file <- paste0(geo_name, "_fires_", year, ".shp")
  
  if (file.exists(shp_file)) {
    fire_pts   <- st_read(shp_file, quiet = TRUE)
    clust_rast <- rasterize(vect(fire_pts), ref_raster,
                            field = "clustsize", fun = "max")
    count_rast <- rasterize(vect(fire_pts), ref_raster,
                            field = "clustsize", fun = "count")
    clust_rast <- resample(clust_rast, ref_raster, method = "bilinear")
    count_rast <- resample(count_rast, ref_raster, method = "bilinear")
  } else {
    clust_rast <- ref_raster * NA
    count_rast <- ref_raster * NA
  }
  
  stack        <- c(doy_rast, frp_rast, clust_rast, count_rast)
  names(stack) <- c("doy", "frp", "clust_size", "fire_count")
  stack        <- subst(stack, NA, 0)
  
  return(stack)
}

# ─────────────────────────────────────────────────────────────
# SECTION 4: PREPARE TRAINING DATA (2013-2015)
# ─────────────────────────────────────────────────────────────

cat("\n[1/5] Building training feature stacks...\n")

train_list <- list()

for (yr in train_years) {
  cat("  Processing year:", yr, "\n")
  feat_stack <- build_feature_stack(yr, burn_mask_aea)
  
  if (!is.null(feat_stack)) {
    feat_df      <- as.data.frame(feat_stack, xy = TRUE, na.rm = FALSE)
    feat_df$year <- yr
    train_list[[length(train_list) + 1]] <- feat_df
  }
  gc()
}

train_features <- bind_rows(train_list)
rm(train_list)
gc()

# ─────────────────────────────────────────────────────────────
# SECTION 5: PREPARE VALIDATION DATA (2016)
# ─────────────────────────────────────────────────────────────

cat("[2/5] Building validation feature stack for", val_year, "...\n")

val_stack   <- build_feature_stack(val_year, burn_mask_aea)
val_df      <- as.data.frame(val_stack, xy = TRUE, na.rm = FALSE)
val_df$year <- val_year

burn_vals        <- as.data.frame(burn_mask_aea, xy = TRUE, na.rm = FALSE)
names(burn_vals)[3] <- "burned"

dnbr_vals        <- as.data.frame(dnbr_aea, xy = TRUE, na.rm = FALSE)
names(dnbr_vals)[3] <- "dnbr"

val_full <- val_df %>%
  left_join(burn_vals, by = c("x", "y")) %>%
  left_join(dnbr_vals, by = c("x", "y")) %>%
  filter(!is.na(burned))

rm(val_stack, val_df, burn_vals, dnbr_vals)
gc()

# ─────────────────────────────────────────────────────────────
# SECTION 6: COMBINE AND CLEAN DATASET
# ─────────────────────────────────────────────────────────────

cat("[3/5] Preparing final dataset...\n")

train_features <- train_features %>%
  mutate(burned = ifelse(frp > 0, 1, 0),
         dnbr   = 0)

all_data <- bind_rows(
  train_features,
  val_full %>% select(x, y, doy, frp, clust_size,
                      fire_count, year, burned, dnbr)
) %>%
  filter(!is.na(doy), !is.na(frp),
         !is.na(clust_size), !is.na(fire_count)) %>%
  mutate(burned = as.factor(burned))

rm(train_features, val_full)
gc()

cat("\nClass distribution:\n")
print(table(all_data$burned))
cat("Burn prevalence:",
    round(mean(all_data$burned == 1) * 100, 2), "%\n\n")

# ─────────────────────────────────────────────────────────────
# SECTION 7: TRAIN / TEST SPLIT
# ─────────────────────────────────────────────────────────────

train_df <- all_data %>% filter(year != val_year)
test_df  <- all_data %>% filter(year == val_year)

set.seed(42)
n_minority <- min(table(train_df$burned))
train_bal  <- train_df %>%
  group_by(burned) %>%
  slice_sample(n = n_minority) %>%
  ungroup()

cat("Balanced training set size:", nrow(train_bal), "\n")
cat("Test set size (2016):",       nrow(test_df),   "\n\n")

rm(all_data, train_df)
gc()

# ─────────────────────────────────────────────────────────────
# SECTION 8: TRAIN RANDOM FOREST MODEL
# ─────────────────────────────────────────────────────────────

cat("[4/5] Training Random Forest model...\n")

rf_model <- ranger(
  formula     = burned ~ doy + frp + clust_size + fire_count + dnbr,
  data        = train_bal,
  num.trees   = 500,
  mtry        = 3,
  importance  = "impurity",
  probability = TRUE,
  seed        = 42
)

cat("OOB Prediction Error:",
    round(rf_model$prediction.error * 100, 2), "%\n\n")

rm(train_bal)
gc()

# ─────────────────────────────────────────────────────────────
# SECTION 9: PREDICT ON 2016 TEST DATA
# ─────────────────────────────────────────────────────────────

cat("[5/5] Predicting and evaluating...\n")

pred_prob  <- predict(rf_model, data = test_df)$predictions[, "1"]
pred_class <- ifelse(pred_prob >= 0.5, 1, 0)

test_df$pred_prob  <- pred_prob
test_df$pred_class <- as.factor(pred_class)

# ─────────────────────────────────────────────────────────────
# SECTION 10: EVALUATION METRICS
# ─────────────────────────────────────────────────────────────

cat("\n===== MODEL EVALUATION =====\n")

cm <- confusionMatrix(
  data      = as.factor(pred_class),
  reference = test_df$burned,
  positive  = "1"
)
print(cm)

roc_obj <- roc(as.numeric(test_df$burned) - 1, pred_prob, quiet = TRUE)
cat("\nAUC:", round(auc(roc_obj), 4), "\n")

importance_df <- data.frame(
  variable   = names(rf_model$variable.importance),
  importance = rf_model$variable.importance
) %>% arrange(desc(importance))

cat("\nVariable Importance:\n")
print(importance_df)

# ─────────────────────────────────────────────────────────────
# SECTION 11: VISUALIZATIONS (displays in RStudio Plots panel)
# ─────────────────────────────────────────────────────────────

# --- Plot 1: Variable Importance ---
p1 <- ggplot(importance_df,
             aes(x = reorder(variable, importance),
                 y = importance, fill = importance)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_gradient(low = "#fee8c8", high = "#d7301f") +
  labs(title    = "Random Forest: Variable Importance",
       subtitle = "Burn Probability in Guinean Lowland Forests",
       x = "Predictor Variable",
       y = "Gini Importance") +
  theme_minimal(base_size = 13)

print(p1)

# --- Plot 2: ROC Curve ---
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p2 <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#d7301f", linewidth = 1.2) +
  geom_abline(linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.7, y = 0.2,
           label = paste0("AUC = ", round(auc(roc_obj), 3)),
           size = 5, color = "#d7301f") +
  labs(title    = "ROC Curve - Burn Probability Classification",
       subtitle = "2016 Validation | EW Guinean Lowland Forests",
       x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_minimal(base_size = 13)

print(p2)

# --- Plot 3: Spatial Burn Probability Map ---
prob_vals <- test_df %>% select(x, y, pred_prob)
prob_rast <- rast(prob_vals, type = "xyz", crs = "ESRI:102022")
prob_rast <- resample(prob_rast, burn_mask_aea, method = "bilinear")

prob_df        <- as.data.frame(prob_rast, xy = TRUE, na.rm = TRUE)
names(prob_df)[3] <- "burn_prob"

p3 <- ggplot(prob_df, aes(x = x, y = y, fill = burn_prob)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("#1a9641", "#ffffbf", "#d7191c"),
    name    = "Burn\nProbability"
  ) +
  labs(title    = "Predicted Burn Probability (2016)",
       subtitle = "EW Guinean Lowland Forests | Random Forest",
       x = NULL, y = NULL) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank())

print(p3)

# --- Plot 4: Confusion Matrix Heatmap ---
cm_df <- as.data.frame(cm$table)
names(cm_df) <- c("Predicted", "Actual", "Count")

p4 <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "#fff7ec", high = "#d7301f") +
  labs(title    = "Confusion Matrix",
       subtitle = "Random Forest | 2016 Burn Classification") +
  theme_minimal(base_size = 13)

print(p4)

# ─────────────────────────────────────────────────────────────
# SECTION 12: SUMMARY TABLE
# ─────────────────────────────────────────────────────────────

summary_table <- data.frame(
  Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity",
             "Precision (PPV)", "F1 Score", "AUC"),
  Value  = c(
    round(cm$overall["Accuracy"],    3),
    round(cm$byClass["Sensitivity"], 3),
    round(cm$byClass["Specificity"], 3),
    round(cm$byClass["Precision"],   3),
    round(cm$byClass["F1"],          3),
    round(auc(roc_obj),              3)
  )
)

cat("\n===== SUMMARY TABLE =====\n")
print(summary_table)

# Also save raster and summary to file
writeRaster(prob_rast, "rf_burn_probability_2016.tif", overwrite = TRUE)
write.csv(summary_table, "rf_model_summary.csv", row.names = FALSE)

cat("\n✓ Done!\n")
cat("  Plots displayed in RStudio Plots panel\n")
cat("  Use arrow buttons in Plots panel to navigate between plots\n")
cat("  Saved: rf_burn_probability_2016.tif\n")
cat("  Saved: rf_model_summary.csv\n")