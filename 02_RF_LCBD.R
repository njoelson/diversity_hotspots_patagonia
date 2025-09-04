# Random forest LCBD model
# Dataset: Veglots and header klimnem 2021-2024 + raster layers
# Author: NJ
# Date: 10.4.2025

# Clean
rm(list = ls())

# Load required libraries
library(randomForest)
library(caret)
library(dplyr)
library(ggplot2)
library(vegan)
library(terra)
library(tmap)
library(RColorBrewer)
library(dismo)
library(spThin)
options(warn = -1)
library(scico)
library(adespatial)

# Load data
env <- read.csv("01_Data/02_Clean/predictors_planar_selected_w_div.csv",
                dec = ".", sep = ",", header = TRUE, row.names = 1)
species <- read.csv("01_Data/02_Clean/species.csv",
                    dec = ".", sep = ",", header = TRUE, row.names = 1)
predictor_masked <- rast("01_Data/02_Clean/predictor_masked_EPSG32719_250m.tif")

# Step 1: Prepare predictors data
env <- na.omit(env)

# Raster class_code as factor + RAT from env
cc <- predictor_masked[["class_code"]] |> as.factor()
ids <- sort(unique(terra::values(cc, mat = FALSE)))
map_df <- env |> dplyr::select(class_code, class) |> dplyr::distinct()
map_df$class_code <- as.integer(as.character(map_df$class_code))
lvl_tab <- data.frame(ID = ids) |> dplyr::left_join(map_df, by = c("ID" = "class_code"))
lvl_tab$class[is.na(lvl_tab$class)] <- as.character(lvl_tab$ID)
lvl_tab <- dplyr::rename(lvl_tab, class_code = class)
levels(cc) <- lvl_tab
predictor_masked[["class_code"]] <- cc

# Tabular class_code as factor with same labels
env$class_code <- factor(env$class, levels = lvl_tab$class_code)
stopifnot(is.factor(env$class_code), is.factor(predictor_masked[["class_code"]]))

# Step 2: Prepare diversity data
species <- species[match(rownames(env), rownames(species)), ]
isTRUE(all(rownames(env) == rownames(species)))
species <- species %>% mutate(across(everything(), as.numeric))

## PA + LCBD
species_pa <- decostand(species, method = "pa")
env$LCBD <- 100 * adespatial::beta.div(vegan::decostand(species, "pa"),
                                       method = "hellinger")$LCBD

# Step 3: Thin the datapoints
env$plot_id <- rownames(env); rownames(env) <- NULL
species_pa$plot_id <- rownames(species_pa); rownames(species_pa) <- NULL
thin_input <- env[, c("longitude", "latitude", "plot_id")]

set.seed(42)
thin_results <- thin(loc.data = thin_input, lat.col = "latitude", long.col = "longitude", # in coords
                     spec.col = "plot_id", thin.par = 0.1, reps = 100, # 0.1 degrees = ~ 11km
                     locs.thinned.list.return = TRUE, write.files = FALSE)
best_run <- which.max(sapply(thin_results, nrow))
thin_indices <- thin_results[[best_run]]

env <- env[rownames(thin_indices), , drop = FALSE]
species_pa <- species_pa[rownames(thin_indices), , drop = FALSE]
LCBD <- env$LCBD

# -------------------- SPLIT (STRATIFIED BY CLASS) --------------------

set.seed(42)
folds_LCBD <- createFolds(env$class_code, k = 6)    #stratified by class
train_indices <- unlist(folds_LCBD[2:6])
test_indices  <- folds_LCBD[[1]]

# Step 3: RFE on TRAIN only (unchanged logic)
set.seed(42)
ctrl <- rfeControl(functions = rfFuncs, method = "repeatedcv",
                   repeats = 10, verbose = FALSE, saveDetails = TRUE,
                   returnResamp = "all")
feat_idx <- c(5, 7, 8, 10:16, 18, 19)
subsets  <- 1:12

RFEresult1 <- rfe(x = env[train_indices, feat_idx, drop = FALSE],
                  y = LCBD[train_indices],
                  sizes = subsets, rfeControl = ctrl)

RFE_RMSE <- RFEresult1$results$RMSE
RFE_R2   <- RFEresult1$results$Rsquared
bestset_RMSE <- which((RFE_RMSE/min(RFE_RMSE)) <= 1.05)[1]
bestset_R2   <- which((RFE_R2/max(RFE_R2)*100) >= 95)[1]
chosenvars1  <- pickVars(RFEresult1$variables, min(bestset_RMSE, bestset_R2))
print(chosenvars1)

# Step 4: Model matrices
joint_data <- cbind(env, LCBD = LCBD)
joint_vars <- c("LCBD", chosenvars1)

LCBD_train <- joint_data[train_indices, joint_vars]
LCBD_test  <- joint_data[test_indices,  joint_vars]
table(LCBD_test$class_code)

# ---------------RANDOM FOREST FIT --------------------
### ensure factor in train/test (unchanged)
if ("class_code" %in% names(LCBD_train)) {
  LCBD_train$class_code <- droplevels(as.factor(LCBD_train$class_code))
  LCBD_test$class_code  <- droplevels(as.factor(LCBD_test$class_code))
}

## imbalance-aware weights for regression (inverse class frequency)
rf_weights <- NULL
if ("class_code" %in% names(LCBD_train)) {
  tab <- table(LCBD_train$class_code)
  inv <- 1 / as.numeric(tab[LCBD_train$class_code])  # per-row inverse freq
  rf_weights <- as.numeric(inv / mean(inv))          # normalize
}

set.seed(42)
rf_model <- randomForest(
  LCBD ~ mean_temp + canopy_height + ndvi_max + temp_range + annual_precip + class_code,
  data       = LCBD_train,
  importance = TRUE,
  ntree      = 500,
  replace    = TRUE,
  weights    = rf_weights
)
print(rf_model)

# -------------------- METRICS --------------------
# 1) OOB (from RF object)
OOB_RMSE <- sqrt(tail(rf_model$mse, 1))
OOB_R2   <- tail(rf_model$rsq, 1)

# 2) Single test
predicted <- predict(rf_model, newdata = LCBD_test)
observed  <- LCBD_test$LCBD
SSE <- sum((observed - predicted)^2)
SST <- sum((observed - mean(observed))^2)
TEST_R2   <- 1 - SSE/SST
TEST_RMSE <- sqrt(mean((observed - predicted)^2))
TEST_MAE  <- mean(abs(observed - predicted))
cal <- lm(observed ~ predicted)
TEST_CAL_SLOPE     <- unname(coef(cal)[2])
TEST_CAL_INTERCEPT <- unname(coef(cal)[1])

# 3) 6-fold CV (stratified by class_code)
fm <- as.formula(paste("LCBD ~", paste(chosenvars1, collapse = " + ")))
set.seed(42)
folds_cv <- createFolds(env$class_code, k = 6)

r2_cv   <- numeric(length(folds_cv))
rmse_cv <- numeric(length(folds_cv))
mae_cv  <- numeric(length(folds_cv))

for (i in seq_along(folds_cv)) {
  te_idx <- folds_cv[[i]]
  tr_idx <- setdiff(seq_len(nrow(joint_data)), te_idx)
  
  tr <- joint_data[tr_idx, c("LCBD", chosenvars1), drop = FALSE]
  te <- joint_data[te_idx,  c("LCBD", chosenvars1), drop = FALSE]
  
  if ("class_code" %in% names(tr)) {
    tr$class_code <- droplevels(as.factor(tr$class_code))
    te$class_code <- droplevels(as.factor(te$class_code))
  }
  
  # same imbalance-aware weights inside each fold
  w <- NULL
  if ("class_code" %in% names(tr)) {
    tabf <- table(tr$class_code)
    invf <- 1 / as.numeric(tabf[tr$class_code])
    w    <- as.numeric(invf / mean(invf))
  }
  
  set.seed(42)
  rf_cv <- randomForest(fm, data = tr, ntree = 500, importance = TRUE,
                        replace = TRUE, weights = w)
  
  yhat <- predict(rf_cv, newdata = te)
  y    <- te$LCBD
  sse  <- sum((y - yhat)^2)
  sst  <- sum((y - mean(y))^2)
  
  r2_cv[i]   <- 1 - sse/sst
  rmse_cv[i] <- sqrt(mean((y - yhat)^2))
  mae_cv[i]  <- mean(abs(y - yhat))
}

# CV summary (mean ± SD) 
CV_R2      <- mean(r2_cv);   CV_R2_SD   <- sd(r2_cv)
CV_RMSE    <- mean(rmse_cv); CV_RMSE_SD <- sd(rmse_cv)
CV_MAE     <- mean(mae_cv);  CV_MAE_SD  <- sd(mae_cv)

# names for plotting
LCBD_cv_r2  <- CV_R2
LCBD_cv_rmse<- CV_RMSE
LCBD_cv_mae <- CV_MAE

cat("\n--- Performance ---\n")
cat(sprintf("OOB (RF):        R² = %.3f   RMSE = %.3f\n", OOB_R2, OOB_RMSE))
cat(sprintf("Test (holdout):  R² = %.3f   RMSE = %.3f   MAE = %.3f\n",
            TEST_R2, TEST_RMSE, TEST_MAE))
cat(sprintf("CV (6-fold, strat. by class): R² = %.3f ± %.3f   RMSE = %.3f ± %.3f   MAE = %.3f ± %.3f\n",
            CV_R2, CV_R2_SD, CV_RMSE, CV_RMSE_SD, CV_MAE, CV_MAE_SD))
cat(sprintf("Calibration (holdout): slope = %.3f, intercept = %.3f\n",
            TEST_CAL_SLOPE, TEST_CAL_INTERCEPT))


# -------------------- IMPORTANCE & PDPs (unchanged) --------------------
varImpPlot(rf_model)

importance(rf_model)

# Save the variable importance and partial plots
png("02_Analysis/03_Figures/varImpPlot_LCBD.png", width = 1600, height = 1200, res = 300)
varImpPlot(rf_model, main = "Variable Importance - LCBD")
dev.off()

# Save multi-panel partial dependence plots for species LCBD

png("02_Analysis/03_Figures/PartialPlots_LCBD_combined1.png", width = 2400, height = 1800, res = 300)

par(mfrow = c(3,2))
partialPlot(x = rf_model, pred.data = LCBD_train, x.var = canopy_height,
            main = "Canopy Height", xlab = "Canopy Height", ylab = "LCBD", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = LCBD_train, x.var = mean_temp, main = "Mean Annual Temperature (°C)", xlab = "Mean Annual Temperature (°C)", ylab = "LCBD", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = LCBD_train, x.var = ndvi_max,
            main = "NDVI Max", xlab = "NDVI Max", ylab = "LCBD", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = LCBD_train, x.var = temp_range, main = "Temperature Range(°C)", xlab = "Temperature Range(°C)", ylab = "LCBD", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = LCBD_train, x.var = annual_precip, main = "Annual Precipitation (mm)", xlab = "Annual Precipitation (mm)", ylab = "LCBD", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = LCBD_train, x.var = class_code,
            main = "Land cover class", xlab = "Land cover class", ylab = "LCBD", las = 1, lwd = 2)
dev.off()

# -------------------- PREDICT RASTER --------------------
rf_prediction <- predict(predictor_masked, rf_model)

# -------------------- MAP -------------------
pal_batlow <- scico::scico(
  n = 255,
  palette = "batlow",
  begin = 0.05,
  end = 0.98,
  alpha = 1.0
)
plot_LCBD_map <- function(r, pal = pal_batlow, legend_title = "Predicted LCBD (%)",
                          oob_r2 = NULL, cv_r2 = NULL, oob_rmse = NULL, cv_rmse = NULL) {
  if (terra::nlyr(r) > 1) r <- r[[1]]
  minmax <- as.numeric(terra::global(r, range, na.rm = TRUE))
  labs   <- format(round(minmax, 2), nsmall = 2)
  v4 <- utils::packageVersion("tmap") >= "4.0.0"
  if (v4) {
    tm <- tm_shape(r) +
      tm_raster(
        col.scale  = tm_scale_continuous(values = pal, limits = minmax, ticks = minmax, labels = labs),
        col.legend = tm_legend(title = legend_title, reverse = FALSE, frame = FALSE, bg.color = NA, bg.alpha = 0)
      ) +
      tm_layout(legend.outside = FALSE, legend.position = c("left","top"),
                legend.text.size = 0.8, legend.title.size = 0.9,
                frame = FALSE, outer.margins = 0, inner.margins = c(0,0,0,0))
  } else {
    tm <- tm_shape(r) + tm_raster(style = "cont", palette = pal, title = legend_title) +
      tm_layout(legend.outside = FALSE, legend.position = c("left","top"),
                legend.frame = FALSE, legend.text.size = 0.8, legend.title.size = 0.9,
                frame = FALSE, outer.margins = 0, inner.margins = c(0,0,0,0))
  }
  line1 <- sprintf("OOB R\u00B2: %.2f   CV R\u00B2: %.2f",
                   oob_r2 %||% NA_real_, cv_r2 %||% NA_real_)
  line2 <- sprintf("OOB RMSE: %.2f   CV RMSE: %.2f",
                   oob_rmse %||% NA_real_, cv_rmse %||% NA_real_)
  tm + tm_credits(paste(line1, line2, sep = "\n"),
                  position = c("RIGHT", "TOP"), size = 1.1)
}
`%||%` <- function(a,b) if (is.null(a)) b else a

# background raster
bg_path <- "01_Data/02_Clean/shp_bg_matched_to_raster.gpkg"
tm_LCBD <-
  if (file.exists(bg_path)) {
    bg <- terra::vect(bg_path)
    tm_shape(bg) + tm_polygons(col = "grey90", border.col = NA) +
      plot_LCBD_map(rf_prediction, pal = pal_batlow,
                    oob_r2 = OOB_R2, cv_r2 = CV_R2,
                    oob_rmse = OOB_RMSE, cv_rmse = CV_RMSE)
  } else {
    plot_LCBD_map(rf_prediction, pal = pal_batlow,
                  oob_r2 = OOB_R2, cv_r2 = CV_R2,
                  oob_rmse = OOB_RMSE, cv_rmse = CV_RMSE)
  }

bb  <- terra::ext(rf_prediction)
asp <- (terra::ymax(bb) - terra::ymin(bb)) / (terra::xmax(bb) - terra::xmin(bb))
tmap_save(tm_LCBD, "02_Analysis/03_Figures/RF_LCBD_model_250m.svg",
          width = 8, height = 8 * asp, units = "in", dpi = 300)
tmap_save(tm_LCBD, "02_Analysis/03_Figures/RF_LCBD_model_250m.png",
          width = 8, height = 8 * asp, units = "in", dpi = 300)

# Save model & raster
saveRDS(rf_model, "02_Analysis/02_Results/rf_LCBD_model.rds")
terra::writeRaster(rf_prediction, "02_Analysis/02_Results/RF_LCBD_250m_full.tif", overwrite = TRUE)
saveRDS(rf_prediction, "02_Analysis/02_Results/LCBD_prediction_250m.rds")

