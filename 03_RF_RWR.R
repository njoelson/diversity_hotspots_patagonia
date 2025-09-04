# Random forest RWR model
# Dataset: Veglots and header klimnem 2021-2024 + raster layers
# Author: NJ
# Date: 26.08.25
rm(list = ls())

# Packages
library(randomForest); library(caret); library(dplyr); library(ggplot2)
library(vegan); library(terra); library(tmap); library(spThin); library(scico)

# Data
env <- read.csv("01_Data/02_Clean/predictors_planar_selected_w_div.csv",
                dec=".", sep=",", header=TRUE, row.names=1)
species <- read.csv("01_Data/02_Clean/species.csv",
                    dec=".", sep=",", header=TRUE, row.names=1)
predictor_masked <- rast("01_Data/02_Clean/predictor_masked_EPSG32719_250m.tif")

# --------------------------- Predictors -----------------------
env <- na.omit(env)

# Raster class_code as factor + RAT from env
cc  <- predictor_masked[["class_code"]] |> as.factor()
ids <- sort(unique(terra::values(cc, mat = FALSE)))
map_df <- env |> dplyr::select(class_code, class) |> dplyr::distinct()
map_df$class_code <- as.integer(as.character(map_df$class_code))
lvl_tab <- data.frame(ID = ids) |> dplyr::left_join(map_df, by = c("ID" = "class_code"))
lvl_tab$class[is.na(lvl_tab$class)] <- as.character(lvl_tab$ID)
lvl_tab <- dplyr::rename(lvl_tab, class_code = class)
levels(cc) <- lvl_tab
predictor_masked[["class_code"]] <- cc

# Tabular class_code as factor with the same labels
env$class_code <- factor(env$class, levels = lvl_tab$class_code)
stopifnot(is.factor(env$class_code), is.factor(predictor_masked[["class_code"]]))

# --------------------------- Species data ----------------------------------
# Align species rows to env, ensure numeric
species <- species[match(rownames(env), rownames(species)), ]
stopifnot(isTRUE(all(rownames(env) == rownames(species))))
species <- species %>% mutate(across(everything(), as.numeric))

# Remove species with total 0 (avoid divide-by-zero)
species <- species[, colSums(species, na.rm = TRUE) > 0, drop = FALSE]

# Presence/absence
species_pa <- vegan::decostand(species, method = "pa")

# RWR: inverse-frequency weighted richness (Williams 1996)
sp_freq    <- colSums(species_pa)          # frequency per species
sp_weight  <- 1 / sp_freq                  # weight = 1 / frequency
sp_weighted <- sweep(species_pa, 2, sp_weight, `*`)
RWR_full   <- rowSums(sp_weighted)
env$RWR    <- RWR_full

# --------------------------- Thinning  -----------------------------
env$plot_id <- rownames(env); rownames(env) <- NULL
species_pa$plot_id <- rownames(species_pa); rownames(species_pa) <- NULL
thin_input <- env[, c("longitude", "latitude", "plot_id")]

set.seed(42)
thin_results <- thin(loc.data = thin_input, lat.col = "latitude", long.col = "longitude",
                     spec.col = "plot_id", thin.par = 0.1, reps = 100,
                     locs.thinned.list.return = TRUE, write.files = FALSE)
best_run     <- which.max(sapply(thin_results, nrow))
thin_indices <- thin_results[[best_run]]

env        <- env[rownames(thin_indices), , drop = FALSE]
species_pa <- species_pa[rownames(thin_indices), , drop = FALSE]
RWR        <- env$RWR  # already thinned with env

# -------------------- Stratified split by class_code -----------------
set.seed(42)
folds_RWR     <- createFolds(env$class_code, k = 6)    # stratified by class
train_indices <- unlist(folds_RWR[2:6])
test_indices  <- folds_RWR[[1]]

# --------------------------- Feature selection (train data only) -----------------
set.seed(42)
ctrl    <- rfeControl(functions = rfFuncs, method = "repeatedcv",
                      repeats = 10, verbose = FALSE, saveDetails = TRUE,
                      returnResamp = "all")
feat_idx <- c(5, 7, 8, 10:16, 18, 19)  # same set of predictors in all models
subsets  <- 1:12

RFE_RWR <- rfe(x = env[train_indices, feat_idx, drop = FALSE],
               y = RWR[train_indices],
               sizes = subsets, rfeControl = ctrl)

RFE_RMSE <- RFE_RWR$results$RMSE
RFE_R2   <- RFE_RWR$results$Rsquared
best_RMSE <- which((RFE_RMSE / min(RFE_RMSE)) <= 1.05)[1]
best_R2   <- which((RFE_R2   / max(RFE_R2)  *100) >= 95)[1]
chosenvars1 <- pickVars(RFE_RWR$variables, min(best_RMSE, best_R2))
print(chosenvars1)

# 7 vars: "ndvi_max"      "mean_temp"     "annual_precip" "temp_range"    "canopy_height" "class_code"    "dist_pw"  
# --------------------------- Model matrices --------------------------------------
joint_data <- cbind(env, RWR = RWR)
joint_vars <- c("RWR", chosenvars1)
rwr_train  <- joint_data[train_indices, joint_vars]
rwr_test   <- joint_data[test_indices,  joint_vars]

# Ensure factor + imbalance-aware weights (inverse class frequency)
if ("class_code" %in% names(rwr_train)) {
  rwr_train$class_code <- droplevels(as.factor(rwr_train$class_code))
  rwr_test$class_code  <- droplevels(as.factor(rwr_test$class_code))
}
rf_weights <- NULL
if ("class_code" %in% names(rwr_train)) {
  tab <- table(rwr_train$class_code)
  inv <- 1 / as.numeric(tab[rwr_train$class_code])
  rf_weights <- as.numeric(inv / mean(inv))
}

# --------------------------- Random Forest fit -----------------------------------
fm <- as.formula(paste("RWR ~", paste(chosenvars1, collapse = " + ")))
set.seed(42)
rf_model <- randomForest(fm, data = rwr_train,
                         importance = TRUE, ntree = 500,
                         replace = TRUE, weights = rf_weights)
print(rf_model)

# var exp = 32.5
#rmse 1.571805 of 0.06778906 9.19133726

# --------------------------- Metrics: OOB, Test, CV ------------------------------
# OOB
OOB_RMSE <- sqrt(tail(rf_model$mse, 1))
OOB_R2   <- tail(rf_model$rsq, 1)

# Single test
pred_test <- predict(rf_model, newdata = rwr_test)
obs_test  <- rwr_test$RWR
SSE <- sum((obs_test - pred_test)^2)
SST <- sum((obs_test - mean(obs_test))^2)
TEST_R2   <- 1 - SSE/SST
TEST_RMSE <- sqrt(mean((obs_test - pred_test)^2))
TEST_MAE  <- mean(abs(obs_test - pred_test))
cal <- lm(obs_test ~ pred_test)
TEST_CAL_SLOPE     <- unname(coef(cal)[2])
TEST_CAL_INTERCEPT <- unname(coef(cal)[1])

# 6-fold CV (stratified by class_code)
set.seed(42)
folds_cv <- createFolds(env$class_code, k = 6)
r2_cv <- rmse_cv <- mae_cv <- numeric(length(folds_cv))

for (i in seq_along(folds_cv)) {
  te_idx <- folds_cv[[i]]
  tr_idx <- setdiff(seq_len(nrow(joint_data)), te_idx)
  tr <- joint_data[tr_idx, c("RWR", chosenvars1), drop = FALSE]
  te <- joint_data[te_idx,  c("RWR", chosenvars1), drop = FALSE]
  if ("class_code" %in% names(tr)) {
    tr$class_code <- droplevels(as.factor(tr$class_code))
    te$class_code <- droplevels(as.factor(te$class_code))
  }
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
  y    <- te$RWR
  sse  <- sum((y - yhat)^2)
  sst  <- sum((y - mean(y))^2)
  r2_cv[i]   <- 1 - sse/sst
  rmse_cv[i] <- sqrt(mean((y - yhat)^2))
  mae_cv[i]  <- mean(abs(y - yhat))
}

CV_R2      <- mean(r2_cv);   CV_R2_SD   <- sd(r2_cv)
CV_RMSE    <- mean(rmse_cv); CV_RMSE_SD <- sd(rmse_cv)
CV_MAE     <- mean(mae_cv);  CV_MAE_SD  <- sd(mae_cv)

cat("\n--- Performance (RWR) ---\n")
cat(sprintf("OOB (RF):        R² = %.3f   RMSE = %.3f\n", OOB_R2, OOB_RMSE))
cat(sprintf("Test (holdout):  R² = %.3f   RMSE = %.3f   MAE = %.3f\n",
            TEST_R2, TEST_RMSE, TEST_MAE))
cat(sprintf("CV (6-fold, strat. by class): R² = %.3f ± %.3f   RMSE = %.3f ± %.3f   MAE = %.3f ± %.3f\n",
            CV_R2, CV_R2_SD, CV_RMSE, CV_RMSE_SD, CV_MAE, CV_MAE_SD))
cat(sprintf("Calibration (holdout): slope = %.3f, intercept = %.3f\n",
            TEST_CAL_SLOPE, TEST_CAL_INTERCEPT))

# --------------------------- Importance & PDPs -----------------------------------
png("02_Analysis/03_Figures/varImpPlot_RWR.png", width = 1600, height = 1200, res = 300)
varImpPlot(rf_model, main = "Variable Importance - RWR")
dev.off()

png("02_Analysis/03_Figures/PartialPlots_rwr_combined_1.png", width = 2400, height = 1800, res = 300)
par(mfrow = c(2,2))
partialPlot(x = rf_model, pred.data = rwr_train, x.var = canopy_height,
            main = "Canopy Height", xlab = "Canopy Height", ylab = "RWR", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rwr_train, x.var = mean_temp, main = "Mean Annual Temperature (°C)", xlab = "Mean Annual Temperature (°C)", ylab = "RWR", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rwr_train, x.var = ndvi_max,
            main = "NDVI Max", xlab = "NDVI Max", ylab = "RWR", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rwr_train, x.var = temp_range, main = "Temperature Range(°C)", xlab = "Temperature Range(°C)", ylab = "RWR", las = 1, lwd = 2)
dev.off()

png("02_Analysis/03_Figures/PartialPlots_rwr_combined_2.png", width = 2400, height = 1800, res = 300)
par(mfrow = c(2,2))
partialPlot(x = rf_model, pred.data = rwr_train, x.var = annual_precip, main = "Annual Precipitation (mm)", xlab = "Annual Precipitation (mm)", ylab = "RWR", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rwr_train, x.var = class_code,
            main = "Land cover class", xlab = "Land cover class", ylab = "RWR", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rwr_train, x.var = dist_pw,
            main = "Distance to permanent water bodies (m)", xlab = "Distance to permanent water bodies (m)", ylab = "RWR", las = 1, lwd = 2)

dev.off()

# --------------------------- Predict to raster -----------------------------------
rwr_pred <- terra::predict(predictor_masked, rf_model)
terra::writeRaster(rwr_pred, "02_Analysis/02_Results/RWR_prediction_250m.tif", overwrite = TRUE)
saveRDS(rwr_pred,  "02_Analysis/02_Results/RWR_prediction_250m.rds")
saveRDS(rf_model,  "02_Analysis/02_Results/rf_RWR_model.rds")

# --------------------------- Map ----------------------------
`%||%` <- function(a,b) if (is.null(a)) b else a

pal_batlow <- scico::scico(n = 255, palette = "batlow", begin = 0.05, end = 0.98, alpha = 1.0)

plot_RWR_map <- function(r, pal = pal_batlow, legend_title = "Predicted RWR",
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
  line1 <- sprintf("OOB R\u00B2: %.2f   CV R\u00B2: %.2f",    oob_r2 %||% NA_real_, cv_r2 %||% NA_real_)
  line2 <- sprintf("OOB RMSE: %.2f   CV RMSE: %.2f", oob_rmse %||% NA_real_, cv_rmse %||% NA_real_)
  tm + tm_credits(paste(line1, line2, sep = "\n"), position = c("RIGHT","TOP"), size = 1.1)
}

# raster for background
bg_path <- "01_Data/02_Clean/shp_bg_matched_to_raster.gpkg"
tm_RWR <-
  if (file.exists(bg_path)) {
    bg <- terra::vect(bg_path)
    if (!terra::same.crs(bg, rwr_pred)) bg <- terra::project(bg, terra::crs(rwr_pred))
    tm_shape(bg) + tm_polygons(col = "grey95", border.col = NA) +
      plot_RWR_map(rwr_pred, pal = pal_batlow,
                   oob_r2 = OOB_R2, cv_r2 = CV_R2,
                   oob_rmse = OOB_RMSE, cv_rmse = CV_RMSE)
  } else {
    plot_RWR_map(rwr_pred, pal = pal_batlow,
                 oob_r2 = OOB_R2, cv_r2 = CV_R2,
                 oob_rmse = OOB_RMSE, cv_rmse = CV_RMSE)
  }

bb  <- terra::ext(rwr_pred)
asp <- (terra::ymax(bb) - terra::ymin(bb)) / (terra::xmax(bb) - terra::xmin(bb))
tmap_save(tm_RWR, "02_Analysis/03_Figures/RF_RWR_model_250m.svg",
          width = 8, height = 8 * asp, units = "in", dpi = 300)
tmap_save(tm_RWR, "02_Analysis/03_Figures/RF_RWR_model_250m.png",
          width = 8, height = 8 * asp, units = "in", dpi = 300)

