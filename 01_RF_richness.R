# Random forest species richness model
# Dataset: Veglots and header klimnem 2021-2024 + raster layers
# Author: NJ
# Date: 10.4.2025

# Clean
rm(list = ls())

# Libraries
library(randomForest); library(caret); library(dplyr); library(ggplot2)
library(vegan); library(terra); library(tmap); library(RColorBrewer)
library(dismo); library(spThin); options(warn = -1); library(scico)

# --------------------------- LOAD DATA ------------------------------------------
env <- read.csv("01_Data/02_Clean/predictors_planar_selected_w_div.csv",
                dec=".", sep=",", header=TRUE, row.names=1)
species <- read.csv("01_Data/02_Clean/species.csv",
                    dec=".", sep=",", header=TRUE, row.names=1)
predictor_masked <- rast("01_Data/02_Clean/predictor_masked_EPSG32719_250m.tif")

# --------------------------- PREPARE PREDICTORS ---------------------------------
env <- na.omit(env)

# class_code factor in raster + RAT from table
cc  <- predictor_masked[["class_code"]] |> as.factor()
ids <- sort(unique(terra::values(cc, mat = FALSE)))
map_df <- env |> dplyr::select(class_code, class) |> dplyr::distinct()
map_df$class_code <- as.integer(as.character(map_df$class_code))
lvl_tab <- data.frame(ID = ids) |> dplyr::left_join(map_df, by = c("ID"="class_code"))
lvl_tab$class[is.na(lvl_tab$class)] <- as.character(lvl_tab$ID)
lvl_tab <- dplyr::rename(lvl_tab, class_code = class)
levels(cc) <- lvl_tab
predictor_masked[["class_code"]] <- cc

# tabular factor with same labels
env$class_code <- factor(env$class, levels = lvl_tab$class_code)
stopifnot(is.factor(env$class_code), is.factor(predictor_masked[["class_code"]]))

# --------------------------- PREPARE DIVERSITY ----------------------------------
species <- species[match(rownames(env), rownames(species)), ]
stopifnot(isTRUE(all(rownames(env) == rownames(species))))
species <- species %>% mutate(across(everything(), as.numeric))
species_pa <- decostand(species, method = "pa")
species_pa$richness <- rowSums(species_pa)

# --------------------------- THINNING -------------------------------------------
env$plot_id <- rownames(env);           rownames(env) <- NULL
species_pa$plot_id <- rownames(species_pa); rownames(species_pa) <- NULL
thin_input <- env[, c("longitude","latitude","plot_id")]

set.seed(42)
thin_results <- thin(loc.data = thin_input, lat.col="latitude", long.col="longitude", # in wgs84
                     spec.col="plot_id", thin.par=0.1, reps=100, # 0.1 degrees = ~ 11km
                     locs.thinned.list.return=TRUE, write.files=FALSE)
best_run     <- which.max(sapply(thin_results, nrow))
thin_indices <- thin_results[[best_run]]

env        <- env[rownames(thin_indices), , drop = FALSE]
species_pa <- species_pa[rownames(thin_indices), , drop = FALSE]
richness   <- species_pa$richness

# --------------------------- SPLIT (STRATIFIED BY CLASS) ------------------------
set.seed(42)
folds_rich   <- createFolds(env$class_code, k = 6)  # stratified by land class levels
train_indices <- unlist(folds_rich[2:6])
test_indices  <- folds_rich[[1]]

# --------------------------- FEATURE SELECTION (TRAIN DATA ONLY) ---------------------
set.seed(42)
ctrl    <- rfeControl(functions = rfFuncs, method="repeatedcv",
                      repeats=10, verbose=FALSE, saveDetails=TRUE,
                      returnResamp="all")
feat_idx <- c(5,7,8,10:16,18,19)
subsets  <- 1:12

RFE_rich <- rfe(x = env[train_indices, feat_idx, drop = FALSE],
                y = richness[train_indices],
                sizes = subsets, rfeControl = ctrl)

RFE_RMSE <- RFE_rich$results$RMSE
RFE_R2   <- RFE_rich$results$Rsquared
best_RMSE <- which((RFE_RMSE/min(RFE_RMSE)) <= 1.05)[1] # set of vars within 5% of min RMSE
best_R2   <- which((RFE_R2/max(RFE_R2)*100) >= 95)[1] # set of vars within 5% of max R2
chosenvars1 <- pickVars(RFE_rich$variables, min(best_RMSE, best_R2))
print(chosenvars1)

# --------------------------- MODEL DATA -----------------------------------------
joint_data <- cbind(env, richness = richness)
joint_vars <- c("richness", chosenvars1)
rich_train <- joint_data[train_indices, joint_vars]
rich_test  <- joint_data[test_indices,  joint_vars]

# ensure class factor in both splits
if ("class_code" %in% names(rich_train)) {
  rich_train$class_code <- droplevels(as.factor(rich_train$class_code))
  rich_test$class_code  <- droplevels(as.factor(rich_test$class_code))
}

# inverse-frequency regression weights by class
rf_weights <- NULL
if ("class_code" %in% names(rich_train)) {
  tab <- table(rich_train$class_code)
  inv <- 1 / as.numeric(tab[rich_train$class_code])
  rf_weights <- as.numeric(inv / mean(inv))
}

# --------------------------- RANDOM FOREST --------------------------------------
fm <- as.formula(paste("richness ~", paste(chosenvars1, collapse = " + ")))

set.seed(42)
rf_model <- randomForest(fm,
                         data       = rich_train,
                         importance = TRUE,
                         ntree      = 500,
                         replace    = TRUE,
                         weights    = rf_weights)
print(rf_model)

# --------------------------- METRICS --------------------------------------------
# 1) OOB
OOB_RMSE <- sqrt(tail(rf_model$mse, 1))
OOB_R2   <- tail(rf_model$rsq, 1)

# 2) Single test
predicted <- predict(rf_model, newdata = rich_test)
observed  <- rich_test$richness
SSE <- sum((observed - predicted)^2)
SST <- sum((observed - mean(observed))^2)
TEST_R2   <- 1 - SSE/SST
TEST_RMSE <- sqrt(mean((observed - predicted)^2))
TEST_MAE  <- mean(abs(observed - predicted))
cal <- lm(observed ~ predicted)
TEST_CAL_SLOPE     <- unname(coef(cal)[2])
TEST_CAL_INTERCEPT <- unname(coef(cal)[1])

# 3) 6-fold CV (stratified by class_code)
set.seed(42)
folds_cv <- createFolds(env$class_code, k = 6)

r2_cv <- rmse_cv <- mae_cv <- numeric(length(folds_cv))
for (i in seq_along(folds_cv)) {
  te <- folds_cv[[i]]
  tr <- setdiff(seq_len(nrow(joint_data)), te)
  dtr <- joint_data[tr, joint_vars, drop = FALSE]
  dte <- joint_data[te, joint_vars, drop = FALSE]
  if ("class_code" %in% names(dtr)) {
    dtr$class_code <- droplevels(as.factor(dtr$class_code))
    dte$class_code <- droplevels(as.factor(dte$class_code))
  }
  w <- NULL
  if ("class_code" %in% names(dtr)) {
    tabf <- table(dtr$class_code)
    invf <- 1 / as.numeric(tabf[dtr$class_code])
    w    <- as.numeric(invf / mean(invf))
  }
  set.seed(42)
  rf_cv <- randomForest(fm, data = dtr, ntree = 500, importance = TRUE,
                        replace = TRUE, weights = w)
  yhat <- predict(rf_cv, newdata = dte)
  y    <- dte$richness
  sse  <- sum((y - yhat)^2)
  sst  <- sum((y - mean(y))^2)
  r2_cv[i]   <- 1 - sse/sst
  rmse_cv[i] <- sqrt(mean((y - yhat)^2))
  mae_cv[i]  <- mean(abs(y - yhat))
}
CV_R2 <- mean(r2_cv);   CV_R2_SD <- sd(r2_cv)
CV_RMSE <- mean(rmse_cv); CV_RMSE_SD <- sd(rmse_cv)
CV_MAE <- mean(mae_cv);   CV_MAE_SD <- sd(mae_cv)

cat("\n--- Performance (Richness) ---\n")
cat(sprintf("OOB (RF):        R² = %.3f   RMSE = %.3f\n", OOB_R2, OOB_RMSE))
cat(sprintf("Test (holdout):  R² = %.3f   RMSE = %.3f   MAE = %.3f\n",
            TEST_R2, TEST_RMSE, TEST_MAE))
cat(sprintf("CV (6-fold, strat.): R² = %.3f ± %.3f   RMSE = %.3f ± %.3f   MAE = %.3f ± %.3f\n",
            CV_R2, CV_R2_SD, CV_RMSE, CV_RMSE_SD, CV_MAE, CV_MAE_SD))
cat(sprintf("Calibration (holdout): slope = %.3f, intercept = %.3f\n",
            TEST_CAL_SLOPE, TEST_CAL_INTERCEPT))

# --------------------------- IMPORTANCE & PDPs ----------------------------------
png("02_Analysis/03_Figures/varImpPlot_richness.png", width=1600, height=1200, res=300)
varImpPlot(rf_model, main="Variable Importance - Richness")
dev.off()

# PDPs only for selected predictors (skip factor class_code)
png("02_Analysis/03_Figures/PartialPlots_richness_combined.png", width = 2400, height = 1800, res = 300)

par(mfrow = c(3,2))

partialPlot(x = rf_model, pred.data = rich_train, x.var = mean_temp, main = "Mean Annual Temperature (°C)", xlab = "Mean Annual Temperature (°C)", ylab = "richness", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rich_train, x.var = annual_precip, main = "Annual Precipitation (mm)", xlab = "Annual Precipitation (mm)", ylab = "richness", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rich_train, x.var = ndvi_max,
            main = "NDVI Max", xlab = "NDVI Max", ylab = "richness", las = 1, lwd = 2)
partialPlot(x = rf_model, pred.data = rich_train, x.var = class_code,
            main = "Land Cover Class", xlab = "Land Cover Class", ylab = "richness", las = 1, lwd = 2)
dev.off()

# --------------------------- SPATIAL PREDICTION ----------------------------------
rf_prediction <- predict(predictor_masked, rf_model)

# --------------------------- MAP  -----------------------
pal_batlow <- scico::scico(n = 255, palette = "batlow", begin = 0.05, end = 0.98, alpha = 1.0)
`%||%` <- function(a,b) if (is.null(a)) b else a

plot_rich_map <- function(r, pal = pal_batlow, legend_title = "Predicted richness",
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
    tm <- tm_shape(r) + tm_raster(style="cont", palette=pal, title=legend_title) +
      tm_layout(legend.outside=FALSE, legend.position=c("left","top"),
                legend.frame=FALSE, legend.text.size=0.8, legend.title.size=0.9,
                frame=FALSE, outer.margins=0, inner.margins=c(0,0,0,0))
  }
  line1 <- sprintf("OOB R\u00B2: %.2f   CV R\u00B2: %.2f", oob_r2, cv_r2)
  line2 <- sprintf("OOB RMSE: %.2f   CV RMSE: %.2f", oob_rmse, cv_rmse)
  tm + tm_credits(paste(line1, line2, sep = "\n"), position = c("RIGHT","TOP"), size = 1.1)
}

# raster for grey background
bg_path <- "01_Data/02_Clean/shp_bg_matched_to_raster.gpkg"

# plot
tm_rich <-
  if (file.exists(bg_path)) {
    bg <- terra::vect(bg_path)
    tm_shape(bg) + tm_polygons(col="grey90", border.col = NA) +
      plot_rich_map(rf_prediction, pal = pal_batlow,
                    oob_r2 = OOB_R2, cv_r2 = CV_R2,
                    oob_rmse = OOB_RMSE, cv_rmse = CV_RMSE)
  } else {
    plot_rich_map(rf_prediction, pal = pal_batlow,
                  oob_r2 = OOB_R2, cv_r2 = CV_R2,
                  oob_rmse = OOB_RMSE, cv_rmse = CV_RMSE)
  }

bb  <- terra::ext(rf_prediction)
asp <- (terra::ymax(bb) - terra::ymin(bb)) / (terra::xmax(bb) - terra::xmin(bb))
tmap_save(tm_rich, "02_Analysis/03_Figures/RF_richness_model_250m.svg",
          width = 8, height = 8 * asp, units = "in", dpi = 300)
tmap_save(tm_rich, "02_Analysis/03_Figures/RF_richness_model_250m.png",
          width = 8, height = 8 * asp, units = "in", dpi = 300)

# --------------------------- SAVE OBJECTS ---------------------------------------
saveRDS(rf_model, "02_Analysis/02_Results/rf_richness_model.rds")
terra::writeRaster(rf_prediction, "02_Analysis/02_Results/RF_richness_250m_full.tif", overwrite = TRUE)
saveRDS(rf_prediction, "02_Analysis/02_Results/richness_prediction_250m.rds")
