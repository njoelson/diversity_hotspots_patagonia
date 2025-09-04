# SnR null model generation and testing against observed
# Dataset: Prediction from RF diversity models
# Author: NJ
# Date: 18.8.25

# ---- Clean ----
rm(list = ls())


# ---- Load required libraries ----
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

# --- Load functions and helpers ----

source("02_Analysis/01_Scripts/sources/Shift&RotFunctions.R")   # provides shiftrotXY() for nulls

# ovelap over quantiles plotting function
overlap_curve <- function(v1, v2, qs = seq(0.05, 0.50, by = 0.05)) {
  sapply(qs, function(q) {
    t1 <- stats::quantile(v1, probs = 1 - q, na.rm = TRUE, type = 7)
    t2 <- stats::quantile(v2, probs = 1 - q, na.rm = TRUE, type = 7)
    A  <- v1 >= t1; B <- v2 >= t2
    inter <- sum(A & B, na.rm = TRUE); uni <- sum(A | B, na.rm = TRUE)
    if (uni == 0) NA_real_ else inter / uni
  })
}

# null distribution rho histograms plotting funtion
save_hist <- function(x, x_obs, title, subtitle, xlab, file) {
  png(file, 1400, 1000, res = 200); print(
    ggplot(data.frame(x), aes(x)) + geom_histogram(bins = 40) +
      geom_vline(xintercept = mean(x, na.rm = TRUE), linetype = 2) +
      geom_vline(xintercept = x_obs, linewidth = 1) +
      labs(title = title, subtitle = subtitle, x = xlab, y = "Count") + theme_minimal()
  ); dev.off()
}
save_curve <- function(df, title, subtitle, xlab, ylab, file) {
  png(file, 1400, 1000, res = 200); print(
    ggplot(df, aes(q, obs)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20) +
      geom_line(aes(y = mean), linetype = 2) +
      geom_line(linewidth = 1) +
      labs(title = title, subtitle = subtitle, x = xlab, y = ylab) + theme_minimal()
  ); dev.off()
}

# ---- Load data for observed model ----
env <- read.csv("01_Data/02_Clean/predictors_planar_selected_w_div.csv",
                dec = ".", 
                sep = ",",
                header = T,
                row.names = 1)
species <- read.csv("01_Data/02_Clean/species.csv",
                    dec = ".", 
                    sep = ",",
                    header = T,
                    row.names = 1)

predictor_masked <- rast("01_Data/02_Clean/predictor_masked_EPSG32719_250m.tif")

# ---- Load observed model prediction rasters ----

lcbd_path     <- "02_Analysis/02_Results/02_Predictions/LCBD_prediction_250m.tif"
richness_path <- "02_Analysis/02_Results/02_Predictions/RF_richness_250m.tif"

# ---- Create output directories ----
out_dir_fig   <- "02_Analysis/03_Figures"
out_dir_res   <- "02_Analysis/02_Results"
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_res, showWarnings = FALSE, recursive = TRUE)

# ---- Step 1: Prepare predictors data ----

# Keep only predictors of interest
env <- na.omit(env)  # remove rows with NAs

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

# tabular factor with same labels, drop unused levels
env$class_code <- droplevels(factor(env$class, levels = lvl_tab$class_code))
predictor_masked[["class_code"]] <- droplevels(predictor_masked[["class_code"]])
stopifnot(is.factor(env$class_code), is.factor(predictor_masked[["class_code"]]))

# ---- Step 2: Prepare diversity data ----

# Match predictors and species matrix
species <- species[match(rownames(env), rownames(species)), ]
# check if they are in the same order
isTRUE(all(rownames(env) == rownames(species)))

# Convert all values to numeric
species <- species %>% 
  mutate(across(everything(), as.numeric))

str(species) #check

# Convert species data to presence/absence
species_pa <- decostand(species, method = "pa") 

# Calculate species richness per plot
species_pa$richness <- rowSums(species_pa)

# ---- Step 3: Thin the datapoints ----
# Set rownames as plot id column
env$plot_id <- rownames(env)
rownames(env) <- NULL
species_pa$plot_id <- rownames(species_pa)
rownames(species_pa) <- NULL

# Load coordinates
thin_input <- env[, c("longitude", "latitude", "plot_id")]

# Run spThin with many reps
set.seed(42)
thin_results <- thin(
  loc.data = thin_input,
  lat.col = "latitude",
  long.col = "longitude",
  spec.col = "plot_id",
  thin.par = 0.1,         # 100 meters minimum distance
  reps = 100,             # 100 repetitions
  locs.thinned.list.return = TRUE,
  write.files = FALSE
)

# Find the thinning that keeps the most points
thin_lengths <- sapply(thin_results, nrow)
best_run <- which.max(thin_lengths)
thin_indices <- thin_results[[best_run]]

# Filter datasets
env_thinned <- env[rownames(thin_indices), , drop = FALSE]
species_thinned <- species_pa[rownames(thin_indices), , drop = FALSE]

# Update your working data
env <- env_thinned
species_pa <- species_thinned

# Prepare data for ggplot
thin_input$dataset <- "Original"
thin_indices$dataset <- "Thinned"

# Drop the plot_id column
thin_input$plot_id <- NULL

# Make the names in thin_indices lower case
names(thin_indices) <- tolower(names(thin_indices))

# Pull richness vector from the thinned dataset
# Check if richness column exists, if not calculate it
if ("richness" %in% names(species_thinned)) {
  richness <- species_thinned$richness
} else {
  # Calculate richness as row sums (excluding plot_id if present)
  species_cols <- setdiff(names(species_thinned), "plot_id")
  richness <- rowSums(species_thinned[, species_cols, drop = FALSE])
}

env$richness <- richness

# ---- Step 3: Feature selection ----


# Select variables from observed model feature selection
chosenvars1 <- c("mean_temp","ndvi_max","annual_precip","class_code")
chosenvars1

# ---- Step 4: Spatial model ----

# Split train and test dataset stratified by class
set.seed(42)
folds_rich   <- createFolds(env$class_code, k = 6)  # same as LCBD/RWR
train_indices <- unlist(folds_rich[2:6])
test_indices  <- folds_rich[[1]]

joint_data <- cbind(env, richness = richness)
joint_vars <- c("richness", chosenvars1)
rich_train <- joint_data[train_indices, joint_vars]
rich_test  <- joint_data[test_indices,  joint_vars]

# ensure class factor in both splits
if ("class_code" %in% names(rich_train)) {
  rich_train$class_code <- droplevels(as.factor(rich_train$class_code))
  rich_test$class_code  <- droplevels(as.factor(rich_test$class_code))
}

# inverse-frequency regression weights by class (like LCBD/RWR)
rf_weights <- NULL
if ("class_code" %in% names(rich_train)) {
  tab <- table(rich_train$class_code)
  inv <- 1 / as.numeric(tab[rich_train$class_code])
  rf_weights <- as.numeric(inv / mean(inv))
}

# --------------------------- RANDOM FOREST (observed) --------------------------------------
fm <- as.formula(paste("richness ~", paste(chosenvars1, collapse = " + ")))

set.seed(42)
rf_model <- randomForest(fm,
                         data       = rich_train,
                         importance = TRUE,
                         ntree      = 500,
                         replace    = TRUE,
                         weights    = rf_weights)

# View model summary
print(rf_model) # 30% var explained 37.8 sq residuals

cat("RMSE:", sqrt(rf_model$mse[500]), "\n") # RMSE: 6.145438 (range 4-45)
cat("Pseudo R²:", rf_model$rsq[500], "\n") # Pseudo R²: 0.29

# ---- Step 5: Prediction and evaluation ----

# Predict
rf_prediction <- predict(predictor_masked, rf_model)

# Performance metrics
predicted <- predict(rf_model, newdata = rich_test)
observed <- rich_test$richness

rmse <- sqrt(mean((predicted - observed)^2))
r2 <- cor(predicted, observed)^2

cat("\nModel performance on test data:\n")
cat("RMSE:", rmse, "\n")
cat("R-squared:", r2, "\n")

# ---- Step 6: S&R nulls + LCBD vs Richness test ----

# Consistent predictors
rf_vars <- c("mean_temp","ndvi_max","annual_precip","class_code")
pred_base <- predictor_masked[[rf_vars]] # subset to selected variables
class_levels <- levels(pred_base[["class_code"]])[[1]]$class_code


pred_src <- pred_base  # Rename 
valid_mask <- terra::app(pred_src, fun = function(...) as.integer(all(!is.na(c(...)))))
names(valid_mask) <- "valid"
ex2     <- terra::ext(pred_src)
X.range <- c(ex2$xmin, ex2$xmax)
Y.range <- c(ex2$ymin, ex2$ymax)
coord   <- as.matrix(env[, c("x","y")])   # planar coords
cat("Total cells:", ncell(pred_src[[1]]), "\n")
cat("Valid cells:", sum(terra::values(valid_mask), na.rm = TRUE), "\n")

# Observed rasters 
r_lcbd <- terra::rast(lcbd_path)
r_rich <- terra::rast(richness_path)

# Check if rasters match predictor_masked geometry
if (!terra::compareGeom(r_lcbd, predictor_masked, stopOnError = FALSE)) {
  cat("Resampling LCBD raster to match predictor_masked geometry...\n")
  r_lcbd <- terra::resample(r_lcbd, predictor_masked, method = "bilinear")
}
if (!terra::compareGeom(r_rich, predictor_masked, stopOnError = FALSE)) {
  cat("Resampling richness raster to match predictor_masked geometry...\n")
  r_rich <- terra::resample(r_rich, predictor_masked, method = "bilinear")
}

# Common mask: areas where all three rasters have data
msk <- !is.na(r_lcbd) & !is.na(r_rich) & !is.na(predictor_masked[[1]])
r_lcbd <- terra::mask(r_lcbd, msk)
r_rich <- terra::mask(r_rich, msk)
cat("Valid pixels after masking:", sum(terra::values(msk), na.rm = TRUE), "\n")

# observed vectors & stats
v_lcbd      <- as.vector(r_lcbd[])
v_rich      <- as.vector(r_rich[])
qs          <- seq(0.05, 0.50, by = 0.05)
obs_rho     <- suppressWarnings(cor(v_lcbd, v_rich, method = "spearman", use = "complete.obs")) # spearman rho, droppin NAs
obs_overlap <- overlap_curve(v_lcbd, v_rich, qs)
obs_overlap
# grid to predict on (same geometry as observed richness raster)
template <- r_rich
grid_xy  <- terra::crds(template, df = TRUE)

# tiny rigid-transform helpers
estT <- function(X,Y){X<-as.matrix(X);Y<-as.matrix(Y);mx<-colMeans(X);my<-colMeans(Y) # computes rotation matrix and shift vector to go from original coords to rotated
Xc<-sweep(X,2,mx);Yc<-sweep(Y,2,my);sv<-svd(t(Xc)%*%Yc);R<-sv$v%*%t(sv$u)
if (det(R)<0){sv$v[,ncol(sv$v)]<- -sv$v[,ncol(sv$v)];R<-sv$v%*%t(sv$u)}
list(R=R,mx=mx,my=my)}
applyT <- function(P,T){sweep(as.matrix(P),2,T$mx,"-")%*%T$R + matrix(T$my,nrow(P),2,byrow=TRUE)} # shifts coords (e.g. prediction raster) with tsr function (equally to rotated points)

# Build S&R nulls and test vs LCBD (observed model)
nsim_snr <- 999
set.seed(42)
null_rho <- numeric(nsim_snr) # vector of ρ from S&R null replicates vs LCBD
null_ovl <- matrix(NA_real_, nrow = nsim_snr, ncol = length(qs)) # a matrix of overlap curves from the nulls, dimension = nsim_snr × length(qs)

for (i in seq_len(nsim_snr)) {
  # 1) find a valid shift+rotation - allow some NAs but ensure enough valid data
  max_attempts <- 500
  attempt <- 0
  repeat {
    attempt <- attempt + 1
    coordR <- shiftrotXY(coord, X.range = X.range, Y.range = Y.range,
                         rotation = TRUE, X.shift = TRUE, Y.shift = TRUE,
                         mirror = "r", verbose = FALSE, MAXtrial = 10000)
    if (is.null(coordR)) {
      if (attempt > max_attempts) {
        warning("Could not find valid shift+rotation after ", max_attempts, " attempts in iteration ", i)
        break
      }
      next
    }
    # Extract and check that we have enough valid data (at least 50% of sites)
    site_test <- terra::extract(pred_src, coordR)
    site_complete <- stats::complete.cases(site_test[, rf_vars, drop = FALSE])
    valid_prop <- sum(site_complete) / length(site_complete)
    if (valid_prop >= 0.50) break  # Reduced from 70% to 50%
    if (attempt > max_attempts) {
      warning("Could not find shift+rotation with >=50% valid data after ", max_attempts, " attempts in iteration ", i)
      break
    }
  }
  
  if (attempt > max_attempts) {
    cat("Skipping iteration", i, "due to insufficient valid data\n")
    next
  }
  
  # 2) fit null RF using SAME train/test split as observed model
  # Extract predictors at transformed coordinates
  site_df <- terra::extract(pred_src, coordR)
  site_df <- site_df[, rf_vars, drop = FALSE]
  if ("class_code" %in% names(site_df)) {
    site_df$class_code <- factor(as.character(site_df$class_code), levels = class_levels)
  }
  
  site_df$richness <- env$richness
  
  # Filter to complete cases and apply SAME train indices as observed model
  complete_rows <- stats::complete.cases(site_df[, rf_vars, drop = FALSE])
  if (sum(complete_rows) < 10) {
    cat("Skipping iteration", i, "- too few complete cases:", sum(complete_rows), "\n")
    next
  }
  
  # Use SAME train indices as observed model 
  train_idx_in_complete <- intersect(train_indices, which(complete_rows))
  if (length(train_idx_in_complete) < 5) {
    cat("Skipping iteration", i, "- too few training cases after filtering NAs\n")
    next
  }
  
  null_train_data <- site_df[train_idx_in_complete, , drop = FALSE]
  
  # Apply SAME weighting strategy as observed model
  null_weights <- NULL
  if ("class_code" %in% names(null_train_data) && nrow(null_train_data) > 0) {
    tab_null <- table(null_train_data$class_code)
    inv_null <- 1 / as.numeric(tab_null[null_train_data$class_code])
    null_weights <- as.numeric(inv_null / mean(inv_null))
  }
  
  set.seed(4242 + i)
  rf_null_i <- randomForest(
    richness ~ mean_temp + ndvi_max + annual_precip + class_code,
    data = null_train_data, ntree = 500, importance = FALSE, weights = null_weights
  )
  
  # 3) apply the SAME transform to the grid, extract predictors, predict
  Tsr <- estT(coord, coordR)
  xy_moved <- applyT(grid_xy, Tsr)
  vals <- terra::extract(pred_src, as.matrix(xy_moved))
  vals <- vals[, rf_vars, drop = FALSE]
  if ("class_code" %in% names(vals)) {
    vals$class_code <- factor(as.character(vals$class_code), levels = class_levels)
  }
  
  ok <- stats::complete.cases(vals)
  pred_full_vec <- rep(NA_real_, nrow(vals))
  if (sum(ok) > 0) {
    pred_full_vec[ok] <- predict(rf_null_i, newdata = vals[ok, ])
  }
  
  # mask to observed support and compute stats
  null_r <- template; terra::values(null_r) <- pred_full_vec
  null_r <- terra::mask(null_r, msk)
  v_null <- as.vector(null_r[])
  
  # Only compute stats if we have enough non-NA values
  if (sum(!is.na(v_null)) > 100) {
    null_rho[i] <- suppressWarnings(cor(v_lcbd, v_null, method = "spearman", use = "complete.obs"))
    null_ovl[i, ] <- overlap_curve(v_lcbd, v_null, qs)
  } else {
    cat("Iteration", i, "- insufficient non-NA predictions:", sum(!is.na(v_null)), "\n")
  }
  
  if (i %% 10 == 0) cat(".. completed", i, "iterations\n")
}
cat("\nS&R nulls done:", nsim_snr, "\n")

# Sanity check: compare original vs rotated points
png("02_Analysis/03_Figures/Original vs SnR points.png")
par(mfrow = c(1,1))
plot(pred_src[[1]], colNA = "grey90")
points(coord[,1],   coord[,2],   pch = 1,  col = "black")  # original
points(coordR[,1],  coordR[,2],  pch = 16, col = "red")    # rotated
legend("bottom", c("original","rotated"), pch = c(1,16), col = c("black","red"))
dev.off()

# ---- Envelopes + p-values ----
p_rho_repulsion <- (1 + sum(null_rho <= obs_rho)) / (nsim_snr + 1)
# Counts how many null ρ are ≤ the observed ρ # Can be flipped to the opposite tail
# Adds 1 to numerator and denominator = rank test with a +1 correction (avoids 0 or 1 p-values)
# Divides by total permutations nsim_snr + 1 → a one-sided permutation p-value

null_mean <- apply(null_ovl, 2, mean, na.rm = TRUE) # the average null overlap at that quantile
null_lo   <- apply(null_ovl, 2, quantile, probs = 0.025, na.rm = TRUE) # the 2.5% lower envelope (pointwise)
null_hi   <- apply(null_ovl, 2, quantile, probs = 0.975, na.rm = TRUE) # the 97.5% upper envelope (pointwise)

sd_q        <- apply(null_ovl, 2, sd, na.rm = TRUE) # null SD overlap at each quantile (column-wise)
z_sum_obs   <- sum((null_mean - obs_overlap) / sd_q, na.rm = TRUE) # z-score at each quantile (how far obs is from null mean in SD units)
z_sum_null  <- apply((null_mean - null_ovl) / matrix(sd_q, nrow = nsim_snr, ncol = length(qs), byrow = TRUE),
                     1, function(z) sum(z, na.rm = TRUE)) # sums those z-scores across the whole curve
p_curve     <- (1 + sum(z_sum_null >= z_sum_obs)) / (nsim_snr + 1) # global p-val for the overlap curve

# ---- Save data outputs for later reuse ----
# Create comprehensive results list
snr_results <- list(
  # Null distributions
  null_rho = null_rho,
  null_overlap_matrix = null_ovl,
  quantiles = qs,
  
  # Observed values
  obs_rho = obs_rho,
  obs_overlap = obs_overlap,
  
  # Summary statistics
  null_rho_mean = mean(null_rho, na.rm = TRUE),
  null_rho_sd = sd(null_rho, na.rm = TRUE),
  null_overlap_mean = null_mean,
  null_overlap_lo = null_lo,
  null_overlap_hi = null_hi,
  
  # P-values and test statistics
  p_rho_repulsion = p_rho_repulsion,
  p_curve_global = p_curve,
  z_sum_obs = z_sum_obs,
  z_sum_null = z_sum_null,
  
  # Metadata
  nsim = nsim_snr,
  analysis_type = "LCBD_vs_Richness",
  date = Sys.Date()
)

# Save as RDS (preserves all R data structures)
saveRDS(snr_results, file.path(out_dir_res, "SnR_nulls_results_richness.rds"))

# Save correlation results as CSV
rho_results <- data.frame(
  iteration = 1:nsim_snr,
  null_rho = null_rho,
  stringsAsFactors = FALSE
)
write.csv(rho_results, file.path(out_dir_res, "SnR_nulls_rho_richness.csv"), row.names = FALSE)

# Save overlap curves as CSV
overlap_results <- data.frame(
  quantile = rep(qs, each = nsim_snr),
  iteration = rep(1:nsim_snr, length(qs)),
  null_overlap = as.vector(null_ovl),
  stringsAsFactors = FALSE
)
write.csv(overlap_results, file.path(out_dir_res, "SnR_nulls_overlap_richness.csv"), row.names = FALSE)

# Save summary statistics as CSV
summary_stats <- data.frame(
  metric = c("obs_rho", "null_rho_mean", "null_rho_sd", "p_rho_repulsion", 
             "p_curve_global", "nsim", "valid_nulls_rho", "valid_nulls_overlap"),
  value = c(obs_rho, mean(null_rho, na.rm = TRUE), sd(null_rho, na.rm = TRUE), 
            p_rho_repulsion, p_curve, nsim_snr, 
            sum(!is.na(null_rho)), sum(apply(null_ovl, 1, function(x) any(!is.na(x))))),
  stringsAsFactors = FALSE
)
write.csv(summary_stats, file.path(out_dir_res, "SnR_nulls_summary_stats_richness.csv"), row.names = FALSE)

# Save observed vs null comparison for each quantile
quantile_comparison <- data.frame(
  quantile = qs,
  obs_overlap = obs_overlap,
  null_mean = null_mean,
  null_sd = apply(null_ovl, 2, sd, na.rm = TRUE),
  null_lo = null_lo,
  null_hi = null_hi,
  z_score = (obs_overlap - null_mean) / apply(null_ovl, 2, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(quantile_comparison, file.path(out_dir_res, "SnR_nulls_quantile_comparison_richness.csv"), row.names = FALSE)

cat("Results saved to:\n")
cat("- RDS:", file.path(out_dir_res, "SnR_nulls_results_richness.rds"), "\n")
cat("- CSV files:", file.path(out_dir_res, "SnR_nulls_*_richness.csv"), "\n")

# ---- Plots + summary ----
save_hist(
  x = null_rho, x_obs = obs_rho,
  title = "S&R nulls — rank correlation (LCBD vs Richness)",
  subtitle = sprintf("obs rho = %.3f; p(repulsion) = %.4f; nsim = %d", obs_rho, p_rho_repulsion, nsim_snr),
  xlab = "Spearman rho",
  file = file.path(out_dir_fig, "SnR_nulls_rho_hist.png")
)
dfA <- data.frame(q = qs, obs = obs_overlap, lo = null_lo, hi = null_hi, mean = null_mean)
save_curve(
  df = dfA,
  title = "Hotspot overlap across quantiles (S&R null envelope)",
  subtitle = sprintf("Global p (lower-than-expected overlap) = %.4f; nsim = %d", p_curve, nsim_snr),
  xlab = "Top-quantile threshold",
  ylab = "Jaccard overlap",
  file = file.path(out_dir_fig, "SnR_nulls_overlap_curve.png")
)

sink(file.path(out_dir_res, "SnR_nulls_summary.txt"))
cat("S&R nulls (in-memory)\n")
cat("n null rasters:", nsim_snr, "\n")
cat(sprintf("Observed Spearman rho (LCBD vs observed richness): %.4f\n", obs_rho))
cat(sprintf("p(repulsion, rho): %.4f\n", p_rho_repulsion))
cat("Quantiles: ", paste(qs, collapse = ", "), "\n")
cat(sprintf("Global p (overlap curve, lower-than-expected): %.4f\n", p_curve), "\n")
sink()

# ---- Final publication-style figures ----

# Histogram of rho
df_rho <- data.frame(null_rho = null_rho)
# compute histogram counts manually
h <- hist(df_rho$null_rho, breaks = 40, plot = FALSE)
ymax <- max(h$counts, na.rm = TRUE)

p1 <- ggplot(df_rho, aes(null_rho)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(xintercept = mean(null_rho, na.rm = TRUE), linetype = 2) +
  geom_vline(xintercept = obs_rho, linewidth = 1) +
  annotate("text", x = min(null_rho, na.rm = TRUE), 
           y = ymax, hjust = 0, vjust = -0.5,
           label = sprintf("p = %.4f\nnsim = %d", p_rho_repulsion, nsim_snr)) +
  labs(x = "Spearman rho", y = "Count") +
  theme_minimal()

ggsave(file.path(out_dir_fig, "SnR_nulls_rho_pub.png"), p1,
       width = 7, height = 5, dpi = 300)


# Overlap curve with CIs + pointwise obs
df_plot <- data.frame(q = qs, obs = obs_overlap,
                      lo = null_lo, hi = null_hi, mean = null_mean)

p2 <- ggplot(df_plot, aes(q, obs)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "grey70") +
  geom_line(aes(y = mean), linetype = 2, color = "black") +
  geom_point(color = "#0072B2", size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.01, color = "#0072B2") +
  annotate("text", x = min(df_plot$q), y = max(df_plot$obs, na.rm = TRUE),
           hjust = 0, vjust = -0.5,
           label = sprintf("global p = %.4f\nnsim = %d", p_curve, nsim_snr)) +
  labs(x = "Top-quantile threshold", y = "Jaccard overlap") +
  theme_minimal()

ggsave(file.path(out_dir_fig, "SnR_nulls_overlap_curve_pub.png"), p2,
       width = 7, height = 5, dpi = 300)
# ---- Comments ----
# Overlap function explanation

# 1. For each quantile q, we compute thresholds (t1, t2) so that only the 
#    top q% of values in v1 and v2 are kept.
#
# 2. These thresholds turn continuous vectors into binary sets (A and B), 
#    where TRUE marks membership in the top q%.
#
# 3. The overlap between these sets is measured with the Jaccard index:
#       J = intersection(A,B) / union(A,B)
#
# 4. Intersection = number of positions where both v1 and v2 are in 
#    their respective top q%.
#
# 5. Union = number of positions where at least one of them is in the top q%.
#
# 6. If union = 0 (i.e., both sets are empty), return NA.
#
# 7. The output is a vector of Jaccard similarities, one per quantile in qs.
#
# 8. Interpretation: each value tells how much the two vectors overlap 
#    in their top q% of values. Higher values = more shared "hotspots."