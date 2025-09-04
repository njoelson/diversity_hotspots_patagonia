# Hotspot delineation
# Author: NJ
# Date: 8.8.2025

# Step 0: Clean environment
rm(list = ls())

library(terra)
library(tmap)
library(scico)


# Step 1: Custom colors
scico(10, palette = "bamako")
chosen_colors_S <- c("grey95", "#E1C76D")  # 0 = non-hotspot, 1 = hotspot
chosen_colors_RWR <- c("grey95", "#988C02")  # 0 = non-hotspot, 1 = hotspot
chosen_colors_LCBD <- c("grey95", "#728202")  # 0 = non-hotspot, 1 = hotspot


# Step 2: Load predictions
pred_rich <- rast("02_Analysis/02_Results/02_Predictions/RF_richness_250m.tif")
pred_rwr  <- rast("02_Analysis/02_Results/02_Predictions/RWR_prediction_250m.tif")
pred_lcbd <- rast("02_Analysis/02_Results/02_Predictions/LCBD_prediction_250m.tif")

# Step 3: Calculate thresholds
threshold_rich <- quantile(values(pred_rich), 0.95, na.rm = TRUE)
threshold_rwr  <- quantile(values(pred_rwr),  0.95, na.rm = TRUE)
threshold_lcbd <- quantile(values(pred_lcbd), 0.95, na.rm = TRUE)

# Step 4: Define classification matrices
m_high <- function(thresh) {
  matrix(c(-Inf, thresh, 0,
           thresh, Inf, 1), ncol = 3, byrow = TRUE)
}
m_low <- function(thresh) {
  matrix(c(-Inf, thresh, 1,
           thresh, Inf, 0), ncol = 3, byrow = TRUE)
}

# Step 5: Create binary hotspot rasters
hotspot_rich <- classify(pred_rich, m_high(threshold_rich), include.lowest = TRUE)
hotspot_rwr  <- classify(pred_rwr,  m_high(threshold_rwr),  include.lowest = TRUE)
hotspot_lcbd <- classify(pred_lcbd, m_high(threshold_lcbd), include.lowest = TRUE)

# Step 6: Export binary hotspots
writeRaster(hotspot_rich, "02_Analysis/02_Results/hotspot_rich_5top.tif", overwrite = TRUE)
writeRaster(hotspot_rwr,  "02_Analysis/02_Results/hotspot_rwr_5top.tif",  overwrite = TRUE)
writeRaster(hotspot_lcbd, "02_Analysis/02_Results/hotspot_lcbd_5top.tif", overwrite = TRUE)

# Step 7: Congruence across all three   
hotspot_stack  <- c(hotspot_rich, hotspot_rwr, hotspot_lcbd)
congruence_all <- sum(hotspot_stack, na.rm = TRUE)  # 0..3
writeRaster(congruence_all, "02_Analysis/02_Results/hotspot_congruence_sum.tif", overwrite = TRUE)

# Step 8: Plotting function 
plot_hotspot <- function(r, title, file, palette) {
  tmap_save(
    tm_shape(r) +
      tm_raster(
        col.scale  = tm_scale_categorical(values = palette, labels = c("Non-hotspot", "Hotspot")),
        col.legend = tm_legend(
          title    = title,
          position = c("left","top"),   # place inside, top-left
          just     = c("left","top"),
          bg.color = NA,                # no legend background
          frame.lwd= 0                  # no legend frame
        )
      ) +
      tm_layout(
        legend.outside = FALSE,         # legend inside the map
        outer.margins  = 0,             # use the whole page
        inner.margins  = c(0.00, 0.12, 0.10, 0.00),  # reserve space left/top
        frame.lwd      = 1.5            # thicker map frame (optional)
      ),
    filename = file,
    width = 8, height = 6
  )
}

# Step 9: Plot all binary hotspots
plot_hotspot(hotspot_rich,"S", "02_Analysis/03_Figures/hotspot_rich.png", chosen_colors_S)
plot_hotspot(hotspot_rwr, "RWR", "02_Analysis/03_Figures/hotspot_rwr.png", chosen_colors_RWR)
plot_hotspot(hotspot_lcbd, "LCBD", "02_Analysis/03_Figures/hotspot_lcbd.png", chosen_colors_LCBD)

# ---- Step 10: Combination map (colors for each intersection) ----
# Bitmask: Rich=1, RWR=2, LCBD=4  → 0..7 encodes all combos
combo <- hotspot_rich*1 + hotspot_rwr*2 + hotspot_lcbd*4

combo_labels <- c(
  "0"="None",
  "1"="S",
  "2"="RWR",
  "3"="S ∩ RWR",
  "4"="LCBD",
  "5"="S ∩ LCBD",
  "6"="RWR ∩ LCBD",
  "7"="S ∩ RWR ∩ LCBD"
)

combo_palette <- c(
  "0"="grey95",   # None
  "1"="#E1C76D",  # S 
  "2"="#988C02",  # RWR
  "3"="tomato",  # S ∩ RWR
  "4"="#728202",  # LCBD
  "5"="white",    # S ∩ LCBD
  "6"="gold2",    # RWR ∩ LCBD
  "7"="white"     # S ∩ RWR ∩ LCBD
)

# Step 2: Factor and drop unused
combo_f <- as.factor(combo)
combo_f <- droplevels(combo_f)

# Step 3: Desired order (legend order)
desired_order <- c("0","1","2","4","3","6","5","7")

# Restrict to levels actually present
# extract current levels table
lev <- levels(combo_f)[[1]]

# keep only those in desired_order (and reorder accordingly)
lev <- lev[match(desired_order, lev$ID), , drop = FALSE]
lev <- lev[!is.na(lev$ID), , drop = FALSE]

# reassign levels in new order
levels(combo_f) <- lev
order_used <- c(1:5,7)

# check
levels(combo_f)[[1]]$ID

# Step 4: Align labels & palette to this order
pal_used <- unname(combo_palette[order_used])
lab_used <- unname(combo_labels[order_used])

# Step 5: Plot
p_combo <-
  tm_shape(combo_f) +
  tm_raster(
    col.scale  = tm_scale_categorical(values = pal_used, labels = lab_used),
    col.legend = tm_legend(
      title    = "Hotspot overlap",
      position = c("left","top"),   # place inside, top-left
      just     = c("left","top"),
      bg.color = NA,                # no legend background
      frame.lwd= 0                  # no legend frame
    )
  ) +
  tm_layout(
    legend.outside = FALSE,         # legend inside the map
    outer.margins  = 0,             # use the whole page
    inner.margins  = c(0.00, 0.12, 0.10, 0.00),  # reserve space left/top
    frame.lwd      = 1.5            # thicker map frame (optional)
  )

tmap_save(
  p_combo,
  filename = "02_Analysis/03_Figures/hotspot_intersections_combo.png",
  width = 8, height = 6, units = "in", dpi = 300
)

# Save GeoTIFF of combo classes
writeRaster(combo, "02_Analysis/02_Results/hotspot_combination_classes.tif", overwrite = TRUE)

# ---- Step 11: Pairwise TSS among hotspot rasters ----
# Helper to compute TSS between two binary rasters (0/1)
tss_pair <- function(pred, obs) {
  vals <- values(c(pred, obs), mat = TRUE)
  if (is.null(vals)) return(NA_real_)
  ok <- stats::complete.cases(vals)
  if (!any(ok)) return(NA_real_)
  v <- vals[ok, , drop = FALSE]
  p <- v[,1]; o <- v[,2]
  TP <- sum(p == 1 & o == 1)
  TN <- sum(p == 0 & o == 0)
  FP <- sum(p == 1 & o == 0)
  FN <- sum(p == 0 & o == 1)
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  if (is.na(sens) || is.na(spec)) return(NA_real_)
  sens + spec - 1
}

# Collect rasters in a named list
hotspots <- list(
  Richness = hotspot_rich,
  RWR      = hotspot_rwr,
  LCBD     = hotspot_lcbd
)

# Compute pairwise TSS matrix (prediction = row, observation = column)
metric_names <- names(hotspots)
tss_mat <- matrix(NA_real_, nrow = length(hotspots), ncol = length(hotspots),
                  dimnames = list(metric_names, metric_names))
for (i in seq_along(hotspots)) {
  for (j in seq_along(hotspots)) {
    tss_mat[i, j] <- tss_pair(hotspots[[i]], hotspots[[j]])
  }
}

# Save pairwise TSS
tss_mat_df <- as.data.frame(tss_mat)
tss_mat_df$Prediction <- rownames(tss_mat_df)
tss_mat_df <- tss_mat_df[, c(ncol(tss_mat_df), 1:(ncol(tss_mat_df)-1))]
write.csv(tss_mat_df, "02_Analysis/02_Results/hotspot_pairwise_TSS_2.csv", row.names = FALSE)

# ---- Step 11b: Pairwise % overlap table (Jaccard) ----  
pair_overlap_stats <- function(A, B, nameA, nameB) {
  vals <- values(c(A, B), mat = TRUE)
  ok <- stats::complete.cases(vals)
  if (!any(ok)) return(data.frame())
  v <- vals[ok, , drop = FALSE]
  a <- v[,1]; b <- v[,2]
  nA  <- sum(a == 1)
  nB  <- sum(b == 1)
  inter <- sum(a == 1 & b == 1)
  uni   <- sum(a == 1 | b == 1)
  jaccard_pct <- if (uni > 0) 100 * inter / uni else NA_real_
  dice_pct    <- if ((nA + nB) > 0) 100 * (2 * inter) / (nA + nB) else NA_real_
  overlapA    <- if (nA > 0) 100 * inter / nA else NA_real_
  overlapB    <- if (nB > 0) 100 * inter / nB else NA_real_
  data.frame(
    Pair        = paste0(nameA, "_x_", nameB),
    A           = nameA,
    B           = nameB,
    Cells_A     = nA,
    Cells_B     = nB,
    Cells_Inter = inter,
    Cells_Union = uni,
    Jaccard_pct = jaccard_pct,
    Dice_pct    = dice_pct,
    OverlapA_pct= overlapA,   # % of A that overlaps B
    OverlapB_pct= overlapB,   # % of B that overlaps A
    stringsAsFactors = FALSE
  )
}

dir.create("02_Analysis/02_Results/overlaps", showWarnings = FALSE, recursive = TRUE)
dir.create("02_Analysis/03_Figures/overlaps", showWarnings = FALSE, recursive = TRUE)

pairs <- combn(metric_names, 2, simplify = FALSE)
overlap_table <- do.call(rbind, lapply(pairs, function(pr) {
  pair_overlap_stats(hotspots[[pr[1]]], hotspots[[pr[2]]], pr[1], pr[2])
}))
write.csv(overlap_table, "02_Analysis/02_Results/hotspot_pairwise_overlap_percent.csv", row.names = FALSE)

# ---- Step 12: Consensus (combined metrics) & TSS vs. each single metric ----
# Build binary consensus rasters for thresholds k = 1..3 (cells with >= k hotspots)
consensus_list <- lapply(1:3, function(k) {
  classify(congruence_all, matrix(c(-Inf, k, 0,
                                    k,   Inf, 1), ncol = 3, byrow = TRUE),
           include.lowest = TRUE)
})
names(consensus_list) <- paste0("Consensus_ge_", 1:3)

# Compute TSS of each consensus (prediction) against each single metric (observation)
consensus_tss <- do.call(rbind, lapply(names(consensus_list), function(cn) {
  cons <- consensus_list[[cn]]
  vals <- sapply(hotspots, function(h) tss_pair(cons, h))
  data.frame(Consensus = cn, Metric = names(vals), TSS = as.numeric(vals), row.names = NULL)
}))
write.csv(consensus_tss, "02_Analysis/02_Results/hotspot_consensus_vs_metric_TSS_2.csv", row.names = FALSE)

# Also save the consensus rasters
for (cn in names(consensus_list)) {
  writeRaster(consensus_list[[cn]],
              filename = file.path("02_Analysis/02_Results", paste0(tolower(cn), ".tif")),
              overwrite = TRUE)
}

# ---- Step 13: Pairwise overlap maps (metric A ∩ metric B) ----
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  rA <- hotspots[[a]]
  rB <- hotspots[[b]]
  
  # 0/1 AND overlap
  overlap <- rA * rB
  overlap <- classify(
    overlap,
    matrix(c(-Inf, 0.5, 0,
             0.5,  Inf, 1),
           ncol = 3, byrow = TRUE),
    include.lowest = TRUE
  )
  
  # choose the hotspot color for the "1" class based on the pair
  # (Richness == S)
  if ((a == "Richness" && b == "RWR") || (a == "RWR" && b == "Richness")) {
    one_color <- "tomato"              # S ∩ RWR
  } else if ((a == "RWR" && b == "LCBD") || (a == "LCBD" && b == "RWR")) {
    one_color <- "gold1"               # RWR ∩ LCBD
  } else if ((a == "Richness" && b == "LCBD") || (a == "LCBD" && b == "Richness")) {
    
    one_color <- unname(combo_palette["5"])  # currently "white"
  } else {
    one_color <- "#000000"              # only if there is an error will appear
  }
  
  # 2-color palette for binary overlap map: non-hotspot (0) + overlap (1)
  ov_palette <- c("grey95", one_color)
  
  out_tif <- file.path("02_Analysis/02_Results/overlaps",
                       paste0("overlap_", tolower(a), "_", tolower(b), ".tif"))
  writeRaster(overlap, out_tif, overwrite = TRUE)
  
  out_png <- file.path("02_Analysis/03_Figures/overlaps",
                       paste0("overlap_", tolower(a), "_", tolower(b), "_custom.png"))
  
  plot_hotspot(
    overlap,
    title   = paste0(a, " ∩ ", b),  # legend title only 
    file    = out_png,
    palette = ov_palette
  )
}

# ---- Step 13: Pairwise overlap maps (metric A ∩ metric B) ----


for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  rA <- hotspots[[a]]
  rB <- hotspots[[b]]
  
  # Binary overlap (0/1)
  overlap <- classify(
    rA * rB,
    matrix(c(-Inf, 0.5, 0,
             0.5,  Inf, 1),
           ncol = 3, byrow = TRUE),
    include.lowest = TRUE
  )
  
  # ----- compute stats: TSS and Jaccard % overlap -----
  vals <- values(c(rA, rB), mat = TRUE)
  ok <- stats::complete.cases(vals)
  v  <- vals[ok, , drop = FALSE]
  a1 <- v[, 1]; b1 <- v[, 2]
  inter <- sum(a1 == 1 & b1 == 1)
  uni   <- sum(a1 == 1 | b1 == 1)
  jaccard_pct <- if (uni > 0) 100 * inter / uni else NA_real_
  tss_val <- tss_pair(rA, rB)
  
  ann_text <- sprintf("TSS = %.2f\nJaccard = %.1f%%",
                      tss_val, jaccard_pct)
  
  # ----- choose overlap color for the '1' class -----
  # Richness == S
  if ((a == "Richness" && b == "RWR") || (a == "RWR" && b == "Richness")) {
    one_color <- "tomato"     # S ∩ RWR
  } else if ((a == "RWR" && b == "LCBD") || (a == "LCBD" && b == "RWR")) {
    one_color <- "gold1"      # RWR ∩ LCBD
  } else if ((a == "Richness" && b == "LCBD") || (a == "LCBD" && b == "Richness")) {
    
    one_color <- "#80476D"    
  } else {
    one_color <- "#000000"    # fallback (shouldn't happen)
  }
  ov_palette <- c("grey95", one_color)
  ov_labels  <- c("None", paste0(a, " ∩ ", b))
  
  # ---------- PLOT -----------------
  p <-
    tm_shape(overlap) +
    tm_raster(
      col.scale  = tm_scale_categorical(values = ov_palette, labels = ov_labels),
      col.legend = tm_legend(
        title    = paste0( a, " ∩ ", b),
        position = c("left", "top"),
        just     = c("left", "top"),
        bg.color = NA,
        frame.lwd= 0
      )
    ) +
    # Big annotation bottom-right (inside map)
    tm_credits(
      text     = ann_text,
      position = c("left", "bottom"),
      just     = c("left", "bottom"),
      size     = 1.1
    ) +
    tm_layout(
      legend.outside = FALSE,
      outer.margins  = 0,
      inner.margins  = c(0.00, 0.12, 0.10, 0.00),
      frame.lwd      = 1.5
    )
  
  # save
  out_tif <- file.path("02_Analysis/02_Results/overlaps",
                       paste0("overlap_", tolower(a), "_", tolower(b), ".tif"))
  writeRaster(overlap, out_tif, overwrite = TRUE)
  
  out_png <- file.path("02_Analysis/03_Figures/overlaps",
                       paste0("overlap_", tolower(a), "_", tolower(b), ".png"))
  tmap_save(p, filename = out_png, width = 8, height = 6, units = "in", dpi = 300)
}

