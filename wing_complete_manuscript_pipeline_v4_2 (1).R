# ============================================================================
# COMPLETE MANUSCRIPT PIPELINE
# Drosophila adult wing: A/P proportion, regional geometry, whole-wing shape,
# fixed-anchor sensitivity, and A/P-associated boundary geometry
#
# Input format:
#   one row = one adult wing
#   px1/py1 ... px112/py112 = DataMartin coordinates
#
# Biological unit:
#   one fly = mean of the left and right wings
#
# Main manuscript figure generated:
#   A. DataMartin A0-A8 scheme + A/P-associated boundary
#   B. Posterior / whole-wing area
#   C. A3 and A4 absolute regional areas
#   D. Size-adjusted whole-wing deformation map
#
# Supplementary analyses generated:
#   - A0-A8 regional statistics with BH-FDR
#   - Global A0-A8 compositional test (CLR + RRPP)
#   - Sliding-semilandmark Procrustes ANOVA
#   - Fixed-anchor-only Procrustes ANOVA sensitivity analysis
#   - Size-adjusted PCA
#   - A/P-associated boundary scalar metrics
#   - Equal-arc-length-resampled, size-adjusted boundary RRPP
#   - Boundary difference profile
#
# Interpretation note:
#   The adult DataMartin A/P-associated geometric boundary is NOT a direct
#   visualization of the larval lineage/compartment boundary.
# ============================================================================

# ----------------------------------------------------------------------------
# 0. PACKAGES
# ----------------------------------------------------------------------------
# Install once if necessary:
# install.packages(c(
#   "readxl", "dplyr", "tidyr", "ggplot2", "ggbeeswarm",
#   "patchwork", "geomorph", "RRPP"
# ))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
library(geomorph)
library(RRPP)

# ----------------------------------------------------------------------------
# 1. USER SETTINGS
# ----------------------------------------------------------------------------
#remove previous datasets
rm(list=ls())

#set working directory
setwd("/Volumes/IODATASSD/2024 Igaki lab/Data/2609/260902 nub p35 miRHG egg30 wing size")

FILE  <- "results.csv"
SHEET <- "Sheet1"

CONTROL <- "nubGFP"
MUTANT  <- "nubp35"      # change to the final genotype name

N_PERM <- 9999
OUTDIR <- "wing_complete_manuscript_results"
dir.create(OUTDIR, showWarnings = FALSE)

# Figure settings
FIG_WIDTH_MM  <- 180
FIG_HEIGHT_MM <- 145
DEFORMATION_MAG <- 1

# Panel D display settings
PANEL_D_CTRL_COLOR   <- "#5A5A5A"
PANEL_D_MUT_COLOR    <- "#D55E00"
PANEL_D_VECTOR_COLOR <- "#1F4E79"
N_MAJOR_VECTORS <- 6   # largest fixed-anchor displacement vectors shown in simplified versions

# User previously requested the wing drawings in Panels A and D to be vertically flipped.
FLIP_WING_VERTICAL <- TRUE

# Optional secondary analyses
RUN_SIZE_ADJUSTED_PCA <- TRUE
RUN_CLR_COMPOSITION <- TRUE
RUN_BOUNDARY_ANALYSIS <- TRUE

# Boundary resampling resolution
N_BOUNDARY_RESAMPLE <- 51

# ----------------------------------------------------------------------------
# 2. FIGURE THEME AND HELPERS
# ----------------------------------------------------------------------------
paper_theme <- theme_classic(base_size = 15, base_family = "Arial") +
  theme(
    axis.title.x = element_text(size = 10, margin = margin(t = 3)),
    axis.title.y = element_text(size = 10, margin = margin(r = 3)),
    axis.text.x  = element_text(size = 6.5, lineheight = 0.9),
    axis.text.y  = element_text(size = 8),
    axis.line    = element_line(linewidth = 0.7),
    axis.ticks   = element_line(linewidth = 0.7),
    legend.title = element_blank(),
    plot.title   = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 7.5)
  )

format_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.0001) return("P < 0.0001")
  if (p < 0.001) return(paste0("P = ", formatC(p, format = "f", digits = 4)))
  if (p < 0.01)  return(paste0("P = ", formatC(p, format = "f", digits = 3)))
  return(paste0("P = ", formatC(p, format = "f", digits = 2)))
}

plot_y <- function(y) {
  if (FLIP_WING_VERTICAL) y else -y
}

rrpp_table <- function(fit) {
  a <- anova(fit)
  if (!is.null(a$table)) {
    tab <- as.data.frame(a$table)
  } else {
    tab <- as.data.frame(a)
  }
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab
}

# ----------------------------------------------------------------------------
# 3. DATAMARTIN TOPOLOGY
# ----------------------------------------------------------------------------
# 15 fixed topological/anatomical anchors
fixed_idx <- c(
  1, 5, 13, 17, 23, 26, 38, 42, 53, 65, 75, 78, 89, 97, 102
)

# Continuous curves; internal points are sliding semilandmarks.
curve_sets <- list(
  c(1, 2, 3, 4, 5),
  c(5, 6, 7, 8, 9, 10, 11, 12, 13),
  c(13, 18, 19, 20, 21, 22, 23),
  c(23, 24, 25, 26),
  c(26, 90, 91, 92, 93, 94, 95, 96, 97),
  c(97, 112, 111, 110, 109, 108, 107, 106, 105, 104, 103, 102),
  c(5, 14, 15, 16, 17),
  c(38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 13),
  c(38, 39, 40, 41, 42),
  c(42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 23),
  c(42, 52, 53),
  c(53, 54, 55, 56, 57, 58, 59, 60, 75),
  c(75, 76, 77, 78),
  c(78, 79, 80, 81, 82, 83, 84, 85, 89),
  c(53, 61, 62, 63, 64, 65),
  c(65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 26),
  c(65, 86, 87, 88, 89),
  c(89, 100, 99, 98, 97),
  c(78, 101, 102)
)

# Fixed-to-fixed graph edges without internal semilandmarks
fixed_edges <- list(
  c(1, 17),
  c(17, 38),
  c(17, 75),
  c(75, 102)
)

make_curve_triples <- function(v) {
  if (length(v) < 3) stop("Each semilandmark curve must contain >=3 points.")
  cbind(
    v[1:(length(v) - 2)],
    v[2:(length(v) - 1)],
    v[3:length(v)]
  )
}

curves_mat <- do.call(rbind, lapply(curve_sets, make_curve_triples))
colnames(curves_mat) <- c("before", "slide", "after")

slider_idx <- sort(unique(curves_mat[, "slide"]))
stopifnot(length(fixed_idx) == 15)
stopifnot(length(slider_idx) == 97)
stopifnot(length(intersect(fixed_idx, slider_idx)) == 0)
stopifnot(setequal(c(fixed_idx, slider_idx), 1:112))

# Adult A/P-associated DataMartin boundary:
# anterior  = A0 + A1 + A2 + A3 + A5
# posterior = A4 + A6 + A7 + A8
boundary_idx <- c(
  17, 75,
  60, 59, 58, 57, 56, 55, 54, 53,
  61, 62, 63, 64, 65,
  66, 67, 68, 69, 70, 71, 72, 73, 74, 26
)

boundary_fixed_nodes <- c(17, 75, 53, 65, 26)

# ----------------------------------------------------------------------------
# 4. READ AND VALIDATE DATA
# ----------------------------------------------------------------------------
dat <- read_xlsx("results.xlsx")

dat <- dat %>%
  filter(Genotype %in% c("nubGFP","nubp35")) %>%
  filter(date != "260902")

required_cols <- c(
  "Genotype", "Sex", "LR", "sample", "filename", "width", "height",
  "whole", "anterior", "posterior", "Whole(px)",
  paste0("Area ", 0:8),
  paste0("px", 1:112), paste0("py", 1:112)
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat <- dat %>%
  mutate(
    Genotype = factor(Genotype, levels = c(CONTROL, MUTANT)),
    LR = factor(tolower(LR), levels = c("left", "right")),
    fly_id = paste(as.character(Genotype), sample, sep = "_")
  )

if (any(is.na(dat$Genotype))) stop("Unexpected genotype label. Check CONTROL/MUTANT.")
if (any(is.na(dat$LR))) stop("LR must contain left/right.")

sex_levels <- unique(na.omit(as.character(dat$Sex)))
if (length(sex_levels) > 1) {
  warning("More than one sex is present. Filter to one sex or model Sex explicitly.")
}

pair_check <- dat %>%
  distinct(fly_id, LR) %>%
  count(fly_id, name = "n_sides")

complete_ids <- pair_check %>%
  filter(n_sides == 2) %>%
  pull(fly_id)

dat_pair <- dat %>%
  filter(fly_id %in% complete_ids) %>%
  arrange(Genotype, sample, LR)

side_check <- dat_pair %>% count(fly_id, LR)
if (any(side_check$n != 1)) stop("Duplicate rows detected for at least one fly/side.")

cat("\nComplete bilateral pairs:\n")
print(dat_pair %>% distinct(fly_id, Genotype) %>% count(Genotype))

# ----------------------------------------------------------------------------
# 5. COORDINATE ARRAY AND L/R ORIENTATION MATCHING
# ----------------------------------------------------------------------------
coords <- array(
  NA_real_,
  dim = c(112, 2, nrow(dat_pair)),
  dimnames = list(paste0("P", 1:112), c("x", "y"), NULL)
)


for (j in seq_len(nrow(dat_pair))) {
  coords[, 1, j] <- as.numeric(unlist(dat_pair[j, paste0("px", 1:112)]))
  coords[, 2, j] <- as.numeric(unlist(dat_pair[j, paste0("py", 1:112)]))

  # Match handedness of left/right wings.
  if (dat_pair$LR[j] == "left") {
    coords[, 2, j] <- dat_pair$height[j] - coords[, 2, j]
  }
}

if (any(!is.finite(coords))) stop("Non-finite coordinate value detected.")

# ----------------------------------------------------------------------------
# 6. PER-FLY AREA DATA
# ----------------------------------------------------------------------------
# Convert DataMartin A0-A8 pixel regions to the physical area unit of 'whole'.
for (i in 0:8) {
  dat_pair[[paste0("A", i, "_area")]] <-
    dat_pair[[paste0("Area ", i)]] / dat_pair[["Whole(px)"]] * dat_pair$whole
}

fly_area <- dat_pair %>%
  group_by(fly_id, Genotype) %>%
  summarise(
    whole = mean(whole, na.rm = TRUE),
    anterior = mean(anterior, na.rm = TRUE),
    posterior = mean(posterior, na.rm = TRUE),
    across(all_of(paste0("A", 0:8, "_area")), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    P_area_fraction = posterior / whole,
    P_area_percent = 100 * P_area_fraction,
    # Secondary/local descriptor; interpret as regional allocation, not a new primary endpoint.
    A34_P_fraction = A4_area / (A3_area + A4_area)
  )

write.csv(fly_area, file.path(OUTDIR, "per_fly_area_data.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 7. PRIMARY SCALAR ENDPOINT: POSTERIOR / WHOLE
# ----------------------------------------------------------------------------
test_Pfrac <- suppressWarnings(
  wilcox.test(P_area_fraction ~ Genotype, data = fly_area, exact = TRUE)
)

capture.output(test_Pfrac, file = file.path(OUTDIR, "posterior_fraction_Wilcoxon.txt"))

cat("\nPosterior / whole-wing area:\n")
print(test_Pfrac)

# ----------------------------------------------------------------------------
# 8. A0-A8 REGIONAL LOCALIZATION + BH-FDR
# ----------------------------------------------------------------------------
region_stats <- lapply(0:8, function(i) {
  nm <- paste0("A", i, "_area")
  x_ctrl <- fly_area[[nm]][fly_area$Genotype == CONTROL]
  x_mut  <- fly_area[[nm]][fly_area$Genotype == MUTANT]
  wt <- suppressWarnings(wilcox.test(x_ctrl, x_mut, exact = TRUE))

  data.frame(
    region = paste0("A", i),
    control_mean = mean(x_ctrl),
    control_sd = sd(x_ctrl),
    mutant_mean = mean(x_mut),
    mutant_sd = sd(x_mut),
    percent_change = 100 * (mean(x_mut) / mean(x_ctrl) - 1),
    W = unname(wt$statistic),
    p_value = wt$p.value
  )
}) %>% bind_rows()

region_stats$FDR <- p.adjust(region_stats$p_value, method = "BH")
write.csv(region_stats, file.path(OUTDIR, "A0_A8_regional_area_statistics.csv"), row.names = FALSE)

cat("\nA0-A8 regional areas:\n")
print(region_stats)

# Secondary local A3/A4 allocation statistic
A34_test <- suppressWarnings(
  wilcox.test(A34_P_fraction ~ Genotype, data = fly_area, exact = TRUE)
)
capture.output(A34_test, file = file.path(OUTDIR, "A3_A4_local_fraction_Wilcoxon.txt"))

# ----------------------------------------------------------------------------
# 9. GLOBAL A0-A8 COMPOSITION (OPTIONAL ROBUSTNESS ANALYSIS)
# ----------------------------------------------------------------------------
if (RUN_CLR_COMPOSITION) {
  area_cols <- paste0("A", 0:8, "_area")
  Pmat <- as.matrix(fly_area[, area_cols])
  Pmat <- Pmat / rowSums(Pmat)

  if (any(Pmat <= 0)) stop("CLR requires strictly positive regional areas.")

  clr_mat <- log(Pmat) - rowMeans(log(Pmat))
  colnames(clr_mat) <- paste0("clr_A", 0:8)

  comp_meta <- data.frame(
    genotype = factor(fly_area$Genotype, levels = c(CONTROL, MUTANT))
  )

  fit_comp <- lm.rrpp(
    clr_mat ~ genotype,
    data = comp_meta,
    iter = N_PERM,
    RRPP = TRUE,
    print.progress = FALSE
  )

  capture.output(
    anova(fit_comp),
    file = file.path(OUTDIR, "A0_A8_global_CLR_RRPP.txt")
  )

  clr_contrast <- data.frame(
    region = paste0("A", 0:8),
    mutant_minus_control_CLR =
      colMeans(clr_mat[comp_meta$genotype == MUTANT, , drop = FALSE]) -
      colMeans(clr_mat[comp_meta$genotype == CONTROL, , drop = FALSE])
  )
  write.csv(
    clr_contrast,
    file.path(OUTDIR, "A0_A8_CLR_genotype_contrast.csv"),
    row.names = FALSE
  )
}

# ----------------------------------------------------------------------------
# 10. SLIDING-SEMILANDMARK GPA
# ----------------------------------------------------------------------------
# Primary dense-morphometric analysis: minimum bending energy sliding.
gpa_all <- gpagen(
  A = coords,
  curves = curves_mat,
  ProcD = FALSE,
  PrinAxes = FALSE,
  Proj = TRUE,
  print.progress = FALSE
)

fly_ids <- unique(dat_pair$fly_id)

fly_mean_shape <- array(
  NA_real_,
  dim = c(112, 2, length(fly_ids)),
  dimnames = list(paste0("P", 1:112), c("x", "y"), fly_ids)
)
fly_logCS <- numeric(length(fly_ids))

for (k in seq_along(fly_ids)) {
  ii <- which(dat_pair$fly_id == fly_ids[k])
  if (length(ii) != 2) stop("Expected exactly two wings for ", fly_ids[k])

  fly_mean_shape[, , k] <- apply(
    gpa_all$coords[, , ii, drop = FALSE],
    c(1, 2), mean
  )

  fly_logCS[k] <- mean(log(gpa_all$Csize[ii]))
}

fly_mean_gpa <- gpagen(
  fly_mean_shape,
  PrinAxes = FALSE,
  Proj = TRUE,
  print.progress = FALSE
)

fly_meta <- dat_pair %>%
  distinct(fly_id, Genotype) %>%
  slice(match(fly_ids, fly_id))

fly_meta$logCS <- fly_logCS
fly_meta$genotype <- factor(fly_meta$Genotype, levels = c(CONTROL, MUTANT))
stopifnot(all(fly_meta$fly_id == fly_ids))

shape_df <- geomorph.data.frame(
  shape = fly_mean_gpa$coords,
  logCS = fly_meta$logCS,
  genotype = fly_meta$genotype
)

# First test whether allometric slopes differ between genotypes.
fit_shape_interaction <- procD.lm(
  shape ~ logCS * genotype,
  data = shape_df,
  iter = N_PERM,
  RRPP = TRUE,
  print.progress = FALSE
)

# Main model: genotype after accounting for centroid size.
fit_shape_additive <- procD.lm(
  shape ~ logCS + genotype,
  data = shape_df,
  iter = N_PERM,
  RRPP = TRUE,
  print.progress = FALSE
)

capture.output(
  {
    cat("SEMILANDMARK ANALYSIS\n")
    cat("Sliding criterion: minimum bending energy\n\n")
    cat("Interaction model\n")
    print(anova(fit_shape_interaction))
    cat("\nAdditive model\n")
    print(anova(fit_shape_additive))
  },
  file = file.path(OUTDIR, "shape_Procrustes_ANOVA_semilandmarks.txt")
)

shape_tab_int <- rrpp_table(fit_shape_interaction)
shape_tab_add <- rrpp_table(fit_shape_additive)

shape_genotype <- shape_tab_add %>% filter(term == "genotype")
shape_size     <- shape_tab_add %>% filter(term == "logCS")
shape_interact <- shape_tab_int %>% filter(term == "logCS:genotype")

if (nrow(shape_genotype) != 1) stop("Could not extract genotype term from shape ANOVA.")

write.csv(
  bind_rows(
    data.frame(model = "additive", shape_tab_add),
    data.frame(model = "interaction", shape_tab_int)
  ),
  file.path(OUTDIR, "shape_Procrustes_ANOVA_table.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 11. FIXED-ANCHOR SENSITIVITY ANALYSIS
# ----------------------------------------------------------------------------
coords_fixed <- coords[fixed_idx, , , drop = FALSE]

gpa_fixed <- gpagen(
  coords_fixed,
  PrinAxes = FALSE,
  Proj = TRUE,
  print.progress = FALSE
)

fixed_mean_shape <- array(
  NA_real_,
  dim = c(length(fixed_idx), 2, length(fly_ids)),
  dimnames = list(paste0("P", fixed_idx), c("x", "y"), fly_ids)
)
fixed_logCS <- numeric(length(fly_ids))

for (k in seq_along(fly_ids)) {
  ii <- which(dat_pair$fly_id == fly_ids[k])
  fixed_mean_shape[, , k] <- apply(
    gpa_fixed$coords[, , ii, drop = FALSE],
    c(1, 2), mean
  )
  fixed_logCS[k] <- mean(log(gpa_fixed$Csize[ii]))
}

fixed_mean_gpa <- gpagen(
  fixed_mean_shape,
  PrinAxes = FALSE,
  Proj = TRUE,
  print.progress = FALSE
)

fixed_df <- geomorph.data.frame(
  shape = fixed_mean_gpa$coords,
  logCS = fixed_logCS,
  genotype = fly_meta$genotype
)

fit_fixed_interaction <- procD.lm(
  shape ~ logCS * genotype,
  data = fixed_df,
  iter = N_PERM,
  RRPP = TRUE,
  print.progress = FALSE
)

fit_fixed_additive <- procD.lm(
  shape ~ logCS + genotype,
  data = fixed_df,
  iter = N_PERM,
  RRPP = TRUE,
  print.progress = FALSE
)

capture.output(
  {
    cat("FIXED-ANCHOR-ONLY SENSITIVITY ANALYSIS\n")
    cat("Fixed nodes: ", paste(fixed_idx, collapse = ", "), "\n\n", sep = "")
    cat("Interaction model\n")
    print(anova(fit_fixed_interaction))
    cat("\nAdditive model\n")
    print(anova(fit_fixed_additive))
  },
  file = file.path(OUTDIR, "shape_Procrustes_ANOVA_fixed_anchors_sensitivity.txt")
)

fixed_tab_add <- rrpp_table(fit_fixed_additive)
fixed_genotype <- fixed_tab_add %>% filter(term == "genotype")

# ----------------------------------------------------------------------------
# 12. SIZE-ADJUSTED PREDICTED WHOLE-WING SHAPES
# ----------------------------------------------------------------------------
Y_shape <- two.d.array(fly_mean_gpa$coords)
model_dat <- data.frame(
  logCS = fly_meta$logCS,
  genotype = fly_meta$genotype
)

X_shape <- model.matrix(~ logCS + genotype, data = model_dat)
B_shape <- qr.solve(X_shape, Y_shape)
common_logCS <- mean(model_dat$logCS)

newdat_shape <- data.frame(
  logCS = c(common_logCS, common_logCS),
  genotype = factor(c(CONTROL, MUTANT), levels = c(CONTROL, MUTANT))
)
Xnew_shape <- model.matrix(~ logCS + genotype, data = newdat_shape)
pred_2d <- Xnew_shape %*% B_shape
pred_arr <- arrayspecs(pred_2d, p = 112, k = 2)

pred_ctrl <- pred_arr[, , 1]
pred_mut  <- pred_arr[, , 2]

normalize_config <- function(M) {
  M <- sweep(M, 2, colMeans(M), "-")
  M / sqrt(sum(M^2))
}

pred_ctrl <- normalize_config(pred_ctrl)
pred_mut  <- normalize_config(pred_mut)
pred_mut_mag <- pred_ctrl + DEFORMATION_MAG * (pred_mut - pred_ctrl)

# ----------------------------------------------------------------------------
# 13. SIZE-ADJUSTED PCA (OPTIONAL)
# ----------------------------------------------------------------------------
if (RUN_SIZE_ADJUSTED_PCA) {
  Ypca <- two.d.array(fly_mean_gpa$coords)
  Xsize <- model.matrix(~ logCS, data = model_dat)
  Bsize <- qr.solve(Xsize, Ypca)
  Yres <- Ypca - Xsize %*% Bsize

  pca_adj <- prcomp(Yres, center = FALSE, scale. = FALSE)
  pc_var <- 100 * (pca_adj$sdev^2) / sum(pca_adj$sdev^2)

  pc_df <- data.frame(
    fly_id = fly_meta$fly_id,
    Genotype = fly_meta$genotype,
    PC1 = pca_adj$x[, 1],
    PC2 = pca_adj$x[, 2]
  )
  write.csv(pc_df, file.path(OUTDIR, "shape_PCA_size_adjusted_scores.csv"), row.names = FALSE)

  pPCA <- ggplot(pc_df, aes(PC1, PC2, shape = Genotype)) +
    geom_point(size = 2.4, alpha = 0.85) +
    labs(
      x = sprintf("Size-adjusted PC1 (%.1f%%)", pc_var[1]),
      y = sprintf("Size-adjusted PC2 (%.1f%%)", pc_var[2])
    ) +
    paper_theme

  ggsave(
    file.path(OUTDIR, "shape_PCA_size_adjusted.pdf"),
    pPCA,
    width = 85, height = 70, units = "mm", device = cairo_pdf
  )
}

# ----------------------------------------------------------------------------
# 14. BOUNDARY ANALYSIS
# ----------------------------------------------------------------------------
# Primary robust boundary model:
#   equal arc-length resampling -> endpoint normalization -> size-adjusted RRPP.
# Scalar metrics:
#   tortuosity, RMS/chord, boundary length/sqrt(wing area).

calc_boundary_metrics <- function(xy, whole_px) {
  dxy <- xy[-1, , drop = FALSE] - xy[-nrow(xy), , drop = FALSE]
  seg <- sqrt(rowSums(dxy^2))
  path_length <- sum(seg)

  A <- xy[1, ]
  B <- xy[nrow(xy), ]
  v <- B - A
  chord <- sqrt(sum(v^2))

  rel <- sweep(xy, 2, A, "-")
  perp <- abs(v[1] * rel[, 2] - v[2] * rel[, 1]) / chord

  # Arc-length weights reduce dependence on nonuniform point density.
  w <- numeric(nrow(xy))
  w[1] <- seg[1] / 2
  w[nrow(xy)] <- seg[length(seg)] / 2
  if (nrow(xy) > 2) {
    w[2:(nrow(xy) - 1)] <- (seg[-length(seg)] + seg[-1]) / 2
  }

  rms_to_chord <- sqrt(sum(w * perp^2) / sum(w))

  data.frame(
    boundary_length_px = path_length,
    chord_length_px = chord,
    tortuosity = path_length / chord,
    rms_to_chord_norm = rms_to_chord / chord,
    boundary_length_sqrtarea = path_length / sqrt(whole_px)
  )
}

resample_by_arc <- function(xy, n = N_BOUNDARY_RESAMPLE) {
  seg <- sqrt(rowSums((xy[-1, , drop = FALSE] - xy[-nrow(xy), , drop = FALSE])^2))
  s <- c(0, cumsum(seg))
  if (max(s) <= 0) stop("Zero-length boundary encountered.")
  s <- s / max(s)
  target <- seq(0, 1, length.out = n)

  cbind(
    approx(s, xy[, 1], xout = target, ties = "ordered")$y,
    approx(s, xy[, 2], xout = target, ties = "ordered")$y
  )
}

normalize_boundary_profile <- function(xy) {
  xy <- sweep(xy, 2, xy[1, ], "-")
  v <- xy[nrow(xy), ]
  chord <- sqrt(sum(v^2))
  theta <- atan2(v[2], v[1])

  ct <- cos(-theta)
  st <- sin(-theta)
  R <- matrix(c(ct, -st, st, ct), nrow = 2, byrow = TRUE)

  out <- xy %*% t(R)
  out / chord
}

if (RUN_BOUNDARY_ANALYSIS) {
  nb <- length(boundary_idx)
  nr <- N_BOUNDARY_RESAMPLE

  boundary_wing_metrics <- vector("list", nrow(dat_pair))
  resampled_profiles <- array(
    NA_real_, dim = c(nr, 2, nrow(dat_pair)),
    dimnames = list(NULL, c("x", "y"), NULL)
  )

  for (j in seq_len(nrow(dat_pair))) {
    xy <- coords[boundary_idx, , j]

    boundary_wing_metrics[[j]] <- cbind(
      data.frame(
        fly_id = dat_pair$fly_id[j],
        Genotype = dat_pair$Genotype[j],
        LR = dat_pair$LR[j]
      ),
      calc_boundary_metrics(xy, as.numeric(dat_pair[["Whole(px)"]][j]))
    )

    xy_resampled <- resample_by_arc(xy, n = nr)
    resampled_profiles[, , j] <- normalize_boundary_profile(xy_resampled)
  }

  boundary_wing_metrics <- bind_rows(boundary_wing_metrics)

  boundary_fly_metrics <- boundary_wing_metrics %>%
    group_by(fly_id, Genotype) %>%
    summarise(
      across(
        c(
          boundary_length_px,
          chord_length_px,
          tortuosity,
          rms_to_chord_norm,
          boundary_length_sqrtarea
        ),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  write.csv(
    boundary_fly_metrics,
    file.path(OUTDIR, "AP_boundary_geometry_per_fly.csv"),
    row.names = FALSE
  )

  # Per-fly mean endpoint-normalized/resampled boundary profiles
  boundary_fly_profiles <- array(
    NA_real_, dim = c(nr, 2, length(fly_ids)),
    dimnames = list(NULL, c("x", "y"), fly_ids)
  )

  for (k in seq_along(fly_ids)) {
    ii <- which(dat_pair$fly_id == fly_ids[k])
    boundary_fly_profiles[, , k] <- apply(
      resampled_profiles[, , ii, drop = FALSE],
      c(1, 2), mean
    )
  }

  # Multivariate response excludes the fixed endpoints.
  Y_boundary <- matrix(
    NA_real_,
    nrow = length(fly_ids),
    ncol = (nr - 2) * 2
  )

  for (k in seq_along(fly_ids)) {
    Y_boundary[k, ] <- as.vector(boundary_fly_profiles[2:(nr - 1), , k])
  }

  boundary_dat <- data.frame(
    logCS_c = fly_meta$logCS - mean(fly_meta$logCS),
    genotype = factor(fly_meta$genotype, levels = c(CONTROL, MUTANT))
  )

  # Interaction check
  fit_boundary_interaction <- lm.rrpp(
    Y_boundary ~ logCS_c * genotype,
    data = boundary_dat,
    iter = N_PERM,
    RRPP = TRUE,
    print.progress = FALSE
  )

  # Main robust boundary model
  fit_boundary_additive <- lm.rrpp(
    Y_boundary ~ logCS_c + genotype,
    data = boundary_dat,
    iter = N_PERM,
    RRPP = TRUE,
    print.progress = FALSE
  )

  capture.output(
    {
      cat("EQUAL ARC-LENGTH RESAMPLED, ENDPOINT-NORMALIZED BOUNDARY\n\n")
      cat("Interaction model\n")
      print(anova(fit_boundary_interaction))
      cat("\nAdditive model\n")
      print(anova(fit_boundary_additive))
    },
    file = file.path(OUTDIR, "AP_boundary_resampled_size_adjusted_RRPP.txt")
  )

  boundary_tab_int <- rrpp_table(fit_boundary_interaction)
  boundary_tab_add <- rrpp_table(fit_boundary_additive)
  boundary_genotype <- boundary_tab_add %>% filter(term == "genotype")
  boundary_size <- boundary_tab_add %>% filter(term == "logCS_c")
  boundary_interact <- boundary_tab_int %>% filter(term == "logCS_c:genotype")

  # Scalar metrics: three prespecified summaries, BH-FDR corrected.
  metric_names <- c(
    "tortuosity",
    "rms_to_chord_norm",
    "boundary_length_sqrtarea"
  )

  boundary_metric_stats <- lapply(metric_names, function(nm) {
    x0 <- boundary_fly_metrics[[nm]][boundary_fly_metrics$Genotype == CONTROL]
    x1 <- boundary_fly_metrics[[nm]][boundary_fly_metrics$Genotype == MUTANT]
    wt <- suppressWarnings(wilcox.test(x0, x1, exact = TRUE))

    data.frame(
      metric = nm,
      control_mean = mean(x0),
      control_sd = sd(x0),
      mutant_mean = mean(x1),
      mutant_sd = sd(x1),
      percent_change = 100 * (mean(x1) / mean(x0) - 1),
      W = unname(wt$statistic),
      p_value = wt$p.value
    )
  }) %>%
    bind_rows() %>%
    mutate(FDR = p.adjust(p_value, method = "BH"))

  write.csv(
    boundary_metric_stats,
    file.path(OUTDIR, "AP_boundary_scalar_metric_statistics.csv"),
    row.names = FALSE
  )

  # Size-adjusted predicted boundary profiles at mean logCS.
  Xb <- model.matrix(~ logCS_c + genotype, data = boundary_dat)
  Bb <- qr.solve(Xb, Y_boundary)

  newb <- data.frame(
    logCS_c = c(0, 0),
    genotype = factor(c(CONTROL, MUTANT), levels = c(CONTROL, MUTANT))
  )
  Xb_new <- model.matrix(~ logCS_c + genotype, data = newb)
  pred_b <- Xb_new %*% Bb

  # Reinsert the fixed endpoints (0,0) and (1,0).
  vec_to_profile <- function(v, n = nr) {
    mid <- matrix(v, nrow = n - 2, ncol = 2)
    rbind(c(0, 0), mid, c(1, 0))
  }

  boundary_pred_ctrl <- vec_to_profile(pred_b[1, ])
  boundary_pred_mut  <- vec_to_profile(pred_b[2, ])

  s_norm <- seq(0, 1, length.out = nr)

  # Difference in transverse displacement (mutant - control)
  difference_df <- data.frame(
    arc_length = s_norm,
    delta_transverse = boundary_pred_mut[, 2] - boundary_pred_ctrl[, 2]
  )

  pBoundaryDiff <- ggplot(difference_df, aes(arc_length, delta_transverse)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.7) +
    geom_line(linewidth = 1.0) +
    labs(
      title = paste0(MUTANT, " - ", CONTROL, " boundary displacement"),
      subtitle = "Size-adjusted difference after equal arc-length resampling",
      x = "Normalized arc length",
      y = expression(Delta ~ "transverse displacement")
    ) +
    paper_theme

  ggsave(
    file.path(OUTDIR, "AP_boundary_difference_profile_size_adjusted.pdf"),
    pBoundaryDiff,
    width = 90, height = 70, units = "mm", device = cairo_pdf
  )

  # Size-adjusted mean overlay
  boundary_overlay <- bind_rows(
    data.frame(
      x = boundary_pred_ctrl[, 1],
      y = plot_y(boundary_pred_ctrl[, 2]),
      Genotype = CONTROL
    ),
    data.frame(
      x = boundary_pred_mut[, 1],
      y = plot_y(boundary_pred_mut[, 2]),
      Genotype = MUTANT
    )
  )

  pBoundaryOverlay <- ggplot(
    boundary_overlay,
    aes(x, y, linetype = Genotype, group = Genotype)
  ) +
    geom_path(linewidth = 0.9) +
    coord_cartesian() +
    labs(
      title = "Adult A/P-associated boundary geometry",
      subtitle = "Equal arc-length resampling; endpoint-normalized; size-adjusted",
      x = "Normalized boundary axis",
      y = "Normalized transverse displacement"
    ) +
    paper_theme +
    theme(legend.position = "top")

  ggsave(
    file.path(OUTDIR, "AP_boundary_mean_profile_size_adjusted.pdf"),
    pBoundaryOverlay,
    width = 105, height = 70, units = "mm", device = cairo_pdf
  )

  # Scalar metric figure
  metrics_long <- boundary_fly_metrics %>%
    select(fly_id, Genotype, all_of(metric_names)) %>%
    pivot_longer(
      cols = all_of(metric_names),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = metric_names,
        labels = c(
          "Tortuosity (length/chord)",
          "RMS deviation/chord",
          "Boundary length/sqrt(wing area)"
        )
      )
    )

  boundary_ann <- boundary_metric_stats %>%
    mutate(
      metric = factor(
        metric,
        levels = metric_names,
        labels = c(
          "Tortuosity (length/chord)",
          "RMS deviation/chord",
          "Boundary length/sqrt(wing area)"
        )
      ),
      label = paste0("BH-FDR ", sapply(FDR, format_p))
    )

  boundary_ann_y <- metrics_long %>%
    group_by(metric) %>%
    summarise(
      ymin = min(value, na.rm = TRUE),
      ymax = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y = ymax + pmax((ymax - ymin) * 0.18, abs(ymax) * 0.03))

  boundary_ann <- left_join(boundary_ann, boundary_ann_y, by = "metric")

  pBoundaryMetrics <- ggplot(
    metrics_long,
    aes(x = Genotype, y = value, fill = Genotype, color = Genotype)
  ) +
    geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.7, alpha = 0.25) +
    geom_quasirandom(width = 0.16, size = 1.8, alpha = 0.85) +
    geom_text(
      data = boundary_ann,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.6,
      vjust = 0
    ) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    # Extra headroom prevents BH-FDR labels from being clipped in free-y facets.
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.30))) +
    coord_cartesian(clip = "off") +
    labs(title = "A/P-associated boundary metrics", x = NULL, y = NULL) +
    paper_theme +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 10, r = 8, b = 6, l = 8)
    )

  ggsave(
    file.path(OUTDIR, "AP_boundary_scalar_metrics.pdf"),
    pBoundaryMetrics,
    width = 180, height = 75, units = "mm", device = cairo_pdf
  )
}

# ----------------------------------------------------------------------------
# 15. FIGURE HELPERS
# ----------------------------------------------------------------------------
curve_df <- function(M, curve_list = curve_sets, prefix = "curve") {
  bind_rows(lapply(seq_along(curve_list), function(i) {
    id <- curve_list[[i]]
    data.frame(
      x = M[id, 1],
      y = plot_y(M[id, 2]),
      group = paste0(prefix, i)
    )
  }))
}

edge_df <- function(M, edge_list = fixed_edges, prefix = "edge") {
  bind_rows(lapply(seq_along(edge_list), function(i) {
    id <- edge_list[[i]]
    data.frame(
      x = M[id, 1],
      y = plot_y(M[id, 2]),
      group = paste0(prefix, i)
    )
  }))
}

path_df <- function(M, id, group_name = "path") {
  data.frame(
    x = M[id, 1],
    y = plot_y(M[id, 2]),
    group = group_name
  )
}

# ----------------------------------------------------------------------------
# 16. PANEL A: DATAMARTIN SCHEMATIC
# ----------------------------------------------------------------------------
qc_candidates <- which(dat_pair$Genotype == CONTROL & dat_pair$LR == "right")
if (length(qc_candidates) == 0) stop("No control right wing available for Panel A.")
qc_med <- median(dat_pair$whole[qc_candidates], na.rm = TRUE)
qc_i <- qc_candidates[which.min(abs(dat_pair$whole[qc_candidates] - qc_med))]
qc <- coords[, , qc_i]

qc_curves <- curve_df(qc)
qc_edges  <- edge_df(qc)
qc_ap     <- path_df(qc, boundary_idx, "AP")

region_label_nodes <- list(
  A0 = c(2, 3, 14, 15, 16),
  A1 = c(7, 8, 30, 31, 33, 34),
  A2 = c(18, 19, 20, 27, 28, 43, 44),
  A3 = c(47, 48, 49, 67, 68, 69, 70, 71, 72),
  A4 = c(74, 90, 91, 92, 94, 95, 96, 98, 99),
  A5 = c(39, 40, 41, 54, 55, 56, 57, 58),
  A6 = c(61, 62, 63, 79, 80, 81, 82, 86, 87),
  A7 = c(84, 85, 99, 104, 105, 106, 107, 108, 109, 110),
  A8 = c(76, 77, 101, 102, 103)
)

region_labels <- bind_rows(lapply(names(region_label_nodes), function(nm) {
  id <- region_label_nodes[[nm]]
  data.frame(
    region = nm,
    x = mean(qc[id, 1]),
    y = mean(plot_y(qc[id, 2]))
  )
}))

pA <- ggplot() +
  geom_path(data = qc_curves, aes(x, y, group = group), linewidth = 0.45) +
  geom_path(data = qc_edges,  aes(x, y, group = group), linewidth = 0.45) +
  geom_path(data = qc_ap, aes(x, y, group = group), linewidth = 1.1) +
  geom_point(
    data = data.frame(x = qc[fixed_idx, 1], y = plot_y(qc[fixed_idx, 2])),
    aes(x, y), shape = 21, fill = "white", size = 1.6, stroke = 0.55
  ) +
  geom_text(
    data = region_labels,
    aes(x, y, label = region),
    fontface = "bold", size = 4
  ) +
  coord_equal() +
  labs(
    title = "DataMartin regional scheme",
    subtitle = paste0(
      "Anterior: A0, A1, A2, A3, A5   |   Posterior: A4, A6, A7, A8\n",
      "Thick line: adult A/P-associated geometric boundary"
    ),
    x = NULL, y = NULL
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 7.2),
    plot.margin = margin(4, 4, 4, 4)
  )

# ----------------------------------------------------------------------------
# 17. PANEL B: POSTERIOR / WHOLE
# ----------------------------------------------------------------------------
yB <- range(fly_area$P_area_percent, na.rm = TRUE)
padB <- max(diff(yB) * 0.20, 0.2)
label_yB <- yB[2] + padB

pB <- ggplot(
  fly_area,
  aes(x = Genotype, y = P_area_percent, fill = Genotype, color = Genotype)
) +
  geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.7, alpha = 0.25) +
  geom_quasirandom(width = 0.16, size = 2.0, alpha = 0.85) +
  annotate(
    "text", x = 1.5, y = label_yB,
    label = format_p(test_Pfrac$p.value), size = 3
  ) +
  coord_cartesian(
    ylim = c(yB[1] - padB * 0.25, label_yB + padB * 0.35),
    clip = "off"
  ) +
  labs(
    title = "Posterior compartment proportion",
    x = NULL,
    y = "Posterior area / whole wing area (%)"
  ) +
  paper_theme +
  theme(legend.position = "none")

# ----------------------------------------------------------------------------
# 18. PANEL C: A3 AND A4 ABSOLUTE AREAS
# ----------------------------------------------------------------------------
a34_long <- fly_area %>%
  select(fly_id, Genotype, A3_area, A4_area) %>%
  pivot_longer(
    cols = c(A3_area, A4_area),
    names_to = "Region",
    values_to = "Area"
  ) %>%
  mutate(
    Region = recode(
      Region,
      A3_area = "A3 (anterior)",
      A4_area = "A4 (posterior)"
    ),
    Region = factor(Region, levels = c("A3 (anterior)", "A4 (posterior)"))
  )

a34_ann <- region_stats %>%
  filter(region %in% c("A3", "A4")) %>%
  mutate(
    Region = recode(region,
                    A3 = "A3 (anterior)",
                    A4 = "A4 (posterior)"),
    Region = factor(Region, levels = c("A3 (anterior)", "A4 (posterior)")),
    label = paste0("BH-FDR ", sapply(FDR, format_p))
  ) %>%
  select(Region, label)

a34_ann_y <- a34_long %>%
  group_by(Region) %>%
  summarise(
    y_max = max(Area, na.rm = TRUE),
    y_min = min(Area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y = y_max + pmax((y_max - y_min) * 0.18, y_max * 0.03))

a34_ann <- left_join(a34_ann, a34_ann_y, by = "Region")

pC <- ggplot(
  a34_long,
  aes(x = Genotype, y = Area, fill = Genotype, color = Genotype)
) +
  geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.7, alpha = 0.25) +
  geom_quasirandom(width = 0.16, size = 2.0, alpha = 0.85) +
  geom_text(
    data = a34_ann,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.7
  ) +
  facet_wrap(~ Region, scales = "free_y", nrow = 1) +
  labs(
    title = "Regional localization",
    x = NULL,
    y = "Regional wing area"
  ) +
  paper_theme +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold")
  )

# ----------------------------------------------------------------------------
# 19. PANEL D: SIZE-ADJUSTED WHOLE-WING DEFORMATION
# ----------------------------------------------------------------------------
ctrl_curves_D <- curve_df(pred_ctrl, prefix = "ctrl")
mut_curves_D  <- curve_df(pred_mut_mag, prefix = "mut")
ctrl_edges_D  <- edge_df(pred_ctrl, prefix = "ctrl_edge")
mut_edges_D   <- edge_df(pred_mut_mag, prefix = "mut_edge")
ctrl_ap_D     <- path_df(pred_ctrl, boundary_idx, "ctrl_AP")
mut_ap_D      <- path_df(pred_mut_mag, boundary_idx, "mut_AP")

anchor_vectors <- data.frame(
  point = fixed_idx,
  x = pred_ctrl[fixed_idx, 1],
  y = plot_y(pred_ctrl[fixed_idx, 2]),
  xend = pred_mut_mag[fixed_idx, 1],
  yend = plot_y(pred_mut_mag[fixed_idx, 2])
) %>%
  mutate(
    displacement = sqrt((xend - x)^2 + (yend - y)^2),
    displacement_rank = rank(-displacement, ties.method = "first")
  )

n_major <- min(N_MAJOR_VECTORS, nrow(anchor_vectors))

major_anchor_vectors <- anchor_vectors %>%
  arrange(desc(displacement)) %>%
  slice_head(n = n_major)

write.csv(
  anchor_vectors %>% arrange(displacement_rank),
  file.path(OUTDIR, "Panel_D_fixed_anchor_displacement_ranking.csv"),
  row.names = FALSE
)

shape_R2 <- shape_genotype$Rsq[1]
shape_F  <- shape_genotype$F[1]
shape_P  <- shape_genotype[["Pr(>F)"]][1]
size_R2  <- shape_size$Rsq[1]
size_F   <- shape_size$F[1]
size_P   <- shape_size[["Pr(>F)"]][1]
int_R2   <- shape_interact$Rsq[1]
int_F    <- shape_interact$F[1]
int_P    <- shape_interact[["Pr(>F)"]][1]

shape_label <- paste0(
  "Procrustes ANOVA: genotype R² = ", formatC(shape_R2, digits = 3, format = "f"),
  ", F = ", formatC(shape_F, digits = 2, format = "f"),
  ", ", format_p(shape_P)
)
shape_label_extra <- paste0(
  "logCS: R² = ", formatC(size_R2, digits = 3, format = "f"),
  ", F = ", formatC(size_F, digits = 2, format = "f"),
  ", ", format_p(size_P),
  "; interaction: R² = ", formatC(int_R2, digits = 3, format = "f"),
  ", F = ", formatC(int_F, digits = 2, format = "f"),
  ", ", format_p(int_P)
)

pD <- ggplot() +
  # Control/reference shape
  geom_path(
    data = ctrl_curves_D, aes(x, y, group = group),
    linewidth = 0.70, color = PANEL_D_CTRL_COLOR
  ) +
  geom_path(
    data = ctrl_edges_D, aes(x, y, group = group),
    linewidth = 0.70, color = PANEL_D_CTRL_COLOR
  ) +
  geom_path(
    data = ctrl_ap_D, aes(x, y, group = group),
    linewidth = 1.10, color = PANEL_D_CTRL_COLOR
  ) +

  # Deformed mutant prediction (magnified for visualization)
  geom_path(
    data = mut_curves_D, aes(x, y, group = group),
    linewidth = 0.95, color = PANEL_D_MUT_COLOR,
    linetype = 2, alpha = 0.95
  ) +
  geom_path(
    data = mut_edges_D, aes(x, y, group = group),
    linewidth = 0.95, color = PANEL_D_MUT_COLOR,
    linetype = 2, alpha = 0.95
  ) +
  geom_path(
    data = mut_ap_D, aes(x, y, group = group),
    linewidth = 1.15, color = PANEL_D_MUT_COLOR,
    linetype = 2, alpha = 0.95
  ) +

  # Displacement vectors at fixed anchors
  geom_segment(
    data = anchor_vectors,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.65,
    color = PANEL_D_VECTOR_COLOR,
    alpha = 0.95,
    arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")
  ) +

  # Control anchor positions
  geom_point(
    data = anchor_vectors,
    aes(x, y), shape = 21, fill = "white",
    color = PANEL_D_CTRL_COLOR, size = 2.0, stroke = 0.65
  ) +

  # Deformed anchor positions
  geom_point(
    data = anchor_vectors,
    aes(xend, yend), shape = 16,
    color = PANEL_D_MUT_COLOR, size = 1.5, alpha = 0.95
  ) +
  coord_equal() +
  labs(
    title = "Size-adjusted wing-shape deformation",
    subtitle = paste0(
      CONTROL, " (gray) → ", MUTANT,
      " (orange dashed); blue arrows = displacement; deformation ×", DEFORMATION_MAG,
      "
", shape_label,
      "
", shape_label_extra
    ),
    x = NULL, y = NULL
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 6.6),
    plot.margin = margin(5, 5, 5, 5)
  )

# ----------------------------------------------------------------------------
# 19B. ALTERNATIVE PANEL D VERSIONS
# ----------------------------------------------------------------------------
# Version 1: show only the largest fixed-anchor displacement vectors.
# This selection is for visualization only; it does not alter the statistical test.
pD_major <- ggplot() +
  geom_path(data = ctrl_curves_D, aes(x, y, group = group),
            linewidth = 0.70, color = PANEL_D_CTRL_COLOR) +
  geom_path(data = ctrl_edges_D, aes(x, y, group = group),
            linewidth = 0.70, color = PANEL_D_CTRL_COLOR) +
  geom_path(data = ctrl_ap_D, aes(x, y, group = group),
            linewidth = 1.10, color = PANEL_D_CTRL_COLOR) +
  geom_path(data = mut_curves_D, aes(x, y, group = group),
            linewidth = 0.95, color = PANEL_D_MUT_COLOR,
            linetype = 2, alpha = 0.95) +
  geom_path(data = mut_edges_D, aes(x, y, group = group),
            linewidth = 0.95, color = PANEL_D_MUT_COLOR,
            linetype = 2, alpha = 0.95) +
  geom_path(data = mut_ap_D, aes(x, y, group = group),
            linewidth = 1.15, color = PANEL_D_MUT_COLOR,
            linetype = 2, alpha = 0.95) +
  geom_segment(
    data = major_anchor_vectors,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.75, color = PANEL_D_VECTOR_COLOR, alpha = 0.98,
    arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")
  ) +
  geom_point(
    data = anchor_vectors,
    aes(x, y), shape = 21, fill = "white",
    color = PANEL_D_CTRL_COLOR, size = 1.8, stroke = 0.55
  ) +
  geom_point(
    data = major_anchor_vectors,
    aes(xend, yend), shape = 16,
    color = PANEL_D_MUT_COLOR, size = 1.7, alpha = 0.98
  ) +
  geom_text(
    data = major_anchor_vectors,
    aes(xend, yend, label = paste0("P", point)),
    color = PANEL_D_VECTOR_COLOR, size = 2.0, vjust = -0.7
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "Size-adjusted wing-shape deformation",
    subtitle = paste0(
      "Top ", nrow(major_anchor_vectors),
      " fixed-anchor displacement vectors; deformation ×", DEFORMATION_MAG,
      "\n", shape_label,
      "\n", shape_label_extra
    ),
    x = NULL, y = NULL
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 6.4),
    plot.margin = margin(7, 8, 6, 6)
  )

# Version 2a: shape overlay only, with no vectors.
pD_overlay <- ggplot() +
  geom_path(data = ctrl_curves_D, aes(x, y, group = group),
            linewidth = 0.80, color = PANEL_D_CTRL_COLOR) +
  geom_path(data = ctrl_edges_D, aes(x, y, group = group),
            linewidth = 0.80, color = PANEL_D_CTRL_COLOR) +
  geom_path(data = ctrl_ap_D, aes(x, y, group = group),
            linewidth = 1.15, color = PANEL_D_CTRL_COLOR) +
  geom_path(data = mut_curves_D, aes(x, y, group = group),
            linewidth = 1.05, color = PANEL_D_MUT_COLOR,
            linetype = 2, alpha = 0.98) +
  geom_path(data = mut_edges_D, aes(x, y, group = group),
            linewidth = 1.05, color = PANEL_D_MUT_COLOR,
            linetype = 2, alpha = 0.98) +
  geom_path(data = mut_ap_D, aes(x, y, group = group),
            linewidth = 1.25, color = PANEL_D_MUT_COLOR,
            linetype = 2, alpha = 0.98) +
  geom_point(data = anchor_vectors,
             aes(x, y), shape = 21, fill = "white",
             color = PANEL_D_CTRL_COLOR, size = 1.7, stroke = 0.55) +
  geom_point(data = anchor_vectors,
             aes(xend, yend), shape = 16,
             color = PANEL_D_MUT_COLOR, size = 1.25, alpha = 0.95) +
  coord_equal() +
  labs(
    title = "Shape overlay",
    subtitle = paste0(CONTROL, " gray; ", MUTANT,
                      " orange dashed; deformation ×", DEFORMATION_MAG),
    x = NULL, y = NULL
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.5, face = "bold"),
    plot.subtitle = element_text(size = 6.8),
    plot.margin = margin(5, 5, 5, 5)
  )

# Version 2b: major displacement vectors, with a faint reference skeleton.
pD_vectors_only <- ggplot() +
  geom_path(data = ctrl_curves_D, aes(x, y, group = group),
            linewidth = 0.45, color = PANEL_D_CTRL_COLOR, alpha = 0.25) +
  geom_path(data = ctrl_edges_D, aes(x, y, group = group),
            linewidth = 0.45, color = PANEL_D_CTRL_COLOR, alpha = 0.25) +
  geom_path(data = ctrl_ap_D, aes(x, y, group = group),
            linewidth = 0.65, color = PANEL_D_CTRL_COLOR, alpha = 0.30) +
  geom_segment(
    data = major_anchor_vectors,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.85, color = PANEL_D_VECTOR_COLOR, alpha = 1,
    arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed")
  ) +
  geom_point(data = major_anchor_vectors,
             aes(x, y), shape = 21, fill = "white",
             color = PANEL_D_VECTOR_COLOR, size = 2.1, stroke = 0.7) +
  geom_point(data = major_anchor_vectors,
             aes(xend, yend), shape = 16,
             color = PANEL_D_MUT_COLOR, size = 1.8) +
  geom_text(data = major_anchor_vectors,
            aes(xend, yend, label = paste0("P", point)),
            color = PANEL_D_VECTOR_COLOR, size = 2.0, vjust = -0.7) +
  coord_equal(clip = "off") +
  labs(
    title = "Major displacement vectors",
    subtitle = paste0("Top ", nrow(major_anchor_vectors),
                      " fixed anchors; deformation ×", DEFORMATION_MAG),
    x = NULL, y = NULL
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.5, face = "bold"),
    plot.subtitle = element_text(size = 6.8),
    plot.margin = margin(5, 8, 5, 5)
  )

pD_two_panel <- (pD_overlay + pD_vectors_only) +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    title = "Size-adjusted wing-shape deformation",
    subtitle = paste0(shape_label, "
", shape_label_extra)
  ) &
  theme(
    plot.title = element_text(size = 9, face = "bold", family = "Arial"),
    plot.subtitle = element_text(size = 6.4, family = "Arial")
  )

# ----------------------------------------------------------------------------
# 20. SAVE MAIN FIGURE PANELS AND COMBINED FIGURE
# ----------------------------------------------------------------------------
ggsave(
  file.path(OUTDIR, "Panel_A_DataMartin_A0_A8_schematic.pdf"),
  pA, width = 100, height = 65, units = "mm", device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Panel_B_posterior_fraction.pdf"),
  pB, width = 75, height = 70, units = "mm", device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Panel_C_A3_A4_absolute_area.pdf"),
  pC, width = 105, height = 70, units = "mm", device = cairo_pdf
)

# Keep the all-vector version for backward compatibility.
ggsave(
  file.path(OUTDIR, "Panel_D_size_adjusted_deformation.pdf"),
  pD, width = 105, height = 72, units = "mm", device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Panel_D_size_adjusted_deformation_all_vectors.pdf"),
  pD, width = 105, height = 72, units = "mm", device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Panel_D_size_adjusted_deformation_major_vectors.pdf"),
  pD_major, width = 105, height = 72, units = "mm", device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Panel_D_size_adjusted_deformation_two_panel.pdf"),
  pD_two_panel, width = 180, height = 78, units = "mm", device = cairo_pdf
)

# Main Figure with all-vector D.
fig_ABCD <- (
  pA + pB + plot_layout(widths = c(1.25, 0.75))
) / (
  pC + pD + plot_layout(widths = c(1.0, 1.15))
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))

ggsave(
  file.path(OUTDIR, "Figure_ABCD_main.pdf"),
  fig_ABCD,
  width = FIG_WIDTH_MM,
  height = FIG_HEIGHT_MM,
  units = "mm",
  device = cairo_pdf
)

# Alternate Main Figure with simplified major-vector D.
fig_ABCD_major <- (
  pA + pB + plot_layout(widths = c(1.25, 0.75))
) / (
  pC + pD_major + plot_layout(widths = c(1.0, 1.15))
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))

ggsave(
  file.path(OUTDIR, "Figure_ABCD_main_D_major_vectors.pdf"),
  fig_ABCD_major,
  width = FIG_WIDTH_MM,
  height = FIG_HEIGHT_MM,
  units = "mm",
  device = cairo_pdf
)

# Alternate Main Figure with D split into shape overlay + vector map.
fig_ABCD_two_panelD <- (
  pA + pB + plot_layout(widths = c(1.25, 0.75))
) / (
  pC + pD_two_panel + plot_layout(widths = c(0.85, 1.55))
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Arial"))

ggsave(
  file.path(OUTDIR, "Figure_ABCD_main_D_two_panel.pdf"),
  fig_ABCD_two_panelD,
  width = 205,
  height = 150,
  units = "mm",
  device = cairo_pdf
)

# ----------------------------------------------------------------------------
# 21. OPTIONAL SUPPLEMENTARY REGIONAL CHANGE PLOT
# ----------------------------------------------------------------------------
region_stats_plot <- region_stats %>%
  mutate(
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  )

pRegionChange <- ggplot(
  region_stats_plot,
  aes(x = region, y = percent_change)
) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_point(size = 2.2) +
  geom_text(aes(label = significance), vjust = -0.8, size = 3) +
  # Add vertical space so significance labels above positive regions are not clipped.
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.28))) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Regional area change (%)"
  ) +
  paper_theme +
  theme(
    plot.margin = margin(t = 10, r = 8, b = 6, l = 8)
  )

ggsave(
  file.path(OUTDIR, "Supplementary_A0_A8_percent_change.pdf"),
  pRegionChange,
  width = 105, height = 78, units = "mm", device = cairo_pdf
)

# ----------------------------------------------------------------------------
# 22. COMPACT STATISTICAL SUMMARY
# ----------------------------------------------------------------------------
summary_rows <- list()

summary_rows[[length(summary_rows) + 1]] <- data.frame(
  analysis = "Posterior / whole-wing area",
  statistic = unname(test_Pfrac$statistic),
  R2 = NA_real_,
  p_value = test_Pfrac$p.value,
  FDR = NA_real_
)

summary_rows[[length(summary_rows) + 1]] <- data.frame(
  analysis = "A3 absolute area",
  statistic = region_stats$W[region_stats$region == "A3"],
  R2 = NA_real_,
  p_value = region_stats$p_value[region_stats$region == "A3"],
  FDR = region_stats$FDR[region_stats$region == "A3"]
)

summary_rows[[length(summary_rows) + 1]] <- data.frame(
  analysis = "A4 absolute area",
  statistic = region_stats$W[region_stats$region == "A4"],
  R2 = NA_real_,
  p_value = region_stats$p_value[region_stats$region == "A4"],
  FDR = region_stats$FDR[region_stats$region == "A4"]
)

summary_rows[[length(summary_rows) + 1]] <- data.frame(
  analysis = "Shape: sliding semilandmarks, genotype after size",
  statistic = shape_genotype$F[1],
  R2 = shape_genotype$Rsq[1],
  p_value = shape_genotype[["Pr(>F)"]][1],
  FDR = NA_real_
)

summary_rows[[length(summary_rows) + 1]] <- data.frame(
  analysis = "Shape: fixed-anchor sensitivity, genotype after size",
  statistic = fixed_genotype$F[1],
  R2 = fixed_genotype$Rsq[1],
  p_value = fixed_genotype[["Pr(>F)"]][1],
  FDR = NA_real_
)

summary_rows[[length(summary_rows) + 1]] <- data.frame(
  analysis = "Shape: size x genotype interaction",
  statistic = shape_interact$F[1],
  R2 = shape_interact$Rsq[1],
  p_value = shape_interact[["Pr(>F)"]][1],
  FDR = NA_real_
)

if (RUN_BOUNDARY_ANALYSIS) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    analysis = "Boundary profile: genotype after size",
    statistic = boundary_genotype$F[1],
    R2 = boundary_genotype$Rsq[1],
    p_value = boundary_genotype[["Pr(>F)"]][1],
    FDR = NA_real_
  )

  for (i in seq_len(nrow(boundary_metric_stats))) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      analysis = paste0("Boundary scalar: ", boundary_metric_stats$metric[i]),
      statistic = boundary_metric_stats$W[i],
      R2 = NA_real_,
      p_value = boundary_metric_stats$p_value[i],
      FDR = boundary_metric_stats$FDR[i]
    )
  }
}

summary_stats <- bind_rows(summary_rows)
write.csv(summary_stats, file.path(OUTDIR, "statistical_summary.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 23. SAVE ANALYSIS METADATA / SESSION INFO
# ----------------------------------------------------------------------------
analysis_metadata <- data.frame(
  setting = c(
    "CONTROL",
    "MUTANT",
    "N_PERM",
    "N_BOUNDARY_RESAMPLE",
    "DEFORMATION_MAG",
    "FLIP_WING_VERTICAL"
  ),
  value = c(
    CONTROL,
    MUTANT,
    as.character(N_PERM),
    as.character(N_BOUNDARY_RESAMPLE),
    as.character(DEFORMATION_MAG),
    as.character(FLIP_WING_VERTICAL)
  )
)
write.csv(analysis_metadata, file.path(OUTDIR, "analysis_metadata.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

cat("\n============================================================\n")
cat("COMPLETE MANUSCRIPT PIPELINE FINISHED\n")
cat("Output directory: ", OUTDIR, "\n", sep = "")
cat("Main figure: ", file.path(OUTDIR, "Figure_ABCD_main.pdf"), "\n", sep = "")
cat("============================================================\n")

