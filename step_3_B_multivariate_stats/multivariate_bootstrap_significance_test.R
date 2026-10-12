# Library loading
library(tidyverse)
library(ggplot2)
library(metanorm)
library(MetabolAnalyze)
library(ggrepel)
library(corrplot)
library(Hmisc)
library(mixOmics)
library(gridGraphics)
library(gridExtra)
library(pheatmap)
library(factoextra)
library(ggrepel)
library(grid)
library(ggtext)
conflicts()

# Plasma (Blood) - Case Control

patient_set <- read.table(
  "patient/scaled_patient.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = " "
)

#############################################################################################################################################################
#############################################################################################################################################################
##           PCA
#############################################################################################################################################################
#############################################################################################################################################################

candidate_outliers <- c("52", "3")
# 
patient_set <- patient_set[!as.character(patient_set[, 1]) %in% candidate_outliers, ]

## 1. Extract sample IDs
sample_id <- patient_set[, 1]

## 2. Extract metabolite signals
signals <- patient_set[, -1]

## 3. Convert all signal columns to numeric
signals <- as.data.frame(lapply(signals, function(x) {
  x <- gsub(",", ".", x)
  as.numeric(x)
}))

rownames(signals) <- sample_id

signals <- as.matrix(signals)

storage.mode(signals) <- "double"

## 4. Prepare X
## rows = samples, columns = metabolites
X <- signals

# ============================================================
# 1. PCA
# ============================================================

pca <- prcomp(X, center = FALSE, scale. = FALSE)

# Variance explained
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

variance_table <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  Variance_explained = var_explained * 100,
  Cumulative_variance = cumsum(var_explained) * 100
)

variance_table[1:9, ]

# ============================================================
# 2. Explained variance plot
# ============================================================

variance_table_10 <- variance_table[1:9, ]

p_var <- ggplot(variance_table_10,
                aes(x = factor(PC, levels = PC))) +
  
  geom_col(
    aes(y = Variance_explained),
    fill = "#4DBBD5",
    width = 0.75,
    alpha = 0.85
  ) +
  
  geom_line(
    aes(
      y = Cumulative_variance,
      colour = "Cumulative variance",
      group = 1
    ),
    linewidth = 1
  ) +
  
  geom_point(
    aes(
      y = Cumulative_variance,
      colour = "Cumulative variance"
    ),
    size = 2.8
  ) +
  
  geom_text(
    aes(
      y = Variance_explained,
      label = paste0(round(Variance_explained, 1), "%")
    ),
    vjust = -1.6,
    size = 4.5
  ) +
  
  scale_colour_manual(
    values = c(
      "Cumulative variance" = "#E64B35"
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      max(variance_table_10$Cumulative_variance) * 1.08
    ),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    title = "PCA explained variance",
    subtitle = "First 9 principal components",
    x = "Principal component",
    y = "Variance Explained and Cumulative variance (%)",
    colour = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    
    plot.subtitle = element_text(
      size = 13
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    legend.position = "top",
    legend.title = element_blank()
  )

print(p_var)
ggsave("variances_explained.png", p_var,  width = 12,height = 8,dpi = 600)

# ============================================================
# 3. PCA scores
# ============================================================

scores <- as.data.frame(pca$x)
scores$SampleID <- rownames(X)

p <- ggplot(scores, aes(PC1, PC2, label = SampleID)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  theme_bw() +
  labs(
    title = "PCA score plot - PC1 vs PC2",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
  )
p
ggsave("PCA_1_2.png", p, width = 12, height = 8, dpi = 300)

p <- ggplot(scores, aes(PC1, PC3, label = SampleID)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  theme_bw() +
  labs(
    title = "PCA score plot - PC1 vs PC3",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC3 (", round(var_explained[3] * 100, 1), "%)")
  )
p
ggsave("PCA_1_3.png", p, width = 12, height = 8, dpi = 300)


p <- ggplot(scores, aes(PC1, PC4, label = SampleID)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  theme_bw() +
  labs(
    title = "PCA score plot - PC1 vs PC4",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC4 (", round(var_explained[4] * 100, 1), "%)")
  )
p
ggsave("PCA_1_4.png", p, width = 12, height = 8, dpi = 300)


p <- fviz_pca_biplot(
  pca,
  axes = c(1, 2),
  repel = TRUE,
  label = "all",
  col.var = "steelblue",
  col.ind = "gray30",
  geom.ind = c("point", "text"), 
) +
  ggtitle("PCA biplot - PC1 vs PC2")
p
ggsave("biplot_PCA_1_2.png", p, width = 18, height = 12, dpi = 300)

################################ TEST

# Individual PCA coordinates
scores_df <- as.data.frame(pca$x[, 1:2])
colnames(scores_df) <- c("PC1", "PC2")

scores_df$SampleID <- rownames(scores_df)

# Variable loadings
loadings_df <- as.data.frame(pca$rotation[, 1:2])
colnames(loadings_df) <- c("PC1", "PC2")

loadings_df$Variable <- rownames(loadings_df)

# ------------------------------------------------------------
# Rescale loadings so arrows fit visibly in the score-plot space
# ------------------------------------------------------------

score_range <- max(
  diff(range(scores_df$PC1, na.rm = TRUE)),
  diff(range(scores_df$PC2, na.rm = TRUE))
)

loading_range <- max(
  abs(loadings_df$PC1),
  abs(loadings_df$PC2)
)

arrow_scale <- 0.45 * score_range / loading_range

loadings_df$PC1_scaled <- loadings_df$PC1 * arrow_scale
loadings_df$PC2_scaled <- loadings_df$PC2 * arrow_scale


loadings_df$nudge_x <- ifelse(
  loadings_df$PC1_scaled >= 0,
  0.035,
  -0.035
)

loadings_df$nudge_y <- ifelse(
  loadings_df$PC2_scaled >= 0,
  0.035,
  -0.035
)

scores_df$nudge_y <- 0.05

scores_df$nudge_y[
  scores_df$SampleID %in% c("53", "94", "101", "141", "44", "100", "14", "113")
] <- -0.15

# ------------------------------------------------------------
# PCA variance explained
# ------------------------------------------------------------

var_explained <- pca$sdev^2 / sum(pca$sdev^2)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot() +
  
  # Reference lines
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    color = "grey75"
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    color = "grey75"
  ) +
  
  # Individuals
  geom_point(
    data = scores_df,
    aes(x = PC1, y = PC2),
    color = "darkred",
    alpha = 0.85,
    size = 3.5
  ) +
  
  # Individual labels
  geom_label_repel(
    data = scores_df,
    aes(x = PC1, y = PC2, label = SampleID),
    nudge_y = scores_df$nudge_y,
    color = "darkred",
    fill = "white",
    fontface = "bold",
    size = 5,
    label.size = 0.25,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    force = 1.2,
    force_pull = 3,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.35
  ) +
  
  # Loading arrows
  geom_segment(
    data = loadings_df,
    aes(
      x = 0,
      y = 0,
      xend = PC1_scaled,
      yend = PC2_scaled
    ),
    color = "steelblue",
    linewidth = 0.4,
    arrow = arrow(
      length = unit(0.25, "cm"),
      type = "closed"
    )
  ) +
  
  # Variable labels, positioned beyond arrow heads
  geom_label_repel(
    data = loadings_df,
    aes(
      x = PC1_scaled,
      y = PC2_scaled,
      label = Variable
    ),
    nudge_x = loadings_df$nudge_x,
    nudge_y = loadings_df$nudge_y,
    color = "steelblue",
    fill = "white",
    alpha = 0.85,
    fontface = "bold",
    size = 4.8,
    label.size = 0,
    box.padding = 0.25,
    point.padding = 0.15,
    force = 1.5,
    force_pull = 4,
    max.overlaps = Inf,
    max.iter = 50000,
    segment.color = "steelblue",
    segment.linewidth = 0.3
  ) +
  
  coord_equal() +
  
  labs(
    title = "PCA biplot: PC1 vs PC2",
    x = paste0(
      "PC1 (",
      round(var_explained[1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_explained[2] * 100, 1),
      "%)"
    )
  ) +
  
  theme_bw(base_size = 18) +
  
  theme(
    plot.title = element_text(
      size = 24,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 16
    ),
    
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(
      t = 15,
      r = 40,
      b = 15,
      l = 15
    )
  )

print(p)

ggsave(
  "biplot_PCA_1_2_manual.png",
  plot = p,
  width = 22,
  height = 14,
  dpi = 600
)

  ################################ END OF TEST

p <- fviz_pca_biplot(
  pca,
  axes = c(1, 3),
  repel = TRUE,
  label = "all",
  col.var = "steelblue",
  col.ind = "gray30",
  geom.ind = c("point", "text")
) +
  ggtitle("PCA biplot - PC1 vs PC3")
p
ggsave("biplot_PCA_1_3.png", p, width = 18, height = 12, dpi = 300)

p <- fviz_pca_biplot(
  pca,
  axes = c(1, 4),
  repel = TRUE,
  label = "all",
  col.var = "steelblue",
  col.ind = "gray30",
  geom.ind = c("point", "text")
) +
  ggtitle("PCA biplot - PC1 vs PC4")
p
ggsave("biplot_PCA_1_4.png", p, width = 18, height = 12, dpi = 300)

# ============================================================
# 5. Extract loadings and contributions
# ============================================================

# Loadings
loadings <- as.data.frame(pca$rotation)
loadings$Metabolite <- rownames(loadings)

# Contributions from factoextra
pca_var <- get_pca_var(pca)

contrib <- as.data.frame(pca_var$contrib)
contrib$Metabolite <- rownames(contrib)

# Combine loadings + contributions
pca_metabolite_table <- loadings %>%
  dplyr::select(Metabolite, dplyr::starts_with("PC")) %>%
  dplyr::rename_with(
    ~ paste0("Loading_", .x),
    dplyr::starts_with("PC")
  ) %>%
  dplyr::left_join(
    contrib %>%
      dplyr::select(Metabolite, dplyr::starts_with("Dim")) %>%
      dplyr::rename_with(
        ~ gsub("Dim", "Contribution_PC", .x),
        dplyr::starts_with("Dim")
      ),
    by = "Metabolite"
  )

# Fix possible Contribution_PC.1 naming issue
colnames(pca_metabolite_table) <- gsub(
  "Contribution_PC\\.",
  "Contribution_PC",
  colnames(pca_metabolite_table)
)

colnames(pca_metabolite_table)



# ============================================================
# 5. Extract loadings and contributions
# ============================================================

# Loadings
loadings <- as.data.frame(pca$rotation)
loadings$Metabolite <- rownames(loadings)

# Contributions from factoextra
pca_var <- get_pca_var(pca)

contrib <- as.data.frame(pca_var$contrib)
contrib$Metabolite <- rownames(contrib)

# Combine loadings + contributions
pca_metabolite_table <- loadings %>%
  dplyr::select(Metabolite, dplyr::starts_with("PC")) %>%
  dplyr::rename_with(
    ~ paste0("Loading_", .x),
    dplyr::starts_with("PC")
  ) %>%
  dplyr::left_join(
    contrib %>%
      dplyr::select(Metabolite, dplyr::starts_with("Dim")) %>%
      dplyr::rename_with(
        ~ gsub("Dim", "Contribution_PC", .x),
        dplyr::starts_with("Dim")
      ),
    by = "Metabolite"
  )

# Fix possible Contribution_PC.1 naming issue
colnames(pca_metabolite_table) <- gsub(
  "Contribution_PC\\.",
  "Contribution_PC",
  colnames(pca_metabolite_table)
)

colnames(pca_metabolite_table)


# ============================================================
# 6. Top contributors for PC1 to PC4
# ============================================================

top_PC1 <- pca_metabolite_table %>%
  dplyr::arrange(dplyr::desc(Contribution_PC1)) %>%
  dplyr::select(
    Metabolite,
    Loading_PC1,
    Contribution_PC1
  )

top_PC2 <- pca_metabolite_table %>%
  dplyr::arrange(dplyr::desc(Contribution_PC2)) %>%
  dplyr::select(
    Metabolite,
    Loading_PC2,
    Contribution_PC2
  )

top_PC3 <- pca_metabolite_table %>%
  dplyr::arrange(dplyr::desc(Contribution_PC3)) %>%
  dplyr::select(
    Metabolite,
    Loading_PC3,
    Contribution_PC3
  )

top_PC4 <- pca_metabolite_table %>%
  dplyr::arrange(dplyr::desc(Contribution_PC4)) %>%
  dplyr::select(
    Metabolite,
    Loading_PC4,
    Contribution_PC4
  )

top_PC1
top_PC2
top_PC3
top_PC4


# ============================================================
# 7. Average contribution threshold
# ============================================================

average_contribution <- 100 / ncol(X)
average_contribution



# ============================================================
# 8. Final PCA driver table: PC1 to PC4
# ============================================================

final_pca_drivers <- pca_metabolite_table %>%
  dplyr::mutate(
    Above_average_PC1 = Contribution_PC1 > average_contribution,
    Above_average_PC2 = Contribution_PC2 > average_contribution,
    Above_average_PC3 = Contribution_PC3 > average_contribution,
    Above_average_PC4 = Contribution_PC4 > average_contribution,
    
    Abs_Loading_PC1 = abs(Loading_PC1),
    Abs_Loading_PC2 = abs(Loading_PC2),
    Abs_Loading_PC3 = abs(Loading_PC3),
    Abs_Loading_PC4 = abs(Loading_PC4),
    
    Weighted_contribution_PC1_PC4 =
      Contribution_PC1 * var_explained[1] +
      Contribution_PC2 * var_explained[2] +
      Contribution_PC3 * var_explained[3] +
      Contribution_PC4 * var_explained[4]
  ) %>%
  dplyr::arrange(dplyr::desc(Weighted_contribution_PC1_PC4)) %>%
  dplyr::select(
    Metabolite,
    
    Loading_PC1,
    Contribution_PC1,
    Above_average_PC1,
    
    Loading_PC2,
    Contribution_PC2,
    Above_average_PC2,
    
    Loading_PC3,
    Contribution_PC3,
    Above_average_PC3,
    
    Loading_PC4,
    Contribution_PC4,
    Above_average_PC4,
    
    Weighted_contribution_PC1_PC4
  )

final_pca_drivers


# ============================================================
# 9. Contribution plots for PC1 to PC4
# ============================================================

top_n <- 20

final_pca_drivers %>%
  dplyr::slice_max(Contribution_PC1, n = top_n) %>%
  ggplot(aes(x = reorder(Metabolite, Contribution_PC1), y = Contribution_PC1)) +
  geom_col() +
  geom_hline(yintercept = average_contribution, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top metabolite contributions to PC1",
    x = "Metabolite",
    y = "Contribution to PC1 (%)"
  )
final_pca_drivers %>%
  dplyr::slice_max(Contribution_PC2, n = top_n) %>%
  ggplot(aes(x = reorder(Metabolite, Contribution_PC2), y = Contribution_PC2)) +
  geom_col() +
  geom_hline(yintercept = average_contribution, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top metabolite contributions to PC2",
    x = "Metabolite",
    y = "Contribution to PC2 (%)"
  )
final_pca_drivers %>%
  dplyr::slice_max(Contribution_PC3, n = top_n) %>%
  ggplot(aes(x = reorder(Metabolite, Contribution_PC3), y = Contribution_PC3)) +
  geom_col() +
  geom_hline(yintercept = average_contribution, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top metabolite contributions to PC3",
    x = "Metabolite",
    y = "Contribution to PC3 (%)"
  )
final_pca_drivers %>%
  dplyr::slice_max(Contribution_PC4, n = top_n) %>%
  ggplot(aes(x = reorder(Metabolite, Contribution_PC4), y = Contribution_PC4)) +
  geom_col() +
  geom_hline(yintercept = average_contribution, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top metabolite contributions to PC4",
    x = "Metabolite",
    y = "Contribution to PC4 (%)"
  )

# ============================================================
# 10. Weighted contribution plot across PC1 to PC4
# ============================================================

final_pca_drivers %>%
  dplyr::slice_max(Weighted_contribution_PC1_PC4, n = top_n) %>%
  ggplot(aes(
    x = reorder(Metabolite, Weighted_contribution_PC1_PC4),
    y = Weighted_contribution_PC1_PC4
  )) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top PCA metabolite drivers across PC1-PC4",
    subtitle = "Weighted by variance explained by each component",
    x = "Metabolite",
    y = "Weighted contribution"
  )


# ============================================================
# 11. Bootstrap significance / stability for weighted PCA drivers
# ============================================================

set.seed(123)

B <- 1000        # number of bootstrap repetitions
ncomps <- 4     # PC1 to PC4

# Expected weighted contribution threshold
# Since each PC has average contribution = 100 / number of metabolites
# and your weighted score uses var_explained as weights
weighted_average_threshold <- average_contribution * sum(var_explained[1:ncomps])

weighted_average_threshold

# Storage matrix
boot_weighted <- matrix(
  NA,
  nrow = B,
  ncol = ncol(X)
)

colnames(boot_weighted) <- colnames(X)


# ============================================================
# Bootstrap loop - safer version
# ============================================================

for (b in 1:B) {
  
  # Resample samples with replacement
  boot_rows <- sample(seq_len(nrow(X)), replace = TRUE)
  X_boot <- X[boot_rows, , drop = FALSE]
  
  # PCA on bootstrap sample
  pca_boot <- prcomp(X_boot, center = FALSE, scale. = FALSE)
  
  # Variance explained in bootstrap PCA
  var_boot <- (pca_boot$sdev^2) / sum(pca_boot$sdev^2)
  
  # Variable contributions
  boot_var <- factoextra::get_pca_var(pca_boot)
  boot_contrib <- as.data.frame(boot_var$contrib)
  
  # Keep only first ncomps columns by position, not by name
  boot_contrib_numeric <- boot_contrib[, 1:ncomps, drop = FALSE]
  
  # Force numeric matrix
  boot_contrib_numeric <- as.matrix(boot_contrib_numeric)
  storage.mode(boot_contrib_numeric) <- "double"
  
  # Match metabolites correctly
  rownames(boot_contrib_numeric) <- rownames(boot_contrib)
  boot_contrib_numeric <- boot_contrib_numeric[colnames(X), , drop = FALSE]
  
  # Weighted contribution across PC1-PC4
  boot_weighted[b, ] <- as.vector(
    boot_contrib_numeric %*% var_boot[1:ncomps]
  )
}

# ============================================================
# 12. Bootstrap summary table
# ============================================================

bootstrap_summary <- data.frame(
  Metabolite = colnames(X),
  
  Boot_mean_weighted = apply(boot_weighted, 2, mean, na.rm = TRUE),
  Boot_median_weighted = apply(boot_weighted, 2, median, na.rm = TRUE),
  
  Boot_CI_low = apply(boot_weighted, 2, quantile, probs = 0.025, na.rm = TRUE),
  Boot_CI_high = apply(boot_weighted, 2, quantile, probs = 0.975, na.rm = TRUE),
  
  # Frequency of being above the weighted average contribution threshold
  Bootstrap_frequency_above_threshold =
    apply(boot_weighted, 2, function(z) mean(z > weighted_average_threshold, na.rm = TRUE))
)
# Bootstrap support p-value:
# low value means the metabolite is consistently above the threshold
bootstrap_summary$Bootstrap_p <- apply(
  boot_weighted,
  2,
  function(z) {
    (1 + sum(z <= weighted_average_threshold, na.rm = TRUE)) / (length(na.omit(z)) + 1)
  }
)

bootstrap_summary$Bootstrap_FDR <- p.adjust(
  bootstrap_summary$Bootstrap_p,
  method = "BH"
)
# ============================================================
# 13. Final PCA driver table with bootstrap significance
# ============================================================

final_pca_drivers_bootstrap <- final_pca_drivers %>%
  dplyr::left_join(bootstrap_summary, by = "Metabolite") %>%
  dplyr::mutate(
    
    Stable_bootstrap_075 =
      Weighted_contribution_PC1_PC4 > weighted_average_threshold &
      Bootstrap_frequency_above_threshold >= 0.75,
    
    Stable_bootstrap_080 =
      Weighted_contribution_PC1_PC4 > weighted_average_threshold &
      Bootstrap_frequency_above_threshold >= 0.80,
    
    Stable_bootstrap_090 =
      Weighted_contribution_PC1_PC4 > weighted_average_threshold &
      Bootstrap_frequency_above_threshold >= 0.90,
    
    Bootstrap_priority = dplyr::case_when(
      Stable_bootstrap_090 == TRUE ~ "Very stable PCA driver",
      Stable_bootstrap_080 == TRUE ~ "Stable PCA driver",
      Stable_bootstrap_075 == TRUE ~ "Moderately stable PCA driver",
      Weighted_contribution_PC1_PC4 > weighted_average_threshold ~ "Above-average PCA driver",
      TRUE ~ "Low PCA contribution"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(Stable_bootstrap_080),
    dplyr::desc(Weighted_contribution_PC1_PC4)
  )

final_pca_drivers_bootstrap


# ============================================================
# 14. Plot weighted PCA drivers with bootstrap significance
# ============================================================

top_n <- 20

final_pca_drivers_bootstrap %>%
  dplyr::slice_max(Weighted_contribution_PC1_PC4, n = top_n) %>%
  ggplot(aes(
    x = reorder(Metabolite, Weighted_contribution_PC1_PC4),
    y = Weighted_contribution_PC1_PC4,
    fill = Stable_bootstrap_075
  )) +
  geom_col() +
  geom_hline(
    yintercept = weighted_average_threshold,
    linetype = "dashed",
    color = "red"
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Bootstrap-stable PCA metabolite drivers across PC1-PC4",
    subtitle = paste0("Bootstrap repetitions = ", B, "; stability threshold = 75%"),
    x = "Metabolite",
    y = "Weighted contribution across PC1-PC4",
    fill = "Bootstrap stable"
  )
p <- final_pca_drivers_bootstrap %>%
  dplyr::slice_max(Bootstrap_frequency_above_threshold, n = top_n) %>%
  ggplot(aes(
    x = reorder(Metabolite, Bootstrap_frequency_above_threshold),
    y = Bootstrap_frequency_above_threshold,
    fill = Stable_bootstrap_075
  )) +
  geom_col() +
  geom_hline(
    yintercept = 0.75,
    linetype = "dashed",
    color = "red"
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Bootstrap stability of PCA metabolite drivers",
    subtitle = "Frequency of weighted contribution above expected threshold",
    x = "Metabolite",
    y = "Bootstrap frequency above threshold",
    fill = "Bootstrap stable"
  )
p

###################################################### TEST

p <- final_pca_drivers_bootstrap %>%
  slice_max(
    Bootstrap_frequency_above_threshold,
    n = top_n
  ) %>%
  mutate(
    Metabolite_label = ifelse(
      Stable_bootstrap_075,
      paste0("**", Metabolite, "**"),
      Metabolite
    ),
    Metabolite_label = reorder(
      Metabolite_label,
      Bootstrap_frequency_above_threshold
    )
  ) %>%
  ggplot(
    aes(
      x = Metabolite_label,
      y = Bootstrap_frequency_above_threshold,
      fill = Stable_bootstrap_075
    )
  ) +
  
  geom_col(
    width = 0.8
  ) +
  
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Bootstrap_frequency_above_threshold
      )
    ),
    hjust = -0.1,
    size = 4
  ) +
  
  geom_vline(
    xintercept = NULL
  ) +
  
  geom_hline(
    yintercept = 0.75,
    linetype = "dashed",
    linewidth = 1,
    colour = "firebrick"
  ) +
  
  annotate(
    "text",
    x = 2,
    y = 0.78,
    label = " ",
    colour = "firebrick",
    fontface = "bold",
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "FALSE" = "grey80",
      "TRUE"  = "#2C7FB8"
    ),
    labels = c(
      "False",
      "True"
    )
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(
      0,
      1.05
    ),
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  
  labs(
    title =
      "Bootstrap stability of PCA CSF metabolites",
    
    
    x = NULL,
    
    y =
      "Bootstrap frequency",
    
    fill =
      "Stable\n(>= 0.75)"
  ) +
  
  theme_bw(
    base_size = 16
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 18
      ),
    
    plot.subtitle =
      element_text(
        size = 13
      ),
    
    axis.text.y = ggtext::element_markdown(),
    
    panel.grid.major.y =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
    
  )

p

###################################################### END OF TEST
ggsave("boostrap_pca_stability.png", p, width = 24, height = 16, dpi = 300)
