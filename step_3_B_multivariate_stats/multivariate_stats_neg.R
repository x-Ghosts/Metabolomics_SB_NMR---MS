#remotes::install_github("HimesGroup/qmtools")
#install.packages(c('SummarizedExperiment', 'heatmaply', 'VIM', 'limma'))

library(tidyverse)
library(qmtools)
library(dplyr)
library(omicsTools)
library(MetabolAnalyze)
library(ggrepel)
library(pheatmap)
library(mixOmics)
library(gridGraphics)
library(gridExtra)


## ==> Data

dataset_neg_blood_selected <- read.csv("dataset_scaled/neg_Waters_scaled_blood.csv")

metadata_neg_blood <- dataset_neg_blood_selected[, 1:5]
signal_neg_metabolites <- dataset_neg_blood_selected[, 6:ncol(dataset_neg_blood_selected)]

sample_names <- colnames(signal_neg_metabolites)
sample_group <- case_when(
  grepl("qc", sample_names, ignore.case = TRUE) ~ "QC",
  grepl("zak", sample_names, ignore.case = TRUE) ~ "Patient",
  grepl("biobank", sample_names, ignore.case = TRUE) ~ "Biobank",
  TRUE ~ NA_character_
)
sample_metadata <- data.frame(
  Sample = sample_names,
  Group = sample_group
)

###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

## Principal Component Analysis - PCA :
pca_data <- signal_neg_metabolites
m_pca <- reduceFeatures(x = pca_data, i = "knn_vsn", method = "pca", ncomp = 10)
summary(m_pca)


m_pca_df <- as.data.frame(m_pca)
m_pca_df$sample <- rownames(m_pca)
m_pca_with_labels <- m_pca_df %>%
  mutate(
    Type_of_sample = case_when(
      grepl("QC", sample, ignore.case = TRUE)       ~ "QC samples",
      grepl("biobank", sample, ignore.case = TRUE)  ~ "Biobank",
      grepl("zak", sample, ignore.case = TRUE)       ~ "Zak - Infectious",
      TRUE ~ "other"
    )
  )

p_pca_plot_1 <-plotReduced(m_pca, group = m_pca_with_labels$Type_of_sample, point_size = 3, title = "PCA — Negative Ionization (Components 1 & 2)", legend = TRUE, comp = c(1,2))
p_pca_plot_1

###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

## ==> Remove the QC:

keep_samples <- sample_metadata$Group != "QC"
mat_stats <- as.matrix(signal_neg_metabolites[, keep_samples])
group <- sample_metadata$Group[keep_samples] # table(group)

###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

#PLSDA V2

rownames(mat_stats) <- metadata_neg_blood$Name

## 2. Extract sample names from columns
sample_names <- colnames(mat_stats)

## 3. Define groups from sample names
group <- case_when(
  grepl("zak", sample_names) ~ "Patient",
  grepl("biobank", sample_names) ~ "Biobank (Control)",
  TRUE ~ NA_character_
)

group <- factor(group)
table(group)

## 4. Keep only samples with valid labels
keep_samples <- !is.na(group)

mat_stats_plsda <- mat_stats[, keep_samples, drop = FALSE]
group_plsda <- droplevels(group[keep_samples])

## 5. Prepare X for mixOmics
## rows = samples, columns = metabolites
X <- t(mat_stats_plsda)
X <- as.matrix(X)
storage.mode(X) <- "double"


plsda_init <- plsda(X, group_plsda, ncomp = 10, scale = TRUE)

set.seed(123)

perf_res <- perf(
  plsda_init,
  validation = "Mfold",
  folds = 3,        # 
  nrepeat = 50,     # increase to 50 later if needed
  progressBar = TRUE
)

plot(perf_res)
perf_res$choice.ncomp

# Graph shows Comp 2

final.plsda.blood <- plsda(X,group_plsda, ncomp = 2)
var_exp <- final.plsda.blood$prop_expl_var$X*100

#################################################################
## PERMUTATION TEST FOR PLS-DA OVERFITTING
#################################################################

set.seed(123)

nperm <- 100      # use 100 first for testing, then 1000 for final
ncomp_final <- 2   # same number of components used in final model

# Observed model performance
obs_perf <- perf(
  final.plsda.blood,
  validation = "Mfold",
  folds = 3,
  nrepeat = 50,
  dist = "max.dist",
  progressBar = TRUE
)

# Observed BER at selected component
obs_BER <- obs_perf$error.rate$BER[paste0("comp", ncomp_final), "max.dist"]

perm_BER <- numeric(nperm)

for (i in seq_len(nperm)) {
  
  # Randomly shuffle class labels
  y_perm <- sample(group_plsda)
  
  # Fit PLS-DA with permuted labels
  perm_model <- plsda(
    X,
    y_perm,
    ncomp = ncomp_final,
    scale = TRUE
  )
  
  # Cross-validate permuted model
  perm_perf <- perf(
    perm_model,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,   # lower repeat for speed; increase if needed
    dist = "max.dist",
    progressBar = TRUE
  )
  
  perm_BER[i] <- perm_perf$error.rate$BER[paste0("comp", ncomp_final), "max.dist"]
}

# Empirical p-value
# Lower BER = better model
perm_pvalue <- mean(perm_BER <= obs_BER)

perm_pvalue
obs_BER
summary(perm_BER)

p <- plotIndiv(
  final.plsda.blood,
  ind.names = FALSE,
  legend = TRUE,
  legend.title = "Groups",
  comp = c(1,2),
  ellipse = TRUE,
  title = 'PLS-DA on blood (neg)',
  X.label = paste0('Comp 1 (', round(var_exp[1],1), '%)'),
  Y.label = paste0('Comp 2 (', round(var_exp[2],1), '%)')
)
p
ggsave("PLS-DA_blood_neg.png",p$graph, width = 9, height = 6, dpi = 300)

#################################################    Variable importance  in projection

vip_scores <- vip(final.plsda.blood)

vip_max <- apply(vip_scores[, 1:2, drop = FALSE], 1, max)

vip_df <- data.frame(
  Metabolite = names(vip_max),
  VIP = as.numeric(vip_max)
)

vip_df <- vip_df %>%
  mutate(Significance = VIP > 1)

p <- ggplot(vip_df, aes(x = reorder(Metabolite, VIP), y = VIP, fill = Significance)) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, color = "black") +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "grey70")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "VIP across Comp1 to Comp2",
    x = "Metabolite",
    y = "VIP Score"
  )

p
ggsave("VIP_Waters_neg.png", p, width = 9, height = 6, dpi = 300)


############################## ROC CURVE at 0.7 ########### Correlation & Loadings

sig_vars <- vip_df %>%
  filter(VIP > 1) %>%
  pull(Metabolite)

var_colors <- ifelse(
  colnames(X) %in% sig_vars,
  "red",
  "grey70"
)

png("CorrCircle_blood_neg.png", width = 9, height = 6, units = "in", res = 300)
plotVar(
  final.plsda.blood,
  comp = c(1,2),
  var.names = TRUE,
  col = list(X = var_colors),
  cex = 3
)
dev.off()




png("loadings_comp1_blood_neg.png", width = 7, height = 4, units = "in", res = 300)
plotLoadings(final.plsda.blood, comp = 1, method = "mean", contrib = "max", title = "Comp1 Contrib.")
dev.off()

png("loadings_comp2_blood_neg.png", width = 7, height = 6, units = "in", res = 300)
plotLoadings(final.plsda.blood, comp = 2, method = "mean", contrib = "max", title = "Comp2 Contrib.")
dev.off()


#################################################################
# Significant metabolites from VIP
sig_vars <- vip_df %>%
  filter(VIP > 1) %>%
  pull(Metabolite)

# Sample group colors
row_cols <- color.mixo(group_plsda)

# Metabolite significance colors
col_cols <- ifelse(colnames(X) %in% sig_vars, "red", "grey70")

png("heatmap_blood_VIP_significant.png", width = 2200, height = 2000, res = 200)

cim(
  X,
  row.sideColors = row_cols,   # sample groups
  col.sideColors = col_cols,   # significant metabolites
  row.names = FALSE,           # hide patient/sample names
  col.names = TRUE             # show metabolite names
)

dev.off()
#################################################################
png("heatmap_Blood_neg_2.png", width = 2000, height = 2000, res = 200)
row_cols <- color.mixo(group_plsda)

cim(
  final.plsda.blood,
  row.sideColors = row_cols
)

dev.off()
