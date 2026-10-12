setwd("DIRNAME/step_1_dataset_filtering_imputations/")

library(tidyverse)
library(qmtools)
library(dplyr)
library(missForest)
library(metabolomicsR)
library(impute)
library(imputeLCMD)


# Loading Sets
dataset_neg_blood_selected <- read.csv("blood/neg Waters wybrane.csv")
dataset_neg_blood_selected <- dataset_neg_blood_selected[-1,]


metadata_neg_blood <- dataset_neg_blood_selected[,1:5]
signal_neg_metabolites <- dataset_neg_blood_selected[,6:ncol(dataset_neg_blood_selected)]  

# We leave the dataset intact, we copy the matrix, and make it in NA for the values that are == 0
signal_neg_metabolites_num <- as.data.frame(
  lapply(signal_neg_metabolites, function(x) {
    if (is.character(x)) type.convert(x, as.is = TRUE) else x
  })
)
signal_neg_metabolites_num[signal_neg_metabolites_num == 0] <- NA

# Fix the problem of Row numbers misaligned
rownames(metadata_neg_blood) <- NULL
rownames(signal_neg_metabolites_num) <- NULL

# Assign metabolite names from metadata to the signal matrix
rownames(signal_neg_metabolites_num) <- metadata_neg_blood$Name

## NO NORMAL DISTRIBUTION IN ALL SAMPLES

# dataset_neg_long <- signal_neg_metabolites_num %>%
#   as.data.frame() %>%
#   mutate(Name = rownames(.)) %>%
#   pivot_longer(
#     cols = -Name,
#     names_to = "Sample",
#     values_to = "signal"
#   )
# dataset_neg_long <- dataset_neg_long %>%
#   mutate(
#     Type_of_sample = case_when(
#       grepl("qc", Sample, ignore.case = TRUE) ~ "qc",
#       grepl("zak", Sample, ignore.case = TRUE) ~ "zak",
#       TRUE ~ "biobank"
#     )
#   )
# 
# p <- ggplot(dataset_neg_long, aes(x = signal, fill = Type_of_sample)) +
#   geom_histogram(position = "identity", bins = 30, alpha = 0.6) +
#   facet_wrap(~ Name, scales = "free", ncol = 4) +
#   theme_minimal() +
#   labs(
#     title = "Raw Signal Distribution",
#     x = "Signal Intensity",
#     y = "Count"
#   )

######################################################################################################################
######################################################################################################################
######################################################################################################################


# Compute Missingness
missingness <- data.frame(
  Metabolite = rownames(signal_neg_metabolites_num),
  Missing_Count = rowSums(is.na(signal_neg_metabolites_num)),
  Total_Samples = ncol(signal_neg_metabolites_num)
) %>%
  mutate(
    Missing_Percent = (Missing_Count / Total_Samples) * 100,
    Present_Percent = 100 - Missing_Percent
  )
keep_metabolites <- missingness %>%
  filter(Present_Percent >= 75) %>% # Keep 75 % of data and up.
  pull(Metabolite)

signal_neg_metabolites_filtered <- signal_neg_metabolites_num[rownames(signal_neg_metabolites_num) %in% keep_metabolites,
                                                              ,
                                                              drop = FALSE
]

p_missingness <- ggplot(
  missingness,
  aes(x = reorder(Metabolite, Missing_Percent), y = Missing_Percent, fill = Missing_Percent)
) +
  geom_col(color = "black") + geom_hline(yintercept = 25, linetype = "dashed", linewidth = 2) +
  scale_fill_gradient(
    low = "#1f77b4",
    high = "#d62728",
    name = "% Missing"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Missingness per Blood Metabolite (Negative Ionization)",
    x = "Metabolite",
    y = "Missingness (%)"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

p_missingness

ggsave("missingness_blood_neg.png", p_missingness, width = 16, height = 12, dpi = 300)

######################################################################################################################
######################################################################################################################
######################################################################################################################

# Base raw matrix
mat <- as.matrix(signal_neg_metabolites_filtered)
storage.mode(mat) <- "double"
mat[mat == 0] <- NA

# Log2-transformed matrix
mat_log <- log2(mat)

# Half-min equivalent on LOG10 scale
impute_halfmin_log2 <- function(xmat) {
  xmat <- as.matrix(xmat)
  storage.mode(xmat) <- "double"
  
  t(apply(xmat, 1, function(x) {
    min_log <- min(x[!is.na(x)], na.rm = TRUE)
    if (!is.finite(min_log)) min_log <- log2(1e-6)
    x[is.na(x)] <- min_log + log(0.5)
    x
  }))
}

# NRMSE
compute_nrmse <- function(true, pred) {
  true <- as.numeric(true)
  pred <- as.numeric(pred)
  sqrt(mean((true - pred)^2, na.rm = TRUE) / var(true, na.rm = TRUE))
}

set.seed(123)

# Mask observed values on log scale
obs_idx <- which(!is.na(mat_log))
mask_idx <- sample(obs_idx, size = floor(0.2 * length(obs_idx)))

mat_log_masked <- mat_log
true_values <- mat_log[mask_idx]
mat_log_masked[mask_idx] <- NA

# Imputation on LOG2 data
mat_knn <- impute::impute.knn(mat_log_masked)$data

mat_rf <- missForest::missForest(as.data.frame(mat_log_masked), maxiter = 20)$ximp
mat_rf <- as.matrix(mat_rf)

mat_qrilc <- imputeLCMD::impute.QRILC(mat_log_masked, tune.sigma = 0.5)[[1]]
mat_qrilc <- as.matrix(mat_qrilc)

mat_halfmin <- impute_halfmin_log2(mat_log_masked)

# NRMSE comparison
nrmse_results <- c(
  halfmin = compute_nrmse(true_values, mat_halfmin[mask_idx]),
  knn     = compute_nrmse(true_values, mat_knn[mask_idx]),
  rf      = compute_nrmse(true_values, mat_rf[mask_idx]),
  qrilc   = compute_nrmse(true_values, mat_qrilc[mask_idx])
)

nrmse_results
rank_results <- rank(nrmse_results, ties.method = "average")

rank_results


###########################################################################################################
###########################################################################################################
###########################################################################################################


### RF for BLOOD:

###########################################################################################################
# FINAL IMPUTATION FOR BLOOD NEG (chosen method = RF)

# Convert filtered data to numeric matrix
mat <- as.matrix(signal_neg_metabolites_filtered)
storage.mode(mat) <- "double"
mat[mat == 0] <- NA

mat_log <- log2(mat)
rf_imputed <- missForest::missForest(as.data.frame(mat_log), maxiter = 20)
mat_imputed <- as.matrix(rf_imputed$ximp)

# Restore row/column names
rownames(mat_imputed) <- rownames(mat_log)
colnames(mat_imputed) <- colnames(mat_log)
mat_imputed


# Recombining the original sets


metadata_neg_blood_filtered <- metadata_neg_blood[
  metadata_neg_blood$Name %in% rownames(mat_imputed),
  ,
  drop = FALSE
]

metadata_neg_blood_filtered <- metadata_neg_blood_filtered[
  match(rownames(mat_imputed), metadata_neg_blood_filtered$Name),
  ,
  drop = FALSE
]

dataset_neg_blood_imputed <- cbind(
  metadata_neg_blood_filtered,
  as.data.frame(mat_imputed)
)

write.csv(
  dataset_neg_blood_imputed,
  "dataset_imputed_filtered/neg_Waters_imputed_blood.csv",
  row.names = FALSE
)
