#remotes::install_github("HimesGroup/qmtools")
#install.packages(c('SummarizedExperiment', 'heatmaply', 'VIM', 'limma'))

library(tidyverse)
library(qmtools)
library(dplyr)
library(omicsTools)
library(MetabolAnalyze)


## ==> Data Transformation to LOG2

dataset_pos_blood_selected <- read.csv("dataset_imputed_filtered/pos_Waters_imputed_blood.csv")

metadata_pos_blood <- dataset_pos_blood_selected[, 1:5]
signal_pos_metabolites <- dataset_pos_blood_selected[, 6:ncol(dataset_pos_blood_selected)]

signal_pos_metabolites_num <- as.data.frame(
  lapply(signal_pos_metabolites, function(x) {
    if (is.character(x)) type.convert(x, as.is = TRUE) else x
  })
)

rownames(signal_pos_metabolites_num) <- metadata_pos_blood$Name

## Current data are LOG2 -> back-transform to raw intensity
signal_mat_log2 <- as.matrix(signal_pos_metabolites_num)
storage.mode(signal_mat_log2) <- "double"

signal_mat_raw <- 2 ^ signal_mat_log2


## Transpose for pqn_normalize():
## rows = samples, columns = metabolites
signal_for_pqn <- as.data.frame(t(signal_mat_raw))

## Add first column as sample IDs (required by pqn_normalize)
signal_for_pqn <- cbind(SampleID = rownames(signal_for_pqn), signal_for_pqn)

signal_pqn <- omicsTools::pqn_normalize(signal_for_pqn)

## Convert output back to data frame
signal_pqn_df <- as.data.frame(signal_pqn)

## Restore sample IDs as row names, remove first column
rownames(signal_pqn_df) <- signal_pqn_df[, 1]
signal_pqn_df <- signal_pqn_df[, -1, drop = FALSE]

## Log2 transform AFTER PQN
signal_pqn_log2 <- log2(signal_pqn_df)

## Return to your usual orientation:
## rows = metabolites, columns = samples
signal_pqn_log2_mat <- t(as.matrix(signal_pqn_log2))






###############################################################################################################"
###############################################################################################################"
###############################################################################################################"
###############################################################################################################"

signal_pqn_log2_mat   # rows = metabolites, cols = samples

############# Centering & Scaling
data_transposed <- t(signal_pqn_log2_mat)
data_centered <- scale(data_transposed, center = TRUE, scale = FALSE)
data_scaled <- scaling(data_centered, type = "pareto")
mat_final <- t(data_scaled)



## Rebuild full dataset
metadata_pos_blood_reordered <- metadata_pos_blood[
  match(rownames(mat_final), metadata_pos_blood$Name),
  ,
  drop = FALSE
]

dataset_pos_blood_scaled <- cbind(
  metadata_pos_blood_reordered,
  as.data.frame(mat_final)
)

write.csv(
  dataset_pos_blood_scaled,
  "dataset_scaled/pos_Waters_scaled_blood.csv",
  row.names = FALSE
)

