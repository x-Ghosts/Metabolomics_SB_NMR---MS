## Summary Results of the Imputations Algorthims
library(tidyverse)

results <- data.frame(
  Dataset = c(
    "Blood (-)", "Blood (+)",
    "Spine Mehl (-)", "Spine Mehl (+)",
    "Spine Waters (-)", "Spine Waters (+)"
  ),
  Patients = c(89, 89, 46, 46, 45, 25),
  Metabolites = c(20, 56, 8, 39, 9, 38),
  HalfMin = c(1.313477, 1.4678345, 1.832166, 1.4189394, 1.443531, 1.2272372),
  KNN     = c(0.624945, 0.5336926, 1.498500, 0.6043642, 1.112500, 0.5565355),
  RF      = c(0.574211, 0.4152741, 1.155859, 0.5595743, 1.088890, 0.5456310),
  QRILC   = c(2.007093, 2.0203112, 1.763094, 1.9404639, 1.819287, 2.0516730)
)
results$Category <- ifelse(grepl("Blood", results$Dataset), "Blood", "Spine")
rank_results <- rank(nrmse_results, ties.method = "average")
rank_results
category_means <- results %>%
  group_by(Category) %>%
  summarise(
    HalfMin = mean(HalfMin, na.rm = TRUE),
    KNN     = mean(KNN, na.rm = TRUE),
    RF      = mean(RF, na.rm = TRUE),
    QRILC   = mean(QRILC, na.rm = TRUE)
  )


category_means



####################################"


# Long format + rank per dataset
results_ranked <- results %>%
  pivot_longer(
    cols = c(HalfMin, KNN, RF, QRILC),
    names_to = "Method",
    values_to = "NRMSE"
  ) %>%
  group_by(Dataset) %>%
  mutate(Rank = rank(NRMSE, ties.method = "average")) %>%
  ungroup()

# Table 1: dataset-level results
dataset_table <- results_ranked %>%
  group_by(Dataset, Patients, Metabolites, Category) %>%
  summarise(
    HalfMin_NRMSE = NRMSE[Method == "HalfMin"],
    HalfMin_Rank  = Rank[Method == "HalfMin"],
    KNN_NRMSE     = NRMSE[Method == "KNN"],
    KNN_Rank      = Rank[Method == "KNN"],
    RF_NRMSE      = NRMSE[Method == "RF"],
    RF_Rank       = Rank[Method == "RF"],
    QRILC_NRMSE   = NRMSE[Method == "QRILC"],
    QRILC_Rank    = Rank[Method == "QRILC"],
    Best_Method   = Method[which.min(NRMSE)],
    Best_NRMSE    = min(NRMSE),
    .groups = "drop"
  )

dataset_table

sum_ranks_table <- results_ranked %>%
  group_by(Method) %>%
  summarise(
    Sum_Rank = sum(Rank),
    Mean_Rank = mean(Rank),
    Mean_NRMSE = mean(NRMSE),
    .groups = "drop"
  ) %>%
  arrange(Sum_Rank)

sum_ranks_table


write.table(dataset_table, "dataset_imputation_results.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sum_ranks_table, "sum_ranks_results.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
