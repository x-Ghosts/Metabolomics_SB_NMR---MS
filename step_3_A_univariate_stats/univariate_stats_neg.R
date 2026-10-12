#remotes::install_github("HimesGroup/qmtools")
#install.packages(c('SummarizedExperiment', 'heatmaply', 'VIM', 'limma'))

setwd("DIRNAME/step_3_A_univariate_stats/")

library(tidyverse)
library(qmtools)
library(dplyr)
library(omicsTools)
library(MetabolAnalyze)
library(ggrepel)
library(pheatmap)

## ==> Data

dataset_neg_blood_selected <- read.csv("dataset_scaled/neg_Waters_scaled_blood.csv")

metadata_neg_blood <- dataset_neg_blood_selected[, 1:5]
signal_neg_metabolites <- dataset_neg_blood_selected[, 6:ncol(dataset_neg_blood_selected)]

sample_names <- colnames(signal_neg_metabolites)
sample_group <- case_when(
  grepl("qc", sample_names, ignore.case = TRUE) ~ "QC",
  grepl("zak", sample_names, ignore.case = TRUE) ~ "Case",
  grepl("biobank", sample_names, ignore.case = TRUE) ~ "Control",
  TRUE ~ NA_character_
)
sample_metadata <- data.frame(
  Sample = sample_names,
  Group = sample_group
)

## ==> Remove the QC:

keep_samples <- sample_metadata$Group != "QC"
mat_stats <- as.matrix(signal_neg_metabolites[, keep_samples])
group <- sample_metadata$Group[keep_samples] # table(group)





pvals <- apply(mat_stats, 1, function(x) {
  t.test(
    x[group == "Case"],
    x[group == "Control"],
    alternative = "two.sided",
    var.equal = FALSE
  )$p.value
})

pvals_adj <- p.adjust(pvals, method = "fdr")

logFC <- apply(mat_stats, 1, function(x) {
  mean(x[group == "Case"], na.rm = TRUE) - mean(x[group == "Control"], na.rm = TRUE)
})

results_univariate <- data.frame(
  Metabolite = metadata_neg_blood$Name,
  logFC = logFC,
  pvalue = pvals,
  FDR = pvals_adj
)

results_univariate <- results_univariate[order(results_univariate$FDR), ]
results_univariate
head(results_univariate)

results_univariate$significance <- ifelse(
  results_univariate$FDR < 0.001, "***",
  ifelse(
    results_univariate$FDR < 0.01, "**",
    ifelse(
      results_univariate$FDR < 0.05 & abs(results_univariate$logFC) >= 1.5, "*",
      ""
    )
  )
)

write.table(results_univariate, "t-test_blood_neg.txt", quote = FALSE, sep = "\t", col.names = T, row.names = F)




################################ Distributions:

dataset_neg_long <- dataset_neg_blood_selected %>%
  pivot_longer(
    cols = 6:ncol(dataset_neg_blood_selected),
    names_to = "sample",
    values_to = "signal"
  ) %>%
  drop_na() %>%
  mutate(
    Type_of_sample = case_when(
      grepl("qc", sample, ignore.case = TRUE) ~ "qc",
      grepl("biobank", sample, ignore.case = TRUE) ~ "biobank",
      grepl("zak", sample, ignore.case = TRUE) ~ "zak",
      TRUE ~ NA_character_
    )
  ) %>%
  relocate(signal, .after = Type_of_sample)

p <- ggplot(dataset_neg_long, aes(x = signal, fill = Type_of_sample)) +
  geom_histogram(
    position = "identity",
    bins = 30,
    alpha = 0.5,
    colour = "black"
  ) +
  facet_wrap(~ Name, scales = "free", ncol = 4) +
  scale_fill_manual(
    values = c(
      "biobank" = "cyan",
      "qc" = "green",
      "zak" = "orange"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Histogram of transformed and scaled metabolite signals (NEG)",
    x = "Signal Intensity",
    y = "Count",
    fill = "Type of Sample"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

ggsave("distributions_blood_neg.png", p, width = 16, height = 10, dpi = 300)



################################ Volcano:


fc_cutoff <- 1.5
FDR_cutoff <- 0.05

volcano_data <- results_univariate %>%
  mutate(
    negLog10FDR = -log10(FDR),
    Status = case_when(
      FDR < FDR_cutoff & logFC >= fc_cutoff ~ "Significant up",
      FDR < FDR_cutoff & logFC <= -fc_cutoff ~ "Significant down",
      TRUE ~ "Not significant"
    )
  )

top_labels <- volcano_data %>%
  filter(FDR < 0.05 & abs(logFC) >= fc_cutoff) %>%
  arrange(FDR)

p_1 <- ggplot(volcano_data, aes(x = logFC, y = negLog10FDR, color = Status)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(
    data = top_labels,
    aes(label = Metabolite),
    size = 4,
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_cutoff), linetype = "dashed") +
  scale_color_manual(values = c(
    "Significant up" = "steelblue",
    "Significant down" = "seagreen4",
    "Not significant" = "goldenrod2"
  )) +
  labs(
    title = "Volcano Plot (Infectious vs Control) -- Blood // Negative",
    x = "log2 Fold Change",
    y = expression(-log[10](FDR)),
    color = "FDR & FC significance"
  ) +
  theme_minimal(base_size = 13)

p_1
ggsave("Volcano_plot_blood_neg.png", p, width = 16, height = 10, dpi = 300)



################################################################################################

sig3 <- results_univariate %>%
  filter(FDR < 0.05, abs(logFC) >= 1.5) %>%
  pull(Metabolite)

sig3

# make sure the stats matrix is numeric first
mat_stats <- as.matrix(mat_stats)
storage.mode(mat_stats) <- "double"

# keep metabolite names
rownames(mat_stats) <- metadata_neg_blood$Name

# z-score per metabolite
mat_heatmap <- t(scale(t(mat_stats)))
mat_heatmap <- mat_heatmap[complete.cases(mat_heatmap), , drop = FALSE]

# force numeric matrix again, just to be safe
mat_heatmap <- as.matrix(mat_heatmap)
storage.mode(mat_heatmap) <- "double"

# column annotation
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(mat_heatmap)

# row annotation
annotation_row <- data.frame(
  Significance = ifelse(rownames(mat_heatmap) %in% sig3, "True", "False")
)
rownames(annotation_row) <- rownames(mat_heatmap)

p <- pheatmap::pheatmap(
  mat_heatmap,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(
    Group = c(Case = "red", Control = "blue"),
    Significance = c("True" = "red", "False" = "grey85")
  ),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  
  fontsize = 20,          # general text
  fontsize_row = 26,      # metabolite labels
  fontsize_col = 18,
  fontsize_number = 18,
  
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  
  color = colorRampPalette(c("seagreen4", "black", "steelblue"))(200),
  main = "Heatmap -- Metabolites in Blood (NEG)"
)

ggsave("Heatmap_plot_blood_neg.png", p, width = 32, height = 24, dpi = 200, bg = "white")


