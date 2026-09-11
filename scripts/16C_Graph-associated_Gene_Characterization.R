# ==============================================================================
# Stage 16C: Characterization of Graph-associated Genes
#
# Purpose:
#   Characterize the Stage 16B graph_test() results without rerunning
#   graph_test() and without modifying the frozen Monocle3 trajectory.
#
# IMPORTANT:
#   Stage 15 = FROZEN
#   Stage 16A = FROZEN
#   Stage 16B = FROZEN
#
# This stage does NOT:
#   - reconstruct the trajectory
#   - change partitions
#   - change the root
#   - recalculate pseudotime
#   - rerun graph_test()
#   - modify the Stage 16B result
#
# Main goals:
#   1. QC of graph_test output
#   2. Characterize Moran's I distribution
#   3. Characterize p-values and q-values
#   4. Rank trajectory-associated genes
#   5. Examine effect-size distribution
#   6. Generate diagnostic plots
#
# No final biological cutoff is defined in this stage.
#
# Repository architecture:
#
# results/
# └── 16_Pseudotime_Gene_Dynamics/
#     ├── figures/
#     │   ├── main/
#     │   └── supplementary/
#     ├── tables/
#     │   ├── main/
#     │   ├── supplementary/
#     │   └── validation/
#     ├── logs/
#     └── README.md
#
# ==============================================================================


# ==============================================================================
# 1. Load Required Libraries
# ==============================================================================

suppressPackageStartupMessages({
  
  library(dplyr)
  library(ggplot2)
  
})


# ==============================================================================
# 2. Define Project Paths
# ==============================================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"

stage16_dir <- file.path(
  project_dir,
  "results",
  "16_Pseudotime_Gene_Dynamics")

figures_dir <- file.path(
  stage16_dir,
  "figures"
)

figures_main_dir <- file.path(
  figures_dir,
  "main"
)

figures_supp_dir <- file.path(
  figures_dir,
  "supplementary"
)

tables_dir <- file.path(
  stage16_dir,
  "tables"
)

tables_main_dir <- file.path(
  tables_dir,
  "main"
)

tables_supp_dir <- file.path(
  tables_dir,
  "supplementary"
)

tables_validation_dir <- file.path(
  tables_dir,
  "validation"
)

logs_dir <- file.path(
  stage16_dir,
  "logs"
)

objects_dir <- file.path(
  project_dir,
  "objects"
)


# ==============================================================================
# 3. Create Required Directories
# ==============================================================================

dir.create(
  figures_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_supp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_supp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 4. Define Stage 16B Input
# ==============================================================================

stage16b_raw_csv <- file.path(
  tables_supp_dir,
  "01_graph_test_RAW.csv"
)


if (!file.exists(stage16b_raw_csv)) {
  
  stop(
    paste0(
      "\nStage 16B RAW graph_test result was not found.\n\n",
      "Expected file:\n",
      stage16b_raw_csv,
      "\n\n",
      "Stage 16C requires the completed Stage 16B result and does not rerun ",
      "graph_test()."
    )
  )
  
}


# ==============================================================================
# 5. Stage 16C Header
# ==============================================================================

cat("\n============================================================\n")
cat("TNBC Single-Cell Analysis Pipeline\n")
cat("Stage 16C: Graph-associated Gene Characterization\n")
cat("============================================================\n\n")

cat(
  "Stage 16B input:\n",
  stage16b_raw_csv,
  "\n\n"
)

cat(
  "Stage 15: FROZEN\n",
  "Stage 16A: FROZEN\n",
  "Stage 16B: FROZEN\n\n",
  sep = ""
)

cat(
  "graph_test() rerun: NOT PERFORMED\n",
  "Trajectory reconstruction: NOT PERFORMED\n",
  "Pseudotime recalculation: NOT PERFORMED\n\n",
  sep = ""
)


# ==============================================================================
# 6. Load Stage 16B RAW graph_test Result
# ==============================================================================

graph_test_df <- read.csv(
  stage16b_raw_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (
  is.null(graph_test_df) ||
  nrow(graph_test_df) == 0
) {
  
  stop(
    "Stage 16B graph_test result is empty."
  )
  
}


# ==============================================================================
# 7. Recover Gene ID Safely
# ==============================================================================

if ("gene_id" %in% colnames(graph_test_df)) {
  
  graph_test_df$gene_id <- as.character(
    graph_test_df$gene_id
  )
  
} else {
  
  possible_gene_id_columns <- intersect(
    c("id", "ID", "gene", "Gene"),
    colnames(graph_test_df)
  )
  
  if (length(possible_gene_id_columns) > 0) {
    
    graph_test_df$gene_id <- as.character(
      graph_test_df[[possible_gene_id_columns[1]]]
    )
    
  } else {
    
    graph_test_df$gene_id <- seq_len(
      nrow(graph_test_df)
    ) |>
      as.character()
    
  }
  
}


# ==============================================================================
# 8. Validate Required Columns
# ==============================================================================

required_columns <- c(
  "status",
  "p_value",
  "morans_test_statistic",
  "morans_I",
  "gene_short_name",
  "q_value"
)


missing_columns <- setdiff(
  required_columns,
  colnames(graph_test_df)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "\nRequired Stage 16B graph_test columns are missing:\n",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "\n\nAvailable columns:\n",
      paste(
        colnames(graph_test_df),
        collapse = ", "
      )
    )
  )
  
}


# ==============================================================================
# 9. Standardize Data Types
# ==============================================================================

graph_test_df$gene_id <- as.character(
  graph_test_df$gene_id
)

graph_test_df$gene_short_name <- as.character(
  graph_test_df$gene_short_name
)

graph_test_df$status <- as.character(
  graph_test_df$status
)


numeric_columns <- c(
  "p_value",
  "morans_test_statistic",
  "morans_I",
  "q_value"
)


for (col in numeric_columns) {
  
  graph_test_df[[col]] <- suppressWarnings(
    as.numeric(
      graph_test_df[[col]]
    )
  )
  
}


# ==============================================================================
# 10. Basic Input Validation
# ==============================================================================

input_validation <- data.frame(
  
  Check = c(
    "Input_file_exists",
    "Input_has_rows",
    "Required_columns_present",
    "Unique_gene_ids",
    "NonNA_Morans_I_available",
    "NonNA_p_values_available",
    "NonNA_q_values_available"
  ),
  
  Result = c(
    file.exists(stage16b_raw_csv),
    nrow(graph_test_df) > 0,
    length(missing_columns) == 0,
    length(unique(graph_test_df$gene_id)) ==
      nrow(graph_test_df),
    sum(!is.na(graph_test_df$morans_I)) > 0,
    sum(!is.na(graph_test_df$p_value)) > 0,
    sum(!is.na(graph_test_df$q_value)) > 0
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  input_validation,
  file.path(
    tables_validation_dir,
    "Stage16C_Input_Validation.csv"
  ),
  row.names = FALSE
)


if (!all(input_validation$Result)) {
  
  stop(
    "Stage 16C input validation failed."
  )
  
}


# ==============================================================================
# 11. Basic Result QC
# ==============================================================================

result_qc <- data.frame(
  
  Metric = c(
    "Total_genes",
    "Genes_with_nonNA_p_value",
    "Genes_with_nonNA_q_value",
    "Genes_with_nonNA_Morans_I",
    "Genes_with_q_lt_0.05",
    "Genes_with_Morans_I_ge_0.10",
    "Genes_with_Morans_I_gt_0",
    "Genes_with_Morans_I_le_0",
    "NA_p_value",
    "NA_q_value",
    "NA_Morans_I"
  ),
  
  Value = c(
    
    nrow(graph_test_df),
    
    sum(
      !is.na(graph_test_df$p_value)
    ),
    
    sum(
      !is.na(graph_test_df$q_value)
    ),
    
    sum(
      !is.na(graph_test_df$morans_I)
    ),
    
    sum(
      graph_test_df$q_value < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      graph_test_df$morans_I >= 0.10,
      na.rm = TRUE
    ),
    
    sum(
      graph_test_df$morans_I > 0,
      na.rm = TRUE
    ),
    
    sum(
      graph_test_df$morans_I <= 0,
      na.rm = TRUE
    ),
    
    sum(
      is.na(graph_test_df$p_value)
    ),
    
    sum(
      is.na(graph_test_df$q_value)
    ),
    
    sum(
      is.na(graph_test_df$morans_I)
    )
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  result_qc,
  file.path(
    tables_main_dir,
    "01_Stage16C_GraphTest_Result_QC.csv"
  ),
  row.names = FALSE
)


cat("Genes loaded from Stage 16B:", nrow(graph_test_df), "\n\n")

print(result_qc)


# ==============================================================================
# 12. Examine graph_test Status
# ==============================================================================

status_summary <- graph_test_df %>%
  
  count(
    status,
    name = "n"
  ) %>%
  
  arrange(
    desc(n)
  )


write.csv(
  status_summary,
  file.path(
    tables_supp_dir,
    "02_Stage16C_GraphTest_Status_Summary.csv"
  ),
  row.names = FALSE
)


cat("\nGraph_test status:\n")
print(status_summary)


# ==============================================================================
# 13. Summary Statistics of Moran's I
# ==============================================================================

morans_values <- graph_test_df$morans_I[
  is.finite(
    graph_test_df$morans_I
  )
]


if (length(morans_values) == 0) {
  
  stop(
    "No finite Moran's I values are available."
  )
  
}


morans_summary <- data.frame(
  
  Statistic = c(
    "Min",
    "1st_Quartile",
    "Median",
    "Mean",
    "3rd_Quartile",
    "95th_Percentile",
    "99th_Percentile",
    "Max"
  ),
  
  Value = c(
    
    min(
      morans_values
    ),
    
    quantile(
      morans_values,
      0.25
    ),
    
    median(
      morans_values
    ),
    
    mean(
      morans_values
    ),
    
    quantile(
      morans_values,
      0.75
    ),
    
    quantile(
      morans_values,
      0.95
    ),
    
    quantile(
      morans_values,
      0.99
    ),
    
    max(
      morans_values
    )
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  morans_summary,
  file.path(
    tables_main_dir,
    "03_Morans_I_Summary.csv"
  ),
  row.names = FALSE
)


print(morans_summary)


# ==============================================================================
# 14. Summary Statistics of P-values
# ==============================================================================

pvalue_values <- graph_test_df$p_value[
  is.finite(
    graph_test_df$p_value
  )
]


pvalue_summary <- data.frame(
  
  Statistic = c(
    "Min",
    "1st_Quartile",
    "Median",
    "Mean",
    "3rd_Quartile",
    "95th_Percentile",
    "99th_Percentile",
    "Max"
  ),
  
  Value = c(
    
    min(
      pvalue_values
    ),
    
    quantile(
      pvalue_values,
      0.25
    ),
    
    median(
      pvalue_values
    ),
    
    mean(
      pvalue_values
    ),
    
    quantile(
      pvalue_values,
      0.75
    ),
    
    quantile(
      pvalue_values,
      0.95
    ),
    
    quantile(
      pvalue_values,
      0.99
    ),
    
    max(
      pvalue_values
    )
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  pvalue_summary,
  file.path(
    tables_main_dir,
    "04_P_value_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Summary Statistics of q-values
# ==============================================================================

qvalue_values <- graph_test_df$q_value[
  is.finite(
    graph_test_df$q_value
  )
]


qvalue_summary <- data.frame(
  
  Statistic = c(
    "Min",
    "1st_Quartile",
    "Median",
    "Mean",
    "3rd_Quartile",
    "95th_Percentile",
    "99th_Percentile",
    "Max"
  ),
  
  Value = c(
    
    min(
      qvalue_values
    ),
    
    quantile(
      qvalue_values,
      0.25
    ),
    
    median(
      qvalue_values
    ),
    
    mean(
      qvalue_values
    ),
    
    quantile(
      qvalue_values,
      0.75
    ),
    
    quantile(
      qvalue_values,
      0.95
    ),
    
    quantile(
      qvalue_values,
      0.99
    ),
    
    max(
      qvalue_values
    )
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  qvalue_summary,
  file.path(
    tables_main_dir,
    "05_Q_value_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Rank All Genes by Moran's I
# ==============================================================================

ranked_genes <- graph_test_df %>%
  
  filter(
    is.finite(morans_I)
  ) %>%
  
  arrange(
    desc(morans_I),
    q_value,
    p_value
  )


write.csv(
  ranked_genes,
  file.path(
    tables_supp_dir,
    "06_All_Genes_Ranked_by_Morans_I.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. Top 100 Genes by Moran's I
# ==============================================================================

top100_morans <- ranked_genes %>%
  
  slice_head(
    n = 100
  )


write.csv(
  top100_morans,
  file.path(
    tables_main_dir,
    "07_Top100_Morans_I_Genes.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. Top 50 Genes by Moran's I
# ==============================================================================

top50_morans <- ranked_genes %>%
  
  slice_head(
    n = 50
  )


write.csv(
  top50_morans,
  file.path(
    tables_main_dir,
    "08_Top50_Morans_I_Genes.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 19. Genes with q < 0.05
# ==============================================================================

q05_genes <- graph_test_df %>%
  
  filter(
    !is.na(q_value),
    q_value < 0.05,
    is.finite(morans_I)
  ) %>%
  
  arrange(
    desc(morans_I),
    q_value,
    p_value
  )


write.csv(
  q05_genes,
  file.path(
    tables_main_dir,
    "09_qvalue_lt_0.05_Genes.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 20. Distribution of Moran's I by Effect-size Bins
#
# Exploratory only.
# These bins are NOT biological cutoffs.
# ==============================================================================

morans_bins <- graph_test_df %>%
  
  filter(
    is.finite(morans_I)
  ) %>%
  
  mutate(
    
    Morans_I_Bin = case_when(
      
      morans_I <= 0 ~
        "<= 0",
      
      morans_I > 0 &
        morans_I < 0.05 ~
        "0–<0.05",
      
      morans_I >= 0.05 &
        morans_I < 0.10 ~
        "0.05–<0.10",
      
      morans_I >= 0.10 &
        morans_I < 0.20 ~
        "0.10–<0.20",
      
      morans_I >= 0.20 &
        morans_I < 0.30 ~
        "0.20–<0.30",
      
      morans_I >= 0.30 ~
        ">=0.30",
      
      TRUE ~
        NA_character_
      
    )
    
  ) %>%
  
  count(
    Morans_I_Bin,
    name = "n"
  )


write.csv(
  morans_bins,
  file.path(
    tables_supp_dir,
    "10_Morans_I_Effect_Size_Bins.csv"
  ),
  row.names = FALSE
)


print(morans_bins)


# ==============================================================================
# 21. Top Genes with Gene Symbols
# ==============================================================================

top100_gene_symbols <- top100_morans %>%
  
  select(
    gene_id,
    gene_short_name,
    status,
    morans_I,
    morans_test_statistic,
    p_value,
    q_value
  )


write.csv(
  top100_gene_symbols,
  file.path(
    tables_main_dir,
    "11_Top100_Genes_Compact.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 22. Moran's I Distribution Plot
# ==============================================================================

morans_plot_df <- graph_test_df %>%
  
  filter(
    is.finite(morans_I)
  )


p_morans <- ggplot(
  morans_plot_df,
  aes(
    x = morans_I
  )
) +
  
  geom_histogram(
    bins = 80
  ) +
  
  geom_vline(
    xintercept = 0.10,
    linetype = "dashed"
  ) +
  
  labs(
    title = "Distribution of Moran's I",
    x = "Moran's I",
    y = "Number of genes"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title =
      element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
  )


ggsave(
  file.path(
    figures_supp_dir,
    "01_Morans_I_Distribution.pdf"
  ),
  p_morans,
  width = 8,
  height = 6
)


# ==============================================================================
# 23. P-value Distribution Plot
# ==============================================================================

pvalue_plot_df <- graph_test_df %>%
  
  filter(
    is.finite(p_value),
    p_value > 0
  )


p_pvalue <- ggplot(
  pvalue_plot_df,
  aes(
    x = -log10(p_value)
  )
) +
  
  geom_histogram(
    bins = 80
  ) +
  
  labs(
    title = "Distribution of graph_test P-values",
    x = expression(-log[10](p)),
    y = "Number of genes"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title =
      element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
  )


ggsave(
  file.path(
    figures_supp_dir,
    "02_P_value_Distribution.pdf"
  ),
  p_pvalue,
  width = 8,
  height = 6
)


# ==============================================================================
# 24. q-value Distribution Plot
# ==============================================================================

qvalue_plot_df <- graph_test_df %>%
  
  filter(
    is.finite(q_value),
    q_value > 0
  )


p_qvalue <- ggplot(
  qvalue_plot_df,
  aes(
    x = -log10(q_value)
  )
) +
  
  geom_histogram(
    bins = 80
  ) +
  
  labs(
    title = "Distribution of graph_test q-values",
    x = expression(-log[10](q)),
    y = "Number of genes"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title =
      element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
  )


ggsave(
  file.path(
    figures_supp_dir,
    "03_Q_value_Distribution.pdf"
  ),
  p_qvalue,
  width = 8,
  height = 6
)


# ==============================================================================
# 25. Moran's I versus -log10(q-value)
# ==============================================================================

plot_df <- graph_test_df %>%
  
  filter(
    is.finite(morans_I),
    is.finite(q_value),
    q_value > 0
  )


p_morans_q <- ggplot(
  plot_df,
  aes(
    x = morans_I,
    y = -log10(q_value)
  )
) +
  
  geom_point(
    alpha = 0.35,
    size = 0.7
  ) +
  
  geom_vline(
    xintercept = 0.10,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  labs(
    title = "Moran's I versus Statistical Significance",
    x = "Moran's I",
    y = expression(-log[10](q))
  ) +
  
  theme_bw() +
  
  theme(
    plot.title =
      element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
  )


ggsave(
  file.path(
    figures_main_dir,
    "04_Morans_I_vs_Q_value.pdf"
  ),
  p_morans_q,
  width = 8,
  height = 6
)


# ==============================================================================
# 26. Top 30 Moran's I Genes
# ==============================================================================

top30_plot_df <- ranked_genes %>%
  
  slice_head(
    n = 30
  ) %>%
  
  mutate(
    
    gene_label = ifelse(
      is.na(gene_short_name) |
        gene_short_name == "",
      gene_id,
      gene_short_name
    )
    
  ) %>%
  
  arrange(
    morans_I
  )


p_top30 <- ggplot(
  top30_plot_df,
  aes(
    x = morans_I,
    y = reorder(
      gene_label,
      morans_I
    )
  )
) +
  
  geom_point(
    size = 2
  ) +
  
  labs(
    title = "Top 30 Genes by Moran's I",
    x = "Moran's I",
    y = "Gene"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title =
      element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
  )


ggsave(
  file.path(
    figures_main_dir,
    "05_Top30_Morans_I_Genes.pdf"
  ),
  p_top30,
  width = 9,
  height = 8
)


# ==============================================================================
# 27. Stage 16C Characterization Summary
# ==============================================================================

stage16c_summary <- data.frame(
  
  Metric = c(
    "Stage16B_Input_Genes",
    "Genes_with_finite_Morans_I",
    "Genes_with_q_lt_0.05",
    "Genes_with_Morans_I_ge_0.10",
    "Genes_with_Morans_I_gt_0",
    "Genes_with_Morans_I_le_0",
    "Top100_Genes_Generated",
    "Top50_Genes_Generated"
  ),
  
  Value = c(
    
    nrow(graph_test_df),
    
    length(morans_values),
    
    sum(
      graph_test_df$q_value < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      graph_test_df$morans_I >= 0.10,
      na.rm = TRUE
    ),
    
    sum(
      graph_test_df$morans_I > 0,
      na.rm = TRUE
    ),
    
    sum(
      graph_test_df$morans_I <= 0,
      na.rm = TRUE
    ),
    
    nrow(top100_morans),
    
    nrow(top50_morans)
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  stage16c_summary,
  file.path(
    tables_main_dir,
    "12_Stage16C_Characterization_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 28. Save Characterization Object
#
# Stored in the central repository objects directory.
# This does NOT modify any frozen Stage 15 or Stage 16B object.
# ==============================================================================

stage16c_results <- list(
  
  graph_test = graph_test_df,
  
  input_validation = input_validation,
  
  result_qc = result_qc,
  
  status_summary = status_summary,
  
  morans_summary = morans_summary,
  
  pvalue_summary = pvalue_summary,
  
  qvalue_summary = qvalue_summary,
  
  ranked_genes = ranked_genes,
  
  top100_morans = top100_morans,
  
  top50_morans = top50_morans,
  
  q05_genes = q05_genes,
  
  morans_bins = morans_bins,
  
  characterization_summary = stage16c_summary
  
)


stage16c_rds <- file.path(
  objects_dir,
  "Stage16C_Characterization_Results.rds"
)


saveRDS(
  stage16c_results,
  stage16c_rds
)


# ==============================================================================
# 29. Stage 16C Completion Log
# ==============================================================================

stage16c_log <- file.path(
  logs_dir,
  "Stage16C_Completion.log"
)


log_lines <- c(
  
  "============================================================",
  
  "TNBC Single-Cell Analysis Pipeline",
  
  "Stage 16C: Graph-associated Gene Characterization",
  
  "============================================================",
  
  "",
  
  paste(
    "Completion time:",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  paste(
    "Stage 16B input:",
    stage16b_raw_csv
  ),
  
  paste(
    "Genes analyzed:",
    nrow(graph_test_df)
  ),
  
  paste(
    "Genes with finite Moran's I:",
    length(morans_values)
  ),
  
  paste(
    "Genes with q < 0.05:",
    sum(
      graph_test_df$q_value < 0.05,
      na.rm = TRUE
    )
  ),
  
  paste(
    "Genes with Moran's I >= 0.10:",
    sum(
      graph_test_df$morans_I >= 0.10,
      na.rm = TRUE
    )
  ),
  
  "",
  
  "Stage 15: FROZEN",
  
  "Stage 16A: FROZEN",
  
  "Stage 16B: FROZEN",
  
  "Stage 16C: COMPLETED",
  
  "",
  
  "Trajectory reconstruction: NOT PERFORMED",
  
  "Root selection: NOT PERFORMED",
  
  "Pseudotime recalculation: NOT PERFORMED",
  
  "graph_test(): NOT RERUN",
  
  "Frozen Stage 16B result: NOT MODIFIED",
  
  "",
  
  "STATUS: STAGE 16C COMPLETED SUCCESSFULLY",
  
  "============================================================"
  
)


writeLines(
  log_lines,
  con = stage16c_log
)


# ==============================================================================
# 30. Final Console Summary
# ==============================================================================

cat("\n============================================================\n")
cat("STAGE 16C COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat(
  "Genes analyzed:",
  nrow(graph_test_df),
  "\n"
)

cat(
  "Genes with finite Moran's I:",
  length(morans_values),
  "\n"
)

cat(
  "Genes with q < 0.05:",
  sum(
    graph_test_df$q_value < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Genes with Moran's I >= 0.10:",
  sum(
    graph_test_df$morans_I >= 0.10,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nStage 16C RDS:\n",
  stage16c_rds,
  "\n"
)

cat(
  "\nStage 16C log:\n",
  stage16c_log,
  "\n"
)

cat("\n============================================================\n")

