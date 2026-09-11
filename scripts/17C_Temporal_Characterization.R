# ============================================================
# STAGE 17C
# TEMPORAL GENE PROGRAM CHARACTERIZATION
#
# Purpose:
#   Characterize and summarize the temporal gene clusters
#   identified in Stage 17B.
#
# Input:
#   Frozen Stage 17B clustering object
#
# Outputs:
#   - Cluster summary tables
#   - Temporal pattern classification
#   - Top genes per cluster
#   - Cluster centroid profiles
#   - Gene-level temporal profiles
#   - Main and supplementary PDF figures
#   - Validation tables
#   - Central Stage 17C RDS object
#
# Important:
#   - Stage 17B is treated as frozen.
#   - No trajectory reconstruction is performed.
#   - No root selection is performed.
#   - No pseudotime recalculation is performed.
#   - No upstream stage is rerun.
# ============================================================


# ============================================================
# 0. INITIALIZATION
# ============================================================

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  width = 120
)

set.seed(12345)


# ============================================================
# 1. PACKAGE LOADING
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Required packages are missing: ",
      paste(missing_packages, collapse = ", "),
      ". Please install them before running Stage 17C."
    )
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"
results_dir <-
  file.path(
    project_dir,
    "results",
    "17_Temporal_Gene_Dynamics"
  )

objects_dir <-
  file.path(
    project_dir,
    "objects"
  )

figures_main_dir <-
  file.path(
    results_dir,
    "figures",
    "main"
  )

figures_supp_dir <-
  file.path(
    results_dir,
    "figures",
    "supplementary"
  )

tables_main_dir <-
  file.path(
    results_dir,
    "tables",
    "main"
  )

tables_supp_dir <-
  file.path(
    results_dir,
    "tables",
    "supplementary"
  )

tables_validation_dir <-
  file.path(
    results_dir,
    "tables",
    "validation"
  )

logs_dir <-
  file.path(
    results_dir,
    "logs"
  )


# Create canonical directories if necessary

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


# ============================================================
# 3. LOGGING
# ============================================================

stage17c_log <-
  file.path(
    logs_dir,
    "Stage17C_Completion.log"
  )

writeLines(
  character(0),
  con = stage17c_log
)

write_log <- function(...) {
  
  txt <-
    paste0(
      ...,
      collapse = ""
    )
  
  cat(
    txt,
    "\n"
  )
  
  write(
    txt,
    file = stage17c_log,
    append = TRUE
  )
}


write_log(
  "============================================================"
)

write_log(
  "STAGE 17C - TEMPORAL GENE PROGRAM CHARACTERIZATION"
)

write_log(
  "============================================================"
)

write_log(
  "Project directory: ",
  project_dir
)

write_log(
  "Results directory: ",
  results_dir
)


# ============================================================
# 4. LOAD FROZEN STAGE 17B OBJECT
# ============================================================

stage17b_rds <-
  file.path(
    objects_dir,
    "Stage17B_Temporal_Clustering_Results.rds"
  )

if (!file.exists(stage17b_rds)) {
  
  stop(
    paste0(
      "Frozen Stage 17B RDS was not found:\n",
      stage17b_rds
    )
  )
}

write_log(
  "Loading frozen Stage 17B object..."
)

stage17b_results <-
  readRDS(
    stage17b_rds
  )

write_log(
  "Stage 17B object loaded successfully."
)


# ============================================================
# 5. VERIFY STAGE 17B OBJECT
# ============================================================

required_fields <- c(
  "stage",
  "analysis_name",
  "synchronized_cells",
  "pseudotime",
  "pseudotime_bins",
  "bin_breaks",
  "bin_counts",
  "n_bins",
  "strong_temporal_genes",
  "strong_temporal_genes_in_expression",
  "mean_expression_matrix",
  "zscore_matrix",
  "k_candidates",
  "silhouette_results",
  "selected_k",
  "selected_mean_silhouette",
  "final_kmeans",
  "gene_clusters",
  "cluster_centroids",
  "cluster_assignments",
  "final_silhouette_width",
  "gene_cluster_correlation",
  "gene_cluster_correlation_table",
  "top_genes_per_cluster",
  "cluster_summary",
  "validation_results",
  "validation_checks",
  "methodology",
  "upstream_frozen",
  "trajectory_reconstruction",
  "root_selection",
  "pseudotime_recalculation"
)

missing_fields <-
  setdiff(
    required_fields,
    names(stage17b_results)
  )

if (length(missing_fields) > 0) {
  
  stop(
    paste0(
      "Stage 17B object is missing required fields:\n",
      paste(
        missing_fields,
        collapse = ", "
      )
    )
  )
}


write_log(
  "Stage 17B schema validation passed."
)


# ============================================================
# 6. VERIFY FROZEN-UPSTREAM STATUS
# ============================================================

if (!isTRUE(stage17b_results$upstream_frozen)) {
  
  stop(
    "Stage 17B does not indicate that upstream analysis is frozen."
  )
}

if (!identical(
  stage17b_results$trajectory_reconstruction,
  FALSE
)) {
  
  stop(
    "Unexpected trajectory_reconstruction status in Stage 17B object."
  )
}

if (!identical(
  stage17b_results$root_selection,
  FALSE
)) {
  
  stop(
    "Unexpected root_selection status in Stage 17B object."
  )
}

if (!identical(
  stage17b_results$pseudotime_recalculation,
  FALSE
)) {
  
  stop(
    "Unexpected pseudotime_recalculation status in Stage 17B object."
  )
}


write_log(
  "Frozen-upstream status verified."
)

write_log(
  "No trajectory reconstruction will be performed."
)

write_log(
  "No root selection will be performed."
)

write_log(
  "No pseudotime recalculation will be performed."
)


# ============================================================
# 7. EXTRACT STAGE 17B OBJECT COMPONENTS
# ============================================================

cluster_assignments <-
  as.data.frame(
    stage17b_results$cluster_assignments,
    stringsAsFactors = FALSE
  )

zscore_matrix <-
  as.matrix(
    stage17b_results$zscore_matrix
  )

mean_expression_matrix <-
  as.matrix(
    stage17b_results$mean_expression_matrix
  )

cluster_centroids <-
  as.matrix(
    stage17b_results$cluster_centroids
  )

gene_cluster_correlation_table <-
  as.data.frame(
    stage17b_results$gene_cluster_correlation_table,
    stringsAsFactors = FALSE
  )

synchronized_cells <-
  as.data.frame(
    stage17b_results$synchronized_cells,
    stringsAsFactors = FALSE
  )

pseudotime <-
  as.numeric(
    stage17b_results$pseudotime
  )

pseudotime_bins <-
  as.integer(
    stage17b_results$pseudotime_bins
  )

bin_breaks <-
  as.numeric(
    stage17b_results$bin_breaks
  )

bin_counts <-
  stage17b_results$bin_counts

n_bins <-
  as.integer(
    stage17b_results$n_bins
  )

strong_temporal_genes <-
  as.character(
    stage17b_results$strong_temporal_genes
  )

selected_k <-
  as.integer(
    stage17b_results$selected_k
  )

selected_mean_silhouette <-
  as.numeric(
    stage17b_results$selected_mean_silhouette
  )


# ============================================================
# 8. STRUCTURAL VALIDATION
# ============================================================

write_log(
  "Running structural validation..."
)


if (nrow(cluster_assignments) != 569) {
  
  stop(
    paste0(
      "Unexpected number of clustered genes: ",
      nrow(cluster_assignments),
      ". Expected 569."
    )
  )
}


if (nrow(zscore_matrix) != 569) {
  
  stop(
    "Unexpected number of rows in zscore_matrix."
  )
}


if (ncol(zscore_matrix) != n_bins) {
  
  stop(
    "zscore_matrix column count does not match n_bins."
  )
}


if (nrow(mean_expression_matrix) != 569) {
  
  stop(
    "Unexpected number of rows in mean_expression_matrix."
  )
}


if (ncol(mean_expression_matrix) != n_bins) {
  
  stop(
    "mean_expression_matrix column count does not match n_bins."
  )
}


if (nrow(cluster_centroids) != selected_k) {
  
  stop(
    "cluster_centroids row count does not match selected_k."
  )
}


if (ncol(cluster_centroids) != n_bins) {
  
  stop(
    "cluster_centroids column count does not match n_bins."
  )
}


if (length(strong_temporal_genes) != 569) {
  
  stop(
    "Unexpected number of Stage 17 strong temporal genes."
  )
}


if (length(pseudotime) != nrow(synchronized_cells)) {
  
  stop(
    "Pseudotime length does not match synchronized-cell count."
  )
}


if (length(pseudotime_bins) != length(pseudotime)) {
  
  stop(
    "Pseudotime-bin vector length does not match pseudotime."
  )
}


if (length(bin_breaks) != n_bins + 1) {
  
  stop(
    "bin_breaks length is inconsistent with n_bins."
  )
}


if (length(bin_counts) != n_bins) {
  
  stop(
    "bin_counts length is inconsistent with n_bins."
  )
}


if (!all(
  cluster_assignments$Cluster %in%
  seq_len(selected_k)
)) {
  
  stop(
    "Invalid cluster labels detected."
  )
}


write_log(
  "Structural validation passed."
)


# ============================================================
# 9. STANDARDIZE GENE ORDER
# ============================================================

gene_names <-
  rownames(
    zscore_matrix
  )

if (is.null(gene_names)) {
  
  stop(
    "zscore_matrix does not contain gene row names."
  )
}


if (anyDuplicated(gene_names) > 0) {
  
  stop(
    "Duplicate gene names detected in zscore_matrix."
  )
}


if (!setequal(
  gene_names,
  cluster_assignments$Gene
)) {
  
  stop(
    "Gene identities in zscore_matrix and cluster_assignments do not match."
  )
}


if (!setequal(
  gene_names,
  gene_cluster_correlation_table$Gene
)) {
  
  stop(
    "Gene identities in zscore_matrix and correlation table do not match."
  )
}


# Reorder cluster assignments to zscore_matrix row order

cluster_assignments <-
  cluster_assignments[
    match(
      gene_names,
      cluster_assignments$Gene
    ),
    ,
    drop = FALSE
  ]


gene_cluster_correlation_table <-
  gene_cluster_correlation_table[
    match(
      gene_names,
      gene_cluster_correlation_table$Gene
    ),
    ,
    drop = FALSE
  ]


# ============================================================
# 10. VERIFY CLUSTER ASSIGNMENTS AGAINST STAGE 17B K-MEANS
# ============================================================

kmeans_clusters <-
  as.integer(
    stage17b_results$final_kmeans$cluster
  )

names(kmeans_clusters) <-
  names(
    stage17b_results$final_kmeans$cluster
  )

kmeans_clusters <-
  kmeans_clusters[
    gene_names
  ]

stored_clusters <-
  as.integer(
    stage17b_results$gene_clusters[
      gene_names
    ]
  )

if (!all(
  kmeans_clusters ==
  stored_clusters
)) {
  
  stop(
    "Stored gene_clusters do not match final_kmeans$cluster."
  )
}


if (!all(
  cluster_assignments$Cluster ==
  stored_clusters
)) {
  
  stop(
    "cluster_assignments do not match Stage 17B gene_clusters."
  )
}


write_log(
  "Cluster-assignment consistency verified."
)


# ============================================================
# 11. CONSTRUCT BIN LABELS AND CENTERS
# ============================================================

bin_labels <-
  paste0(
    "Bin_",
    seq_len(n_bins)
  )

bin_centers <-
  (
    bin_breaks[-length(bin_breaks)] +
      bin_breaks[-1]
  ) / 2


# ============================================================
# 12. VERIFY CENTROIDS AGAINST GENE PROFILES
# ============================================================

recomputed_centroids <-
  matrix(
    NA_real_,
    nrow = selected_k,
    ncol = n_bins
  )

rownames(recomputed_centroids) <-
  rownames(
    cluster_centroids
  )

colnames(recomputed_centroids) <-
  colnames(
    cluster_centroids
  )


for (k in seq_len(selected_k)) {
  
  cluster_genes <-
    gene_names[
      stored_clusters == k
    ]
  
  if (length(cluster_genes) == 0) {
    
    stop(
      paste0(
        "Cluster ",
        k,
        " contains zero genes."
      )
    )
    
  }
  
  recomputed_centroids[k, ] <-
    colMeans(
      zscore_matrix[
        cluster_genes,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )
}


centroid_difference <-
  max(
    abs(
      recomputed_centroids -
        cluster_centroids
    ),
    na.rm = TRUE
  )


if (!is.finite(centroid_difference)) {
  
  stop(
    "Centroid validation produced a non-finite difference."
  )
}


if (centroid_difference > 1e-8) {
  
  write_log(
    "WARNING: Recomputed centroids differ slightly from stored centroids."
  )
  
  write_log(
    "Maximum absolute centroid difference: ",
    format(
      centroid_difference,
      scientific = TRUE
    )
  )
  
} else {
  
  write_log(
    "Stored cluster centroids validated against z-score profiles."
  )
}


# ============================================================
# 13. CLUSTER SIZE SUMMARY
# ============================================================

cluster_size_summary <-
  cluster_assignments %>%
  count(
    Cluster,
    name = "Gene_Count"
  ) %>%
  arrange(
    Cluster
  ) %>%
  mutate(
    Percentage =
      100 *
      Gene_Count /
      sum(Gene_Count)
  )


# ============================================================
# 14. DETERMINE TEMPORAL PATTERN OF EACH CLUSTER
# ============================================================

classify_temporal_pattern <- function(profile) {
  
  n <- length(profile)
  
  if (n < 4) {
    
    return(
      "Insufficient_Bins"
    )
  }
  
  early_n <-
    max(
      2,
      floor(
        n * 0.20
      )
    )
  
  late_n <-
    max(
      2,
      floor(
        n * 0.20
      )
    )
  
  early_profile <-
    profile[
      seq_len(
        early_n
      )
    ]
  
  late_profile <-
    profile[
      (n - late_n + 1):n
    ]
  
  early_mean <-
    mean(
      early_profile,
      na.rm = TRUE
    )
  
  late_mean <-
    mean(
      late_profile,
      na.rm = TRUE
    )
  
  early_max_index <-
    which.max(
      profile
    )
  
  early_min_index <-
    which.min(
      profile
    )
  
  overall_slope <-
    coef(
      lm(
        profile ~ seq_along(profile)
      )
    )[2]
  
  early_peak <-
    early_max_index <=
    ceiling(
      n * 0.30
    )
  
  late_peak <-
    early_max_index >=
    ceiling(
      n * 0.70
    )
  
  early_trough <-
    early_min_index <=
    ceiling(
      n * 0.30
    )
  
  late_trough <-
    early_min_index >=
    ceiling(
      n * 0.70
    )
  
  if (
    early_peak &&
    early_mean > late_mean &&
    overall_slope < 0
  ) {
    
    return(
      "Early_Peaking_Decreasing"
    )
    
  }
  
  if (
    late_peak &&
    late_mean > early_mean &&
    overall_slope > 0
  ) {
    
    return(
      "Late_Increasing"
    )
    
  }
  
  if (
    early_trough &&
    late_mean > early_mean &&
    overall_slope > 0
  ) {
    
    return(
      "Early_Low_Late_Increasing"
    )
    
  }
  
  if (
    late_trough &&
    early_mean > late_mean &&
    overall_slope < 0
  ) {
    
    return(
      "Early_High_Late_Decreasing"
    )
    
  }
  
  if (
    abs(
      early_mean -
      late_mean
    ) < 0.20 &&
    abs(
      overall_slope
    ) < 0.03
  ) {
    
    return(
      "Relatively_Stable"
    )
    
  }
  
  if (
    abs(
      overall_slope
    ) < 0.03
  ) {
    
    return(
      "Non_Monotonic_or_Stable"
    )
    
  }
  
  if (
    overall_slope > 0
  ) {
    
    return(
      "Overall_Increasing"
    )
    
  }
  
  return(
    "Overall_Decreasing"
  )
}


# ============================================================
# 15. CHARACTERIZE CLUSTER CENTROIDS
# ============================================================

cluster_pattern_list <-
  lapply(
    seq_len(selected_k),
    function(k) {
      
      profile <-
        as.numeric(
          cluster_centroids[k, ]
        )
      
      peak_bin <-
        which.max(
          profile
        )
      
      trough_bin <-
        which.min(
          profile
        )
      
      early_mean <-
        mean(
          profile[
            1:4
          ],
          na.rm = TRUE
        )
      
      late_mean <-
        mean(
          profile[
            (n_bins - 3):n_bins
          ],
          na.rm = TRUE
        )
      
      data.frame(
        Cluster = k,
        Gene_Count =
          sum(
            stored_clusters == k
          ),
        Pattern =
          classify_temporal_pattern(
            profile
          ),
        Peak_Bin =
          peak_bin,
        Peak_Value =
          profile[
            peak_bin
          ],
        Trough_Bin =
          trough_bin,
        Trough_Value =
          profile[
            trough_bin
          ],
        Early_Mean =
          early_mean,
        Late_Mean =
          late_mean,
        Early_Late_Difference =
          late_mean -
          early_mean,
        Overall_Slope =
          coef(
            lm(
              profile ~
                seq_along(profile)
            )
          )[2],
        stringsAsFactors = FALSE
      )
    }
  )


cluster_pattern_summary <-
  bind_rows(
    cluster_pattern_list
  )


# ============================================================
# 16. CREATE CLUSTER CENTROID TABLE
# ============================================================

cluster_centroid_table <-
  as.data.frame(
    cluster_centroids,
    stringsAsFactors = FALSE
  )

colnames(
  cluster_centroid_table
) <-
  bin_labels

cluster_centroid_table <-
  data.frame(
    Cluster =
      seq_len(
        nrow(
          cluster_centroid_table
        )
      ),
    cluster_centroid_table,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


# ============================================================
# 17. CREATE LONG-FORM CENTROID TABLE
# ============================================================

cluster_centroid_long <-
  cluster_centroid_table %>%
  pivot_longer(
    cols =
      all_of(
        bin_labels
      ),
    names_to =
      "Bin",
    values_to =
      "Centroid_Z"
  ) %>%
  mutate(
    Bin_Number =
      as.integer(
        sub(
          "Bin_",
          "",
          Bin
        )
      ),
    Pseudotime_Bin_Center =
      bin_centers[
        Bin_Number
      ]
  ) %>%
  arrange(
    Cluster,
    Bin_Number
  )


# ============================================================
# 18. CREATE GENE PROFILE TABLE
# ============================================================

gene_zscore_table <-
  as.data.frame(
    zscore_matrix,
    stringsAsFactors = FALSE
  )

gene_zscore_table$Gene <-
  rownames(
    zscore_matrix
  )

gene_zscore_table <-
  gene_zscore_table %>%
  dplyr::select(
    Gene,
    dplyr::everything()
  )

colnames(
  gene_zscore_table
)[
  -1
] <-
  bin_labels


# ============================================================
# 19. CREATE LONG-FORM GENE PROFILE TABLE
# ============================================================

gene_zscore_long <-
  gene_zscore_table %>%
  pivot_longer(
    cols =
      all_of(
        bin_labels
      ),
    names_to =
      "Bin",
    values_to =
      "Z_Score"
  ) %>%
  mutate(
    Bin_Number =
      as.integer(
        sub(
          "Bin_",
          "",
          Bin
        )
      ),
    Pseudotime_Bin_Center =
      bin_centers[
        Bin_Number
      ]
  )


# ============================================================
# 20. ADD CLUSTER INFORMATION TO GENE PROFILES
# ============================================================

gene_zscore_long <-
  gene_zscore_long %>%
  left_join(
    cluster_assignments %>%
      dplyr::select(
        Gene,
        Cluster,
        Silhouette_Width
      ),
    by = "Gene"
  )


# ============================================================
# 21. CALCULATE GENE-LEVEL TEMPORAL FEATURES
# ============================================================

gene_temporal_features <-
  gene_zscore_long %>%
  group_by(
    Gene,
    Cluster,
    Silhouette_Width
  ) %>%
  summarise(
    
    Peak_Bin =
      which.max(
        Z_Score
      ),
    
    Peak_Z_Score =
      max(
        Z_Score,
        na.rm = TRUE
      ),
    
    Trough_Bin =
      which.min(
        Z_Score
      ),
    
    Trough_Z_Score =
      min(
        Z_Score,
        na.rm = TRUE
      ),
    
    Early_Mean =
      mean(
        Z_Score[
          Bin_Number <= 4
        ],
        na.rm = TRUE
      ),
    
    Late_Mean =
      mean(
        Z_Score[
          Bin_Number > n_bins - 4
        ],
        na.rm = TRUE
      ),
    
    Early_Late_Difference =
      Late_Mean -
      Early_Mean,
    
    Overall_Slope =
      coef(
        lm(
          Z_Score ~
            Bin_Number
        )
      )[2],
    
    .groups = "drop"
  )


# ============================================================
# 22. ADD TEMPORAL PATTERN LABELS TO GENE FEATURES
# ============================================================

gene_temporal_features <-
  gene_temporal_features %>%
  left_join(
    gene_zscore_long %>%
      group_by(
        Gene
      ) %>%
      summarise(
        Temporal_Pattern =
          classify_temporal_pattern(
            Z_Score
          ),
        .groups = "drop"
      ),
    by = "Gene"
  )


# ============================================================
# 23. CREATE CLUSTER-LEVEL SUMMARY
# ============================================================

cluster_summary_final <-
  cluster_pattern_summary %>%
  left_join(
    cluster_size_summary,
    by = "Cluster",
    suffix = c(
      "",
      "_from_assignments"
    )
  ) %>%
  dplyr::select(
    Cluster,
    Gene_Count,
    Percentage,
    Pattern,
    Peak_Bin,
    Peak_Value,
    Trough_Bin,
    Trough_Value,
    Early_Mean,
    Late_Mean,
    Early_Late_Difference,
    Overall_Slope
  )


# ============================================================
# 24. TOP GENES PER CLUSTER
# ============================================================

top_genes_final <-
  gene_cluster_correlation_table %>%
  arrange(
    Cluster,
    desc(
      Centroid_Pearson_Correlation
    )
  ) %>%
  group_by(
    Cluster
  ) %>%
  mutate(
    Rank_Within_Cluster =
      row_number()
  ) %>%
  ungroup()


top20_genes_final <-
  top_genes_final %>%
  filter(
    Rank_Within_Cluster <= 20
  )


# ============================================================
# 25. MERGE TEMPORAL FEATURES WITH GENE RANKING
# ============================================================

gene_characterization_table <-
  top_genes_final %>%
  left_join(
    gene_temporal_features %>%
      dplyr::select(
        Gene,
        Peak_Bin,
        Peak_Z_Score,
        Trough_Bin,
        Trough_Z_Score,
        Early_Mean,
        Late_Mean,
        Early_Late_Difference,
        Overall_Slope,
        Temporal_Pattern
      ),
    by = "Gene"
  ) %>%
  arrange(
    Cluster,
    Rank_Within_Cluster
  )


# ============================================================
# 26. CREATE BIN SUMMARY
# ============================================================

bin_summary <-
  data.frame(
    Bin_Number =
      seq_len(
        n_bins
      ),
    Bin =
      bin_labels,
    Lower_Bound =
      bin_breaks[
        -length(bin_breaks)
      ],
    Upper_Bound =
      bin_breaks[
        -1
      ],
    Center =
      bin_centers,
    Cell_Count =
      as.integer(
        bin_counts
      ),
    stringsAsFactors = FALSE
  )


# ============================================================
# 27. MAIN FIGURE 1
#     CLUSTER CENTROID PROFILES
# ============================================================

p_centroids <-
  ggplot(
    cluster_centroid_long,
    aes(
      x = Bin_Number,
      y = Centroid_Z,
      group = factor(Cluster),
      linetype = factor(Cluster)
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 1.8
  ) +
  
  scale_x_continuous(
    breaks =
      seq_len(
        n_bins
      )
  ) +
  
  labs(
    title =
      "Temporal profiles of Stage 17B gene clusters",
    subtitle =
      "Cluster centroids across 20 pseudotime bins",
    x =
      "Pseudotime bin",
    y =
      "Mean gene-wise Z-score",
    linetype =
      "Cluster"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    legend.position =
      "right"
  )


ggsave(
  filename =
    file.path(
      figures_main_dir,
      "Stage17C_Temporal_Cluster_Centroid_Profiles.pdf"
    ),
  plot =
    p_centroids,
  width =
    9,
  height =
    6,
  device =
    cairo_pdf
)


# ============================================================
# 28. MAIN FIGURE 2
#     CENTROID HEATMAP
# ============================================================

centroid_heatmap_data <-
  cluster_centroid_long %>%
  mutate(
    Cluster =
      factor(
        Cluster,
        levels =
          rev(
            seq_len(
              selected_k
            )
          )
      )
  )


p_centroid_heatmap <-
  ggplot(
    centroid_heatmap_data,
    aes(
      x = Bin_Number,
      y = Cluster,
      fill = Centroid_Z
    )
  ) +
  
  geom_tile() +
  
  scale_x_continuous(
    breaks =
      seq_len(
        n_bins
      )
  ) +
  
  labs(
    title =
      "Stage 17B temporal cluster centroids",
    subtitle =
      "Gene-wise Z-score profiles summarized by cluster",
    x =
      "Pseudotime bin",
    y =
      "Cluster",
    fill =
      "Centroid Z-score"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


ggsave(
  filename =
    file.path(
      figures_main_dir,
      "Stage17C_Cluster_Centroid_Heatmap.pdf"
    ),
  plot =
    p_centroid_heatmap,
  width =
    9,
  height =
    4.5,
  device =
    cairo_pdf
)


# ============================================================
# 29. SUPPLEMENTARY FIGURE
#     INDIVIDUAL GENE PROFILES BY CLUSTER
# ============================================================

p_gene_profiles <-
  ggplot(
    gene_zscore_long,
    aes(
      x = Bin_Number,
      y = Z_Score,
      group = Gene
    )
  ) +
  
  geom_line(
    linewidth = 0.25,
    alpha = 0.20
  ) +
  
  facet_wrap(
    ~ Cluster,
    ncol = 1
  ) +
  
  scale_x_continuous(
    breaks =
      seq_len(
        n_bins
      )
  ) +
  
  labs(
    title =
      "Individual temporal profiles of strong temporal genes",
    subtitle =
      "569 Stage 17 strong temporal genes grouped by Stage 17B cluster",
    x =
      "Pseudotime bin",
    y =
      "Gene-wise Z-score"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


ggsave(
  filename =
    file.path(
      figures_supp_dir,
      "Stage17C_Individual_Gene_Temporal_Profiles.pdf"
    ),
  plot =
    p_gene_profiles,
  width =
    9,
  height =
    8,
  device =
    cairo_pdf
)


# ============================================================
# 30. SUPPLEMENTARY FIGURE
#     CLUSTER SILHOUETTE DISTRIBUTION
# ============================================================

p_silhouette <-
  ggplot(
    cluster_assignments,
    aes(
      x =
        factor(
          Cluster
        ),
      y =
        Silhouette_Width
    )
  ) +
  
  geom_boxplot(
    width = 0.65
  ) +
  
  geom_jitter(
    width = 0.10,
    height = 0,
    alpha = 0.30,
    size = 0.8
  ) +
  
  labs(
    title =
      "Silhouette width distribution by temporal gene cluster",
    x =
      "Cluster",
    y =
      "Silhouette width"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


ggsave(
  filename =
    file.path(
      figures_supp_dir,
      "Stage17C_Silhouette_Width_By_Cluster.pdf"
    ),
  plot =
    p_silhouette,
  width =
    7,
  height =
    5.5,
  device =
    cairo_pdf
)


# ============================================================
# 31. SUPPLEMENTARY FIGURE
#     BIN CELL COUNTS
# ============================================================

p_bin_counts <-
  ggplot(
    bin_summary,
    aes(
      x = Bin_Number,
      y = Cell_Count
    )
  ) +
  
  geom_col() +
  
  scale_x_continuous(
    breaks =
      seq_len(
        n_bins
      )
  ) +
  
  labs(
    title =
      "Cell distribution across pseudotime bins",
    x =
      "Pseudotime bin",
    y =
      "Cell count"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


ggsave(
  filename =
    file.path(
      figures_supp_dir,
      "Stage17C_Pseudotime_Bin_Cell_Counts.pdf"
    ),
  plot =
    p_bin_counts,
  width =
    9,
  height =
    5.5,
  device =
    cairo_pdf
)


# ============================================================
# 32. EXPORT MAIN TABLES
# ============================================================

write.csv(
  cluster_summary_final,
  file =
    file.path(
      tables_main_dir,
      "Stage17C_Cluster_Temporal_Summary.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  cluster_centroid_table,
  file =
    file.path(
      tables_main_dir,
      "Stage17C_Cluster_Centroids.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  top20_genes_final,
  file =
    file.path(
      tables_main_dir,
      "Stage17C_Top20_Genes_Per_Cluster.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  gene_characterization_table,
  file =
    file.path(
      tables_main_dir,
      "Stage17C_Gene_Temporal_Characterization.csv"
    ),
  row.names =
    FALSE
)


# ============================================================
# 33. EXPORT SUPPLEMENTARY TABLES
# ============================================================

write.csv(
  cluster_centroid_long,
  file =
    file.path(
      tables_supp_dir,
      "Stage17C_Cluster_Centroids_Long.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  gene_zscore_table,
  file =
    file.path(
      tables_supp_dir,
      "Stage17C_Gene_ZScore_Profiles.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  gene_zscore_long,
  file =
    file.path(
      tables_supp_dir,
      "Stage17C_Gene_ZScore_Profiles_Long.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  gene_temporal_features,
  file =
    file.path(
      tables_supp_dir,
      "Stage17C_Gene_Temporal_Features.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  bin_summary,
  file =
    file.path(
      tables_supp_dir,
      "Stage17C_Pseudotime_Bin_Summary.csv"
    ),
  row.names =
    FALSE
)


write.csv(
  gene_cluster_correlation_table,
  file =
    file.path(
      tables_supp_dir,
      "Stage17C_Gene_Cluster_Correlation_Ranking.csv"
    ),
  row.names =
    FALSE
)


# ============================================================
# 34. VALIDATION TABLE
# ============================================================

validation_checks <-
  data.frame(
    
    Check =
      c(
        "Stage17B object exists",
        "Stage17B schema complete",
        "Stage17B upstream frozen",
        "Trajectory reconstruction disabled",
        "Root selection disabled",
        "Pseudotime recalculation disabled",
        "Expected strong temporal gene count = 569",
        "Z-score matrix has 569 genes",
        "Z-score matrix has 20 bins",
        "Mean expression matrix has 569 genes",
        "Mean expression matrix has 20 bins",
        "Cluster centroid matrix has 2 clusters",
        "Cluster centroid matrix has 20 bins",
        "All 569 genes assigned to a cluster",
        "No duplicate gene identifiers",
        "Stored clusters match k-means clusters",
        "Cluster assignments match stored clusters",
        "Centroids successfully validated",
        "Pseudotime cells = 6735",
        "Pseudotime bins = 20"
      ),
    
    Passed =
      c(
        file.exists(
          stage17b_rds
        ),
        
        length(
          missing_fields
        ) == 0,
        
        isTRUE(
          stage17b_results$upstream_frozen
        ),
        
        identical(
          stage17b_results$trajectory_reconstruction,
          FALSE
        ),
        
        identical(
          stage17b_results$root_selection,
          FALSE
        ),
        
        identical(
          stage17b_results$pseudotime_recalculation,
          FALSE
        ),
        
        length(
          strong_temporal_genes
        ) == 569,
        
        nrow(
          zscore_matrix
        ) == 569,
        
        ncol(
          zscore_matrix
        ) == 20,
        
        nrow(
          mean_expression_matrix
        ) == 569,
        
        ncol(
          mean_expression_matrix
        ) == 20,
        
        nrow(
          cluster_centroids
        ) == 2,
        
        ncol(
          cluster_centroids
        ) == 20,
        
        nrow(
          cluster_assignments
        ) == 569,
        
        anyDuplicated(
          gene_names
        ) == 0,
        
        all(
          kmeans_clusters ==
            stored_clusters
        ),
        
        all(
          cluster_assignments$Cluster ==
            stored_clusters
        ),
        
        is.finite(
          centroid_difference
        ) &&
          centroid_difference <= 1e-8,
        
        length(
          pseudotime
        ) == 6735,
        
        n_bins == 20
      ),
    
    stringsAsFactors =
      FALSE
  )


write.csv(
  validation_checks,
  file =
    file.path(
      tables_validation_dir,
      "Stage17C_Final_Validation.csv"
    ),
  row.names =
    FALSE
)


# ============================================================
# 35. VALIDATION SUMMARY
# ============================================================

validation_summary <-
  data.frame(
    
    Metric =
      c(
        "Stage17B strong temporal genes",
        "Strong temporal genes in expression",
        "Synchronized cells",
        "Pseudotime bins",
        "Genes characterized",
        "Selected K",
        "Selected mean silhouette",
        "Cluster 1 gene count",
        "Cluster 2 gene count",
        "Maximum centroid validation difference",
        "Validation checks passed",
        "Validation checks total"
      ),
    
    Value =
      c(
        length(
          strong_temporal_genes
        ),
        
        length(
          stage17b_results$strong_temporal_genes_in_expression
        ),
        
        nrow(
          synchronized_cells
        ),
        
        n_bins,
        
        nrow(
          gene_characterization_table
        ),
        
        selected_k,
        
        selected_mean_silhouette,
        
        sum(
          stored_clusters == 1
        ),
        
        sum(
          stored_clusters == 2
        ),
        
        centroid_difference,
        
        sum(
          validation_checks$Passed
        ),
        
        nrow(
          validation_checks
        )
      ),
    
    stringsAsFactors =
      FALSE
  )


write.csv(
  validation_summary,
  file =
    file.path(
      tables_validation_dir,
      "Stage17C_Validation_Summary.csv"
    ),
  row.names =
    FALSE
)


# ============================================================
# 36. FINAL VALIDATION STATUS
# ============================================================

all_validation_passed <-
  all(
    validation_checks$Passed
  )


if (!all_validation_passed) {
  
  failed_checks <-
    validation_checks$Check[
      !validation_checks$Passed
    ]
  
  stop(
    paste0(
      "Stage 17C validation failed:\n",
      paste(
        failed_checks,
        collapse = "\n"
      )
    )
  )
}


# ============================================================
# 37. BUILD CENTRAL STAGE 17C RESULT OBJECT
# ============================================================

stage17c_results <-
  list(
    
    stage =
      "17C",
    
    analysis_name =
      "Temporal Gene Program Characterization",
    
    source_stage17b_rds =
      stage17b_rds,
    
    source_stage17_rds =
      stage17b_results$source_stage17_rds,
    
    synchronized_cells =
      synchronized_cells,
    
    pseudotime =
      pseudotime,
    
    pseudotime_bins =
      pseudotime_bins,
    
    bin_breaks =
      bin_breaks,
    
    bin_centers =
      bin_centers,
    
    bin_counts =
      bin_counts,
    
    n_bins =
      n_bins,
    
    strong_temporal_genes =
      strong_temporal_genes,
    
    zscore_matrix =
      zscore_matrix,
    
    mean_expression_matrix =
      mean_expression_matrix,
    
    cluster_centroids =
      cluster_centroids,
    
    cluster_centroid_table =
      cluster_centroid_table,
    
    cluster_centroid_long =
      cluster_centroid_long,
    
    cluster_assignments =
      cluster_assignments,
    
    cluster_summary =
      cluster_summary_final,
    
    cluster_pattern_summary =
      cluster_pattern_summary,
    
    gene_cluster_correlation_table =
      gene_cluster_correlation_table,
    
    top20_genes_per_cluster =
      top20_genes_final,
    
    gene_characterization_table =
      gene_characterization_table,
    
    gene_temporal_features =
      gene_temporal_features,
    
    gene_zscore_profiles =
      gene_zscore_table,
    
    gene_zscore_profiles_long =
      gene_zscore_long,
    
    bin_summary =
      bin_summary,
    
    selected_k =
      selected_k,
    
    selected_mean_silhouette =
      selected_mean_silhouette,
    
    centroid_validation_difference =
      centroid_difference,
    
    validation_checks =
      validation_checks,
    
    validation_summary =
      validation_summary,
    
    methodology =
      list(
        
        source =
          "Frozen Stage 17B temporal gene clustering object",
        
        pseudotime_source =
          "Frozen Stage 17B synchronized_cells$Monocle3_pseudotime",
        
        pseudotime_bins =
          "20 equal-width bins",
        
        temporal_profile =
          "Mean expression per pseudotime bin followed by gene-wise Z-score normalization",
        
        clustering =
          "Stage 17B k-means clustering of gene-wise Z-scored temporal profiles",
        
        cluster_characterization =
          "Temporal characterization of frozen Stage 17B cluster centroids and gene profiles",
        
        pattern_classification =
          "Rule-based descriptive classification using peak/trough location, early-vs-late means, and overall profile slope",
        
        gene_ranking =
          "Frozen Stage 17B Pearson correlation with cluster centroid",
        
        biological_interpretation =
          "Temporal programs are interpreted as gene-expression programs along inferred pseudotime and are not evidence of cell-type conversion",
        
        trajectory_reconstruction =
          FALSE,
        
        root_selection =
          FALSE,
        
        pseudotime_recalculation =
          FALSE
      ),
    
    upstream_frozen =
      TRUE,
    
    trajectory_reconstruction =
      FALSE,
    
    root_selection =
      FALSE,
    
    pseudotime_recalculation =
      FALSE
  )


# ============================================================
# 38. SAVE CENTRAL RDS
# ============================================================

stage17c_rds <-
  file.path(
    objects_dir,
    "Stage17C_Temporal_Characterization_Results.rds"
  )

saveRDS(
  stage17c_results,
  file =
    stage17c_rds
)


# ============================================================
# 39. FINAL LOGGING
# ============================================================

write_log(
  ""
)

write_log(
  "============================================================"
)

write_log(
  "STAGE 17C COMPLETED SUCCESSFULLY"
)

write_log(
  "============================================================"
)

write_log(
  "Stage 17B strong temporal genes: ",
  length(
    strong_temporal_genes
  )
)

write_log(
  "Genes characterized: ",
  nrow(
    gene_characterization_table
  )
)

write_log(
  "Final synchronized cells: ",
  nrow(
    synchronized_cells
  )
)

write_log(
  "Pseudotime bins: ",
  n_bins
)

write_log(
  "Selected K: ",
  selected_k
)

write_log(
  "Mean silhouette: ",
  format(
    selected_mean_silhouette,
    digits = 8
  )
)

write_log(
  "Cluster 1 genes: ",
  sum(
    stored_clusters == 1
  )
)

write_log(
  "Cluster 2 genes: ",
  sum(
    stored_clusters == 2
  )
)

write_log(
  "Cluster 1 temporal pattern: ",
  cluster_pattern_summary$Pattern[
    cluster_pattern_summary$Cluster == 1
  ]
)

write_log(
  "Cluster 2 temporal pattern: ",
  cluster_pattern_summary$Pattern[
    cluster_pattern_summary$Cluster == 2
  ]
)

write_log(
  "Maximum centroid validation difference: ",
  format(
    centroid_difference,
    scientific = TRUE
  )
)

write_log(
  "Validation checks passed: ",
  sum(
    validation_checks$Passed
  ),
  "/",
  nrow(
    validation_checks
  )
)

write_log(
  ""
)

write_log(
  "Central RDS:"
)

write_log(
  stage17c_rds
)

write_log(
  ""
)

write_log(
  "Completion log:"
)

write_log(
  stage17c_log
)

write_log(
  ""
)

write_log(
  "Upstream Stage 17B remained frozen."
)

write_log(
  "No trajectory reconstruction was performed."
)

write_log(
  "No root selection was performed."
)

write_log(
  "No pseudotime recalculation was performed."
)

write_log(
  "============================================================"
)


# ============================================================
# 40. FINAL CONSOLE SUMMARY
# ============================================================

cat(
  "\n\n"
)

cat(
  "STAGE 17C COMPLETED SUCCESSFULLY\n"
)

cat(
  "--------------------------------------------\n"
)

cat(
  "Strong temporal genes: ",
  length(
    strong_temporal_genes
  ),
  "\n"
)

cat(
  "Genes characterized: ",
  nrow(
    gene_characterization_table
  ),
  "\n"
)

cat(
  "Synchronized cells: ",
  nrow(
    synchronized_cells
  ),
  "\n"
)

cat(
  "Pseudotime bins: ",
  n_bins,
  "\n"
)

cat(
  "Selected K: ",
  selected_k,
  "\n"
)

cat(
  "Mean silhouette: ",
  format(
    selected_mean_silhouette,
    digits = 8
  ),
  "\n"
)

cat(
  "Cluster 1 genes: ",
  sum(
    stored_clusters == 1
  ),
  "\n"
)

cat(
  "Cluster 2 genes: ",
  sum(
    stored_clusters == 2
  ),
  "\n"
)

cat(
  "Cluster 1 pattern: ",
  cluster_pattern_summary$Pattern[
    cluster_pattern_summary$Cluster == 1
  ],
  "\n"
)

cat(
  "Cluster 2 pattern: ",
  cluster_pattern_summary$Pattern[
    cluster_pattern_summary$Cluster == 2
  ],
  "\n"
)

cat(
  "Validation: ",
  sum(
    validation_checks$Passed
  ),
  "/",
  nrow(
    validation_checks
  ),
  " passed\n"
)

cat(
  "\nCentral RDS:\n",
  stage17c_rds,
  "\n"
)

cat(
  "\nCompletion log:\n",
  stage17c_log,
  "\n"
)

cat(
  "\nStage 17B remained frozen.\n"
)

cat(
  "No trajectory reconstruction, root selection, or pseudotime recalculation was performed.\n"
)

cat(
  "\n"
)

