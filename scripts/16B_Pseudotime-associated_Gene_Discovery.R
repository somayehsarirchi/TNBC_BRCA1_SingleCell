# ==============================================================================
# STAGE 16B
# PSEUDOTIME-ASSOCIATED GENE DISCOVERY
#
# Purpose:
#   Identify genes whose expression is associated with the frozen Monocle3
#   trajectory using graph-based gene testing.
#
# IMPORTANT:
#   Stage 15 is FROZEN.
#   Stage 16A is FROZEN.
#
#   This stage does NOT:
#     - reconstruct the trajectory
#     - change partitions
#     - select a new root
#     - recalculate pseudotime
#     - modify the frozen Stage 15 CDS
#
#   Analysis population:
#     Frozen target partition = 2
#     Finite pseudotime cells = 6735
#
# Main computational method:
#   Monocle3 graph_test()
#
# Primary significance threshold:
#   q_value < 0.05
#
# Exploratory stronger association set:
#   Moran's I >= 0.10
#
# IMPORTANT FOR REPOSITORY USE:
#   If validated Stage 16B outputs already exist, this script should NOT
#   be rerun merely for repository organization.
#
# ==============================================================================


# ==============================================================================
# 1. Load Required Libraries
# ==============================================================================

suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
})


# ==============================================================================
# 2. Define Project Paths
# ==============================================================================

project_dir <- "C:/Users/asus/Desktop/TNBC_BRCA1_SingleCell"

objects_dir <- file.path(
  project_dir,
  "objects"
)

stage15_cds_dir <- file.path(
  objects_dir,
  "TNBC_Monocle3_CDS_Final"
)

stage16_dir <- file.path(
  project_dir,
  "results",
  "16_Pseudotime_Gene_Dynamics"
)

stage16_main_dir <- file.path(
  stage16_dir,
  "tables",
  "main"
)

stage16_supp_dir <- file.path(
  stage16_dir,
  "tables",
  "supplementary"
)

stage16_validation_dir <- file.path(
  stage16_dir,
  "tables",
  "validation"
)

stage16_log_dir <- file.path(
  stage16_dir,
  "logs"
)


# ==============================================================================
# 3. Create Stage 16 Output Directories
# ==============================================================================

dir.create(
  stage16_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  stage16_supp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  stage16_validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  stage16_log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 4. Validate Frozen Stage 15 CDS
# ==============================================================================

if (!dir.exists(stage15_cds_dir)) {
  
  stop(
    paste0(
      "Frozen Stage 15 CDS directory not found:\n",
      stage15_cds_dir
    )
  )
}


# ==============================================================================
# 5. Load Frozen Stage 15 Monocle3 CDS
# ==============================================================================

cds <- load_monocle_objects(
  stage15_cds_dir
)


if (!inherits(cds, "cell_data_set")) {
  
  stop(
    "Loaded object is not a valid Monocle3 cell_data_set."
  )
}


# ==============================================================================
# 6. Validate Basic CDS Dimensions
# ==============================================================================

if (nrow(cds) == 0 || ncol(cds) == 0) {
  
  stop(
    "Frozen Stage 15 CDS contains zero genes or zero cells."
  )
}


# ==============================================================================
# 7. Validate Required Metadata
# ==============================================================================

required_metadata <- c(
  "CNV_Status"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(colData(cds))
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "Required metadata missing from frozen Stage 15 CDS: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


# ==============================================================================
# 8. Retrieve Frozen Pseudotime
#
# IMPORTANT:
#   pseudotime(cds) retrieves the pseudotime already stored in the frozen
#   Monocle3 object.
#
#   No trajectory reconstruction or pseudotime recalculation is performed.
# ==============================================================================

frozen_pseudotime <- as.numeric(
  pseudotime(cds)
)

if (length(frozen_pseudotime) != ncol(cds)) {
  
  stop(
    paste0(
      "Frozen pseudotime length does not match CDS cell count. ",
      "Expected ",
      ncol(cds),
      " values but found ",
      length(frozen_pseudotime),
      "."
    )
  )
}


# ==============================================================================
# 9. Restore Explicit Pseudotime Metadata in Memory if Necessary
#
# This does NOT recalculate pseudotime.
#
# It only exposes the pseudotime already stored internally by Monocle3.
# The frozen CDS on disk is not overwritten.
# ==============================================================================

if (!"Monocle3_pseudotime" %in% colnames(colData(cds))) {
  
  cds$Monocle3_pseudotime <- frozen_pseudotime
  
} else {
  
  metadata_pseudotime <- as.numeric(
    colData(cds)$Monocle3_pseudotime
  )
  
  if (
    !isTRUE(
      all.equal(
        metadata_pseudotime,
        frozen_pseudotime,
        check.attributes = FALSE
      )
    )
  ) {
    
    stop(
      paste0(
        "Existing Monocle3_pseudotime metadata does not match ",
        "the pseudotime stored internally in the frozen CDS."
      )
    )
  }
}


# ==============================================================================
# 10. Retrieve Frozen Partitions
# ==============================================================================

cds_partition <- as.character(
  partitions(cds)
)

if (length(cds_partition) != ncol(cds)) {
  
  stop(
    "Partition vector length does not match CDS cell count."
  )
}


# ==============================================================================
# 11. Define Frozen Analysis Population
# ==============================================================================

target_partition <- "2"

expected_analysis_cells <- 6735

analysis_cells <- colnames(cds)[
  cds_partition == target_partition &
    is.finite(cds$Monocle3_pseudotime)
]

n_analysis_cells <- length(
  analysis_cells
)


# ==============================================================================
# 12. Strict Frozen Population Validation
# ==============================================================================

if (n_analysis_cells != expected_analysis_cells) {
  
  stop(
    paste0(
      "Frozen Stage 16A validation failed.\n",
      "Expected ",
      expected_analysis_cells,
      " cells in partition ",
      target_partition,
      ", found ",
      n_analysis_cells,
      "."
    )
  )
}


# ==============================================================================
# 13. Prepare Analysis Metadata
# ==============================================================================

analysis_metadata <- as.data.frame(
  colData(cds)[analysis_cells, , drop = FALSE]
)

analysis_metadata$Pseudotime <- cds$Monocle3_pseudotime[
  match(
    analysis_cells,
    colnames(cds)
  )
]

analysis_metadata$Partition <- cds_partition[
  match(
    analysis_cells,
    colnames(cds)
  )
]


# ==============================================================================
# 14. Validate Frozen Pseudotime
# ==============================================================================

if (
  !all(
    is.finite(
      analysis_metadata$Pseudotime
    )
  )
) {
  
  stop(
    "Non-finite pseudotime detected inside the frozen analysis population."
  )
}


# ==============================================================================
# 15. Validate Frozen Root
# ==============================================================================

expected_root <- "Y_48"

umap_graph <- principal_graph(cds)[["UMAP"]]

if (is.null(umap_graph)) {
  
  stop(
    "UMAP principal graph is missing from the frozen Stage 15 CDS."
  )
}

graph_nodes <- igraph::V(
  umap_graph
)$name

if (!(expected_root %in% graph_nodes)) {
  
  stop(
    paste0(
      "Frozen root node not found in principal graph: ",
      expected_root
    )
  )
}


# ==============================================================================
# 16. Prepare Frozen Input Snapshot
# ==============================================================================

input_snapshot <- data.frame(
  
  Parameter = c(
    "Stage",
    "Frozen_Stage15",
    "Frozen_Stage16A",
    "Total_CDS_Genes",
    "Total_CDS_Cells",
    "Target_Partition",
    "Analysis_Cells",
    "Root_Node",
    "Normal_like_Cells",
    "Tumor_like_Cells"
  ),
  
  Value = c(
    "Stage 16B",
    "YES",
    "YES",
    nrow(cds),
    ncol(cds),
    target_partition,
    n_analysis_cells,
    expected_root,
    sum(
      analysis_metadata$CNV_Status == "Normal_like",
      na.rm = TRUE
    ),
    sum(
      analysis_metadata$CNV_Status == "Tumor_like",
      na.rm = TRUE
    )
  ),
  
  stringsAsFactors = FALSE
)


# ==============================================================================
# 17. Save Input Validation Table
# ==============================================================================

write.csv(
  input_snapshot,
  file.path(
    stage16_validation_dir,
    "Stage16B_Input_Snapshot.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. Save Analysis Cell IDs
# ==============================================================================

write.csv(
  data.frame(
    Cell = analysis_cells,
    stringsAsFactors = FALSE
  ),
  file.path(
    stage16_supp_dir,
    "Stage16B_Analysis_Cells.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 19. Create Separate Analysis CDS
#
# IMPORTANT:
#   The original frozen Stage 15 CDS is not modified.
# ==============================================================================

cds_analysis <- cds[
  ,
  analysis_cells
]


# ==============================================================================
# 20. Validate Analysis CDS
# ==============================================================================

if (
  ncol(cds_analysis) != expected_analysis_cells
) {
  
  stop(
    paste0(
      "Analysis CDS does not contain exactly ",
      expected_analysis_cells,
      " frozen cells."
    )
  )
}

analysis_pt <- as.numeric(
  colData(cds_analysis)$Monocle3_pseudotime
)

if (!all(is.finite(analysis_pt))) {
  
  stop(
    "Analysis CDS contains non-finite pseudotime values."
  )
}


# ==============================================================================
# 21. Run Monocle3 graph_test()
#
# MAIN COMPUTATIONAL STEP
#
# This tests genes for spatially structured expression on the existing
# principal graph.
#
# IMPORTANT:
#   graph_test() does not reconstruct the trajectory.
# ==============================================================================

cat("\n============================================================\n")
cat("STAGE 16B — MONOCLE3 graph_test()\n")
cat("============================================================\n\n")

cat(
  "Frozen partition:",
  target_partition,
  "\n"
)

cat(
  "Frozen root:",
  expected_root,
  "\n"
)

cat(
  "Analysis cells:",
  n_analysis_cells,
  "\n\n"
)

cat(
  "Running graph_test(). This may take some time...\n\n"
)


graph_test_result <- graph_test(
  cds_analysis,
  neighbor_graph = "principal_graph",
  cores = 1
)


# ==============================================================================
# 22. Validate graph_test() Output
# ==============================================================================

if (
  is.null(graph_test_result) ||
  nrow(graph_test_result) == 0
) {
  
  stop(
    "graph_test() returned no results."
  )
}


graph_test_df <- as.data.frame(
  graph_test_result
)

graph_test_df$gene_id <- rownames(
  graph_test_df
)


# ==============================================================================
# 23. Validate Required graph_test Columns
# ==============================================================================

required_graph_test_columns <- c(
  "q_value",
  "morans_I"
)

missing_graph_test_columns <- setdiff(
  required_graph_test_columns,
  colnames(graph_test_df)
)

if (length(missing_graph_test_columns) > 0) {
  
  stop(
    paste0(
      "Required graph_test columns missing: ",
      paste(
        missing_graph_test_columns,
        collapse = ", "
      )
    )
  )
}


graph_test_df$q_value <- as.numeric(
  graph_test_df$q_value
)

graph_test_df$morans_I <- as.numeric(
  graph_test_df$morans_I
)


# ==============================================================================
# 24. Save RAW graph_test Results
#
# IMPORTANT:
#   RAW results are preserved before any filtering.
# ==============================================================================

write.csv(
  graph_test_df,
  file.path(
    stage16_supp_dir,
    "01_graph_test_RAW.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 25. Identify Significant Trajectory-associated Genes
#
# Primary threshold:
#   q_value < 0.05
#
# Ranking:
#   Moran's I
# ==============================================================================

trajectory_genes <- graph_test_df %>%
  
  filter(
    !is.na(q_value),
    q_value < 0.05,
    !is.na(morans_I)
  ) %>%
  
  arrange(
    desc(morans_I),
    q_value
  )


# ==============================================================================
# 26. Identify Strong Trajectory-associated Genes
#
# Exploratory threshold:
#   q_value < 0.05
#   Moran's I >= 0.10
#
# This is a ranking/filtering threshold and not a biological cutoff.
# ==============================================================================

strong_trajectory_genes <- trajectory_genes %>%
  
  filter(
    morans_I >= 0.10
  ) %>%
  
  arrange(
    desc(morans_I),
    q_value
  )


# ==============================================================================
# 27. Top 100 Trajectory-associated Genes
# ==============================================================================

top100_trajectory_genes <- trajectory_genes %>%
  
  slice_head(
    n = 100
  )


# ==============================================================================
# 28. Top 50 Trajectory-associated Genes
# ==============================================================================

top50_trajectory_genes <- trajectory_genes %>%
  
  slice_head(
    n = 50
  )


# ==============================================================================
# 29. Save Main Gene-level Results
# ==============================================================================

write.csv(
  trajectory_genes,
  file.path(
    stage16_main_dir,
    "02_Significant_Trajectory_Genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  strong_trajectory_genes,
  file.path(
    stage16_main_dir,
    "03_Strong_Trajectory_Genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  top100_trajectory_genes,
  file.path(
    stage16_main_dir,
    "04_Top100_Trajectory_Genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  top50_trajectory_genes,
  file.path(
    stage16_main_dir,
    "05_Top50_Trajectory_Genes.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 30. Stage 16B Summary Statistics
# ==============================================================================

trajectory_summary <- data.frame(
  
  Parameter = c(
    "Total_genes_tested",
    "Significant_qvalue_lt_0.05",
    "Strong_MoransI_ge_0.10",
    "Top100_available",
    "Top50_available"
  ),
  
  Value = c(
    nrow(graph_test_df),
    nrow(trajectory_genes),
    nrow(strong_trajectory_genes),
    min(
      100,
      nrow(trajectory_genes)
    ),
    min(
      50,
      nrow(trajectory_genes)
    )
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  trajectory_summary,
  file.path(
    stage16_main_dir,
    "06_Stage16B_Gene_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 31. Validate Generated Outputs
# ==============================================================================

expected_outputs <- c(
  
  file.path(
    stage16_validation_dir,
    "Stage16B_Input_Snapshot.csv"
  ),
  
  file.path(
    stage16_supp_dir,
    "Stage16B_Analysis_Cells.csv"
  ),
  
  file.path(
    stage16_supp_dir,
    "01_graph_test_RAW.csv"
  ),
  
  file.path(
    stage16_main_dir,
    "02_Significant_Trajectory_Genes.csv"
  ),
  
  file.path(
    stage16_main_dir,
    "03_Strong_Trajectory_Genes.csv"
  ),
  
  file.path(
    stage16_main_dir,
    "04_Top100_Trajectory_Genes.csv"
  ),
  
  file.path(
    stage16_main_dir,
    "05_Top50_Trajectory_Genes.csv"
  ),
  
  file.path(
    stage16_main_dir,
    "06_Stage16B_Gene_Summary.csv"
  )
)


output_exists <- file.exists(
  expected_outputs
)

if (!all(output_exists)) {
  
  missing_outputs <- expected_outputs[
    !output_exists
  ]
  
  stop(
    paste0(
      "Stage 16B output validation failed. Missing files:\n",
      paste(
        missing_outputs,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# 32. Completion Log
# ==============================================================================

stage16b_log <- file.path(
  stage16_log_dir,
  "Stage16B_Completion.log"
)

log_lines <- c(
  
  "============================================================",
  "TNBC Single-Cell Analysis Pipeline",
  "Stage 16B: Pseudotime-associated Gene Discovery",
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
    "Frozen Stage 15 CDS:",
    stage15_cds_dir
  ),
  
  "",
  
  paste(
    "Total genes:",
    nrow(cds)
  ),
  
  paste(
    "Total cells:",
    ncol(cds)
  ),
  
  paste(
    "Frozen target partition:",
    target_partition
  ),
  
  paste(
    "Frozen root:",
    expected_root
  ),
  
  paste(
    "Analysis cells:",
    n_analysis_cells
  ),
  
  paste(
    "Genes tested:",
    nrow(graph_test_df)
  ),
  
  paste(
    "Significant genes q < 0.05:",
    nrow(trajectory_genes)
  ),
  
  paste(
    "Strong genes Moran's I >= 0.10:",
    nrow(strong_trajectory_genes)
  ),
  
  "",
  
  "Trajectory reconstruction: NOT PERFORMED",
  "Partition modification: NOT PERFORMED",
  "Root selection: NOT PERFORMED",
  "Pseudotime recalculation: NOT PERFORMED",
  "Frozen Stage 15 CDS: NOT MODIFIED",
  
  "",
  
  "Stage 15: FROZEN",
  "Stage 16A: FROZEN",
  "Stage 16B: COMPLETED",
  
  "",
  
  "STATUS: STAGE 16B COMPLETED SUCCESSFULLY",
  
  "============================================================"
)


writeLines(
  log_lines,
  con = stage16b_log
)


# ==============================================================================
# 33. Final Console Summary
# ==============================================================================

cat("\n============================================================\n")
cat("STAGE 16B COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat(
  "Frozen partition:",
  target_partition,
  "\n"
)

cat(
  "Frozen root:",
  expected_root,
  "\n"
)

cat(
  "Analysis cells:",
  n_analysis_cells,
  "\n"
)

cat(
  "Genes tested:",
  nrow(graph_test_df),
  "\n"
)

cat(
  "Significant genes (q < 0.05):",
  nrow(trajectory_genes),
  "\n"
)

cat(
  "Strong genes (Moran's I >= 0.10):",
  nrow(strong_trajectory_genes),
  "\n"
)

cat(
  "\nStage 16B output directory:\n",
  stage16_dir,
  "\n"
)

cat("\n============================================================\n")

