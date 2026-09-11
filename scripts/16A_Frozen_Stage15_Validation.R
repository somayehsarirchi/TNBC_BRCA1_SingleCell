# ==============================================================================
# STAGE 16A
# VALIDATION OF FROZEN STAGE 15 INPUT
#
# Purpose:
#   Validate the frozen Monocle3 trajectory object before downstream
#   pseudotime-associated gene dynamics analysis.
#
# IMPORTANT:
#   Stage 15 is frozen.
#
#   This stage does NOT:
#     - reconstruct the trajectory
#     - recluster cells
#     - learn a new principal graph
#     - select a new root
#     - recalculate pseudotime
#
#   If the explicit Monocle3_pseudotime metadata column is absent, it is
#   restored from the pseudotime already stored internally in the frozen
#   Monocle3 CDS. This is metadata normalization only and does not alter
#   the trajectory or pseudotime calculation.
#
# Output:
#   Validation tables and analysis-cell metadata are written to the
#   Stage 16 repository structure.
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

project_dir <- "YOUR_PROJECT_DIRECTORY"

objects_dir <- file.path(
  project_dir,
  "objects")

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
# 4. Validate Frozen Stage 15 CDS Location
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


# ==============================================================================
# 6. Basic CDS Integrity Validation
# ==============================================================================

if (!inherits(cds, "cell_data_set")) {
  
  stop(
    "Loaded object is not a valid Monocle3 cell_data_set."
  )
}

if (nrow(cds) == 0) {
  
  stop(
    "Frozen Stage 15 CDS contains zero genes."
  )
}

if (ncol(cds) == 0) {
  
  stop(
    "Frozen Stage 15 CDS contains zero cells."
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
# 8. Retrieve Frozen Monocle3 Pseudotime
#
# IMPORTANT:
#   pseudotime(cds) retrieves the pseudotime already stored in the
#   Monocle3 object. No trajectory reconstruction or pseudotime
#   recalculation is performed.
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
# 9. Restore Explicit Pseudotime Metadata if Necessary
#
# This does NOT recalculate pseudotime.
#
# It only exposes the already-stored Monocle3 pseudotime as an explicit
# colData column for downstream reproducible analysis.
# ==============================================================================

if (!"Monocle3_pseudotime" %in% colnames(colData(cds))) {
  
  cds$Monocle3_pseudotime <- frozen_pseudotime
  
} else {
  
  metadata_pseudotime <- as.numeric(
    colData(cds)$Monocle3_pseudotime
  )
  
  if (length(metadata_pseudotime) != length(frozen_pseudotime)) {
    
    stop(
      "Existing Monocle3_pseudotime metadata has an invalid length."
    )
  }
  
  pseudotime_match <- all.equal(
    metadata_pseudotime,
    frozen_pseudotime,
    check.attributes = FALSE
  )
  
  if (!isTRUE(pseudotime_match)) {
    
    stop(
      paste0(
        "Existing Monocle3_pseudotime metadata does not match ",
        "the pseudotime stored internally in the frozen Monocle3 CDS."
      )
    )
  }
}


# ==============================================================================
# 10. Retrieve Pseudotime from Explicit Metadata
# ==============================================================================

cds_pseudotime <- as.numeric(
  colData(cds)$Monocle3_pseudotime
)

if (all(is.na(cds_pseudotime))) {
  
  stop(
    "Monocle3 pseudotime is completely missing."
  )
}


# ==============================================================================
# 11. Validate Frozen Partitions
# ==============================================================================

cds_partition <- as.character(
  partitions(cds)
)

if (length(cds_partition) != ncol(cds)) {
  
  stop(
    paste0(
      "Partition vector length does not match CDS cell count. ",
      "Expected ",
      ncol(cds),
      " values but found ",
      length(cds_partition),
      "."
    )
  )
}


# ==============================================================================
# 12. Define Frozen Analysis Population
#
# Stage 15 frozen parameters:
#   Selected partition = 2
#   Expected finite pseudotime cells = 6735
# ==============================================================================

target_partition <- "2"

expected_analysis_cells <- 6735

analysis_cells <- colnames(cds)[
  cds_partition == target_partition &
    is.finite(cds_pseudotime)
]


# ==============================================================================
# 13. Validate Frozen Analysis Population
# ==============================================================================

n_analysis_cells <- length(
  analysis_cells
)

if (n_analysis_cells != expected_analysis_cells) {
  
  stop(
    paste0(
      "Frozen Stage 15 validation failed.\n",
      "Expected ",
      expected_analysis_cells,
      " finite cells in partition ",
      target_partition,
      ", found ",
      n_analysis_cells,
      "."
    )
  )
}


# ==============================================================================
# 14. Prepare Analysis Metadata
# ==============================================================================

analysis_metadata <- as.data.frame(
  colData(cds)[analysis_cells, , drop = FALSE]
)

analysis_metadata$Pseudotime <- cds_pseudotime[
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
# 15. Validate CNV Composition
# ==============================================================================

analysis_cnv_counts <- table(
  factor(
    analysis_metadata$CNV_Status,
    levels = c(
      "Normal_like",
      "Tumor_like"
    )
  )
)

analysis_normal_like <- unname(
  analysis_cnv_counts["Normal_like"]
)

analysis_tumor_like <- unname(
  analysis_cnv_counts["Tumor_like"]
)

if (is.na(analysis_normal_like)) {
  analysis_normal_like <- 0
}

if (is.na(analysis_tumor_like)) {
  analysis_tumor_like <- 0
}


# ==============================================================================
# 16. Validate Frozen Root Node
#
# Stage 15 frozen root:
#   Y_48
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
      "Expected frozen root node ",
      expected_root,
      " was not found in the principal graph."
    )
  )
}


# ==============================================================================
# 17. Retrieve Cell-to-Graph Projection
# ==============================================================================

principal_graph_aux_umap <- principal_graph_aux(
  cds
)[["UMAP"]]

if (is.null(principal_graph_aux_umap)) {
  
  stop(
    "UMAP principal graph auxiliary data are missing."
  )
}

if (
  is.null(
    principal_graph_aux_umap$pr_graph_cell_proj_closest_vertex
  )
) {
  
  stop(
    "Cell-to-principal-graph projection information is missing."
  )
}


# ==============================================================================
# 18. Map Cells to Actual Principal Graph Vertex Names
# ==============================================================================

pr_graph_proj <- principal_graph_aux_umap$pr_graph_cell_proj_closest_vertex

pr_graph_proj_df <- as.data.frame(
  pr_graph_proj
)

rownames(pr_graph_proj_df) <- colnames(cds)

colnames(pr_graph_proj_df) <- "closest_vertex_index"

real_vertex_names <- igraph::V(
  umap_graph
)$name

closest_vertex_name <- real_vertex_names[
  pr_graph_proj_df$closest_vertex_index
]

if (any(is.na(closest_vertex_name))) {
  
  stop(
    "Some cells could not be mapped to principal graph vertices."
  )
}


# ==============================================================================
# 19. Identify Cells Assigned to Frozen Root
# ==============================================================================

root_cells <- rownames(
  pr_graph_proj_df[
    closest_vertex_name == expected_root,
    ,
    drop = FALSE
  ]
)

if (length(root_cells) == 0) {
  
  stop(
    paste0(
      "No cells were mapped to the frozen root node: ",
      expected_root
    )
  )
}


# ==============================================================================
# 20. Validate Frozen Root Composition
#
# Stage 15 frozen result:
#   Root cells = 59
#   Normal_like = 59
#   Tumor_like = 0
# ==============================================================================

expected_root_cell_count <- 59

expected_root_normal_like <- 59

expected_root_tumor_like <- 0

root_cnv_counts <- table(
  factor(
    colData(cds)[root_cells, "CNV_Status"],
    levels = c(
      "Normal_like",
      "Tumor_like"
    )
  )
)

root_normal_like <- unname(
  root_cnv_counts["Normal_like"]
)

root_tumor_like <- unname(
  root_cnv_counts["Tumor_like"]
)

if (is.na(root_normal_like)) {
  root_normal_like <- 0
}

if (is.na(root_tumor_like)) {
  root_tumor_like <- 0
}

if (length(root_cells) != expected_root_cell_count) {
  
  stop(
    paste0(
      "Frozen root-cell validation failed. ",
      "Expected ",
      expected_root_cell_count,
      " cells at root ",
      expected_root,
      ", found ",
      length(root_cells),
      "."
    )
  )
}

if (root_normal_like != expected_root_normal_like) {
  
  stop(
    paste0(
      "Frozen root CNV validation failed. ",
      "Expected ",
      expected_root_normal_like,
      " Normal_like root cells, found ",
      root_normal_like,
      "."
    )
  )
}

if (root_tumor_like != expected_root_tumor_like) {
  
  stop(
    paste0(
      "Frozen root CNV validation failed. ",
      "Expected ",
      expected_root_tumor_like,
      " Tumor_like root cells, found ",
      root_tumor_like,
      "."
    )
  )
}


# ==============================================================================
# 21. Frozen Input Summary
# ==============================================================================

stage16_input_summary <- data.frame(
  
  Parameter = c(
    "Stage",
    "Frozen_Stage15",
    "Total_CDS_Cells",
    "Total_CDS_Genes",
    "Target_Partition",
    "Finite_Pseudotime_Cells",
    "Root_Node",
    "Root_Cell_Count",
    "Root_Normal_like_Cells",
    "Root_Tumor_like_Cells",
    "Analysis_Normal_like_Cells",
    "Analysis_Tumor_like_Cells"
  ),
  
  Value = c(
    "Stage 16A",
    "YES",
    ncol(cds),
    nrow(cds),
    target_partition,
    n_analysis_cells,
    expected_root,
    length(root_cells),
    root_normal_like,
    root_tumor_like,
    analysis_normal_like,
    analysis_tumor_like
  ),
  
  stringsAsFactors = FALSE
)


# ==============================================================================
# 22. Save Main Validation Table
# ==============================================================================

stage16_input_summary_file <- file.path(
  stage16_main_dir,
  "Stage16A_Frozen_Input_Validation.csv"
)

write.csv(
  stage16_input_summary,
  stage16_input_summary_file,
  row.names = FALSE
)


# ==============================================================================
# 23. Save Analysis Cell IDs
# ==============================================================================

analysis_cells_file <- file.path(
  stage16_supp_dir,
  "Stage16A_Analysis_Cells.csv"
)

write.csv(
  data.frame(
    Cell = analysis_cells,
    stringsAsFactors = FALSE
  ),
  analysis_cells_file,
  row.names = FALSE
)


# ==============================================================================
# 24. Save Analysis Cell Metadata
# ==============================================================================

analysis_metadata_file <- file.path(
  stage16_supp_dir,
  "Stage16A_Analysis_Cell_Metadata.csv"
)

write.csv(
  analysis_metadata,
  analysis_metadata_file,
  row.names = TRUE
)


# ==============================================================================
# 25. Validate Generated Outputs
# ==============================================================================

expected_outputs <- c(
  stage16_input_summary_file,
  analysis_cells_file,
  analysis_metadata_file
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
      "Stage 16A output validation failed. Missing files:\n",
      paste(
        missing_outputs,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# 26. Stage 16A Completion Log
# ==============================================================================

stage16a_log <- file.path(
  stage16_log_dir,
  "Stage16A_Completion.log"
)

stage16a_log_lines <- c(
  "============================================================",
  "TNBC Single-Cell Analysis Pipeline",
  "Stage 16A: Frozen Stage 15 Input Validation",
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
  "",
  paste(
    "Total cells:",
    ncol(cds)
  ),
  "",
  paste(
    "Frozen target partition:",
    target_partition
  ),
  "",
  paste(
    "Finite pseudotime cells:",
    n_analysis_cells
  ),
  "",
  paste(
    "Frozen root node:",
    expected_root
  ),
  "",
  paste(
    "Root cells:",
    length(root_cells)
  ),
  "",
  paste(
    "Root Normal-like cells:",
    root_normal_like
  ),
  "",
  paste(
    "Root Tumor-like cells:",
    root_tumor_like
  ),
  "",
  paste(
    "Analysis Normal-like cells:",
    analysis_normal_like
  ),
  "",
  paste(
    "Analysis Tumor-like cells:",
    analysis_tumor_like
  ),
  "",
  paste(
    "Monocle3_pseudotime metadata present:",
    "YES"
  ),
  "",
  "Trajectory reconstruction: NOT PERFORMED",
  "Root selection: NOT PERFORMED",
  "Pseudotime recalculation: NOT PERFORMED",
  "Stage 15 remained frozen.",
  "",
  "STATUS: STAGE 16A VALIDATED SUCCESSFULLY",
  "============================================================"
)

writeLines(
  stage16a_log_lines,
  con = stage16a_log
)


# ==============================================================================
# 27. Final Status
# ==============================================================================

message(
  "Stage 16A completed successfully. ",
  "Frozen Stage 15 input validated; no trajectory reconstruction ",
  "or pseudotime recalculation was performed."
)