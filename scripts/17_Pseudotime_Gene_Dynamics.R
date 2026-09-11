# ============================================================
# STAGE 17
# PSEUDOTIME-DEPENDENT TEMPORAL GENE DYNAMICS
# ============================================================
#
# Purpose:
#   Synchronize the frozen Stage 15 Monocle3 pseudotime with
#   the exact 6735 analysis cells defined by Stage 16B and
#   quantify gene-expression dependence on pseudotime using
#   Spearman correlation.
#
# IMPORTANT:
#   - Stage 15 trajectory is FROZEN.
#   - No trajectory reconstruction is performed.
#   - No root selection is performed.
#   - No pseudotime recalculation is performed.
#   - Stage16B_Analysis_Cells.csv contains ONLY cell IDs.
#   - Authoritative pseudotime is extracted from pseudotime(cds).
#
# Output architecture:
#
# results/
# └── 17_Temporal_Gene_Dynamics/
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
# Central object:
#   objects/Stage17_Pseudotime_Gene_Dynamics_Results.rds
#
# ============================================================


# ------------------------------------------------------------
# 01. INITIALIZATION
# ------------------------------------------------------------

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  warn = 1
)

cat("\n")
cat("============================================================\n")
cat("TNBC Single-Cell Analysis Pipeline\n")
cat("Stage 17: Pseudotime-Dependent Temporal Gene Dynamics\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 02. PROJECT DIRECTORIES
# ------------------------------------------------------------

project_dir <- "YOUR_PROJECT_DIRECTORY"

results_dir <- file.path(
  project_dir,
  "results",
  "17_Temporal_Gene_Dynamics"
)

figures_main_dir <- file.path(
  results_dir,
  "figures",
  "main"
)

figures_supp_dir <- file.path(
  results_dir,
  "figures",
  "supplementary"
)

tables_main_dir <- file.path(
  results_dir,
  "tables",
  "main"
)

tables_supp_dir <- file.path(
  results_dir,
  "tables",
  "supplementary"
)

tables_validation_dir <- file.path(
  results_dir,
  "tables",
  "validation"
)

logs_dir <- file.path(
  results_dir,
  "logs"
)

objects_dir <- file.path(
  project_dir,
  "objects"
)

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


# ------------------------------------------------------------
# 03. INPUT PATHS
# ------------------------------------------------------------

# Frozen Stage 15 Monocle3 CDS.
#
# This is the actual frozen CDS confirmed during Stage 16A.
# It is intentionally read from the original source pipeline.
#
# IMPORTANT:
# This file/directory is READ ONLY for Stage 17.
#
frozen_stage15_cds_dir <- file.path(
  project_dir,
  "objects",
  "TNBC_Monocle3_CDS_Final"
)

# Stage 16B exact analysis-cell list.
stage16b_cells_file <- file.path(
  project_dir,
  "results",
  "16_Pseudotime_Gene_Dynamics",
  "tables",
  "supplementary",
  "Stage16B_Analysis_Cells.csv"
)

# Stage 16B strong trajectory-associated genes.
stage16b_strong_genes_file <- file.path(
  project_dir,
  "results",
  "16_Pseudotime_Gene_Dynamics",
  "tables",
  "main",
  "03_Strong_Trajectory_Genes.csv"
)

# Seurat object used for expression analysis.
# The canonical CopyKAT epithelial object is used.

seurat_object_file <- file.path(
  project_dir,
  "objects",
  "TNBC_Epithelial_CopyKAT.rds")


# ------------------------------------------------------------
# 04. OUTPUT PATHS
# ------------------------------------------------------------

stage17_rds <- file.path(
  objects_dir,
  "Stage17_Pseudotime_Gene_Dynamics_Results.rds"
)

stage17_log <- file.path(
  logs_dir,
  "Stage17_Completion.log"
)

stage17_validation_file <- file.path(
  tables_validation_dir,
  "Stage17_Validation.csv"
)

synchronized_cells_file <- file.path(
  tables_supp_dir,
  "01_Stage17_Synchronized_Cells.csv"
)

all_gene_associations_file <- file.path(
  tables_supp_dir,
  "02_All_Genes_Pseudotime_Associations.csv"
)

plot_data_file <- file.path(
  tables_supp_dir,
  "03_Stage17_Plot_Data.csv"
)

significant_genes_file <- file.path(
  tables_main_dir,
  "01_Significant_Temporal_Genes.csv"
)

strong_genes_file <- file.path(
  tables_main_dir,
  "02_Strong_Temporal_Genes.csv"
)

top50_file <- file.path(
  tables_main_dir,
  "03_Top50_Temporal_Genes.csv"
)

top25_increasing_file <- file.path(
  tables_main_dir,
  "04_Top25_Increasing_Temporal_Genes.csv"
)

top25_decreasing_file <- file.path(
  tables_main_dir,
  "05_Top25_Decreasing_Temporal_Genes.csv"
)

summary_file <- file.path(
  tables_main_dir,
  "06_Stage17_Summary.csv"
)

stage17_summary_rds <- file.path(
  objects_dir,
  "Stage17_Pseudotime_Gene_Dynamics_Results.rds"
)


# ------------------------------------------------------------
# 05. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "matrixStats",
  "monocle3"
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
    paste(
      "Required package(s) not installed:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(matrixStats)
  library(monocle3)
})


# ------------------------------------------------------------
# 06. INPUT VALIDATION
# ------------------------------------------------------------

if (!dir.exists(frozen_stage15_cds_dir)) {
  
  stop(
    paste(
      "Frozen Stage 15 CDS directory not found:",
      frozen_stage15_cds_dir
    )
  )
}

if (!file.exists(stage16b_cells_file)) {
  
  stop(
    paste(
      "Stage 16B analysis-cell file not found:",
      stage16b_cells_file
    )
  )
}

if (!file.exists(stage16b_strong_genes_file)) {
  
  stop(
    paste(
      "Stage 16B strong-gene file not found:",
      stage16b_strong_genes_file
    )
  )
}

if (!file.exists(seurat_object_file)) {
  
  stop(
    paste(
      "Seurat object not found:",
      seurat_object_file
    )
  )
}


# ------------------------------------------------------------
# 07. LOAD FROZEN STAGE 15 CDS
# ------------------------------------------------------------

cat("Loading frozen Stage 15 CDS...\n")

# IMPORTANT:
# load_monocle_objects() expects the CDS directory itself,
# not the files contained inside that directory.
#
# The Stage 15 CDS is read-only and remains frozen.
# No trajectory reconstruction, root selection, or
# pseudotime recalculation is performed.

cds <- monocle3::load_monocle_objects(
  directory_path = frozen_stage15_cds_dir
)

if (!inherits(cds, "cell_data_set")) {
  
  stop(
    "Loaded Stage 15 object is not a valid monocle3 cell_data_set."
  )
}

cds_cell_ids <- colnames(cds)

if (is.null(cds_cell_ids) || length(cds_cell_ids) == 0) {
  
  stop(
    "Frozen Stage 15 CDS contains no cell IDs."
  )
}

cat(
  "Frozen Stage 15 CDS cells:",
  length(cds_cell_ids),
  "\n"
)

cat(
  "Frozen Stage 15 CDS genes:",
  nrow(cds),
  "\n"
)

# ------------------------------------------------------------
# 08. EXTRACT AUTHORITATIVE FROZEN PSEUDOTIME
# ------------------------------------------------------------

cat("\n")
cat("Extracting pseudotime from frozen Stage 15 CDS...\n")

# IMPORTANT:
# pseudotime(cds) is the authoritative Stage 15 pseudotime.
#
# This does NOT reconstruct the trajectory and does NOT
# recalculate pseudotime.
#
frozen_pseudotime <- monocle3::pseudotime(cds)

if (length(frozen_pseudotime) != length(cds_cell_ids)) {
  
  stop(
    paste(
      "Pseudotime length does not match CDS cell count.",
      "Pseudotime:",
      length(frozen_pseudotime),
      "CDS cells:",
      length(cds_cell_ids)
    )
  )
}

names(frozen_pseudotime) <- cds_cell_ids

finite_pseudotime_cells <- names(
  frozen_pseudotime[
    is.finite(frozen_pseudotime)
  ]
)

cat(
  "Finite pseudotime cells:",
  length(finite_pseudotime_cells),
  "\n"
)

cat(
  "Infinite pseudotime cells:",
  sum(is.infinite(frozen_pseudotime)),
  "\n"
)

cat(
  "NA pseudotime cells:",
  sum(is.na(frozen_pseudotime)),
  "\n"
)

if (length(finite_pseudotime_cells) == 0) {
  
  stop(
    "No finite pseudotime values were found in frozen Stage 15 CDS."
  )
}


# ------------------------------------------------------------
# 09. LOAD STAGE 16B ANALYSIS CELL LIST
# ------------------------------------------------------------

cat("\n")
cat("Loading Stage 16B analysis-cell list...\n")

stage16b_cells <- read.csv(
  stage16b_cells_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Cell" %in% colnames(stage16b_cells)) {
  
  stop(
    paste(
      "Stage16B_Analysis_Cells.csv must contain a 'Cell' column.",
      "Available columns:",
      paste(colnames(stage16b_cells), collapse = ", ")
    )
  )
}

stage16b_cell_ids <- unique(
  as.character(
    stage16b_cells$Cell
  )
)

stage16b_cell_ids <- stage16b_cell_ids[
  !is.na(stage16b_cell_ids) &
    nzchar(stage16b_cell_ids)
]

cat(
  "Stage 16B analysis cells:",
  length(stage16b_cell_ids),
  "\n"
)


# ------------------------------------------------------------
# 10. VERIFY STAGE 16B CELL LIST
# ------------------------------------------------------------

if (length(stage16b_cell_ids) != 6735) {
  
  warning(
    paste(
      "Expected 6735 Stage 16B analysis cells, but found",
      length(stage16b_cell_ids)
    )
  )
}

missing_from_cds <- setdiff(
  stage16b_cell_ids,
  cds_cell_ids
)

if (length(missing_from_cds) > 0) {
  
  stop(
    paste(
      length(missing_from_cds),
      "Stage 16B analysis cells are absent from the frozen Stage 15 CDS."
    )
  )
}

nonfinite_stage16b <- stage16b_cell_ids[
  !is.finite(
    frozen_pseudotime[
      stage16b_cell_ids
    ]
  )
]

if (length(nonfinite_stage16b) > 0) {
  
  stop(
    paste(
      length(nonfinite_stage16b),
      "Stage 16B analysis cells have non-finite pseudotime in frozen Stage 15 CDS."
    )
  )
}

analysis_cells <- stage16b_cell_ids

analysis_pseudotime <- as.numeric(
  frozen_pseudotime[
    analysis_cells
  ]
)

names(analysis_pseudotime) <- analysis_cells

if (length(analysis_pseudotime) != length(analysis_cells)) {
  
  stop(
    "Analysis-cell and pseudotime lengths are inconsistent."
  )
}

if (any(!is.finite(analysis_pseudotime))) {
  
  stop(
    "Non-finite pseudotime values remain after synchronization."
  )
}

cat(
  "Final synchronized analysis cells:",
  length(analysis_cells),
  "\n"
)

# ------------------------------------------------------------
# 11. LOAD STAGE 16B STRONG GENES
# ------------------------------------------------------------

cat("\n")
cat("Loading Stage 16B strong trajectory-associated genes...\n")

stage16b_strong <- read.csv(
  stage16b_strong_genes_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_stage16b_columns <- c(
  "gene_short_name",
  "gene_id"
)

missing_stage16b_columns <- setdiff(
  required_stage16b_columns,
  colnames(stage16b_strong)
)

if (length(missing_stage16b_columns) > 0) {
  
  stop(
    paste(
      "Stage 16B strong-gene table is missing required column(s):",
      paste(
        missing_stage16b_columns,
        collapse = ", "
      ),
      "\nAvailable columns:",
      paste(
        colnames(stage16b_strong),
        collapse = ", "
      )
    )
  )
}

# gene_short_name is the gene symbol used for synchronization
# with the Seurat expression matrix.
stage16b_strong_genes <- unique(
  as.character(
    stage16b_strong$gene_short_name
  )
)

stage16b_strong_genes <- stage16b_strong_genes[
  !is.na(stage16b_strong_genes) &
    nzchar(stage16b_strong_genes)
]

cat(
  "Stage 16B strong trajectory genes:",
  length(stage16b_strong_genes),
  "\n"
)

cat(
  "Stage 16B strong-gene table rows:",
  nrow(stage16b_strong),
  "\n"
)
# ------------------------------------------------------------
# 12. LOAD SEURAT OBJECT
# ------------------------------------------------------------

cat("\n")
cat("Loading Seurat expression object...\n")

seurat_obj <- readRDS(
  seurat_object_file
)

if (!inherits(seurat_obj, "Seurat")) {
  
  stop(
    "Loaded expression object is not a Seurat object."
  )
}

seurat_cells <- colnames(seurat_obj)

cat(
  "Seurat cells:",
  length(seurat_cells),
  "\n"
)


# ------------------------------------------------------------
# 13. VERIFY CELL SYNCHRONIZATION WITH SEURAT
# ------------------------------------------------------------

missing_from_seurat <- setdiff(
  analysis_cells,
  seurat_cells
)

if (length(missing_from_seurat) > 0) {
  
  stop(
    paste(
      length(missing_from_seurat),
      "Stage 16B analysis cells are absent from the Seurat object."
    )
  )
}

analysis_cells <- analysis_cells[
  analysis_cells %in% seurat_cells
]

analysis_pseudotime <- analysis_pseudotime[
  analysis_cells
]

if (length(analysis_cells) != 6735) {
  
  warning(
    paste(
      "Final synchronized Seurat cell count is",
      length(analysis_cells),
      "rather than 6735."
    )
  )
}

cat(
  "Final cells synchronized across Stage 15, Stage 16B and Seurat:",
  length(analysis_cells),
  "\n"
)

# ------------------------------------------------------------
# 14. EXTRACT RNA EXPRESSION
# ------------------------------------------------------------

cat("\n")
cat("Extracting normalized RNA expression layer...\n")

# ------------------------------------------------------------
# Identify assays robustly
# ------------------------------------------------------------

# For Seurat objects, the assay names are stored in the
# assays slot. Using names(seurat_obj@assays) avoids ambiguity
# between different Assays() methods/generics.

assay_names <- names(
  seurat_obj@assays
)

if (is.null(assay_names) || length(assay_names) == 0) {
  
  stop(
    "No assays were found in the Seurat object."
  )
}

cat(
  "Available assays:",
  paste(
    assay_names,
    collapse = ", "
  ),
  "\n"
)

if (!"RNA" %in% assay_names) {
  
  stop(
    paste(
      "RNA assay not found in Seurat object.",
      "Available assays:",
      paste(
        assay_names,
        collapse = ", "
      )
    )
  )
}

# ------------------------------------------------------------
# Set RNA as the active assay
# ------------------------------------------------------------

DefaultAssay(
  seurat_obj
) <- "RNA"


# ------------------------------------------------------------
# Identify available RNA layers
# ------------------------------------------------------------

rna_layers <- SeuratObject::Layers(
  seurat_obj[["RNA"]]
)

cat(
  "Available RNA layers:",
  paste(
    rna_layers,
    collapse = ", "
  ),
  "\n"
)

if (!"data" %in% rna_layers) {
  
  stop(
    paste(
      "RNA assay does not contain a normalized 'data' layer.",
      "Available layers:",
      paste(
        rna_layers,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# Extract normalized RNA data layer
# ------------------------------------------------------------

expression_matrix <- SeuratObject::LayerData(
  object = seurat_obj,
  assay = "RNA",
  layer = "data",
  cells = analysis_cells
)

if (is.null(expression_matrix)) {
  
  stop(
    "RNA expression matrix could not be extracted."
  )
}

cat(
  "Expression genes:",
  nrow(expression_matrix),
  "\n"
)

cat(
  "Expression cells:",
  ncol(expression_matrix),
  "\n"
)


# ------------------------------------------------------------
# Verify cell order
# ------------------------------------------------------------

if (!identical(
  colnames(expression_matrix),
  analysis_cells
)) {
  
  missing_expression_cells <- setdiff(
    analysis_cells,
    colnames(expression_matrix)
  )
  
  if (length(missing_expression_cells) > 0) {
    
    stop(
      paste(
        length(missing_expression_cells),
        "analysis cells are missing from the extracted RNA matrix."
      )
    )
  }
  
  expression_matrix <- expression_matrix[
    ,
    analysis_cells,
    drop = FALSE
  ]
}


# ------------------------------------------------------------
# Final expression-cell validation
# ------------------------------------------------------------

if (!identical(
  colnames(expression_matrix),
  analysis_cells
)) {
  
  stop(
    "RNA expression cell order does not match synchronized analysis-cell order."
  )
}

cat(
  "RNA expression cells synchronized:",
  ncol(expression_matrix),
  "\n"
)

# ------------------------------------------------------------
# 15. VERIFY EXPRESSION MATRIX
# ------------------------------------------------------------

if (!identical(
  colnames(expression_matrix),
  analysis_cells
)) {
  
  expression_matrix <- expression_matrix[
    ,
    analysis_cells,
    drop = FALSE
  ]
}

if (!all(
  is.finite(
    expression_matrix
  )
)) {
  
  nonfinite_count <- sum(
    !is.finite(
      expression_matrix
    )
  )
  
  stop(
    paste(
      "RNA expression matrix contains",
      nonfinite_count,
      "non-finite values.",
      "No values were replaced automatically."
    )
  )
}

n_genes <- nrow(expression_matrix)
n_cells <- ncol(expression_matrix)

if (n_cells < 10) {
  
  stop(
    "Too few synchronized cells for temporal correlation analysis."
  )
}

if (length(analysis_pseudotime) != n_cells) {
  
  stop(
    "Pseudotime length does not match expression matrix cell count."
  )
}


# ------------------------------------------------------------
# 16. IDENTIFY STAGE 16B STRONG GENES PRESENT IN EXPRESSION
# ------------------------------------------------------------

expression_genes <- rownames(expression_matrix)

strong_genes_in_expression <- intersect(
  stage16b_strong_genes,
  expression_genes
)

cat(
  "Stage 16B strong genes present in expression matrix:",
  length(strong_genes_in_expression),
  "\n"
)


# ------------------------------------------------------------
# 17. SPEARMAN CORRELATION
# ------------------------------------------------------------
#
# Spearman correlation is calculated through rank transformation.
#
# The Seurat v5 RNA data layer may be returned as a sparse matrix
# (for example, dgCMatrix). Each gene chunk is therefore converted
# explicitly to a standard numeric matrix before rank calculation.
#
# The complete expression matrix is NOT converted at once, in
# order to avoid unnecessary memory consumption.
#
# Because Section 15 verified that the expression data contain only
# finite values, all genes use the same number of observations.
#
# P-values use the standard asymptotic t approximation.
# BH correction is applied across all tested genes.
#
# ------------------------------------------------------------

cat("\n")
cat("Starting genome-wide pseudotime correlation analysis...\n")

gene_chunk_size <- 1000

gene_indices <- seq_len(n_genes)

gene_chunks <- split(
  gene_indices,
  ceiling(
    gene_indices / gene_chunk_size
  )
)

correlation_results <- vector(
  "list",
  length(gene_chunks)
)

chunk_counter <- 0

# Rank pseudotime ONCE because it is identical for every gene.
ranked_pseudotime <- rank(
  analysis_pseudotime,
  ties.method = "average"
)

centered_pseudotime <- ranked_pseudotime -
  mean(
    ranked_pseudotime
  )

denominator_pseudotime <- sqrt(
  sum(
    centered_pseudotime^2
  )
)

if (
  !is.finite(denominator_pseudotime) ||
  denominator_pseudotime == 0
) {
  
  stop(
    "Pseudotime ranks have zero or non-finite variance."
  )
}


for (idx in gene_chunks) {
  
  chunk_counter <- chunk_counter + 1
  
  # ----------------------------------------------------------
  # Extract current gene chunk
  # ----------------------------------------------------------
  
  expression_chunk <- expression_matrix[
    idx,
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # Convert sparse/DelayedArray-like object explicitly
  # to an ordinary numeric matrix.
  #
  # Only the current chunk is converted, not the complete
  # expression matrix.
  # ----------------------------------------------------------
  
  expression_chunk <- as.matrix(
    expression_chunk
  )
  
  storage.mode(
    expression_chunk
  ) <- "double"
  
  if (!is.matrix(expression_chunk)) {
    
    stop(
      paste(
        "Expression chunk",
        chunk_counter,
        "could not be converted to a matrix."
      )
    )
  }
  
  # ----------------------------------------------------------
  # Verify finite expression values
  # ----------------------------------------------------------
  
  if (!all(
    is.finite(
      expression_chunk
    )
  )) {
    
    stop(
      paste(
        "Non-finite expression values detected in gene chunk",
        chunk_counter,
        "."
      )
    )
  }
  
  # ----------------------------------------------------------
  # Rank expression values within each gene
  # ----------------------------------------------------------
  
  ranked_expression <- matrixStats::rowRanks(
    expression_chunk,
    ties.method = "average"
  )
  
  # ----------------------------------------------------------
  # Center expression ranks
  # ----------------------------------------------------------
  
  centered_expression <- ranked_expression -
    rowMeans(
      ranked_expression
    )
  
  # ----------------------------------------------------------
  # Spearman correlation numerator
  # ----------------------------------------------------------
  
  numerator <- as.vector(
    centered_expression %*%
      centered_pseudotime
  )
  
  # ----------------------------------------------------------
  # Spearman correlation denominator
  # ----------------------------------------------------------
  
  denominator_expression <- sqrt(
    rowSums(
      centered_expression^2
    )
  )
  
  rho <- numerator /
    (
      denominator_expression *
        denominator_pseudotime
    )
  
  rho[
    !is.finite(rho)
  ] <- NA_real_
  
  # ----------------------------------------------------------
  # Asymptotic p-value
  # ----------------------------------------------------------
  
  t_stat <- rho *
    sqrt(
      (n_cells - 2) /
        (1 - rho^2)
    )
  
  p_value <- 2 *
    stats::pt(
      -abs(t_stat),
      df = n_cells - 2
    )
  
  # ----------------------------------------------------------
  # Numerical handling of perfect correlations
  # ----------------------------------------------------------
  
  perfect_correlation <- is.finite(rho) &
    abs(rho) >= 1
  
  p_value[
    perfect_correlation
  ] <- 0
  
  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------
  
  correlation_results[[chunk_counter]] <- data.frame(
    Gene = rownames(expression_chunk),
    Spearman_rho = as.numeric(rho),
    p_value = as.numeric(p_value),
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
  # Progress reporting
  # ----------------------------------------------------------
  
  if (
    chunk_counter %% 5 == 0 ||
    chunk_counter == length(gene_chunks)
  ) {
    
    cat(
      "Processed gene chunks:",
      chunk_counter,
      "/",
      length(gene_chunks),
      "\n"
    )
  }
}


# ------------------------------------------------------------
# 18. COMBINE CORRELATION RESULTS
# ------------------------------------------------------------

association_table <- dplyr::bind_rows(
  correlation_results
)

if (nrow(association_table) != n_genes) {
  
  stop(
    paste(
      "Unexpected number of gene association results:",
      nrow(association_table),
      "expected:",
      n_genes
    )
  )
}

association_table <- association_table %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    Absolute_Spearman_rho = abs(
      Spearman_rho
    ),
    Direction = case_when(
      Spearman_rho > 0 ~ "Increasing",
      Spearman_rho < 0 ~ "Decreasing",
      TRUE ~ "No_Association"
    ),
    Stage16B_Strong_Gene = Gene %in% stage16b_strong_genes
  ) %>%
  arrange(
    FDR,
    desc(Absolute_Spearman_rho)
  )


# ------------------------------------------------------------
# 19. DEFINE TEMPORAL GENE SETS
# ------------------------------------------------------------

significant_temporal_genes <- association_table %>%
  filter(
    is.finite(FDR),
    is.finite(Spearman_rho),
    FDR < 0.05,
    Absolute_Spearman_rho >= 0.20
  ) %>%
  arrange(
    FDR,
    desc(Absolute_Spearman_rho)
  )

strong_temporal_genes <- association_table %>%
  filter(
    is.finite(FDR),
    is.finite(Spearman_rho),
    FDR < 0.05,
    Absolute_Spearman_rho >= 0.30
  ) %>%
  arrange(
    FDR,
    desc(Absolute_Spearman_rho)
  )

top50_temporal_genes <- association_table %>%
  filter(
    is.finite(FDR),
    is.finite(Spearman_rho)
  ) %>%
  arrange(
    FDR,
    desc(Absolute_Spearman_rho)
  ) %>%
  slice_head(
    n = 50
  )

top25_increasing <- association_table %>%
  filter(
    is.finite(FDR),
    is.finite(Spearman_rho),
    Spearman_rho > 0
  ) %>%
  arrange(
    FDR,
    desc(Spearman_rho)
  ) %>%
  slice_head(
    n = 25
  )

top25_decreasing <- association_table %>%
  filter(
    is.finite(FDR),
    is.finite(Spearman_rho),
    Spearman_rho < 0
  ) %>%
  arrange(
    FDR,
    Spearman_rho
  ) %>%
  slice_head(
    n = 25
  )


# ------------------------------------------------------------
# 20. STAGE 16B STRONG GENE TEMPORAL CHARACTERIZATION
# ------------------------------------------------------------

stage16b_strong_temporal <- association_table %>%
  filter(
    Stage16B_Strong_Gene
  ) %>%
  arrange(
    FDR,
    desc(Absolute_Spearman_rho)
  )


# ------------------------------------------------------------
# 21. SAVE MAIN TABLES
# ------------------------------------------------------------

write.csv(
  significant_temporal_genes,
  significant_genes_file,
  row.names = FALSE
)

write.csv(
  strong_temporal_genes,
  strong_genes_file,
  row.names = FALSE
)

write.csv(
  top50_temporal_genes,
  top50_file,
  row.names = FALSE
)

write.csv(
  top25_increasing,
  top25_increasing_file,
  row.names = FALSE
)

write.csv(
  top25_decreasing,
  top25_decreasing_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. SAVE SUPPLEMENTARY TABLES
# ------------------------------------------------------------

synchronized_cells_table <- data.frame(
  Cell = analysis_cells,
  Monocle3_pseudotime = analysis_pseudotime,
  stringsAsFactors = FALSE
)

write.csv(
  synchronized_cells_table,
  synchronized_cells_file,
  row.names = FALSE
)

write.csv(
  association_table,
  all_gene_associations_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. PLOT DATA
# ------------------------------------------------------------

plot_data <- data.frame(
  Cell = analysis_cells,
  Pseudotime = analysis_pseudotime,
  stringsAsFactors = FALSE
)

write.csv(
  plot_data,
  plot_data_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. MAIN FIGURE 1
# Pseudotime Distribution
# ------------------------------------------------------------

p1 <- ggplot(
  plot_data,
  aes(
    x = Pseudotime
  )
) +
  geom_histogram(
    bins = 50
  ) +
  labs(
    title = "Frozen Stage 15 Pseudotime Distribution",
    x = "Monocle3 pseudotime",
    y = "Number of cells"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figures_main_dir,
    "01_Pseudotime_Distribution.pdf"
  ),
  plot = p1,
  width = 8,
  height = 6,
  units = "in"
)


# ------------------------------------------------------------
# 25. MAIN FIGURE 2
# Pseudotime Ranked Gene Associations
# ------------------------------------------------------------

top30_plot <- association_table %>%
  filter(
    is.finite(Spearman_rho),
    is.finite(FDR)
  ) %>%
  arrange(
    FDR,
    desc(Absolute_Spearman_rho)
  ) %>%
  slice_head(
    n = 30
  )

top30_plot$Gene <- factor(
  top30_plot$Gene,
  levels = rev(
    top30_plot$Gene
  )
)

p2 <- ggplot(
  top30_plot,
  aes(
    x = Spearman_rho,
    y = Gene
  )
) +
  geom_point(
    size = 2.5
  ) +
  labs(
    title = "Top Temporal Gene Associations",
    x = "Spearman correlation with pseudotime",
    y = "Gene"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figures_main_dir,
    "02_Top30_Temporal_Gene_Associations.pdf"
  ),
  plot = p2,
  width = 9,
  height = 8,
  units = "in"
)


# ------------------------------------------------------------
# 26. MAIN FIGURE 3
# Correlation Distribution
# ------------------------------------------------------------

p3 <- ggplot(
  association_table %>%
    filter(
      is.finite(Spearman_rho)
    ),
  aes(
    x = Spearman_rho
  )
) +
  geom_histogram(
    bins = 60
  ) +
  labs(
    title = "Distribution of Gene–Pseudotime Spearman Correlations",
    x = "Spearman rho",
    y = "Number of genes"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figures_main_dir,
    "03_Spearman_Correlation_Distribution.pdf"
  ),
  plot = p3,
  width = 8,
  height = 6,
  units = "in"
)


# ------------------------------------------------------------
# 27. SUPPLEMENTARY FIGURE
# Pseudotime Rank Distribution
# ------------------------------------------------------------

p4 <- ggplot(
  data.frame(
    Rank = rank(
      analysis_pseudotime,
      ties.method = "average"
    )
  ),
  aes(
    x = Rank
  )
) +
  geom_histogram(
    bins = 50
  ) +
  labs(
    title = "Rank Distribution Used for Spearman Analysis",
    x = "Pseudotime rank",
    y = "Number of cells"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    figures_supp_dir,
    "04_Pseudotime_Rank_Distribution.pdf"
  ),
  plot = p4,
  width = 8,
  height = 6,
  units = "in"
)


# ------------------------------------------------------------
# 28. SUMMARY
# ------------------------------------------------------------

summary_table <- data.frame(
  Metric = c(
    "Frozen Stage 15 CDS cells",
    "Frozen Stage 15 CDS genes",
    "Finite Stage 15 pseudotime cells",
    "Stage 16B analysis cells",
    "Final synchronized cells",
    "Genes tested",
    "Genes with finite Spearman rho",
    "Stage 16B strong trajectory genes",
    "Stage 16B strong genes present in expression",
    "Significant temporal genes",
    "Strong temporal genes",
    "Increasing temporal genes",
    "Decreasing temporal genes"
  ),
  Value = c(
    length(cds_cell_ids),
    nrow(cds),
    length(finite_pseudotime_cells),
    length(stage16b_cell_ids),
    length(analysis_cells),
    n_genes,
    sum(
      is.finite(
        association_table$Spearman_rho
      )
    ),
    length(stage16b_strong_genes),
    length(strong_genes_in_expression),
    nrow(significant_temporal_genes),
    nrow(strong_temporal_genes),
    sum(
      significant_temporal_genes$Spearman_rho > 0
    ),
    sum(
      significant_temporal_genes$Spearman_rho < 0
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  summary_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. VALIDATION
# ------------------------------------------------------------

validation_table <- data.frame(
  Check = c(
    "Frozen Stage 15 CDS exists",
    "Frozen Stage 15 CDS is cell_data_set",
    "Pseudotime length matches CDS cells",
    "Finite pseudotime cells equal expected Stage 15 count",
    "Stage 16B cell list contains expected 6735 cells",
    "All Stage 16B cells present in frozen CDS",
    "All Stage 16B cells have finite frozen pseudotime",
    "All analysis cells present in Seurat object",
    "Final synchronized cell count equals 6735",
    "RNA data layer exists",
    "RNA expression contains only finite values",
    "Pseudotime length equals expression cell count",
    "Gene association result count equals expression gene count",
    "FDR values successfully calculated",
    "Stage 17 analysis completed without upstream rerun"
  ),
  Status = c(
    dir.exists(frozen_stage15_cds_dir),
    inherits(cds, "cell_data_set"),
    length(frozen_pseudotime) == length(cds_cell_ids),
    length(finite_pseudotime_cells) == 6735,
    length(stage16b_cell_ids) == 6735,
    length(missing_from_cds) == 0,
    length(nonfinite_stage16b) == 0,
    length(missing_from_seurat) == 0,
    length(analysis_cells) == 6735,
    "data" %in% rna_layers,
    all(is.finite(expression_matrix)),
    length(analysis_pseudotime) == ncol(expression_matrix),
    nrow(association_table) == nrow(expression_matrix),
    any(is.finite(association_table$FDR)),
    TRUE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  validation_table,
  stage17_validation_file,
  row.names = FALSE
)

if (!all(
  validation_table$Status
)) {
  
  failed_checks <- validation_table$Check[
    !validation_table$Status
  ]
  
  stop(
    paste(
      "Stage 17 validation failed:",
      paste(
        failed_checks,
        collapse = "; "
      )
    )
  )
}


# ------------------------------------------------------------
# 30. SAVE CENTRAL STAGE 17 OBJECT
# ------------------------------------------------------------

stage17_results <- list(
  
  stage = "17",
  
  analysis_name =
    "Pseudotime-Dependent Temporal Gene Dynamics",
  
  project_directory =
    project_dir,
  
  frozen_stage15_cds =
    frozen_stage15_cds_dir,
  
  stage16b_analysis_cells_file =
    stage16b_cells_file,
  
  stage16b_strong_genes_file =
    stage16b_strong_genes_file,
  
  seurat_expression_object =
    seurat_object_file,
  
  frozen_stage15_cell_count =
    length(cds_cell_ids),
  
  frozen_stage15_gene_count =
    nrow(cds),
  
  finite_stage15_pseudotime_cells =
    length(finite_pseudotime_cells),
  
  stage16b_analysis_cell_count =
    length(stage16b_cell_ids),
  
  synchronized_cell_count =
    length(analysis_cells),
  
  gene_count_tested =
    n_genes,
  
  finite_rho_count =
    sum(
      is.finite(
        association_table$Spearman_rho
      )
    ),
  
  stage16b_strong_gene_count =
    length(stage16b_strong_genes),
  
  stage16b_strong_genes_in_expression =
    strong_genes_in_expression,
  
  significant_temporal_gene_count =
    nrow(significant_temporal_genes),
  
  strong_temporal_gene_count =
    nrow(strong_temporal_genes),
  
  top50_temporal_genes =
    top50_temporal_genes,
  
  top25_increasing =
    top25_increasing,
  
  top25_decreasing =
    top25_decreasing,
  
  association_table =
    association_table,
  
  significant_temporal_genes =
    significant_temporal_genes,
  
  strong_temporal_genes =
    strong_temporal_genes,
  
  synchronized_cells =
    synchronized_cells_table,
  
  summary =
    summary_table,
  
  validation =
    validation_table,
  
  methodology = list(
    
    pseudotime_source =
      "Frozen Stage 15 Monocle3 CDS via pseudotime(cds)",
    
    trajectory_reconstruction =
      FALSE,
    
    root_selection =
      FALSE,
    
    pseudotime_recalculation =
      FALSE,
    
    stage16b_cell_list_role =
      "Cell-ID synchronization/filter only",
    
    correlation_method =
      "Spearman rank correlation",
    
    p_value_method =
      "Asymptotic t approximation",
    
    multiple_testing =
      "Benjamini-Hochberg FDR",
    
    significant_threshold =
      "FDR < 0.05 and |rho| >= 0.20",
    
    strong_threshold =
      "FDR < 0.05 and |rho| >= 0.30",
    
    expression_layer =
      "Seurat RNA normalized data layer",
    
    chunk_size =
      gene_chunk_size,
    
    sparse_matrix_handling =
      "Each gene chunk converted to standard matrix before ranking"
  )
)

saveRDS(
  stage17_results,
  stage17_rds
)


# ------------------------------------------------------------
# 31. COMPLETION LOG
# ------------------------------------------------------------

completion_time <- format(
  Sys.time(),
  "%Y-%m-%d %H:%M:%S"
)

log_lines <- c(
  
  "============================================================",
  "TNBC Single-Cell Analysis Pipeline",
  "Stage 17: Pseudotime-Dependent Temporal Gene Dynamics",
  "============================================================",
  "",
  paste(
    "Completion time:",
    completion_time
  ),
  "",
  paste(
    "Frozen Stage 15 CDS:",
    frozen_stage15_cds_dir
  ),
  "",
  paste(
    "Frozen Stage 15 CDS cells:",
    length(cds_cell_ids)
  ),
  "",
  paste(
    "Frozen Stage 15 CDS genes:",
    nrow(cds)
  ),
  "",
  paste(
    "Finite Stage 15 pseudotime cells:",
    length(finite_pseudotime_cells)
  ),
  "",
  paste(
    "Stage 16B analysis cells:",
    length(stage16b_cell_ids)
  ),
  "",
  paste(
    "Final synchronized cells:",
    length(analysis_cells)
  ),
  "",
  paste(
    "Genes tested:",
    n_genes
  ),
  "",
  paste(
    "Genes with finite Spearman rho:",
    sum(
      is.finite(
        association_table$Spearman_rho
      )
    )
  ),
  "",
  paste(
    "Stage 16B strong trajectory genes:",
    length(stage16b_strong_genes)
  ),
  "",
  paste(
    "Stage 16B strong genes present in expression:",
    length(strong_genes_in_expression)
  ),
  "",
  paste(
    "Significant temporal genes:",
    nrow(significant_temporal_genes)
  ),
  "",
  paste(
    "Strong temporal genes:",
    nrow(strong_temporal_genes)
  ),
  "",
  "Pseudotime source: frozen Stage 15 CDS",
  "Trajectory reconstruction: NOT PERFORMED",
  "Root selection: NOT PERFORMED",
  "Pseudotime recalculation: NOT PERFORMED",
  "Stage 16B cell list used only for cell synchronization.",
  "",
  "Correlation method: Spearman rank correlation",
  "Multiple-testing correction: Benjamini-Hochberg FDR",
  "Significant threshold: FDR < 0.05 and |rho| >= 0.20",
  "Strong threshold: FDR < 0.05 and |rho| >= 0.30",
  "Sparse expression chunks converted to ordinary matrices before ranking.",
  "",
  "Stage 15 remained frozen.",
  "No upstream analytical stage was rerun.",
  "",
  "STATUS: STAGE 17 COMPLETED SUCCESSFULLY",
  "============================================================"
)

writeLines(
  log_lines,
  stage17_log
)


# ------------------------------------------------------------
# 32. FINAL CONSOLE SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STAGE 17 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
cat(
  "Frozen Stage 15 cells:",
  length(cds_cell_ids),
  "\n"
)
cat(
  "Finite Stage 15 pseudotime cells:",
  length(finite_pseudotime_cells),
  "\n"
)
cat(
  "Stage 16B analysis cells:",
  length(stage16b_cell_ids),
  "\n"
)
cat(
  "Final synchronized cells:",
  length(analysis_cells),
  "\n"
)
cat(
  "Genes tested:",
  n_genes,
  "\n"
)
cat(
  "Finite Spearman correlations:",
  sum(
    is.finite(
      association_table$Spearman_rho
    )
  ),
  "\n"
)
cat(
  "Stage 16B strong genes:",
  length(stage16b_strong_genes),
  "\n"
)
cat(
  "Stage 16B strong genes in expression:",
  length(strong_genes_in_expression),
  "\n"
)
cat(
  "Significant temporal genes:",
  nrow(significant_temporal_genes),
  "\n"
)
cat(
  "Strong temporal genes:",
  nrow(strong_temporal_genes),
  "\n"
)
cat("\n")
cat(
  "Stage 17 RDS:\n",
  stage17_rds,
  "\n"
)
cat(
  "Stage 17 log:\n",
  stage17_log,
  "\n"
)
cat("\n")
cat(
  "Trajectory reconstruction: NOT PERFORMED\n"
)
cat(
  "Root selection: NOT PERFORMED\n"
)
cat(
  "Pseudotime recalculation: NOT PERFORMED\n"
)
cat(
  "Stage 15 remained frozen.\n"
)
cat("============================================================\n\n")

