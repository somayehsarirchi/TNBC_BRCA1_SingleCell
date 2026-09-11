# ==============================================================================
# STAGE 17B
# TEMPORAL GENE PROFILE CLUSTERING
# ==============================================================================
#
# Purpose:
#   Cluster the strong temporal genes identified in Stage 17 according to
#   their expression profiles across pseudotime.
#
# Frozen upstream input:
#   Stage 17 canonical RDS:
#   objects/Stage17_Pseudotime_Gene_Dynamics_Results.rds
#
# Stage 17 provides:
#   - strong_temporal_genes
#   - synchronized_cells
#       * Cell
#       * Monocle3_pseudotime
#   - seurat_expression_object
#
# Method:
#   1. Load frozen Stage 17 results
#   2. Load the exact Seurat expression object recorded by Stage 17
#   3. Synchronize 6735 finite-pseudotime cells
#   4. Extract the RNA/data layer
#   5. Restrict analysis to the 569 Stage 17 strong temporal genes
#   6. Divide pseudotime into 20 equal-width bins
#   7. Calculate mean expression per gene per pseudotime bin
#   8. Calculate gene-wise Z-scores across pseudotime bins
#   9. Evaluate K = 2:8 using mean silhouette width
#  10. Select K with the highest mean silhouette
#  11. Rank genes within clusters by Pearson correlation with cluster centroid
#  12. Export tables, figures, validation results, log and central RDS
#
# Important:
#   Stage 15, Stage 16 and Stage 17 remain frozen.
#   No trajectory reconstruction is performed.
#   No root selection is performed.
#   No pseudotime recalculation is performed.
#
# ==============================================================================


# ==============================================================================
# 1. PROJECT PATHS
# ==============================================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"

stage17_dir <-
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
    stage17_dir,
    "figures",
    "main"
  )

figures_supp_dir <-
  file.path(
    stage17_dir,
    "figures",
    "supplementary"
  )

tables_main_dir <-
  file.path(
    stage17_dir,
    "tables",
    "main"
  )

tables_supp_dir <-
  file.path(
    stage17_dir,
    "tables",
    "supplementary"
  )

tables_validation_dir <-
  file.path(
    stage17_dir,
    "tables",
    "validation"
  )

logs_dir <-
  file.path(
    stage17_dir,
    "logs"
  )


# ==============================================================================
# 2. CREATE / VERIFY OUTPUT DIRECTORIES
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


# ==============================================================================
# 3. ROBUST LOGGING
# ==============================================================================

stage17b_log <-
  file.path(
    logs_dir,
    "Stage17B_Completion.log"
  )

writeLines(
  character(0),
  con = stage17b_log
)

write_log <-
  function(...) {
    
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
      file = stage17b_log,
      append = TRUE
    )
  }


write_log(
  "======================================================================"
)

write_log(
  "STAGE 17B - TEMPORAL GENE PROFILE CLUSTERING"
)

write_log(
  "======================================================================"
)

write_log(
  "Project directory: ",
  project_dir
)

write_log(
  "Stage 17 directory: ",
  stage17_dir
)


# ==============================================================================
# 4. REQUIRED PACKAGES
# ==============================================================================

required_packages <-
  c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "dplyr",
    "tibble",
    "tidyr",
    "ggplot2",
    "cluster"
  )


missing_packages <-
  required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]


if (
  length(missing_packages) > 0
) {
  
  stop(
    paste(
      "The following required packages are missing:",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


suppressPackageStartupMessages(
  library(
    Seurat
  )
)

suppressPackageStartupMessages(
  library(
    SeuratObject
  )
)

suppressPackageStartupMessages(
  library(
    Matrix
  )
)

suppressPackageStartupMessages(
  library(
    dplyr
  )
)

suppressPackageStartupMessages(
  library(
    tibble
  )
)

suppressPackageStartupMessages(
  library(
    tidyr
  )
)

suppressPackageStartupMessages(
  library(
    ggplot2
  )
)

suppressPackageStartupMessages(
  library(
    cluster
  )
)


write_log(
  "Required packages loaded successfully."
)


# ==============================================================================
# 5. CANONICAL STAGE 17 INPUT
# ==============================================================================

stage17_rds <-
  file.path(
    objects_dir,
    "Stage17_Pseudotime_Gene_Dynamics_Results.rds"
  )


if (
  !file.exists(stage17_rds)
) {
  
  stop(
    paste(
      "Canonical Stage 17 RDS not found:",
      stage17_rds
    )
  )
}


write_log(
  "Canonical Stage 17 RDS found."
)

write_log(
  "Stage 17 RDS: ",
  stage17_rds
)


# ==============================================================================
# 6. LOAD FROZEN STAGE 17 RESULTS
# ==============================================================================

stage17 <-
  readRDS(
    stage17_rds
  )


write_log(
  "Stage 17 RDS loaded successfully."
)


write_log(
  "Stage 17 object class: ",
  paste(
    class(stage17),
    collapse = ", "
  )
)


# ==============================================================================
# 7. VERIFY EXPECTED STAGE 17 STRUCTURE
# ==============================================================================

required_stage17_fields <-
  c(
    "strong_temporal_genes",
    "synchronized_cells",
    "seurat_expression_object"
  )


missing_stage17_fields <-
  setdiff(
    required_stage17_fields,
    names(stage17)
  )


if (
  length(missing_stage17_fields) > 0
) {
  
  stop(
    paste(
      "Required Stage 17 fields are missing:",
      paste(
        missing_stage17_fields,
        collapse = ", "
      )
    )
  )
}


write_log(
  "Required Stage 17 fields verified."
)


# ==============================================================================
# 8. RECOVER STRONG TEMPORAL GENES
# ==============================================================================

strong_temporal_genes <-
  stage17$strong_temporal_genes


if (
  !is.data.frame(strong_temporal_genes)
) {
  
  stop(
    "stage17$strong_temporal_genes is not a data.frame."
  )
}


write_log(
  "Stage 17 strong temporal gene table loaded."
)

write_log(
  "Strong temporal gene table dimensions: ",
  nrow(strong_temporal_genes),
  " rows x ",
  ncol(strong_temporal_genes),
  " columns."
)

write_log(
  "Strong temporal gene columns: ",
  paste(
    colnames(strong_temporal_genes),
    collapse = ", "
  )
)


# ==============================================================================
# 9. IDENTIFY GENE COLUMN IN STRONG TEMPORAL GENE TABLE
# ==============================================================================

gene_column_candidates <-
  c(
    "Gene",
    "gene",
    "gene_short_name",
    "gene_symbol",
    "GeneSymbol"
  )


gene_column_matches <-
  gene_column_candidates[
    gene_column_candidates %in%
      colnames(strong_temporal_genes)
  ]


if (
  length(gene_column_matches) == 0
) {
  
  stop(
    paste(
      "No recognized gene-symbol column was found.",
      "Available columns:",
      paste(
        colnames(strong_temporal_genes),
        collapse = ", "
      )
    )
  )
}


gene_column <-
  gene_column_matches[1]


write_log(
  "Strong temporal gene column selected: ",
  gene_column
)


strong_temporal_genes$Gene <-
  as.character(
    strong_temporal_genes[[gene_column]]
  )


strong_temporal_genes <-
  strong_temporal_genes[
    !is.na(
      strong_temporal_genes$Gene
    ) &
      nzchar(
        strong_temporal_genes$Gene
      ),
    ,
    drop = FALSE
  ]


strong_temporal_genes <-
  strong_temporal_genes[
    !duplicated(
      strong_temporal_genes$Gene
    ),
    ,
    drop = FALSE
  ]


strong_genes <-
  strong_temporal_genes$Gene


write_log(
  "Unique strong temporal genes: ",
  length(strong_genes)
)


# ==============================================================================
# 10. VERIFY EXPECTED STAGE 17 STRONG-GENE COUNT
# ==============================================================================

expected_strong_gene_count <-
  569


if (
  length(strong_genes) !=
  expected_strong_gene_count
) {
  
  write_log(
    "WARNING: Expected 569 strong temporal genes, but found ",
    length(strong_genes),
    "."
  )
  
} else {
  
  write_log(
    "Verified Stage 17 strong temporal gene count: 569."
  )
}


# ==============================================================================
# 11. RECOVER SYNCHRONIZED CELLS AND PSEUDOTIME
# ==============================================================================

synchronized_cells <-
  stage17$synchronized_cells


if (
  !is.data.frame(synchronized_cells)
) {
  
  stop(
    "stage17$synchronized_cells is not a data.frame."
  )
}


required_sync_columns <-
  c(
    "Cell",
    "Monocle3_pseudotime"
  )


missing_sync_columns <-
  setdiff(
    required_sync_columns,
    colnames(synchronized_cells)
  )


if (
  length(missing_sync_columns) > 0
) {
  
  stop(
    paste(
      "Required synchronized-cell columns are missing:",
      paste(
        missing_sync_columns,
        collapse = ", "
      )
    )
  )
}


write_log(
  "Stage 17 synchronized-cell table verified."
)

write_log(
  "Synchronized-cell table dimensions: ",
  nrow(synchronized_cells),
  " rows x ",
  ncol(synchronized_cells),
  " columns."
)


# ==============================================================================
# 12. CLEAN AND VALIDATE SYNCHRONIZED CELL TABLE
# ==============================================================================

synchronized_cells$Cell <-
  as.character(
    synchronized_cells$Cell
  )

synchronized_cells$Monocle3_pseudotime <-
  as.numeric(
    synchronized_cells$Monocle3_pseudotime
  )


synchronized_cells <-
  synchronized_cells[
    !is.na(
      synchronized_cells$Cell
    ) &
      nzchar(
        synchronized_cells$Cell
      ),
    ,
    drop = FALSE
  ]


if (
  anyDuplicated(
    synchronized_cells$Cell
  ) > 0
) {
  
  stop(
    "Duplicate cell IDs detected in Stage 17 synchronized_cells."
  )
}


finite_pseudotime_mask <-
  is.finite(
    synchronized_cells$Monocle3_pseudotime
  )


finite_synchronized_cells <-
  synchronized_cells[
    finite_pseudotime_mask,
    ,
    drop = FALSE
  ]


write_log(
  "Finite Stage 17 pseudotime cells: ",
  nrow(finite_synchronized_cells)
)


if (
  nrow(finite_synchronized_cells) !=
  6735
) {
  
  write_log(
    "WARNING: Expected 6735 finite pseudotime cells, but found ",
    nrow(finite_synchronized_cells),
    "."
  )
  
} else {
  
  write_log(
    "Verified 6735 finite Stage 17 pseudotime cells."
  )
}


# ==============================================================================
# 13. LOAD THE EXACT SEURAT OBJECT RECORDED BY STAGE 17
# ==============================================================================

seurat_expression_file <-
  as.character(
    stage17$seurat_expression_object[1]
  )


if (
  is.na(seurat_expression_file) ||
  !nzchar(seurat_expression_file)
) {
  
  stop(
    "stage17$seurat_expression_object does not contain a valid file path."
  )
}


write_log(
  "Seurat expression object recorded by Stage 17:"
)

write_log(
  seurat_expression_file
)


if (
  !file.exists(seurat_expression_file)
) {
  
  stop(
    paste(
      "The Seurat expression object recorded by Stage 17 was not found:",
      seurat_expression_file
    )
  )
}


# ==============================================================================
# 14. LOAD SEURAT OBJECT
# ==============================================================================

seurat_obj <-
  readRDS(
    seurat_expression_file
  )


write_log(
  "Seurat expression object loaded successfully."
)

write_log(
  "Seurat object class: ",
  paste(
    class(seurat_obj),
    collapse = ", "
  )
)

write_log(
  "Total Seurat cells: ",
  ncol(seurat_obj)
)


# ==============================================================================
# 15. VERIFY RNA ASSAY
# ==============================================================================

assay_names <-
  names(
    seurat_obj@assays
  )


write_log(
  "Available Seurat assays: ",
  paste(
    assay_names,
    collapse = ", "
  )
)


if (
  !"RNA" %in% assay_names
) {
  
  stop(
    "RNA assay is not present in the Seurat object."
  )
}


rna_layers <-
  SeuratObject::Layers(
    seurat_obj[["RNA"]]
  )


write_log(
  "RNA assay layers: ",
  paste(
    rna_layers,
    collapse = ", "
  )
)


if (
  !"data" %in% rna_layers
) {
  
  stop(
    "RNA data layer is not present in the Seurat object."
  )
}


# ==============================================================================
# 16. SYNCHRONIZE STAGE 17 CELLS WITH SEURAT
# ==============================================================================

seurat_cells <-
  colnames(
    seurat_obj
  )


common_cells <-
  intersect(
    finite_synchronized_cells$Cell,
    seurat_cells
  )


if (
  length(common_cells) == 0
) {
  
  stop(
    "No Stage 17 synchronized cells were found in the Seurat object."
  )
}


analysis_cells <-
  finite_synchronized_cells[
    match(
      common_cells,
      finite_synchronized_cells$Cell
    ),
    ,
    drop = FALSE
  ]


analysis_cells <-
  analysis_cells[
    order(
      analysis_cells$Monocle3_pseudotime
    ),
    ,
    drop = FALSE
  ]


write_log(
  "Common Stage 17 / Seurat cells: ",
  nrow(analysis_cells)
)


if (
  nrow(analysis_cells) !=
  6735
) {
  
  write_log(
    "WARNING: Final synchronized cell count differs from expected 6735."
  )
  
} else {
  
  write_log(
    "Verified final synchronized cell count: 6735."
  )
}


# ==============================================================================
# 17. DEFINE PSEUDOTIME VECTOR
# ==============================================================================

analysis_pseudotime <-
  analysis_cells$Monocle3_pseudotime


names(analysis_pseudotime) <-
  analysis_cells$Cell


if (
  any(
    !is.finite(
      analysis_pseudotime
    )
  )
) {
  
  stop(
    "Non-finite pseudotime values remain after synchronization."
  )
}


pseudotime_min <-
  min(
    analysis_pseudotime
  )

pseudotime_max <-
  max(
    analysis_pseudotime
  )


write_log(
  "Pseudotime minimum: ",
  format(
    pseudotime_min,
    digits = 8
  )
)

write_log(
  "Pseudotime maximum: ",
  format(
    pseudotime_max,
    digits = 8
  )
)


# ==============================================================================
# 18. EXTRACT RNA DATA LAYER
# ==============================================================================

write_log(
  "Extracting RNA/data layer for synchronized cells."
)


rna_data <-
  SeuratObject::LayerData(
    seurat_obj,
    assay = "RNA",
    layer = "data",
    cells = analysis_cells$Cell
  )


write_log(
  "RNA/data dimensions: ",
  nrow(rna_data),
  " genes x ",
  ncol(rna_data),
  " cells."
)


# ==============================================================================
# 19. VERIFY CELL ORDER
# ==============================================================================

if (
  !identical(
    colnames(rna_data),
    analysis_cells$Cell
  )
) {
  
  common_expression_cells <-
    intersect(
      analysis_cells$Cell,
      colnames(rna_data)
    )
  
  
  if (
    length(common_expression_cells) !=
    nrow(analysis_cells)
  ) {
    
    stop(
      "RNA expression layer does not contain all synchronized cells."
    )
  }
  
  
  analysis_cells <-
    analysis_cells[
      match(
        common_expression_cells,
        analysis_cells$Cell
      ),
      ,
      drop = FALSE
    ]
  
  
  analysis_pseudotime <-
    analysis_cells$Monocle3_pseudotime
  
  names(analysis_pseudotime) <-
    analysis_cells$Cell
  
  
  rna_data <-
    rna_data[
      ,
      common_expression_cells,
      drop = FALSE
    ]
}


if (
  !identical(
    colnames(rna_data),
    analysis_cells$Cell
  )
) {
  
  stop(
    "Final RNA expression cell order does not match pseudotime cell order."
  )
}


write_log(
  "RNA expression cell order synchronized successfully."
)


# ==============================================================================
# 20. MATCH STRONG TEMPORAL GENES TO EXPRESSION
# ==============================================================================

expression_genes <-
  rownames(
    rna_data
  )


strong_genes_in_expression <-
  intersect(
    strong_genes,
    expression_genes
  )


write_log(
  "Stage 17 strong temporal genes: ",
  length(strong_genes)
)

write_log(
  "Strong temporal genes found in RNA expression: ",
  length(strong_genes_in_expression)
)


if (
  length(strong_genes_in_expression) == 0
) {
  
  stop(
    "None of the Stage 17 strong temporal genes are present in the RNA expression layer."
  )
}


if (
  length(strong_genes_in_expression) !=
  length(strong_genes)
) {
  
  missing_strong_genes <-
    setdiff(
      strong_genes,
      strong_genes_in_expression
    )
  
  
  write_log(
    "WARNING: ",
    length(missing_strong_genes),
    " Stage 17 strong temporal genes were not found in RNA expression."
  )
  
} else {
  
  write_log(
    "Verified all 569 Stage 17 strong temporal genes are present in expression."
  )
}


# Preserve Stage 17 gene order.

strong_genes_in_expression <-
  strong_genes[
    strong_genes %in%
      strong_genes_in_expression
  ]


rna_data_strong <-
  rna_data[
    strong_genes_in_expression,
    analysis_cells$Cell,
    drop = FALSE
  ]


# ==============================================================================
# 21. VERIFY EXPRESSION MATRIX
# ==============================================================================

if (
  nrow(rna_data_strong) == 0 ||
  ncol(rna_data_strong) == 0
) {
  
  stop(
    "The strong temporal gene expression matrix is empty."
  )
}


write_log(
  "Strong temporal gene expression matrix: ",
  nrow(rna_data_strong),
  " genes x ",
  ncol(rna_data_strong),
  " cells."
)


# ==============================================================================
# 22. DEFINE 20 EQUAL-WIDTH PSEUDOTIME BINS
# ==============================================================================

n_bins <-
  20


if (
  pseudotime_max <=
  pseudotime_min
) {
  
  stop(
    "Pseudotime range is not positive; cannot construct pseudotime bins."
  )
}


bin_breaks <-
  seq(
    from = pseudotime_min,
    to = pseudotime_max,
    length.out = n_bins + 1
  )


pseudotime_bins <-
  cut(
    analysis_pseudotime,
    breaks = bin_breaks,
    include.lowest = TRUE,
    labels = FALSE
  )


if (
  any(
    is.na(
      pseudotime_bins
    )
  )
) {
  
  stop(
    "NA pseudotime-bin assignments detected."
  )
}


bin_counts <-
  table(
    factor(
      pseudotime_bins,
      levels = seq_len(n_bins)
    )
  )


write_log(
  "Pseudotime bins: ",
  n_bins
)

write_log(
  "Pseudotime bin cell counts: ",
  paste(
    as.integer(bin_counts),
    collapse = ", "
  )
)


if (
  any(
    bin_counts == 0
  )
) {
  
  write_log(
    "WARNING: One or more equal-width pseudotime bins contain zero cells."
  )
}


# ==============================================================================
# 23. CALCULATE MEAN EXPRESSION PER PSEUDOTIME BIN
# ==============================================================================

write_log(
  "Calculating mean expression per pseudotime bin."
)


mean_expression_matrix <-
  matrix(
    NA_real_,
    nrow = nrow(rna_data_strong),
    ncol = n_bins,
    dimnames = list(
      rownames(rna_data_strong),
      paste0(
        "Bin_",
        seq_len(n_bins)
      )
    )
  )


for (
  b in seq_len(n_bins)
) {
  
  cell_indices <-
    which(
      pseudotime_bins == b
    )
  
  
  if (
    length(cell_indices) == 0
  ) {
    
    next
  }
  
  
  mean_expression_matrix[
    ,
    b
  ] <-
    Matrix::rowMeans(
      rna_data_strong[
        ,
        cell_indices,
        drop = FALSE
      ]
    )
}


# ==============================================================================
# 24. VALIDATE MEAN EXPRESSION MATRIX
# ==============================================================================

if (
  any(
    !is.finite(
      mean_expression_matrix
    )
  )
) {
  
  stop(
    "Non-finite values detected in mean expression matrix."
  )
}


write_log(
  "Mean expression matrix generated successfully."
)

write_log(
  "Mean expression dimensions: ",
  nrow(mean_expression_matrix),
  " genes x ",
  ncol(mean_expression_matrix),
  " pseudotime bins."
)


# ==============================================================================
# 25. GENE-WISE Z-SCORE ACROSS PSEUDOTIME BINS
# ==============================================================================

write_log(
  "Calculating gene-wise Z-scores across pseudotime bins."
)


row_means <-
  rowMeans(
    mean_expression_matrix
  )


row_sds <-
  apply(
    mean_expression_matrix,
    1,
    sd
  )


zero_variance_mask <-
  !is.finite(row_sds) |
  row_sds == 0


zero_variance_genes <-
  rownames(
    mean_expression_matrix
  )[
    zero_variance_mask
  ]


write_log(
  "Zero-variance genes across pseudotime bins: ",
  length(zero_variance_genes)
)


if (
  all(
    zero_variance_mask
  )
) {
  
  stop(
    "All strong temporal genes have zero variance across pseudotime bins."
  )
}


zscore_matrix <-
  sweep(
    mean_expression_matrix,
    1,
    row_means,
    FUN = "-"
  )


zscore_matrix <-
  sweep(
    zscore_matrix,
    1,
    row_sds,
    FUN = "/"
  )


zscore_matrix <-
  zscore_matrix[
    !zero_variance_mask,
    ,
    drop = FALSE
  ]


if (
  any(
    !is.finite(
      zscore_matrix
    )
  )
) {
  
  stop(
    "Non-finite values detected in the final Z-score matrix."
  )
}


write_log(
  "Final Z-score matrix: ",
  nrow(zscore_matrix),
  " genes x ",
  ncol(zscore_matrix),
  " bins."
)


# ==============================================================================
# 26. SAVE Z-SCORE PROFILE TABLE
# ==============================================================================

zscore_profile_table <-
  as.data.frame(
    zscore_matrix,
    check.names = FALSE
  ) %>%
  tibble::rownames_to_column(
    var = "Gene"
  )


write.csv(
  zscore_profile_table,
  file =
    file.path(
      tables_supp_dir,
      "Stage17B_Gene_ZScore_Profiles.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 27. GENE DISTANCE MATRIX
# ==============================================================================

write_log(
  "Calculating Euclidean gene distance matrix."
)


gene_distance <-
  dist(
    zscore_matrix,
    method = "euclidean"
  )


# ==============================================================================
# 28. K-MEANS CANDIDATE EVALUATION
# ==============================================================================

set.seed(
  170217
)


k_candidates <-
  2:8


silhouette_results <-
  data.frame(
    K = integer(),
    Mean_Silhouette = numeric(),
    stringsAsFactors = FALSE
  )


kmeans_results <-
  list()


write_log(
  "Evaluating K candidates: ",
  paste(
    k_candidates,
    collapse = ", "
  )
)


for (
  k in k_candidates
) {
  
  write_log(
    "Evaluating K = ",
    k
  )
  
  
  km <-
    kmeans(
      zscore_matrix,
      centers = k,
      nstart = 50,
      iter.max = 100
    )
  
  
  sil <-
    cluster::silhouette(
      km$cluster,
      gene_distance
    )
  
  
  mean_silhouette <-
    mean(
      sil[, "sil_width"]
    )
  
  
  silhouette_results <-
    rbind(
      silhouette_results,
      data.frame(
        K =
          k,
        Mean_Silhouette =
          mean_silhouette,
        stringsAsFactors = FALSE
      )
    )
  
  
  kmeans_results[[as.character(k)]] <-
    list(
      model = km,
      silhouette = sil
    )
}


# ==============================================================================
# 29. SELECT OPTIMAL K
# ==============================================================================

best_row <-
  silhouette_results[
    which.max(
      silhouette_results$Mean_Silhouette
    ),
    ,
    drop = FALSE
  ]


selected_k <-
  best_row$K[1]


selected_mean_silhouette <-
  best_row$Mean_Silhouette[1]


write_log(
  "Selected K: ",
  selected_k
)

write_log(
  "Selected mean silhouette: ",
  format(
    selected_mean_silhouette,
    digits = 8
  )
)


# ==============================================================================
# 30. SAVE K-SELECTION TABLE
# ==============================================================================

write.csv(
  silhouette_results,
  file =
    file.path(
      tables_validation_dir,
      "Stage17B_K_Silhouette_Evaluation.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 31. FINAL K-MEANS MODEL
# ==============================================================================

set.seed(
  170217
)


final_kmeans <-
  kmeans(
    zscore_matrix,
    centers = selected_k,
    nstart = 100,
    iter.max = 200
  )


gene_clusters <-
  final_kmeans$cluster


names(gene_clusters) <-
  rownames(
    zscore_matrix
  )


# ==============================================================================
# 32. FINAL SILHOUETTE VALUES
# ==============================================================================

final_silhouette <-
  cluster::silhouette(
    gene_clusters,
    gene_distance
  )


final_silhouette_width <-
  final_silhouette[
    ,
    "sil_width"
  ]


names(final_silhouette_width) <-
  rownames(
    zscore_matrix
  )


# ==============================================================================
# 33. CLUSTER ASSIGNMENTS
# ==============================================================================

cluster_assignments <-
  data.frame(
    Gene =
      names(gene_clusters),
    Cluster =
      as.integer(
        gene_clusters
      ),
    Silhouette_Width =
      as.numeric(
        final_silhouette_width[
          names(gene_clusters)
        ]
      ),
    stringsAsFactors = FALSE
  )


cluster_assignments <-
  cluster_assignments[
    order(
      cluster_assignments$Cluster,
      -cluster_assignments$Silhouette_Width
    ),
    ,
    drop = FALSE
  ]


write.csv(
  cluster_assignments,
  file =
    file.path(
      tables_main_dir,
      "Stage17B_Temporal_Gene_Cluster_Assignments.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 34. CLUSTER CENTROIDS
# ==============================================================================

cluster_centroids <-
  matrix(
    NA_real_,
    nrow = selected_k,
    ncol = n_bins,
    dimnames = list(
      paste0(
        "Cluster_",
        seq_len(selected_k)
      ),
      colnames(zscore_matrix)
    )
  )


for (
  k in seq_len(selected_k)
) {
  
  genes_k <-
    names(
      gene_clusters[
        gene_clusters == k
      ]
    )
  
  
  cluster_centroids[
    k,
  ] <-
    colMeans(
      zscore_matrix[
        genes_k,
        ,
        drop = FALSE
      ]
    )
}


write.csv(
  cluster_centroids,
  file =
    file.path(
      tables_main_dir,
      "Stage17B_Cluster_Centroids.csv"
    ),
  row.names = TRUE
)


# ==============================================================================
# 35. RANK GENES BY PEARSON CORRELATION WITH CLUSTER CENTROID
# ==============================================================================

write_log(
  "Ranking genes by Pearson correlation with cluster centroids."
)


gene_cluster_correlation <-
  numeric(
    nrow(
      zscore_matrix
    )
  )


names(gene_cluster_correlation) <-
  rownames(
    zscore_matrix
  )


for (
  gene in rownames(zscore_matrix)
) {
  
  k <-
    gene_clusters[
      gene
    ]
  
  
  centroid <-
    cluster_centroids[
      k,
    ]
  
  
  gene_profile <-
    zscore_matrix[
      gene,
    ]
  
  
  gene_cluster_correlation[
    gene
  ] <-
    cor(
      gene_profile,
      centroid,
      method = "pearson"
    )
}


gene_cluster_correlation_table <-
  data.frame(
    Gene =
      names(gene_cluster_correlation),
    Cluster =
      as.integer(
        gene_clusters[
          names(gene_cluster_correlation)
        ]
      ),
    Centroid_Pearson_Correlation =
      as.numeric(
        gene_cluster_correlation
      ),
    Silhouette_Width =
      as.numeric(
        final_silhouette_width[
          names(gene_cluster_correlation)
        ]
      ),
    stringsAsFactors = FALSE
  )


gene_cluster_correlation_table <-
  gene_cluster_correlation_table[
    order(
      gene_cluster_correlation_table$Cluster,
      -gene_cluster_correlation_table$Centroid_Pearson_Correlation
    ),
    ,
    drop = FALSE
  ]


gene_cluster_correlation_table$Rank_Within_Cluster <-
  ave(
    gene_cluster_correlation_table$Centroid_Pearson_Correlation,
    gene_cluster_correlation_table$Cluster,
    FUN = function(x) {
      
      rank(
        -x,
        ties.method = "first"
      )
      
    }
  )


gene_cluster_correlation_table <-
  gene_cluster_correlation_table[
    order(
      gene_cluster_correlation_table$Cluster,
      gene_cluster_correlation_table$Rank_Within_Cluster
    ),
    ,
    drop = FALSE
  ]


write.csv(
  gene_cluster_correlation_table,
  file =
    file.path(
      tables_main_dir,
      "Stage17B_Temporal_Gene_Cluster_Ranking.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 36. TOP 20 GENES PER CLUSTER
# ==============================================================================

top_n_per_cluster <-
  20


top_genes_per_cluster <-
  gene_cluster_correlation_table %>%
  dplyr::group_by(
    Cluster
  ) %>%
  dplyr::arrange(
    Rank_Within_Cluster,
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = top_n_per_cluster
  ) %>%
  dplyr::ungroup()


write.csv(
  top_genes_per_cluster,
  file =
    file.path(
      tables_main_dir,
      "Stage17B_Top20_Genes_Per_Cluster.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 37. LONG-FORM Z-SCORE PROFILE TABLE
# ==============================================================================

zscore_profile_long <-
  zscore_profile_table %>%
  tidyr::pivot_longer(
    cols = -Gene,
    names_to = "Pseudotime_Bin",
    values_to = "Z_Score"
  ) %>%
  dplyr::left_join(
    cluster_assignments %>%
      dplyr::select(
        Gene,
        Cluster
      ),
    by = "Gene"
  )


write.csv(
  zscore_profile_long,
  file =
    file.path(
      tables_supp_dir,
      "Stage17B_Gene_ZScore_Profiles_Long.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 38. CLUSTER SUMMARY
# ==============================================================================

cluster_summary <-
  cluster_assignments %>%
  dplyr::left_join(
    gene_cluster_correlation_table %>%
      dplyr::select(
        Gene,
        Centroid_Pearson_Correlation
      ),
    by = "Gene"
  ) %>%
  dplyr::group_by(
    Cluster
  ) %>%
  dplyr::summarise(
    N_Genes =
      dplyr::n(),
    Mean_Silhouette =
      mean(
        Silhouette_Width
      ),
    Median_Silhouette =
      median(
        Silhouette_Width
      ),
    Mean_Centroid_Pearson =
      mean(
        Centroid_Pearson_Correlation
      ),
    Median_Centroid_Pearson =
      median(
        Centroid_Pearson_Correlation
      ),
    .groups = "drop"
  )


write.csv(
  cluster_summary,
  file =
    file.path(
      tables_main_dir,
      "Stage17B_Cluster_Summary.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 39. CLUSTER PROFILE PLOT
# ==============================================================================

cluster_centroid_long <-
  as.data.frame(
    cluster_centroids,
    check.names = FALSE
  ) %>%
  tibble::rownames_to_column(
    var = "Cluster"
  ) %>%
  tidyr::pivot_longer(
    cols = -Cluster,
    names_to = "Pseudotime_Bin",
    values_to = "Mean_Z_Score"
  )


cluster_profile_plot <-
  ggplot(
    cluster_centroid_long,
    aes(
      x = Pseudotime_Bin,
      y = Mean_Z_Score,
      group = Cluster
    )
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~ Cluster,
    scales = "free_y"
  ) +
  labs(
    title =
      paste0(
        "Stage 17B: Temporal Gene Cluster Profiles (K = ",
        selected_k,
        ")"
      ),
    x = "Pseudotime bin",
    y = "Mean gene Z-score"
  ) +
  theme_bw() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    strip.text =
      element_text(
        face = "bold"
      )
  )


ggsave(
  filename =
    file.path(
      figures_main_dir,
      "Stage17B_Temporal_Cluster_Profiles.pdf"
    ),
  plot =
    cluster_profile_plot,
  width = 11,
  height = 8.5,
  units = "in"
)


# ==============================================================================
# 40. SILHOUETTE K-SELECTION PLOT
# ==============================================================================

silhouette_plot <-
  ggplot(
    silhouette_results,
    aes(
      x = K,
      y = Mean_Silhouette
    )
  ) +
  geom_line() +
  geom_point(
    size = 2.5
  ) +
  geom_point(
    data = best_row,
    size = 4
  ) +
  labs(
    title =
      "Stage 17B: K Selection by Mean Silhouette Width",
    x = "Number of clusters (K)",
    y = "Mean silhouette width"
  ) +
  scale_x_continuous(
    breaks = k_candidates
  ) +
  theme_bw()


ggsave(
  filename =
    file.path(
      figures_supp_dir,
      "Stage17B_K_Silhouette_Evaluation.pdf"
    ),
  plot =
    silhouette_plot,
  width = 8,
  height = 6,
  units = "in"
)


# ==============================================================================
# 41. SILHOUETTE DISTRIBUTION PLOT
# ==============================================================================

silhouette_distribution <-
  data.frame(
    Gene =
      names(final_silhouette_width),
    Cluster =
      as.factor(
        gene_clusters[
          names(final_silhouette_width)
        ]
      ),
    Silhouette_Width =
      as.numeric(
        final_silhouette_width
      ),
    stringsAsFactors = FALSE
  )


silhouette_distribution_plot <-
  ggplot(
    silhouette_distribution,
    aes(
      x = Cluster,
      y = Silhouette_Width
    )
  ) +
  geom_boxplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title =
      paste0(
        "Stage 17B: Silhouette Width by Cluster (K = ",
        selected_k,
        ")"
      ),
    x = "Cluster",
    y = "Silhouette width"
  ) +
  theme_bw()


ggsave(
  filename =
    file.path(
      figures_supp_dir,
      "Stage17B_Silhouette_Width_By_Cluster.pdf"
    ),
  plot =
    silhouette_distribution_plot,
  width = 8,
  height = 6,
  units = "in"
)


# ==============================================================================
# 42. OPTIONAL HEATMAP
# ==============================================================================

if (
  requireNamespace(
    "pheatmap",
    quietly = TRUE
  )
) {
  
  heatmap_file <-
    file.path(
      figures_main_dir,
      "Stage17B_Temporal_Gene_Clusters_Heatmap.pdf"
    )
  
  
  ordered_genes <-
    rownames(
      zscore_matrix
    )[
      order(
        gene_clusters
      )
    ]
  
  
  pheatmap::pheatmap(
    zscore_matrix[
      ordered_genes,
      ,
      drop = FALSE
    ],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    main =
      paste0(
        "Stage 17B: Strong Temporal Gene Clusters (K = ",
        selected_k,
        ")"
      ),
    filename =
      heatmap_file,
    width = 9,
    height = 12
  )
  
  
  write_log(
    "Heatmap generated: ",
    heatmap_file
  )
  
} else {
  
  write_log(
    "pheatmap not installed; heatmap generation skipped."
  )
}


# ==============================================================================
# 43. FINAL NUMERICAL VALIDATION
# ==============================================================================

validation_results <-
  data.frame(
    Metric =
      c(
        "Stage17_Strong_Temporal_Genes",
        "Strong_Genes_In_Expression",
        "Final_Synchronized_Cells",
        "Pseudotime_Bins",
        "Genes_After_ZeroVariance_Filter",
        "Selected_K",
        "Selected_Mean_Silhouette",
        "Minimum_Silhouette",
        "Median_Silhouette",
        "Maximum_Silhouette",
        "Mean_Centroid_Correlation",
        "Minimum_Centroid_Correlation"
      ),
    Value =
      c(
        length(strong_genes),
        length(strong_genes_in_expression),
        nrow(analysis_cells),
        n_bins,
        nrow(zscore_matrix),
        selected_k,
        selected_mean_silhouette,
        min(
          final_silhouette_width
        ),
        median(
          final_silhouette_width
        ),
        max(
          final_silhouette_width
        ),
        mean(
          gene_cluster_correlation
        ),
        min(
          gene_cluster_correlation
        )
      ),
    stringsAsFactors = FALSE
  )


write.csv(
  validation_results,
  file =
    file.path(
      tables_validation_dir,
      "Stage17B_Final_Validation.csv"
    ),
  row.names = FALSE
)


# ==============================================================================
# 44. STRUCTURAL VALIDATION
# ==============================================================================

validation_checks <-
  data.frame(
    Check =
      c(
        "Stage17 strong temporal gene table available",
        "Stage17 strong temporal genes > 0",
        "Expected strong temporal gene count = 569",
        "Synchronized-cell table available",
        "Expected synchronized cell count = 6735",
        "All synchronized pseudotime values finite",
        "RNA assay available",
        "RNA data layer available",
        "Strong temporal genes found in expression",
        "Mean expression matrix finite",
        "Z-score matrix finite",
        "Selected K within candidate range",
        "All analyzed genes assigned to a cluster",
        "Cluster centroid dimensions valid",
        "All cluster assignments unique per gene"
      ),
    Passed =
      c(
        is.data.frame(
          strong_temporal_genes
        ),
        length(strong_genes) > 0,
        length(strong_genes) ==
          expected_strong_gene_count,
        is.data.frame(
          synchronized_cells
        ),
        nrow(analysis_cells) ==
          6735,
        all(
          is.finite(
            analysis_pseudotime
          )
        ),
        "RNA" %in% assay_names,
        "data" %in% rna_layers,
        length(strong_genes_in_expression) > 0,
        all(
          is.finite(
            mean_expression_matrix
          )
        ),
        all(
          is.finite(
            zscore_matrix
          )
        ),
        selected_k %in%
          k_candidates,
        length(gene_clusters) ==
          nrow(zscore_matrix),
        all(
          dim(cluster_centroids) ==
            c(
              selected_k,
              n_bins
            )
        ),
        !anyDuplicated(
          names(gene_clusters)
        )
      ),
    stringsAsFactors = FALSE
  )


write.csv(
  validation_checks,
  file =
    file.path(
      tables_validation_dir,
      "Stage17B_Structural_Validation.csv"
    ),
  row.names = FALSE
)


if (
  !all(
    validation_checks$Passed
  )
) {
  
  failed_checks <-
    validation_checks$Check[
      !validation_checks$Passed
    ]
  
  
  stop(
    paste(
      "Stage 17B structural validation failed:",
      paste(
        failed_checks,
        collapse = "; "
      )
    )
  )
}


write_log(
  "All structural validation checks passed."
)


# ==============================================================================
# 45. SAVE CENTRAL STAGE 17B OBJECT
# ==============================================================================

stage17b_results <-
  list(
    
    stage =
      "17B",
    
    analysis_name =
      "Temporal Gene Profile Clustering",
    
    source_stage17_rds =
      stage17_rds,
    
    source_seurat_expression_object =
      seurat_expression_file,
    
    synchronized_cells =
      analysis_cells,
    
    pseudotime =
      analysis_pseudotime,
    
    pseudotime_bins =
      pseudotime_bins,
    
    bin_breaks =
      bin_breaks,
    
    bin_counts =
      bin_counts,
    
    n_bins =
      n_bins,
    
    strong_temporal_genes =
      strong_genes,
    
    strong_temporal_genes_in_expression =
      strong_genes_in_expression,
    
    zero_variance_genes =
      zero_variance_genes,
    
    mean_expression_matrix =
      mean_expression_matrix,
    
    zscore_matrix =
      zscore_matrix,
    
    k_candidates =
      k_candidates,
    
    silhouette_results =
      silhouette_results,
    
    selected_k =
      selected_k,
    
    selected_mean_silhouette =
      selected_mean_silhouette,
    
    final_kmeans =
      final_kmeans,
    
    gene_clusters =
      gene_clusters,
    
    cluster_centroids =
      cluster_centroids,
    
    cluster_assignments =
      cluster_assignments,
    
    final_silhouette_width =
      final_silhouette_width,
    
    gene_cluster_correlation =
      gene_cluster_correlation,
    
    gene_cluster_correlation_table =
      gene_cluster_correlation_table,
    
    top_genes_per_cluster =
      top_genes_per_cluster,
    
    cluster_summary =
      cluster_summary,
    
    validation_results =
      validation_results,
    
    validation_checks =
      validation_checks,
    
    methodology =
      list(
        pseudotime_source =
          "Frozen Stage 17 synchronized_cells$Monocle3_pseudotime",
        
        pseudotime_bins =
          "20 equal-width bins",
        
        expression_summary =
          "Mean RNA/data-layer expression per pseudotime bin",
        
        profile_normalization =
          "Gene-wise Z-score across pseudotime bins",
        
        clustering =
          "K-means on gene-wise Z-scored temporal profiles",
        
        distance =
          "Euclidean",
        
        k_candidates =
          "2:8",
        
        k_selection =
          "Maximum mean silhouette width",
        
        gene_ranking =
          "Pearson correlation with cluster centroid",
        
        top_genes =
          "Top 20 genes per cluster"
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


stage17b_rds <-
  file.path(
    objects_dir,
    "Stage17B_Temporal_Clustering_Results.rds"
  )


saveRDS(
  stage17b_results,
  file =
    stage17b_rds,
  compress = TRUE
)


write_log(
  "Central Stage 17B RDS saved: ",
  stage17b_rds
)


# ==============================================================================
# 46. OUTPUT INVENTORY
# ==============================================================================

output_files <-
  c(
    
    file.path(
      tables_main_dir,
      "Stage17B_Temporal_Gene_Cluster_Assignments.csv"
    ),
    
    file.path(
      tables_main_dir,
      "Stage17B_Cluster_Centroids.csv"
    ),
    
    file.path(
      tables_main_dir,
      "Stage17B_Temporal_Gene_Cluster_Ranking.csv"
    ),
    
    file.path(
      tables_main_dir,
      "Stage17B_Top20_Genes_Per_Cluster.csv"
    ),
    
    file.path(
      tables_main_dir,
      "Stage17B_Cluster_Summary.csv"
    ),
    
    file.path(
      tables_supp_dir,
      "Stage17B_Gene_ZScore_Profiles.csv"
    ),
    
    file.path(
      tables_supp_dir,
      "Stage17B_Gene_ZScore_Profiles_Long.csv"
    ),
    
    file.path(
      tables_validation_dir,
      "Stage17B_K_Silhouette_Evaluation.csv"
    ),
    
    file.path(
      tables_validation_dir,
      "Stage17B_Final_Validation.csv"
    ),
    
    file.path(
      tables_validation_dir,
      "Stage17B_Structural_Validation.csv"
    ),
    
    file.path(
      figures_main_dir,
      "Stage17B_Temporal_Cluster_Profiles.pdf"
    ),
    
    file.path(
      figures_supp_dir,
      "Stage17B_K_Silhouette_Evaluation.pdf"
    ),
    
    file.path(
      figures_supp_dir,
      "Stage17B_Silhouette_Width_By_Cluster.pdf"
    ),
    
    stage17b_rds,
    
    stage17b_log
  )


output_exists <-
  file.exists(
    output_files
  )


write_log(
  "Verified output files: ",
  sum(output_exists),
  " / ",
  length(output_files)
)


# ==============================================================================
# 47. FINAL SUMMARY
# ==============================================================================

write_log(
  "======================================================================"
)

write_log(
  "STAGE 17B COMPLETED SUCCESSFULLY"
)

write_log(
  "======================================================================"
)

write_log(
  "Stage 17 strong temporal genes: ",
  length(strong_genes)
)

write_log(
  "Strong temporal genes in expression: ",
  length(strong_genes_in_expression)
)

write_log(
  "Final synchronized cells: ",
  nrow(analysis_cells)
)

write_log(
  "Pseudotime bins: ",
  n_bins
)

write_log(
  "Genes clustered: ",
  nrow(zscore_matrix)
)

write_log(
  "Candidate K values: ",
  paste(
    k_candidates,
    collapse = ", "
  )
)

write_log(
  "Selected K: ",
  selected_k
)

write_log(
  "Selected mean silhouette: ",
  format(
    selected_mean_silhouette,
    digits = 8
  )
)

write_log(
  "Minimum silhouette: ",
  format(
    min(
      final_silhouette_width
    ),
    digits = 8
  )
)

write_log(
  "Median silhouette: ",
  format(
    median(
      final_silhouette_width
    ),
    digits = 8
  )
)

write_log(
  "Maximum silhouette: ",
  format(
    max(
      final_silhouette_width
    ),
    digits = 8
  )
)

write_log(
  "Trajectory reconstruction: NOT PERFORMED"
)

write_log(
  "Root selection: NOT PERFORMED"
)

write_log(
  "Pseudotime recalculation: NOT PERFORMED"
)

write_log(
  "Stage 15 remained frozen."
)

write_log(
  "Stage 16 remained frozen."
)

write_log(
  "Stage 17 remained frozen."
)

write_log(
  "Central object: ",
  stage17b_rds
)

write_log(
  "Completion log: ",
  stage17b_log
)

write_log(
  "======================================================================"
)


# ==============================================================================
# 48. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "STAGE 17B COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  "Stage 17 strong temporal genes: ",
  length(strong_genes),
  "\n",
  "Strong temporal genes in expression: ",
  length(strong_genes_in_expression),
  "\n",
  "Final synchronized cells: ",
  nrow(analysis_cells),
  "\n",
  "Pseudotime bins: ",
  n_bins,
  "\n",
  "Genes clustered: ",
  nrow(zscore_matrix),
  "\n",
  "Selected K: ",
  selected_k,
  "\n",
  "Mean silhouette: ",
  format(
    selected_mean_silhouette,
    digits = 8
  ),
  "\n",
  "Minimum silhouette: ",
  format(
    min(
      final_silhouette_width
    ),
    digits = 8
  ),
  "\n",
  "Median silhouette: ",
  format(
    median(
      final_silhouette_width
    ),
    digits = 8
  ),
  "\n",
  "Maximum silhouette: ",
  format(
    max(
      final_silhouette_width
    ),
    digits = 8
  ),
  "\n",
  "\n",
  "Central RDS:\n",
  stage17b_rds,
  "\n",
  "\n",
  "Completion log:\n",
  stage17b_log,
  "\n",
  "============================================================\n"
)


# ==============================================================================
# END OF STAGE 17B
# ==============================================================================