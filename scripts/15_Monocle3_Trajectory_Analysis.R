# ============================================================
# STAGE 15
# MONOCLE3 TRAJECTORY INFERENCE AND PSEUDOTIME
# TNBC BRCA1 Single-Cell RNA-seq Project
#
# FINAL CLEAN-RERUN VERSION
#
# Purpose:
#   1. Load the final integrated epithelial + CopyKAT Seurat object
#   2. Construct a Monocle3 cell_data_set
#   3. Transfer the validated Seurat UMAP
#   4. Infer partition-aware trajectories
#   5. Identify a CNV-mixed trajectory partition
#   6. Select a Normal-like root node
#   7. Calculate and validate pseudotime
#   8. Transfer pseudotime back to Seurat
#   9. Export publication-ready figures and tables
#  10. Preserve computational objects in the central objects/
#      directory
#
# IMPORTANT:
#   - The original/source pipeline is NOT modified.
#   - Stages 1-14 are NOT rerun by this script.
#   - The final validated Seurat object from previous stages is used.
#   - This script does not create numbered 01-08 validation files.
#   - This script does not create PNG files.
#   - This script does not use sink(type="message", split=TRUE).
#   - Refactoring this script does not imply re-analysis of earlier stages.
# ============================================================


# ============================================================
# 0. LIBRARIES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(igraph)
})


# ============================================================
# 1. REPRODUCIBILITY
# ============================================================

set.seed(12345)


# ============================================================
# 2. PROJECT PATHS
# ============================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"

objects_dir <- file.path(
  project_dir,
  "objects"
)

results_dir <- file.path(
  project_dir,
  "results",
  "15_Monocle3_Trajectory"
)

figures_dir <- file.path(
  results_dir,
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
  results_dir,
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
  results_dir,
  "logs"
)


# ============================================================
# 3. CREATE REQUIRED DIRECTORIES
# ============================================================

required_dirs <- c(
  project_dir,
  objects_dir,
  results_dir,
  figures_dir,
  figures_main_dir,
  figures_supp_dir,
  tables_dir,
  tables_main_dir,
  tables_supp_dir,
  tables_validation_dir,
  logs_dir
)

for (d in required_dirs) {
  
  if (!dir.exists(d)) {
    
    dir.create(
      d,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
  }
}


# ============================================================
# 4. LOGGING
# ============================================================
#
# IMPORTANT:
#   Only standard output is logged.
#   No message sink with split=TRUE is used because R does not
#   allow splitting the message connection.
# ============================================================

log_file <- file.path(
  logs_dir,
  "Stage15_Completion.log"
)

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_message <- function(...) {
  
  txt <- paste0(
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    " | ",
    paste0(..., collapse = "")
  )
  
  cat(
    txt,
    "\n",
    file = log_file,
    append = TRUE
  )
  
  message(txt)
}

write_validation_csv <- function(
    x,
    filename
) {
  
  write.csv(
    x,
    file = file.path(
      tables_validation_dir,
      filename
    ),
    row.names = FALSE
  )
}

log_message(
  "============================================"
)

log_message(
  "STAGE 15 STARTED"
)

log_message(
  "Monocle3 Trajectory and Pseudotime Analysis"
)

log_message(
  "Final clean-rerun version"
)

log_message(
  "============================================"
)


# ============================================================
# 5. INPUT / OUTPUT OBJECT PATHS
# ============================================================

seurat_path <- file.path(
  objects_dir,
  "TNBC_Epithelial_CopyKAT.rds"
)

trajectory_seurat_path <- file.path(
  objects_dir,
  "TNBC_Epithelial_Trajectory_Added.rds"
)

monocle3_object_path <- file.path(
  objects_dir,
  "TNBC_Monocle3_CDS_Final"
)

validation_rds_path <- file.path(
  objects_dir,
  "Stage15_Pseudotime_Validation_Results.rds"
)


# ============================================================
# 6. VALIDATE INPUT
# ============================================================

if (!file.exists(seurat_path)) {
  
  stop(
    paste0(
      "Required input object not found:\n",
      seurat_path
    )
  )
}

log_message(
  "Input Seurat object found: ",
  seurat_path
)


# ============================================================
# 7. LOAD SEURAT OBJECT
# ============================================================

log_message(
  "Loading Seurat object..."
)

seurat_obj <- readRDS(
  seurat_path
)

if (!inherits(
  seurat_obj,
  "Seurat"
)) {
  
  stop(
    "Input object is not a valid Seurat object."
  )
}

log_message(
  "Seurat object loaded successfully."
)

log_message(
  "Cells: ",
  ncol(seurat_obj)
)

log_message(
  "Features: ",
  nrow(seurat_obj)
)


# ============================================================
# 8. VALIDATE REQUIRED METADATA
# ============================================================

required_metadata <- c(
  "CNV_Status"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(
    seurat_obj@meta.data
  )
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "Required metadata missing from Seurat object: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 9. NORMALIZE CNV STATUS
# ============================================================

cnv_status <- seurat_obj@meta.data$CNV_Status

if (is.list(cnv_status)) {
  
  cnv_status <- vapply(
    cnv_status,
    function(x) {
      
      if (length(x) == 0) {
        
        NA_character_
        
      } else {
        
        as.character(
          x[[1]]
        )
        
      }
      
    },
    character(1)
  )
  
} else {
  
  cnv_status <- as.character(
    cnv_status
  )
}

if (length(cnv_status) != ncol(seurat_obj)) {
  
  stop(
    "CNV_Status length does not match Seurat cell count."
  )
}

seurat_obj@meta.data$CNV_Status <- cnv_status

log_message(
  "CNV_Status normalized successfully."
)


# ============================================================
# 10. SET RNA ASSAY
# ============================================================

if (!"RNA" %in% Assays(seurat_obj)) {
  
  stop(
    "RNA assay was not found in the Seurat object."
  )
}

DefaultAssay(seurat_obj) <- "RNA"

log_message(
  "Default assay set to RNA."
)


# ============================================================
# 11. VALIDATE UMAP
# ============================================================

if (!"umap" %in% Reductions(seurat_obj)) {
  
  stop(
    "Required Seurat UMAP reduction was not found."
  )
}

seurat_umap <- Embeddings(
  seurat_obj,
  reduction = "umap"
)

if (
  is.null(seurat_umap) ||
  nrow(seurat_umap) == 0
) {
  
  stop(
    "Seurat UMAP embeddings are empty."
  )
}

log_message(
  "Seurat UMAP found: ",
  nrow(seurat_umap),
  " cells x ",
  ncol(seurat_umap),
  " dimensions."
)


# ============================================================
# 12. EXTRACT RNA COUNTS
# ============================================================

log_message(
  "Extracting RNA counts..."
)

rna_assay <- seurat_obj[["RNA"]]

if (inherits(
  rna_assay,
  "Assay5"
)) {
  
  available_layers <- Layers(
    rna_assay
  )
  
  if (!"counts" %in% available_layers) {
    
    stop(
      paste0(
        "RNA Assay5 does not contain a counts layer.\n",
        "Available layers: ",
        paste(
          available_layers,
          collapse = ", "
        )
      )
    )
  }
  
  counts_matrix <- LayerData(
    rna_assay,
    layer = "counts"
  )
  
} else {
  
  counts_matrix <- GetAssayData(
    seurat_obj,
    assay = "RNA",
    slot = "counts"
  )
}


# ============================================================
# 13. VALIDATE COUNTS MATRIX
# ============================================================

if (is.null(counts_matrix)) {
  
  stop(
    "RNA counts matrix could not be extracted."
  )
}

if (
  nrow(counts_matrix) == 0 ||
  ncol(counts_matrix) == 0
) {
  
  stop(
    "RNA counts matrix is empty."
  )
}

log_message(
  "RNA counts dimensions: ",
  nrow(counts_matrix),
  " genes x ",
  ncol(counts_matrix),
  " cells."
)


# ============================================================
# 14. ENSURE CELL ALIGNMENT
# ============================================================

seurat_cells <- colnames(
  seurat_obj
)

count_cells <- colnames(
  counts_matrix
)

if (!identical(
  seurat_cells,
  count_cells
)) {
  
  common_cells <- intersect(
    seurat_cells,
    count_cells
  )
  
  if (length(common_cells) == 0) {
    
    stop(
      "No overlapping cell barcodes between Seurat object and RNA counts."
    )
  }
  
  if (!all(
    seurat_cells %in% count_cells
  )) {
    
    stop(
      "RNA counts do not contain all Seurat cells."
    )
  }
  
  log_message(
    "Reordering RNA counts to match Seurat cell order."
  )
  
  counts_matrix <- counts_matrix[
    ,
    seurat_cells,
    drop = FALSE
  ]
}


# ============================================================
# 15. GENE METADATA
# ============================================================

gene_ids <- rownames(
  counts_matrix
)

gene_metadata <- data.frame(
  gene_short_name = gene_ids,
  row.names = gene_ids,
  stringsAsFactors = FALSE
)


# ============================================================
# 16. CELL METADATA
# ============================================================

cell_metadata <- seurat_obj@meta.data

cell_metadata <- cell_metadata[
  seurat_cells,
  ,
  drop = FALSE
]

if (!identical(
  rownames(cell_metadata),
  colnames(counts_matrix)
)) {
  
  stop(
    "Cell metadata and counts matrix are not aligned."
  )
}


# ============================================================
# 17. CREATE MONOCLE3 CDS
# ============================================================

log_message(
  "Constructing Monocle3 cell_data_set..."
)

cds <- new_cell_data_set(
  counts_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

rm(counts_matrix)

gc()

log_message(
  "Monocle3 cell_data_set created."
)


# ============================================================
# 18. PREPROCESS CDS
# ============================================================

log_message(
  "Running Monocle3 preprocessing..."
)

cds <- preprocess_cds(
  cds,
  num_dim = 50,
  method = "PCA"
)

log_message(
  "Monocle3 preprocessing completed."
)


# ============================================================
# 19. TRANSFER VALIDATED SEURAT UMAP
# ============================================================

log_message(
  "Transferring Seurat UMAP into Monocle3..."
)

umap_cells <- rownames(
  seurat_umap
)

cds_cells <- colnames(
  cds
)

if (!identical(
  umap_cells,
  cds_cells
)) {
  
  if (!setequal(
    umap_cells,
    cds_cells
  )) {
    
    stop(
      "Seurat UMAP cells and Monocle3 CDS cells do not match."
    )
  }
  
  seurat_umap <- seurat_umap[
    cds_cells,
    ,
    drop = FALSE
  ]
}

reducedDims(cds)$UMAP <- as.matrix(
  seurat_umap
)

log_message(
  "Seurat UMAP transferred successfully."
)


# ============================================================
# 20. MONOCLE3 CLUSTERING
# ============================================================

log_message(
  "Running Monocle3 clustering..."
)

cds <- cluster_cells(
  cds,
  reduction_method = "UMAP"
)

log_message(
  "Monocle3 clustering completed."
)


# ============================================================
# 21. LEARN PARTITION-AWARE GRAPH
# ============================================================

log_message(
  "Learning partition-aware principal graph..."
)

cds <- learn_graph(
  cds,
  use_partition = TRUE,
  close_loop = FALSE
)

log_message(
  "Principal graph learning completed."
)


# ============================================================
# 22. EXTRACT PARTITIONS
# ============================================================

cds_partitions <- partitions(
  cds
)

if (is.null(cds_partitions)) {
  
  stop(
    "Monocle3 partitions could not be extracted."
  )
}

cds_partitions <- as.character(
  cds_partitions
)

names(cds_partitions) <- colnames(
  cds
)

if (any(
  is.na(cds_partitions)
)) {
  
  stop(
    "Some cells have NA Monocle3 partitions."
  )
}

log_message(
  "Number of Monocle3 partitions: ",
  length(
    unique(
      cds_partitions
    )
  )
)

colData(cds)$Monocle3_Partition <-
  cds_partitions


# ============================================================
# 23. EXTRACT PRINCIPAL GRAPH PROJECTION
# ============================================================

log_message(
  "Extracting principal graph projection..."
)

principal_graph_aux_data <- principal_graph_aux(
  cds
)[["UMAP"]]

if (is.null(
  principal_graph_aux_data
)) {
  
  stop(
    "UMAP principal graph auxiliary data could not be extracted."
  )
}

pr_graph_proj <-
  principal_graph_aux_data[["pr_graph_cell_proj_closest_vertex"]]

if (is.null(pr_graph_proj)) {
  
  stop(
    "Closest principal graph vertex information not found."
  )
}

pr_graph_proj <- as.matrix(
  pr_graph_proj
)

if (
  nrow(pr_graph_proj) != ncol(cds)
) {
  
  stop(
    "Principal graph projection does not contain all cells."
  )
}

pr_graph_proj_df <- data.frame(
  cell = rownames(pr_graph_proj),
  closest_vertex = as.character(
    pr_graph_proj[, 1]
  ),
  stringsAsFactors = FALSE
)

rownames(pr_graph_proj_df) <-
  pr_graph_proj_df$cell


# ============================================================
# 24. EXTRACT ACTUAL GRAPH VERTEX NAMES
# ============================================================

principal_graph_list <- principal_graph(
  cds
)

if (!"UMAP" %in% names(
  principal_graph_list
)) {
  
  stop(
    "UMAP principal graph was not found."
  )
}

graph_obj <- principal_graph_list[["UMAP"]]

real_vertex_names <- V(
  graph_obj
)$name

if (length(real_vertex_names) == 0) {
  
  stop(
    "No principal graph vertices were identified."
  )
}

log_message(
  "Principal graph vertices: ",
  length(real_vertex_names)
)


# ============================================================
# 25. ROBUST VERTEX NAME RESOLUTION
# ============================================================

projection_vertices <- unique(
  pr_graph_proj_df$closest_vertex
)

if (all(
  projection_vertices %in%
  real_vertex_names
)) {
  
  resolved_vertices <-
    pr_graph_proj_df$closest_vertex
  
} else {
  
  numeric_projection <- suppressWarnings(
    as.integer(
      pr_graph_proj_df$closest_vertex
    )
  )
  
  if (
    all(
      !is.na(numeric_projection) &
      numeric_projection >= 1 &
      numeric_projection <=
      length(real_vertex_names)
    )
  ) {
    
    resolved_vertices <-
      real_vertex_names[
        numeric_projection
      ]
    
  } else {
    
    stop(
      paste0(
        "Unable to resolve principal graph vertex identifiers.\n",
        "Example projected values: ",
        paste(
          head(
            projection_vertices,
            10
          ),
          collapse = ", "
        )
      )
    )
  }
}

pr_graph_proj_df$vertex <-
  resolved_vertices


# ============================================================
# 26. GRAPH PROJECTION VALIDATION
# ============================================================

projection_validation <- data.frame(
  Metric = c(
    "Total_CDS_Cells",
    "Projected_Cells",
    "Missing_Projection",
    "Unique_Graph_Vertices_Used"
  ),
  Value = c(
    ncol(cds),
    nrow(pr_graph_proj_df),
    sum(
      !colnames(cds) %in%
        rownames(pr_graph_proj_df)
    ),
    length(
      unique(
        pr_graph_proj_df$vertex
      )
    )
  )
)

write_validation_csv(
  projection_validation,
  "Step15_Graph_Projection_Validation.csv"
)

if (!all(
  colnames(cds) %in%
  rownames(pr_graph_proj_df)
)) {
  
  stop(
    "Not all CDS cells have graph projection information."
  )
}


# ============================================================
# 27. BUILD CELL-LEVEL TRAJECTORY TABLE
# ============================================================

cell_trajectory <- data.frame(
  cell = colnames(cds),
  Partition = cds_partitions[
    colnames(cds)
  ],
  CNV_Status = as.character(
    colData(cds)$CNV_Status
  ),
  Graph_Vertex = pr_graph_proj_df[
    colnames(cds),
    "vertex"
  ],
  stringsAsFactors = FALSE
)

rownames(cell_trajectory) <-
  cell_trajectory$cell


# ============================================================
# 28. PARTITION-LEVEL CNV COMPOSITION
# ============================================================

partition_cnv_stats <-
  cell_trajectory %>%
  dplyr::filter(
    !is.na(CNV_Status)
  ) %>%
  dplyr::group_by(
    Partition
  ) %>%
  dplyr::summarise(
    Total_Cells = dplyr::n(),
    
    Normal_like = sum(
      CNV_Status == "Normal_like"
    ),
    
    Tumor_like = sum(
      CNV_Status == "Tumor_like"
    ),
    
    Other_CNV = sum(
      !CNV_Status %in%
        c(
          "Normal_like",
          "Tumor_like"
        )
    ),
    
    Normal_Fraction = ifelse(
      Total_Cells > 0,
      Normal_like / Total_Cells,
      NA_real_
    ),
    
    Tumor_Fraction = ifelse(
      Total_Cells > 0,
      Tumor_like / Total_Cells,
      NA_real_
    ),
    
    Mixed_Cells =
      Normal_like + Tumor_like,
    
    Mixed_Fraction = ifelse(
      Total_Cells > 0,
      2 * pmin(
        Normal_Fraction,
        Tumor_Fraction
      ),
      0
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Mixed_Score =
      Mixed_Cells * Mixed_Fraction
  ) %>%
  dplyr::arrange(
    desc(Mixed_Score)
  )

write.csv(
  partition_cnv_stats,
  file = file.path(
    tables_supp_dir,
    "Partition_CNV_Statistics.csv"
  ),
  row.names = FALSE
)

log_message(
  "Partition CNV statistics calculated."
)


# ============================================================
# 29. IDENTIFY ELIGIBLE MIXED PARTITIONS
# ============================================================

minimum_cells_per_state <- 50

partition_candidates <-
  partition_cnv_stats %>%
  dplyr::filter(
    Normal_like >= minimum_cells_per_state,
    Tumor_like >= minimum_cells_per_state
  ) %>%
  dplyr::arrange(
    desc(Mixed_Score),
    desc(Mixed_Cells)
  )

write.csv(
  partition_candidates,
  file = file.path(
    tables_supp_dir,
    "Trajectory_Partition_Candidates.csv"
  ),
  row.names = FALSE
)

if (nrow(partition_candidates) == 0) {
  
  stop(
    paste0(
      "No trajectory partition satisfies the minimum requirement of ",
      minimum_cells_per_state,
      " Normal_like and ",
      minimum_cells_per_state,
      " Tumor_like cells."
    )
  )
}

selected_partition <-
  partition_candidates$Partition[1]

selected_partition_stats <-
  partition_candidates[
    1,
    ,
    drop = FALSE
  ]

log_message(
  "Selected trajectory partition: ",
  selected_partition
)

log_message(
  "Normal_like cells: ",
  selected_partition_stats$Normal_like
)

log_message(
  "Tumor_like cells: ",
  selected_partition_stats$Tumor_like
)


# ============================================================
# 30. GRAPH NODE-LEVEL CNV COMPOSITION
# ============================================================

selected_partition_cells <-
  cell_trajectory %>%
  dplyr::filter(
    Partition == selected_partition
  )

node_stats <-
  selected_partition_cells %>%
  dplyr::group_by(
    Graph_Vertex
  ) %>%
  dplyr::summarise(
    Total_Cells = dplyr::n(),
    
    Normal_like = sum(
      CNV_Status == "Normal_like",
      na.rm = TRUE
    ),
    
    Tumor_like = sum(
      CNV_Status == "Tumor_like",
      na.rm = TRUE
    ),
    
    Normal_Fraction = ifelse(
      Total_Cells > 0,
      Normal_like / Total_Cells,
      0
    ),
    
    Tumor_Fraction = ifelse(
      Total_Cells > 0,
      Tumor_like / Total_Cells,
      0
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Mixed_Cells =
      Normal_like + Tumor_like,
    
    Mixed_Fraction =
      ifelse(
        Mixed_Cells > 0,
        2 * pmin(
          Normal_Fraction,
          Tumor_Fraction
        ),
        0
      )
  )

if (nrow(node_stats) == 0) {
  
  stop(
    "No graph nodes were found in the selected partition."
  )
}


# ============================================================
# 31. IDENTIFY LEAF NODES
# ============================================================

selected_graph <- graph_obj

selected_graph_vertices <- V(
  selected_graph
)$name

graph_degree <- degree(
  selected_graph,
  mode = "all"
)

leaf_vertex_names <-
  selected_graph_vertices[
    graph_degree == 1
  ]

leaf_node_stats <-
  node_stats %>%
  dplyr::filter(
    Graph_Vertex %in%
      leaf_vertex_names
  ) %>%
  dplyr::mutate(
    Root_Score =
      Normal_Fraction *
      log1p(
        Normal_like
      )
  ) %>%
  dplyr::arrange(
    desc(Root_Score),
    desc(Normal_like)
  )

write.csv(
  leaf_node_stats,
  file = file.path(
    tables_supp_dir,
    "Root_Node_Candidates.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 32. SELECT NORMAL-LIKE ROOT NODE
# ============================================================

if (nrow(leaf_node_stats) > 0) {
  
  best_root_node <-
    leaf_node_stats$Graph_Vertex[1]
  
  selected_root_stats <-
    leaf_node_stats[
      1,
      ,
      drop = FALSE
    ]
  
  root_selection_method <-
    "Normal_like_leaf_node"
  
} else {
  
  node_stats_fallback <-
    node_stats %>%
    dplyr::mutate(
      Root_Score =
        Normal_Fraction *
        log1p(
          Normal_like
        )
    ) %>%
    dplyr::arrange(
      desc(Root_Score),
      desc(Normal_like)
    )
  
  if (nrow(node_stats_fallback) == 0) {
    
    stop(
      "No valid graph node available for root selection."
    )
  }
  
  best_root_node <-
    node_stats_fallback$Graph_Vertex[1]
  
  selected_root_stats <-
    node_stats_fallback[
      1,
      ,
      drop = FALSE
    ]
  
  root_selection_method <-
    "Normal_like_graph_node_fallback"
}


# ============================================================
# 33. VALIDATE ROOT NODE
# ============================================================

if (
  length(best_root_node) != 1 ||
  is.na(best_root_node) ||
  !best_root_node %in%
  selected_graph_vertices
) {
  
  stop(
    "Selected root node is not present in the principal graph."
  )
}

if (
  selected_root_stats$Normal_like[1] <= 0
) {
  
  stop(
    "Selected root node contains no Normal_like cells."
  )
}

selected_root_table <- data.frame(
  Selected_Partition =
    selected_partition,
  
  Selected_Root_Node =
    best_root_node,
  
  Root_Selection_Method =
    root_selection_method,
  
  Root_Total_Cells =
    selected_root_stats$Total_Cells[1],
  
  Root_Normal_like =
    selected_root_stats$Normal_like[1],
  
  Root_Tumor_like =
    selected_root_stats$Tumor_like[1],
  
  Root_Normal_Fraction =
    selected_root_stats$Normal_Fraction[1],
  
  Root_Tumor_Fraction =
    selected_root_stats$Tumor_Fraction[1],
  
  Root_Score =
    selected_root_stats$Root_Score[1],
  
  stringsAsFactors = FALSE
)

write.csv(
  selected_root_table,
  file = file.path(
    tables_supp_dir,
    "Selected_Trajectory_Root.csv"
  ),
  row.names = FALSE
)

log_message(
  "Selected root node: ",
  best_root_node
)

log_message(
  "Root selection method: ",
  root_selection_method
)


# ============================================================
# 34. ORDER CELLS AND CALCULATE PSEUDOTIME
# ============================================================

log_message(
  "Ordering cells from selected Normal-like root..."
)

cds <- order_cells(
  cds,
  root_pr_nodes = best_root_node
)

log_message(
  "Cell ordering completed."
)


# ============================================================
# 35. EXTRACT PSEUDOTIME
# ============================================================

cds_pseudotime <- pseudotime(
  cds
)

names(cds_pseudotime) <-
  colnames(cds)

if (
  length(cds_pseudotime) !=
  ncol(cds)
) {
  
  stop(
    "Pseudotime vector length does not match CDS cell count."
  )
}


# ============================================================
# 36. GLOBAL PSEUDOTIME VALIDATION
# ============================================================

finite_pt <- is.finite(
  cds_pseudotime
)

infinite_pt <- is.infinite(
  cds_pseudotime
)

na_pt <- is.na(
  cds_pseudotime
)

n_finite <- sum(
  finite_pt
)

n_infinite <- sum(
  infinite_pt
)

n_na <- sum(
  na_pt
)

finite_fraction <-
  n_finite /
  length(cds_pseudotime)

global_pseudotime_validation <- data.frame(
  Metric = c(
    "Total_Cells",
    "Finite_Pseudotime",
    "Infinite_Pseudotime",
    "NA_Pseudotime",
    "Finite_Fraction"
  ),
  
  Value = c(
    length(cds_pseudotime),
    n_finite,
    n_infinite,
    n_na,
    finite_fraction
  )
)

write_validation_csv(
  global_pseudotime_validation,
  "Pseudotime_Global_Validation.csv"
)

if (n_finite == 0) {
  
  stop(
    "No finite pseudotime values were generated."
  )
}


# ============================================================
# 37. CELL-LEVEL PSEUDOTIME METADATA
# ============================================================

pseudotime_metadata <- data.frame(
  Cell = colnames(cds),
  
  Partition = cds_partitions[
    colnames(cds)
  ],
  
  CNV_Status = as.character(
    colData(cds)$CNV_Status
  ),
  
  Graph_Vertex = pr_graph_proj_df[
    colnames(cds),
    "vertex"
  ],
  
  Pseudotime = as.numeric(
    cds_pseudotime[
      colnames(cds)
    ]
  ),
  
  Pseudotime_Finite = is.finite(
    cds_pseudotime[
      colnames(cds)
    ]
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 38. SELECTED-PARTITION PSEUDOTIME VALIDATION
# ============================================================

selected_partition_pt <-
  pseudotime_metadata %>%
  dplyr::filter(
    Partition == selected_partition
  )

selected_partition_finite <-
  selected_partition_pt %>%
  dplyr::filter(
    Pseudotime_Finite
  )

selected_partition_validation <- data.frame(
  Metric = c(
    "Selected_Partition_Total_Cells",
    "Selected_Partition_Finite_Pseudotime",
    "Selected_Partition_Infinite_Pseudotime",
    "Selected_Partition_NA_Pseudotime",
    "Finite_Fraction"
  ),
  
  Value = c(
    nrow(selected_partition_pt),
    
    sum(
      is.finite(
        selected_partition_pt$Pseudotime
      )
    ),
    
    sum(
      is.infinite(
        selected_partition_pt$Pseudotime
      )
    ),
    
    sum(
      is.na(
        selected_partition_pt$Pseudotime
      )
    ),
    
    mean(
      selected_partition_pt$Pseudotime_Finite
    )
  )
)

write_validation_csv(
  selected_partition_validation,
  "Pseudotime_Selected_Partition_Validation.csv"
)


# ============================================================
# 39. CNV-STATE FINITE PSEUDOTIME VALIDATION
# ============================================================

state_finite_counts <-
  selected_partition_finite %>%
  dplyr::filter(
    CNV_Status %in%
      c(
        "Normal_like",
        "Tumor_like"
      )
  ) %>%
  dplyr::count(
    CNV_Status,
    name = "Finite_Pseudotime_Cells"
  )

required_states <- data.frame(
  CNV_Status = c(
    "Normal_like",
    "Tumor_like"
  ),
  stringsAsFactors = FALSE
)

state_finite_counts <-
  required_states %>%
  dplyr::left_join(
    state_finite_counts,
    by = "CNV_Status"
  ) %>%
  dplyr::mutate(
    Finite_Pseudotime_Cells =
      ifelse(
        is.na(
          Finite_Pseudotime_Cells
        ),
        0,
        Finite_Pseudotime_Cells
      )
  )

write_validation_csv(
  state_finite_counts,
  "Pseudotime_Finite_CNV_State_Counts.csv"
)


# ============================================================
# 40. REQUIRE BOTH CNV STATES
# ============================================================

normal_finite_n <-
  state_finite_counts$Finite_Pseudotime_Cells[
    state_finite_counts$CNV_Status ==
      "Normal_like"
  ]

tumor_finite_n <-
  state_finite_counts$Finite_Pseudotime_Cells[
    state_finite_counts$CNV_Status ==
      "Tumor_like"
  ]

if (
  length(normal_finite_n) != 1 ||
  length(tumor_finite_n) != 1
) {
  
  stop(
    "Could not resolve finite pseudotime counts for both CNV states."
  )
}

if (
  normal_finite_n < 1 ||
  tumor_finite_n < 1
) {
  
  stop(
    paste0(
      "Selected trajectory does not contain finite pseudotime cells ",
      "from both Normal_like and Tumor_like states.\n",
      "Normal_like finite cells: ",
      normal_finite_n,
      "\nTumor_like finite cells: ",
      tumor_finite_n
    )
  )
}


# ============================================================
# 41. PSEUDOTIME SUMMARY BY CNV STATUS
# ============================================================

pseudotime_by_cnv <-
  selected_partition_finite %>%
  dplyr::filter(
    CNV_Status %in%
      c(
        "Normal_like",
        "Tumor_like"
      )
  ) %>%
  dplyr::group_by(
    CNV_Status
  ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    
    Min = min(
      Pseudotime,
      na.rm = TRUE
    ),
    
    Q1 = quantile(
      Pseudotime,
      0.25,
      na.rm = TRUE
    ),
    
    Median = median(
      Pseudotime,
      na.rm = TRUE
    ),
    
    Mean = mean(
      Pseudotime,
      na.rm = TRUE
    ),
    
    Q3 = quantile(
      Pseudotime,
      0.75,
      na.rm = TRUE
    ),
    
    Max = max(
      Pseudotime,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

write.csv(
  pseudotime_by_cnv,
  file = file.path(
    tables_main_dir,
    "Pseudotime_Summary_by_CNV_Status.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 42. WILCOXON TEST
# ============================================================

wilcox_data <-
  selected_partition_finite %>%
  dplyr::filter(
    CNV_Status %in%
      c(
        "Normal_like",
        "Tumor_like"
      )
  )

wilcox_result <- wilcox.test(
  Pseudotime ~ CNV_Status,
  data = wilcox_data,
  exact = FALSE
)

wilcox_summary <- data.frame(
  Test =
    "Wilcoxon rank-sum test",
  
  Normal_like_N =
    sum(
      wilcox_data$CNV_Status ==
        "Normal_like"
    ),
  
  Tumor_like_N =
    sum(
      wilcox_data$CNV_Status ==
        "Tumor_like"
    ),
  
  W =
    unname(
      wilcox_result$statistic
    ),
  
  P_value =
    wilcox_result$p.value,
  
  stringsAsFactors = FALSE
)

write_validation_csv(
  wilcox_summary,
  "Pseudotime_Wilcoxon_Normal_vs_Tumor.csv"
)


# ============================================================
# 43. PSEUDOTIME FINITENESS BY PARTITION
# ============================================================

pseudotime_partition_summary <-
  pseudotime_metadata %>%
  dplyr::group_by(
    Partition
  ) %>%
  dplyr::summarise(
    Total_Cells = dplyr::n(),
    
    Finite_Pseudotime = sum(
      Pseudotime_Finite
    ),
    
    Infinite_Pseudotime = sum(
      is.infinite(Pseudotime)
    ),
    
    NA_Pseudotime = sum(
      is.na(Pseudotime)
    ),
    
    Finite_Fraction = mean(
      Pseudotime_Finite
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    Partition
  )

write_validation_csv(
  pseudotime_partition_summary,
  "Pseudotime_Finiteness_by_Partition.csv"
)


# ============================================================
# 44. PSEUDOTIME DECILE TRANSITION ANALYSIS
# ============================================================

transition_data <-
  selected_partition_finite %>%
  dplyr::filter(
    CNV_Status %in%
      c(
        "Normal_like",
        "Tumor_like"
      )
  ) %>%
  dplyr::mutate(
    Pseudotime_Decile =
      dplyr::ntile(
        Pseudotime,
        10
      )
  )

transition_summary <-
  transition_data %>%
  dplyr::group_by(
    Pseudotime_Decile,
    CNV_Status
  ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    Pseudotime_Decile
  ) %>%
  dplyr::mutate(
    Fraction = N / sum(N)
  ) %>%
  dplyr::ungroup()

write_validation_csv(
  transition_summary,
  "Pseudotime_Decile_CNV_Transition.csv"
)


# ============================================================
# 45. EXPORT COMPLETE PSEUDOTIME METADATA
# ============================================================

write.csv(
  pseudotime_metadata,
  file = file.path(
    tables_main_dir,
    "Monocle3_Pseudotime_Metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 46. TRAJECTORY FIGURE PREPARATION
# ============================================================

log_message(
  "Preparing trajectory plotting object..."
)

cds_plot <- cds

cds_plot_metadata <-
  as.data.frame(
    SummarizedExperiment::colData(
      cds_plot
    )
  )

metadata_to_remove <- intersect(
  c(
    "sample",
    "sample_name"
  ),
  colnames(cds_plot_metadata)
)

if (length(metadata_to_remove) > 0) {
  
  cds_plot_metadata[
    metadata_to_remove
  ] <- NULL
}

SummarizedExperiment::colData(
  cds_plot
) <-
  S4Vectors::DataFrame(
    cds_plot_metadata
  )


# ============================================================
# 47. TRAJECTORY FIGURE: CNV STATUS
# ============================================================

log_message(
  "Generating trajectory CNV-status figure..."
)

trajectory_cnv_plot <-
  plot_cells(
    cds_plot,
    reduction_method = "UMAP",
    color_cells_by = "CNV_Status",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    show_trajectory_graph = TRUE
  ) +
  ggtitle(
    "Monocle3 Trajectory by CNV Status"
  )

print(
  trajectory_cnv_plot
)

ggsave(
  filename = file.path(
    figures_main_dir,
    "Trajectory_CNV_Status_Final.pdf"
  ),
  plot = trajectory_cnv_plot,
  width = 10,
  height = 8,
  units = "in"
)

log_message(
  "Trajectory CNV-status PDF saved."
)


# ============================================================
# 48. TRAJECTORY FIGURE: PSEUDOTIME
# ============================================================

log_message(
  "Generating trajectory pseudotime figure..."
)

trajectory_pseudotime_plot <-
  plot_cells(
    cds_plot,
    reduction_method = "UMAP",
    color_cells_by = "pseudotime",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    show_trajectory_graph = TRUE
  ) +
  ggtitle(
    "Monocle3 Trajectory by Pseudotime"
  )

print(
  trajectory_pseudotime_plot
)

ggsave(
  filename = file.path(
    figures_main_dir,
    "Trajectory_Pseudotime_Final.pdf"
  ),
  plot = trajectory_pseudotime_plot,
  width = 10,
  height = 8,
  units = "in"
)

log_message(
  "Trajectory pseudotime PDF saved."
)


# ============================================================
# 49. TRANSFER PSEUDOTIME BACK TO SEURAT
# ============================================================

log_message(
  "Transferring pseudotime back to Seurat metadata..."
)

seurat_obj$Monocle3_Pseudotime <-
  as.numeric(
    cds_pseudotime[
      colnames(seurat_obj)
    ]
  )

seurat_obj$Monocle3_Partition <-
  as.character(
    cds_partitions[
      colnames(seurat_obj)
    ]
  )

seurat_obj$Monocle3_Trajectory_Root <-
  best_root_node

seurat_obj$Monocle3_Trajectory_Partition <-
  cds_partitions[
    colnames(seurat_obj)
  ] ==
  selected_partition


# ============================================================
# 50. VALIDATE TRANSFERRED PSEUDOTIME
# ============================================================

if (
  !"Monocle3_Pseudotime" %in%
  colnames(
    seurat_obj@meta.data
  )
) {
  
  stop(
    "Monocle3 pseudotime was not added to Seurat metadata."
  )
}

if (
  length(
    seurat_obj$Monocle3_Pseudotime
  ) != ncol(seurat_obj)
) {
  
  stop(
    "Transferred pseudotime length does not match Seurat cell count."
  )
}

transfer_validation <- data.frame(
  Metric = c(
    "Seurat_Cells",
    "Transferred_Pseudotime_Values",
    "Transferred_Finite",
    "Transferred_Inf",
    "Transferred_NA"
  ),
  
  Value = c(
    ncol(seurat_obj),
    
    length(
      seurat_obj$Monocle3_Pseudotime
    ),
    
    sum(
      is.finite(
        seurat_obj$Monocle3_Pseudotime
      )
    ),
    
    sum(
      is.infinite(
        seurat_obj$Monocle3_Pseudotime
      )
    ),
    
    sum(
      is.na(
        seurat_obj$Monocle3_Pseudotime
      )
    )
  )
)

write_validation_csv(
  transfer_validation,
  "Pseudotime_Seurat_Transfer_Validation.csv"
)


# ============================================================
# 51. SAVE UPDATED SEURAT OBJECT
# ============================================================

log_message(
  "Saving updated Seurat trajectory object..."
)

saveRDS(
  seurat_obj,
  file = trajectory_seurat_path
)

if (!file.exists(
  trajectory_seurat_path
)) {
  
  stop(
    "Updated Seurat object was not successfully saved."
  )
}

log_message(
  "Updated Seurat object saved: ",
  trajectory_seurat_path
)


# ============================================================
# 52. SAVE COMPLETE MONOCLE3 OBJECT
# ============================================================

log_message(
  "Saving complete Monocle3 CDS..."
)

if (dir.exists(
  monocle3_object_path
)) {
  
  unlink(
    monocle3_object_path,
    recursive = TRUE,
    force = TRUE
  )
}

save_monocle_objects(
  cds,
  monocle3_object_path
)

if (!dir.exists(
  monocle3_object_path
)) {
  
  stop(
    "Monocle3 object directory was not created."
  )
}

log_message(
  "Complete Monocle3 CDS saved: ",
  monocle3_object_path
)


# ============================================================
# 53. SAVE VALIDATION RESULTS AS CENTRAL RDS
# ============================================================

validation_results <- list(
  
  selected_partition =
    selected_partition,
  
  selected_partition_stats =
    selected_partition_stats,
  
  selected_root_node =
    best_root_node,
  
  selected_root_stats =
    selected_root_stats,
  
  root_selection_method =
    root_selection_method,
  
  partition_cnv_stats =
    partition_cnv_stats,
  
  partition_candidates =
    partition_candidates,
  
  node_stats =
    node_stats,
  
  root_candidates =
    leaf_node_stats,
  
  global_pseudotime_validation =
    global_pseudotime_validation,
  
  selected_partition_validation =
    selected_partition_validation,
  
  state_finite_counts =
    state_finite_counts,
  
  pseudotime_by_cnv =
    pseudotime_by_cnv,
  
  wilcox_summary =
    wilcox_summary,
  
  pseudotime_partition_summary =
    pseudotime_partition_summary,
  
  transition_summary =
    transition_summary
)

saveRDS(
  validation_results,
  file = validation_rds_path
)

if (!file.exists(
  validation_rds_path
)) {
  
  stop(
    "Validation results RDS was not successfully saved."
  )
}

log_message(
  "Validation results RDS saved: ",
  validation_rds_path
)


# ============================================================
# 54. STAGE 15 SUMMARY
# ============================================================

stage15_summary <- data.frame(
  
  Metric = c(
    
    "Input_Cells",
    
    "Input_Genes",
    
    "Monocle3_Partitions",
    
    "Selected_Partition",
    
    "Selected_Partition_Total_Cells",
    
    "Selected_Partition_Normal_like",
    
    "Selected_Partition_Tumor_like",
    
    "Selected_Root_Node",
    
    "Root_Selection_Method",
    
    "Global_Finite_Pseudotime",
    
    "Global_Infinite_Pseudotime",
    
    "Global_NA_Pseudotime",
    
    "Selected_Normal_like_Finite_Pseudotime",
    
    "Selected_Tumor_like_Finite_Pseudotime",
    
    "Wilcoxon_P_value"
    
  ),
  
  Value = c(
    
    ncol(seurat_obj),
    
    nrow(seurat_obj),
    
    length(
      unique(
        cds_partitions
      )
    ),
    
    selected_partition,
    
    selected_partition_stats$Total_Cells[1],
    
    selected_partition_stats$Normal_like[1],
    
    selected_partition_stats$Tumor_like[1],
    
    best_root_node,
    
    root_selection_method,
    
    n_finite,
    
    n_infinite,
    
    n_na,
    
    normal_finite_n,
    
    tumor_finite_n,
    
    wilcox_result$p.value
    
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  stage15_summary,
  file = file.path(
    tables_main_dir,
    "Stage15_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 55. OUTPUT FILE VALIDATION
# ============================================================
#
# Canonical Stage 15 output inventory.
#
# No numbered 01-08 validation files.
# No duplicate Partition_CNV_Statistics files.
# No PNG files.
# ============================================================

expected_output_paths <- c(
  
  # ----------------------------------------------------------
  # Main figures
  # ----------------------------------------------------------
  
  file.path(
    figures_main_dir,
    "Trajectory_CNV_Status_Final.pdf"
  ),
  
  file.path(
    figures_main_dir,
    "Trajectory_Pseudotime_Final.pdf"
  ),
  
  # ----------------------------------------------------------
  # Main tables
  # ----------------------------------------------------------
  
  file.path(
    tables_main_dir,
    "Monocle3_Pseudotime_Metadata.csv"
  ),
  
  file.path(
    tables_main_dir,
    "Pseudotime_Summary_by_CNV_Status.csv"
  ),
  
  file.path(
    tables_main_dir,
    "Stage15_Summary.csv"
  ),
  
  # ----------------------------------------------------------
  # Supplementary tables
  # ----------------------------------------------------------
  
  file.path(
    tables_supp_dir,
    "Partition_CNV_Statistics.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "Trajectory_Partition_Candidates.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "Root_Node_Candidates.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "Selected_Trajectory_Root.csv"
  ),
  
  # ----------------------------------------------------------
  # Validation tables
  # ----------------------------------------------------------
  
  file.path(
    tables_validation_dir,
    "Step15_Graph_Projection_Validation.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Pseudotime_Global_Validation.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Pseudotime_Selected_Partition_Validation.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Pseudotime_Finite_CNV_State_Counts.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Pseudotime_Wilcoxon_Normal_vs_Tumor.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Pseudotime_Finiteness_by_Partition.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Pseudotime_Decile_CNV_Transition.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Pseudotime_Seurat_Transfer_Validation.csv"
  )
  
)

output_exists <- file.exists(
  expected_output_paths
)

output_validation <- data.frame(
  
  Expected_Output =
    basename(expected_output_paths),
  
  Full_Path =
    expected_output_paths,
  
  Exists =
    output_exists,
  
  stringsAsFactors = FALSE
)

write.csv(
  output_validation,
  file = file.path(
    tables_validation_dir,
    "Step15_Output_Validation.csv"
  ),
  row.names = FALSE
)

missing_outputs <-
  expected_output_paths[
    !output_exists
  ]

if (length(missing_outputs) > 0) {
  
  log_message(
    "Output validation failed."
  )
  
  print(
    basename(missing_outputs)
  )
  
  stop(
    paste0(
      "Stage 15 output validation failed. Missing outputs: ",
      paste(
        basename(missing_outputs),
        collapse = ", "
      )
    )
  )
}

log_message(
  "All expected Stage 15 output files are present."
)


# ============================================================
# 56. FINAL BIOLOGICAL / COMPUTATIONAL VALIDATION
# ============================================================

final_checks <- c(
  
  "Input_Seurat_Object",
  
  "RNA_Counts_Available",
  
  "UMAP_Available",
  
  "Monocle3_CDS_Created",
  
  "Partition_Inference",
  
  "Selected_Mixed_Partition",
  
  "Normal_like_Root_Selected",
  
  "Finite_Pseudotime_Exists",
  
  "Finite_Normal_like_Pseudotime_Exists",
  
  "Finite_Tumor_like_Pseudotime_Exists",
  
  "Wilcoxon_Test_Completed",
  
  "Seurat_Pseudotime_Transfer",
  
  "Updated_Seurat_Object_Saved",
  
  "Monocle3_CDS_Saved",
  
  "Output_Files_Validated"
  
)

final_status <- c(
  
  file.exists(
    seurat_path
  ),
  
  exists(
    "cds"
  ) &&
    nrow(
      SingleCellExperiment::counts(cds)
    ) > 0,
  
  "umap" %in%
    Reductions(seurat_obj),
  
  TRUE,
  
  length(
    unique(
      cds_partitions
    )
  ) > 0,
  
  nrow(
    partition_candidates
  ) > 0,
  
  !is.na(
    best_root_node
  ),
  
  n_finite > 0,
  
  normal_finite_n > 0,
  
  tumor_finite_n > 0,
  
  TRUE,
  
  "Monocle3_Pseudotime" %in%
    colnames(
      seurat_obj@meta.data
    ),
  
  file.exists(
    trajectory_seurat_path
  ),
  
  dir.exists(
    monocle3_object_path
  ),
  
  all(
    output_validation$Exists
  )
  
)

if (
  length(final_checks) !=
  length(final_status)
) {
  
  stop(
    "Internal validation error: final validation checks and statuses have different lengths."
  )
}

final_validation_table <- data.frame(
  
  Check =
    final_checks,
  
  Status =
    final_status,
  
  stringsAsFactors = FALSE
)

write.csv(
  final_validation_table,
  file = file.path(
    tables_validation_dir,
    "Step15_Final_Validation.csv"
  ),
  row.names = FALSE
)

failed_checks <-
  final_validation_table$Check[
    !final_validation_table$Status
  ]

if (length(failed_checks) > 0) {
  
  log_message(
    "Final validation failed."
  )
  
  stop(
    paste0(
      "Stage 15 final validation failed:\n",
      paste(
        failed_checks,
        collapse = "\n"
      )
    )
  )
}

log_message(
  "Final validation: PASSED."
)


# ============================================================
# 57. FINAL LOG
# ============================================================

log_message(
  "============================================"
)

log_message(
  "STAGE 15 COMPLETED SUCCESSFULLY"
)

log_message(
  "============================================"
)

log_message(
  "Selected partition: ",
  selected_partition
)

log_message(
  "Selected root node: ",
  best_root_node
)

log_message(
  "Root selection method: ",
  root_selection_method
)

log_message(
  "Finite pseudotime cells: ",
  n_finite
)

log_message(
  "Infinite pseudotime cells: ",
  n_infinite
)

log_message(
  "NA pseudotime cells: ",
  n_na
)

log_message(
  "Finite Normal_like cells in selected partition: ",
  normal_finite_n
)

log_message(
  "Finite Tumor_like cells in selected partition: ",
  tumor_finite_n
)

log_message(
  "Wilcoxon P value: ",
  format(
    wilcox_result$p.value,
    scientific = TRUE
  )
)

log_message(
  "Updated Seurat object: ",
  trajectory_seurat_path
)

log_message(
  "Monocle3 CDS: ",
  monocle3_object_path
)

log_message(
  "Validation results RDS: ",
  validation_rds_path
)

log_message(
  "Stage 15 result directory: ",
  results_dir
)

log_message(
  "Final validation: PASSED"
)

log_message(
  "============================================"
)


# ============================================================
# 58. CONSOLE SUMMARY
# ============================================================

cat("\n")

cat(
  "============================================\n"
)

cat(
  "STAGE 15 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================\n"
)

cat(
  "Selected partition: ",
  selected_partition,
  "\n",
  sep = ""
)

cat(
  "Selected root node: ",
  best_root_node,
  "\n",
  sep = ""
)

cat(
  "Root selection method: ",
  root_selection_method,
  "\n",
  sep = ""
)

cat(
  "Finite pseudotime cells: ",
  n_finite,
  "\n",
  sep = ""
)

cat(
  "Infinite pseudotime cells: ",
  n_infinite,
  "\n",
  sep = ""
)

cat(
  "NA pseudotime cells: ",
  n_na,
  "\n",
  sep = ""
)

cat(
  "Finite Normal_like cells: ",
  normal_finite_n,
  "\n",
  sep = ""
)

cat(
  "Finite Tumor_like cells: ",
  tumor_finite_n,
  "\n",
  sep = ""
)

cat(
  "Wilcoxon P value: ",
  format(
    wilcox_result$p.value,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Final validation: PASSED\n"
)

cat(
  "Results: ",
  results_dir,
  "\n",
  sep = ""
)

cat(
  "============================================\n"
)

