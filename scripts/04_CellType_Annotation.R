############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 4: Cell-type annotation
#
# This script:
# - Loads the clustered Seurat object
# - Identifies cluster marker genes
# - Performs automated cell-type annotation using SingleR
# - Adds SingleR predictions to Seurat metadata
# - Performs marker-based biological refinement
# - Generates annotated UMAP plots
# - Saves the annotated Seurat object
#
# NOTE:
# This script represents the original Stage 4 analytical workflow.
# The current version updates project paths and output organization
# to match the canonical repository structure.
#
# No new analysis is intended during repository reorganization.
############################################################


############################################################
# Load libraries
############################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(celldex)
library(SingleR)
library(SingleCellExperiment)


############################################################
# Project directories
############################################################

# Canonical project repository
project_dir <- "YOUR_PROJECT_DIRECTORY"

# Central object directory
objects_dir <- file.path(
  project_dir,
  "objects"
)

# Metadata directory
metadata_dir <- file.path(
  project_dir,
  "metadata"
)

# Canonical Stage 4 results directory
stage_results_dir <- file.path(
  project_dir,
  "results",
  "04_CellType_Annotation"
)

# Stage 4 figure directory
figures_main_dir <- file.path(
  stage_results_dir,
  "figures",
  "main"
)

# Stage 4 table directories
tables_main_dir <- file.path(
  stage_results_dir,
  "tables",
  "main"
)

tables_supplementary_dir <- file.path(
  stage_results_dir,
  "tables",
  "supplementary"
)

# Stage 4 log directory
logs_dir <- file.path(
  stage_results_dir,
  "logs"
)


############################################################
# Create required directories if needed
############################################################

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_supplementary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Load clustered Seurat object
############################################################

srobj <- readRDS(
  file.path(
    objects_dir,
    "TNBC_PCA_Clustering_UMAP.rds"
  )
)


############################################################
# Cluster marker identification
############################################################

Idents(srobj) <- "seurat_clusters"

markers <- FindAllMarkers(
  object = srobj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


############################################################
# Save all cluster markers
#
# Supplementary table
############################################################

write.csv(
  markers,
  file = file.path(
    tables_supplementary_dir,
    "TNBC_cluster_markers_all.csv"
  ),
  row.names = FALSE
)


############################################################
# Filter strong marker genes
############################################################

strong_markers <- markers %>%
  filter(avg_log2FC >= 2)


############################################################
# Save strong cluster markers
#
# Supplementary table
############################################################

write.csv(
  strong_markers,
  file = file.path(
    tables_supplementary_dir,
    "TNBC_cluster_markers_logFC2.csv"
  ),
  row.names = FALSE
)


############################################################
# SingleR reference
############################################################

DefaultAssay(srobj) <- "RNA"


############################################################
# Join RNA layers before conversion
#
# Seurat v5 compatibility
############################################################

srobj <- JoinLayers(
  srobj,
  assay = "RNA"
)


############################################################
# Load reference dataset
############################################################

myref <- celldex::HumanPrimaryCellAtlasData(
  ensembl = FALSE
)


############################################################
# Convert Seurat object to SingleCellExperiment
############################################################

sce <- as.SingleCellExperiment(
  srobj,
  assay = "RNA"
)


############################################################
# Run SingleR
############################################################

mylabels <- SingleR(
  test = sce,
  ref = myref,
  labels = myref$label.main
)


############################################################
# Add SingleR labels to Seurat metadata
############################################################

srobj$SingleR_CellType <- mylabels$labels


############################################################
# Save SingleR predictions
#
# Main table
############################################################

write.csv(
  data.frame(
    Cell = rownames(mylabels),
    SingleR_Label = mylabels$labels
  ),
  file = file.path(
    tables_main_dir,
    "SingleR_CellType_Annotations.csv"
  ),
  row.names = FALSE
)


############################################################
# Cell-type composition
#
# Supplementary table
############################################################

celltype_table <- table(
  srobj$SingleR_CellType
)


write.csv(
  as.data.frame(celltype_table),
  file = file.path(
    tables_supplementary_dir,
    "SingleR_CellType_Abundance.csv"
  ),
  row.names = FALSE
)


############################################################
# Load manual cluster annotations
############################################################

annotation_table <- read.csv(
  file.path(
    metadata_dir,
    "cluster_annotations.csv"
  ),
  stringsAsFactors = FALSE
)


############################################################
# Create named vector:
# names = cluster IDs
############################################################

manual_annotations <-
  setNames(
    annotation_table$Manual_CellType,
    annotation_table$cluster
  )


############################################################
# Add manual annotations to Seurat metadata
############################################################

srobj$Manual_CellType <-
  as.vector(
    manual_annotations[
      as.character(srobj$seurat_clusters)
    ]
  )


############################################################
# Check annotation distribution
############################################################

table(
  srobj$Manual_CellType
)


############################################################
# UMAP: Manual cell-type annotations
############################################################

pdf(
  file = file.path(
    figures_main_dir,
    "UMAP_Manual_CellTypes_Split.pdf"
  ),
  width = 14,
  height = 8
)

DimPlot(
  srobj,
  reduction = "umap",
  group.by = "Manual_CellType",
  label = TRUE,
  split.by = "Class",
  repel = TRUE
) +
  NoLegend()

dev.off()


############################################################
# Save annotated Seurat object
#
# Central project object directory:
# objects/
############################################################

saveRDS(
  srobj,
  file = file.path(
    objects_dir,
    "TNBC_Annotated.rds"
  ),
  compress = TRUE
)

message(
  "Annotated object saved: TNBC_Annotated.rds"
)


############################################################
# Summary log
############################################################

log_file <- file.path(
  logs_dir,
  "Step4_CellType_Annotation.log"
)

log_con <- file(
  log_file,
  open = "wt"
)


writeLines(
  c(
    "TNBC Single-cell RNA-seq Analysis Pipeline",
    "Step 4: Cell-type annotation",
    paste("Date:", Sys.time()),
    "",
    "Summary",
    "-------",
    paste(
      "Number of cells:",
      ncol(srobj)
    ),
    paste(
      "Number of clusters:",
      length(levels(srobj$seurat_clusters))
    ),
    paste(
      "Number of SingleR cell types:",
      length(unique(srobj$SingleR_CellType))
    ),
    paste(
      "Number of manual cell types:",
      length(unique(srobj$Manual_CellType))
    ),
    "",
    "Cells per manual cell type",
    "--------------------------",
    capture.output(
      print(
        table(srobj$Manual_CellType)
      )
    ),
    "",
    "Cells per biological group",
    "--------------------------",
    capture.output(
      print(
        table(srobj$Class)
      )
    ),
    "",
    "Objects and outputs generated",
    "-----------------------------",
    file.path(
      tables_supplementary_dir,
      "TNBC_cluster_markers_all.csv"
    ),
    file.path(
      tables_supplementary_dir,
      "TNBC_cluster_markers_logFC2.csv"
    ),
    file.path(
      tables_main_dir,
      "SingleR_CellType_Annotations.csv"
    ),
    file.path(
      tables_supplementary_dir,
      "SingleR_CellType_Abundance.csv"
    ),
    file.path(
      figures_main_dir,
      "UMAP_Manual_CellTypes_Split.pdf"
    ),
    file.path(
      objects_dir,
      "TNBC_Annotated.rds"
    ),
    "",
    "Step 4 completed successfully."
  ),
  con = log_con
)

close(log_con)


############################################################
# Completion message
############################################################

message(
  "Step 4 completed successfully."
)

message(
  "Log saved to: ",
  log_file
)