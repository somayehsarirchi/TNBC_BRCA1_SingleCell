############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 3: Cell-cycle scoring, PCA, clustering and UMAP
#
# This script:
# - Loads the integrated Seurat object
# - Scores cell-cycle phases
# - Performs diagnostic PCA using cell-cycle genes
# - Regresses cell-cycle effects
# - Performs PCA on the integrated assay
# - Determines the number of PCs
# - Performs graph-based clustering
# - Computes UMAP
# - Generates diagnostic plots
# - Saves the processed Seurat object
#
# NOTE:
# This script represents the original Stage 3 analytical workflow.
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

# Canonical Stage 3 results directory
stage_results_dir <- file.path(
  project_dir,
  "results",
  "03_PCA_Clustering_UMAP"
)

# Stage 3 figure directory
figures_main_dir <- file.path(
  stage_results_dir,
  "figures",
  "main"
)

# Stage 3 log directory
logs_dir <- file.path(
  stage_results_dir,
  "logs"
)

# Create required directories if needed
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
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Load integrated Seurat object
############################################################

srobj <- readRDS(
  file.path(
    objects_dir,
    "TNBC_Integrated.rds"
  )
)


############################################################
# Metadata
############################################################

# The biological group was already stored during Step 1.
# Here we create a standardized metadata field for downstream
# analysis.

srobj$Class <- srobj$group


############################################################
# Check group composition
############################################################

print(
  table(srobj$Class)
)


############################################################
# Cell-cycle scoring
############################################################

DefaultAssay(srobj) <- "RNA"

data("cc.genes.updated.2019")

s.genes <- cc.genes.updated.2019$s.genes

g2m.genes <- cc.genes.updated.2019$g2m.genes


srobj <- CellCycleScoring(
  object = srobj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)


# Cell-cycle difference score
srobj$CC.diff <-
  srobj$S.Score - srobj$G2M.Score


############################################################
# Diagnostic PCA
############################################################

srobj <- RunPCA(
  object = srobj,
  assay = "RNA",
  features = c(
    s.genes,
    g2m.genes
  ),
  verbose = FALSE
)


############################################################
# Cell-cycle diagnostic plot
############################################################

pdf(
  file.path(
    figures_main_dir,
    "CellCycle_Diagnostic_PCA.pdf"
  )
)

DimPlot(
  srobj,
  reduction = "pca",
  group.by = "Phase"
)

dev.off()


############################################################
# Main dimensionality reduction
############################################################

DefaultAssay(srobj) <- "integrated"


srobj <- ScaleData(
  object = srobj,
  vars.to.regress = "CC.diff",
  verbose = FALSE
)


srobj <- RunPCA(
  object = srobj,
  npcs = 50,
  verbose = FALSE
)


############################################################
# PCA diagnostic plot
############################################################

pdf(
  file = file.path(
    figures_main_dir,
    "PCA_ElbowPlot.pdf"
  )
)

ElbowPlot(
  srobj,
  ndims = 50
)

dev.off()


############################################################
# Clustering
############################################################

srobj <- FindNeighbors(
  object = srobj,
  dims = 1:30
)

srobj <- FindClusters(
  object = srobj,
  resolution = 0.5
)


############################################################
# UMAP
############################################################

srobj <- RunUMAP(
  object = srobj,
  dims = 1:30
)


############################################################
# UMAP visualization
############################################################

# UMAP by biological group — split
pdf(
  file = file.path(
    figures_main_dir,
    "UMAP_by_Biological_Group_Split.pdf"
  )
)

DimPlot(
  srobj,
  reduction = "umap",
  label = TRUE,
  split.by = "Class"
) +
  NoLegend()

dev.off()


# UMAP by cell-cycle phase, split by biological group
pdf(
  file = file.path(
    figures_main_dir,
    "UMAP_by_CellCycle_Phase.pdf"
  )
)

DimPlot(
  srobj,
  reduction = "umap",
  group.by = "Phase",
  split.by = "Class"
) +
  NoLegend()

dev.off()


# UMAP by biological group
pdf(
  file = file.path(
    figures_main_dir,
    "UMAP_by_Biological_Group.pdf"
  )
)

DimPlot(
  srobj,
  reduction = "umap",
  group.by = "Class"
)

dev.off()


############################################################
# Save processed Seurat object
#
# Central project object directory:
# objects/
############################################################

saveRDS(
  srobj,
  file = file.path(
    objects_dir,
    "TNBC_PCA_Clustering_UMAP.rds"
  ),
  compress = TRUE
)


############################################################
# Summary log
############################################################

log_file <- file.path(
  logs_dir,
  "Step3_PCA_Clustering_UMAP.log"
)

log_con <- file(
  log_file,
  open = "wt"
)


writeLines(
  c(
    "TNBC Single-cell RNA-seq Analysis Pipeline",
    "Step 3: PCA, clustering and UMAP",
    paste("Date:", Sys.time()),
    "",
    "Summary",
    "-------",
    paste(
      "Number of cells:",
      ncol(srobj)
    ),
    paste(
      "Number of features (integrated assay):",
      nrow(srobj[["integrated"]])
    ),
    paste(
      "Number of features (RNA assay):",
      nrow(srobj[["RNA"]])
    ),
    paste(
      "Number of clusters:",
      length(levels(Idents(srobj)))
    ),
    paste(
      "Cluster IDs:",
      paste(
        levels(Idents(srobj)),
        collapse = ", "
      )
    ),
    "",
    "Cells per cluster",
    "-----------------",
    capture.output(
      print(
        table(Idents(srobj))
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
    "Step 3 completed successfully."
  ),
  con = log_con
)

close(log_con)


############################################################
# Completion message
############################################################

message(
  "Step 3 completed successfully."
)

message(
  "Log saved to: ",
  log_file
)