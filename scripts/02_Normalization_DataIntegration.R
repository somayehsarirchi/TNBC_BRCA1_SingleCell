############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 2: Normalization and Data Integration
#
# This script:
# - Loads QC-filtered group-level Seurat objects
# - Normalizes RNA expression data
# - Identifies highly variable features
# - Selects integration features
# - Finds integration anchors
# - Integrates TotalCell and BRCA1 tumour datasets
# - Saves integration objects
#
# NOTE:
# This script represents the original Stage 2 analytical workflow.
# The current version only updates project paths and output
# organization to match the canonical repository structure.
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

# Canonical Stage 2 results directory
stage_results_dir <- file.path(
  project_dir,
  "results",
  "02_Normalization_Integration"
)

# Stage 2 log directory
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
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Load QC-filtered Seurat objects
############################################################

TotalCell <- readRDS(
  file.path(
    objects_dir,
    "TotalCell_QC_filtered.rds"
  )
)

BRCA1tumour <- readRDS(
  file.path(
    objects_dir,
    "BRCA1tumour_QC_filtered.rds"
  )
)


############################################################
# Check loaded objects
############################################################

message(
  "TotalCell dimensions: ",
  paste(dim(TotalCell), collapse = " x ")
)

message(
  "BRCA1 tumour dimensions: ",
  paste(dim(BRCA1tumour), collapse = " x ")
)


############################################################
# Normalize RNA expression data
############################################################

TotalCell <- NormalizeData(
  TotalCell,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

BRCA1tumour <- NormalizeData(
  BRCA1tumour,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


############################################################
# Identify highly variable features
############################################################

TotalCell <- FindVariableFeatures(
  TotalCell,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000
)

BRCA1tumour <- FindVariableFeatures(
  BRCA1tumour,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000
)


############################################################
# Select features for integration
############################################################

AnchorFeatures <- SelectIntegrationFeatures(
  object.list = list(
    TotalCell,
    BRCA1tumour
  ),
  nfeatures = 2000
)


############################################################
# Find integration anchors
############################################################

anchors <- FindIntegrationAnchors(
  object.list = list(
    TotalCell,
    BRCA1tumour
  ),
  anchor.features = AnchorFeatures
)


############################################################
# Save integration anchors
#
# Central object directory:
# objects/
############################################################

saveRDS(
  anchors,
  file = file.path(
    objects_dir,
    "TNBC_Integration_Anchors.rds"
  ),
  compress = TRUE
)


############################################################
# Integrate datasets
############################################################

srobj <- IntegrateData(
  anchorSet = anchors
)


############################################################
# Save integrated Seurat object
#
# Central object directory:
# objects/
############################################################

saveRDS(
  srobj,
  file = file.path(
    objects_dir,
    "TNBC_Integrated.rds"
  ),
  compress = TRUE
)


############################################################
# Summary log
############################################################

log_file <- file.path(
  logs_dir,
  "Step2_Integration.log"
)

log_con <- file(
  log_file,
  open = "wt"
)

writeLines(
  c(
    "TNBC Single-cell RNA-seq Analysis Pipeline",
    "Step 2: Normalization and data integration",
    paste("Date:", Sys.time()),
    "",
    "Summary",
    "-------",
    paste(
      "TotalCell cells:",
      ncol(TotalCell)
    ),
    paste(
      "BRCA1 tumour cells:",
      ncol(BRCA1tumour)
    ),
    paste(
      "Integrated cells:",
      ncol(srobj)
    ),
    paste(
      "Integrated features (integrated assay):",
      nrow(srobj[["integrated"]])
    ),
    paste(
      "RNA features:",
      nrow(srobj[["RNA"]])
    ),
    paste(
      "Number of integration features:",
      length(AnchorFeatures)
    ),
    "",
    "Objects saved",
    "-------------",
    file.path(
      objects_dir,
      "TNBC_Integration_Anchors.rds"
    ),
    file.path(
      objects_dir,
      "TNBC_Integrated.rds"
    ),
    "",
    "Step 2 completed successfully."
  ),
  con = log_con
)

close(log_con)

message(
  "Step 2 completed successfully."
)

message(
  "Log saved to: ",
  log_file
)