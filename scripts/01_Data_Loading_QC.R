############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 1: Data Loading and Quality Control
#
# Dataset:
# Single-cell RNA-seq of TNBC tumors and BRCA1-mutant tumors
#
# This script:
# - Loads 10X Genomics matrices
# - Creates Seurat objects
# - Merges samples
# - Performs basic QC filtering
#
# NOTE:
# This script is preserved as the Stage 1 analytical source code.
# Current repository organization reflects the finalized project
# structure. No new analysis is performed during repository
# reorganization.
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

# Raw 10X input data
data_dir <- file.path(
  project_dir,
  "data",
  "raw_10X"
)

# Central object directory
objects_dir <- file.path(
  project_dir,
  "objects"
)

# Canonical Stage 1 results directory
stage_results_dir <- file.path(
  project_dir,
  "results",
  "01_Data_Loading_QC"
)

# Stage 1 subdirectories
figures_main_dir <- file.path(
  stage_results_dir,
  "figures",
  "main"
)

logs_dir <- file.path(
  stage_results_dir,
  "logs"
)

# Create canonical directories if required
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
# Stage 1 log
############################################################

log_file <- file.path(
  logs_dir,
  "Step1_QC.log"
)

sink(
  log_file,
  append = TRUE,
  split = TRUE
)

cat("\n")
cat("============================================================\n")
cat("STAGE 1 — DATA LOADING AND QUALITY CONTROL\n")
cat("============================================================\n")
cat("Project directory:\n")
cat(project_dir, "\n\n")


############################################################
# Define samples
############################################################

samples <- c(
  GSM4909281 = "Total1",
  GSM4909282 = "Total2",
  GSM4909283 = "Total3",
  GSM4909284 = "Total4",
  GSM4909285 = "BRCA1_1",
  GSM4909286 = "BRCA1_2",
  GSM4909287 = "BRCA1_3",
  GSM4909288 = "BRCA1_4"
)


############################################################
# Create Seurat objects from 10X matrices
############################################################

obj_list <- lapply(
  names(samples),
  function(sample_id) {
    
    sample_path <- file.path(
      data_dir,
      sample_id,
      "10X"
    )
    
    counts <- Read10X(sample_path)
    
    obj <- CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = 3,
      min.features = 200
    )
    
    # Store sample metadata
    obj$sample <- sample_id
    
    # Define biological groups
    obj$group <- ifelse(
      grepl("GSM490928[1-4]", sample_id),
      "TotalCell",
      "BRCA1_tumour"
    )
    
    return(obj)
  }
)

names(obj_list) <- names(samples)


############################################################
# Merge all samples
############################################################

tnbc_obj <- merge(
  obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = samples
)


############################################################
# Join RNA layers
#
# Seurat v5 compatibility
############################################################

tnbc_obj <- JoinLayers(
  tnbc_obj,
  assay = "RNA"
)


############################################################
# Quality control metrics
############################################################

tnbc_obj[["percent.mt"]] <-
  PercentageFeatureSet(
    tnbc_obj,
    pattern = "^MT-"
  )


############################################################
# QC visualization
#
# Canonical Stage 1 figure locations:
#
# results/01_Data_Loading_QC/figures/main/
#
# Current finalized outputs:
# - QC_violin_BRCA1_tumour.pdf
# - QC_violin_TotalCell.pdf
#
# The analytical content of the QC plots is unchanged.
############################################################

# BRCA1_tumour QC plot
pdf(
  file.path(
    figures_main_dir,
    "QC_violin_BRCA1_tumour.pdf"
  )
)

VlnPlot(
  subset(
    tnbc_obj,
    subset = group == "BRCA1_tumour"
  ),
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

dev.off()


# TotalCell QC plot
pdf(
  file.path(
    figures_main_dir,
    "QC_violin_TotalCell.pdf"
  )
)

VlnPlot(
  subset(
    tnbc_obj,
    subset = group == "TotalCell"
  ),
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

dev.off()


############################################################
# Filtering low-quality cells
############################################################

tnbc_obj <- subset(
  tnbc_obj,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    nCount_RNA < 40000 &
    percent.mt < 20
)


############################################################
# Save filtered Seurat object
#
# Central project object directory:
# objects/
############################################################

saveRDS(
  tnbc_obj,
  file = file.path(
    objects_dir,
    "TNBC_QC_filtered.rds"
  )
)


############################################################
# Completion log
############################################################

cat("\n")
cat("STAGE 1 COMPLETED\n")
cat("============================================================\n")
cat("Samples loaded:", length(samples), "\n")
cat("Final cells:", ncol(tnbc_obj), "\n")
cat("Final genes:", nrow(tnbc_obj), "\n")
cat("Filtered Seurat object:\n")
cat(
  file.path(
    objects_dir,
    "TNBC_QC_filtered.rds"
  ),
  "\n"
)

cat("\n")
cat("Canonical Stage 1 results:\n")
cat(stage_results_dir, "\n")
cat("============================================================\n")

sink()