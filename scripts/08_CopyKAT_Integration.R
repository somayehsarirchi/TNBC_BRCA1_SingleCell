#!/usr/bin/env Rscript

############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 8: Integrate CopyKAT predictions into Seurat
#
# This script:
# - Loads the epithelial Seurat object from Step 5
# - Imports combined CopyKAT predictions from Step 7
# - Checks prediction coverage against epithelial cells
# - Merges CopyKAT annotations into Seurat metadata
# - Labels epithelial cells as tumor-like or normal-like
# - Visualizes CNV status on UMAP
# - Summarizes CNV status by sample
# - Summarizes CNV status by epithelial subtype
# - Saves the CNV-annotated epithelial object
# - Records session information and a summary log
#
# NOTE:
# This script is retained for reproducibility and repository
# documentation. Do not rerun during the current refactoring pass.
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

project_dir <- "YOUR_PROJECT_DIRECTORY"

objects_dir <- file.path(project_dir, "objects")
results_dir <- file.path(project_dir, "results")

stage_results_dir <- file.path(
  results_dir,
  "08_CopyKAT_Integration"
)

figures_main_dir <- file.path(
  stage_results_dir,
  "figures",
  "main"
)

tables_main_dir <- file.path(
  stage_results_dir,
  "tables",
  "main"
)

tables_supp_dir <- file.path(
  stage_results_dir,
  "tables",
  "supplementary"
)

logs_dir <- file.path(
  stage_results_dir,
  "logs"
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
  tables_supp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Record starting session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "sessionInfo_start.txt")
)


############################################################
# Load epithelial object from Step 5
############################################################

epithelial_obj <- readRDS(
  file.path(
    objects_dir,
    "TNBC_Epithelial_Annotated.rds"
  )
)


############################################################
# Load combined CopyKAT predictions from Step 7
############################################################

copykat_pred_file <- file.path(
  results_dir,
  "07_CopyKAT_Batch_Inference",
  "copykat_final",
  "copykat_combined_predictions.csv"
)

copykat_pred <- read.csv(
  copykat_pred_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# Inspect prediction columns
############################################################

print(colnames(copykat_pred))


############################################################
# Validate required cell identifier column
############################################################

if (!"cell.names" %in% colnames(copykat_pred)) {
  stop(
    "Could not find the required 'cell.names' column ",
    "in the CopyKAT prediction file."
  )
}


############################################################
# Detect CopyKAT prediction column automatically
############################################################

pred_candidates <- c(
  "copykat.pred",
  "prediction",
  "pred",
  "copykat_prediction"
)

pred_col <- pred_candidates[
  pred_candidates %in% colnames(copykat_pred)
][1]

if (is.na(pred_col)) {
  stop(
    "Could not find the CopyKAT prediction column."
  )
}


############################################################
# Standardize CopyKAT prediction metadata
############################################################

copykat_meta <- copykat_pred %>%
  dplyr::select(
    cell.names,
    !!pred_col
  ) %>%
  dplyr::rename(
    CopyKAT_Prediction = !!pred_col
  )


############################################################
# Check prediction identifiers
############################################################

if (anyDuplicated(copykat_meta$cell.names) > 0) {
  stop(
    "Duplicate cell identifiers were detected in the ",
    "CopyKAT prediction table."
  )
}


############################################################
# Compare CopyKAT predictions with epithelial cells
############################################################

epithelial_cells <- colnames(epithelial_obj)
copykat_cells <- copykat_meta$cell.names

missing_predictions <- setdiff(
  epithelial_cells,
  copykat_cells
)

extra_predictions <- setdiff(
  copykat_cells,
  epithelial_cells
)


############################################################
# Save missing and extra prediction records
############################################################

write.csv(
  data.frame(
    cell.names = missing_predictions,
    stringsAsFactors = FALSE
  ),
  file.path(
    tables_supp_dir,
    "copykat_missing_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    cell.names = extra_predictions,
    stringsAsFactors = FALSE
  ),
  file.path(
    tables_supp_dir,
    "copykat_extra_predictions.csv"
  ),
  row.names = FALSE
)


############################################################
# Merge CopyKAT metadata with Seurat object
############################################################

copykat_meta <- copykat_meta[
  match(
    epithelial_cells,
    copykat_meta$cell.names
  ),
]


############################################################
# Verify cell-order matching
############################################################

if (!identical(
  epithelial_cells,
  copykat_meta$cell.names
)) {
  stop(
    "CopyKAT cell identifiers could not be aligned ",
    "with the epithelial Seurat object."
  )
}


############################################################
# Add CopyKAT prediction to Seurat metadata
############################################################

epithelial_obj$CopyKAT_Prediction <-
  copykat_meta$CopyKAT_Prediction


############################################################
# Create simplified CNV status
############################################################

epithelial_obj$CNV_Status <- dplyr::case_when(
  
  epithelial_obj$CopyKAT_Prediction %in%
    c("aneuploid", "Aneuploid") ~ "Tumor_like",
  
  epithelial_obj$CopyKAT_Prediction %in%
    c("diploid", "Diploid") ~ "Normal_like",
  
  TRUE ~ "Unknown"
)


############################################################
# Define CNV status factor levels
############################################################

epithelial_obj$CNV_Status <- factor(
  epithelial_obj$CNV_Status,
  levels = c(
    "Tumor_like",
    "Normal_like",
    "Unknown"
  )
)


############################################################
# Check annotation summary
############################################################

print(
  table(
    epithelial_obj$CNV_Status,
    useNA = "ifany"
  )
)

print(
  table(
    epithelial_obj$Class,
    epithelial_obj$CNV_Status,
    useNA = "ifany"
  )
)

print(
  prop.table(
    table(
      epithelial_obj$Class,
      epithelial_obj$CNV_Status,
      useNA = "ifany"
    ),
    margin = 1
  )
)


############################################################
# UMAP: CNV status
############################################################

pdf(
  file.path(
    figures_main_dir,
    "CopyKAT_CNV_Status_UMAP.pdf"
  ),
  width = 12,
  height = 8
)

print(
  DimPlot(
    epithelial_obj,
    reduction = "umap",
    group.by = "CNV_Status",
    label = TRUE,
    repel = TRUE,
    split.by = "Class"
  ) +
    NoLegend()
)

dev.off()


############################################################
# UMAP: epithelial cells by CNV status
############################################################

pdf(
  file.path(
    figures_main_dir,
    "UMAP_Epithelial_CNV_Status.pdf"
  ),
  width = 12,
  height = 8
)

print(
  DimPlot(
    epithelial_obj,
    reduction = "umap",
    group.by = "CNV_Status",
    label = TRUE,
    repel = TRUE
  ) +
    NoLegend()
)

dev.off()


############################################################
# Cell counts by sample and CNV status
############################################################

cnv_table <- table(
  epithelial_obj$orig.ident,
  epithelial_obj$CNV_Status
)

write.csv(
  as.data.frame(cnv_table),
  file.path(
    tables_main_dir,
    "CopyKAT_CNV_by_Sample.csv"
  ),
  row.names = FALSE
)


############################################################
# Proportion of CNV status by sample
############################################################

cnv_prop <- prop.table(
  cnv_table,
  margin = 1
)

write.csv(
  as.data.frame(cnv_prop),
  file.path(
    tables_main_dir,
    "CopyKAT_CNV_Proportions_by_Sample.csv"
  ),
  row.names = FALSE
)


############################################################
# CNV status by epithelial subtype
############################################################

cnv_subtype_table <- table(
  epithelial_obj$Epithelial_Subtype,
  epithelial_obj$CNV_Status
)

write.csv(
  as.data.frame(cnv_subtype_table),
  file.path(
    tables_supp_dir,
    "CNV_by_EpithelialSubtype.csv"
  ),
  row.names = FALSE
)


############################################################
# CNV status proportions by epithelial subtype
############################################################

cnv_subtype_prop <- prop.table(
  cnv_subtype_table,
  margin = 1
)

write.csv(
  as.data.frame(cnv_subtype_prop),
  file.path(
    tables_supp_dir,
    "CNV_by_EpithelialSubtype_Percent.csv"
  ),
  row.names = FALSE
)


############################################################
# Save CNV-annotated epithelial object
############################################################

saveRDS(
  epithelial_obj,
  file = file.path(
    objects_dir,
    "TNBC_Epithelial_CopyKAT.rds"
  ),
  compress = TRUE
)


############################################################
# Record ending session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "sessionInfo_end.txt")
)


############################################################
# Summary log
############################################################

log_file <- file.path(
  logs_dir,
  "Step8_CopyKAT_Integration.log"
)

log_con <- file(
  log_file,
  open = "wt"
)

writeLines(
  c(
    "TNBC Single-cell RNA-seq Analysis Pipeline",
    "Step 8: Integrate CopyKAT predictions into Seurat",
    paste("Date:", Sys.time()),
    "",
    "Input epithelial object",
    "-----------------------",
    file.path(
      objects_dir,
      "TNBC_Epithelial_Annotated.rds"
    ),
    "",
    "Input CopyKAT predictions",
    "-------------------------",
    copykat_pred_file,
    "",
    "Summary",
    "-------",
    paste(
      "Total epithelial cells:",
      ncol(epithelial_obj)
    ),
    paste(
      "CopyKAT prediction records:",
      nrow(copykat_pred)
    ),
    paste(
      "Epithelial cells without CopyKAT predictions:",
      length(missing_predictions)
    ),
    paste(
      "CopyKAT predictions not present in epithelial object:",
      length(extra_predictions)
    ),
    "",
    "CNV status counts",
    "-----------------",
    capture.output(
      print(
        table(
          epithelial_obj$CNV_Status,
          useNA = "ifany"
        )
      )
    ),
    "",
    "Outputs generated",
    "-----------------",
    file.path(
      figures_main_dir,
      "CopyKAT_CNV_Status_UMAP.pdf"
    ),
    file.path(
      figures_main_dir,
      "UMAP_Epithelial_CNV_Status.pdf"
    ),
    file.path(
      tables_main_dir,
      "CopyKAT_CNV_by_Sample.csv"
    ),
    file.path(
      tables_main_dir,
      "CopyKAT_CNV_Proportions_by_Sample.csv"
    ),
    file.path(
      tables_supp_dir,
      "CNV_by_EpithelialSubtype.csv"
    ),
    file.path(
      tables_supp_dir,
      "CNV_by_EpithelialSubtype_Percent.csv"
    ),
    file.path(
      tables_supp_dir,
      "copykat_missing_predictions.csv"
    ),
    file.path(
      tables_supp_dir,
      "copykat_extra_predictions.csv"
    ),
    file.path(
      objects_dir,
      "TNBC_Epithelial_CopyKAT.rds"
    ),
    "",
    "Step 8 completed successfully."
  ),
  con = log_con
)

close(log_con)


############################################################
# Completion message
############################################################

message(
  "Step 8 completed successfully."
)

message(
  "Log saved to: ",
  log_file
)