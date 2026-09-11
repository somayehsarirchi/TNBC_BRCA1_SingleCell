############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 6: CopyKAT preparation
#
# This script:
# - Loads the annotated whole-dataset Seurat object
# - Loads the purified epithelial Seurat object
# - Selects immune and stromal cells as known normal reference cells
# - Combines purified epithelial and normal cells
# - Extracts raw RNA counts
# - Saves CopyKAT input files for Linux execution
############################################################


############################################################
# Load libraries
############################################################

library(Seurat)
library(dplyr)


############################################################
# Project directories
############################################################

project_dir <- "YOUR_PROJECT_DIRECTORY"

objects_dir <- file.path(
  project_dir,
  "objects"
)

stage_results_dir <- file.path(
  project_dir,
  "results",
  "06_CopyKAT_Preparation"
)

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

logs_dir <- file.path(
  stage_results_dir,
  "logs"
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
# Load Seurat objects
############################################################

srobj <- readRDS(
  file.path(
    objects_dir,
    "TNBC_Annotated.rds"
  )
)

# Use the purified epithelial object generated in Step 5
epithelial_obj <- readRDS(
  file.path(
    objects_dir,
    "TNBC_Epithelial_Pure.rds"
  )
)


############################################################
# Join RNA layers (Seurat v5 compatibility)
############################################################

srobj <- JoinLayers(
  srobj,
  assay = "RNA"
)


############################################################
# Define known normal reference populations
############################################################

normal_celltypes <- c(
  "NK_cell-T_cells",
  "B_cell",
  "Macrophage",
  "DC-Monocyte",
  "Fibroblasts-Tissue_stem_cells",
  "Endothelial_cells"
)

normal_cells <- subset(
  srobj,
  subset = Manual_CellType %in% normal_celltypes
)

normal_barcodes <- colnames(normal_cells)


############################################################
# Combine purified epithelial cells and known normal cells
############################################################

cells_for_copykat <- unique(
  c(
    colnames(epithelial_obj),
    normal_barcodes
  )
)

copykat_obj <- subset(
  srobj,
  cells = cells_for_copykat
)


############################################################
# Verify that both experimental groups are present
############################################################

group_table <- table(
  copykat_obj$Class
)

print(group_table)

if (!all(
  c(
    "TotalCell",
    "BRCA1_tumour"
  ) %in% names(group_table)
)) {
  stop(
    "One or both experimental groups are missing from the CopyKAT input object."
  )
}


############################################################
# Extract raw RNA counts
############################################################

raw_matrix <- LayerData(
  copykat_obj,
  assay = "RNA",
  layer = "counts"
)


############################################################
# Save CopyKAT input object
############################################################

copykat_inputs <- list(
  matrix = raw_matrix,
  normal_cells = intersect(
    normal_barcodes,
    colnames(raw_matrix)
  )
)

saveRDS(
  copykat_inputs,
  file = file.path(
    objects_dir,
    "copykat_inputs.rds"
  ),
  compress = TRUE
)


############################################################
# Save normal cell barcodes
############################################################

write.csv(
  data.frame(
    barcode = copykat_inputs$normal_cells
  ),
  file = file.path(
    tables_supplementary_dir,
    "CopyKAT_Normal_Cell_Barcodes.csv"
  ),
  row.names = FALSE
)


############################################################
# Save matrix summary
############################################################

write.csv(
  data.frame(
    genes = nrow(raw_matrix),
    cells = ncol(raw_matrix),
    epithelial_cells = ncol(epithelial_obj),
    normal_cells = length(
      copykat_inputs$normal_cells
    ),
    totalcell_cells = sum(
      copykat_obj$Class == "TotalCell"
    ),
    brca1_cells = sum(
      copykat_obj$Class == "BRCA1_tumour"
    )
  ),
  file = file.path(
    tables_main_dir,
    "CopyKAT_Input_Summary.csv"
  ),
  row.names = FALSE
)


############################################################
# Summary log
############################################################

log_file <- file.path(
  logs_dir,
  "Step6_CopyKAT_Preparation.log"
)

log_con <- file(
  log_file,
  open = "wt"
)

writeLines(
  c(
    "TNBC Single-cell RNA-seq Analysis Pipeline",
    "Step 6: CopyKAT preparation",
    paste("Date:", Sys.time()),
    "",
    "Summary",
    "-------",
    paste(
      "Purified epithelial cells:",
      ncol(epithelial_obj)
    ),
    paste(
      "Known normal cells:",
      length(normal_barcodes)
    ),
    paste(
      "Total cells for CopyKAT:",
      ncol(raw_matrix)
    ),
    paste(
      "Genes:",
      nrow(raw_matrix)
    ),
    paste(
      "TotalCell cells:",
      sum(
        copykat_obj$Class == "TotalCell"
      )
    ),
    paste(
      "BRCA1_tumour cells:",
      sum(
        copykat_obj$Class == "BRCA1_tumour"
      )
    ),
    paste(
      "Normal cells included:",
      length(
        copykat_inputs$normal_cells
      )
    ),
    "",
    "Normal cell types used",
    "----------------------",
    paste(
      normal_celltypes,
      collapse = ", "
    ),
    "",
    "Output files",
    "------------",
    file.path(
      objects_dir,
      "copykat_inputs.rds"
    ),
    file.path(
      tables_supplementary_dir,
      "CopyKAT_Normal_Cell_Barcodes.csv"
    ),
    file.path(
      tables_main_dir,
      "CopyKAT_Input_Summary.csv"
    ),
    "",
    "Step 6 completed successfully."
  ),
  con = log_con
)

close(log_con)


############################################################
# Completion messages
############################################################

message(
  "Step 6 completed successfully."
)

message(
  "CopyKAT input saved to: ",
  file.path(
    objects_dir,
    "copykat_inputs.rds"
  )
)

message(
  "Log saved to: ",
  log_file
)