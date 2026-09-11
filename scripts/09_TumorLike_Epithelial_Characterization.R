#!/usr/bin/env Rscript

############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 9: Tumor-like epithelial compartment characterization
#
# This script:
# - Loads the CopyKAT-annotated epithelial object
# - Extracts Tumor-like epithelial cells
# - Summarizes CNV status across experimental groups
# - Summarizes tumor-like epithelial subtypes by group
# - Generates tumor-like epithelial UMAP and barplot figures
# - Saves the tumor-like epithelial object for downstream analyses
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

objects_dir <- file.path(
  project_dir,
  "objects"
)

results_dir <- file.path(
  project_dir,
  "results"
)

stage_results_dir <- file.path(
  results_dir,
  "09_TumorLike_Epithelial_Characterization"
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
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Load CopyKAT-annotated epithelial object
############################################################

epithelial_obj <- readRDS(
  file.path(
    objects_dir,
    "TNBC_Epithelial_CopyKAT.rds"
  )
)


############################################################
# Verify required metadata
############################################################

required_cols <- c(
  "Class",
  "Epithelial_Subtype",
  "CNV_Status"
)

missing_cols <- setdiff(
  required_cols,
  colnames(epithelial_obj@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}


############################################################
# Extract Tumor-like epithelial cells
############################################################

epithelial_obj_tumor <- subset(
  epithelial_obj,
  subset = CNV_Status == "Tumor_like"
)


############################################################
# Drop unused epithelial subtype levels
############################################################

epithelial_obj_tumor$Epithelial_Subtype <- droplevels(
  epithelial_obj_tumor$Epithelial_Subtype
)


############################################################
# Tumor-like epithelial subtype percentage by group
############################################################

tumor_subtype_percent <- round(
  prop.table(
    table(
      epithelial_obj_tumor$Epithelial_Subtype,
      epithelial_obj_tumor$Class
    ),
    margin = 2
  ) * 100,
  1
)

print(tumor_subtype_percent)

write.csv(
  as.data.frame.matrix(tumor_subtype_percent),
  file.path(
    tables_main_dir,
    "TumorLike_EpithelialSubtype_Percent_by_Group.csv"
  ),
  row.names = TRUE
)


############################################################
# Save tumor-like epithelial object
############################################################

saveRDS(
  epithelial_obj_tumor,
  file = file.path(
    objects_dir,
    "TNBC_Epithelial_TumorLike.rds"
  ),
  compress = TRUE
)


############################################################
# CNV status by experimental group
############################################################

cnv_by_group <- as.data.frame.matrix(
  table(
    epithelial_obj$Class,
    epithelial_obj$CNV_Status
  )
)

write.csv(
  cnv_by_group,
  file.path(
    tables_main_dir,
    "CNV_by_Group.csv"
  ),
  row.names = TRUE
)


############################################################
# CNV status percentage by experimental group
############################################################

cnv_by_group_percent <- round(
  prop.table(
    table(
      epithelial_obj$Class,
      epithelial_obj$CNV_Status
    ),
    margin = 1
  ) * 100,
  1
)

write.csv(
  as.data.frame.matrix(cnv_by_group_percent),
  file.path(
    tables_main_dir,
    "CNV_by_Group_Percent.csv"
  ),
  row.names = TRUE
)


############################################################
# Tumor-like epithelial subtype counts by group
############################################################

tumor_subtype_by_group <- as.data.frame.matrix(
  table(
    epithelial_obj_tumor$Epithelial_Subtype,
    epithelial_obj_tumor$Class
  )
)

write.csv(
  tumor_subtype_by_group,
  file.path(
    tables_main_dir,
    "TumorLike_EpithelialSubtype_by_Group.csv"
  ),
  row.names = TRUE
)


############################################################
# UMAP: tumor-like epithelial subtypes
############################################################

p1 <- DimPlot(
  epithelial_obj_tumor,
  reduction = "umap",
  group.by = "Epithelial_Subtype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  NoLegend() +
  ggtitle(
    "Tumor-like epithelial subtypes"
  )

ggsave(
  file.path(
    figures_main_dir,
    "UMAP_TumorLike_Epithelial_Subtypes.pdf"
  ),
  p1,
  width = 10,
  height = 7
)


############################################################
# Barplot: tumor-like epithelial subtypes by group
############################################################

subtype_order <- epithelial_obj_tumor@meta.data %>%
  count(
    Epithelial_Subtype,
    sort = TRUE
  ) %>%
  pull(Epithelial_Subtype)


df_subtype <- epithelial_obj_tumor@meta.data %>%
  count(
    Epithelial_Subtype,
    Class
  ) %>%
  group_by(
    Epithelial_Subtype
  ) %>%
  mutate(
    percent = 100 * n / sum(n)
  )


df_subtype$Epithelial_Subtype <- factor(
  df_subtype$Epithelial_Subtype,
  levels = subtype_order
)


p2 <- ggplot(
  df_subtype,
  aes(
    x = Epithelial_Subtype,
    y = percent,
    fill = Class
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.2
  ) +
  theme_classic() +
  labs(
    title = "Tumor-like epithelial subtype composition",
    x = "Epithelial subtype",
    y = "Percentage",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    )
  )

ggsave(
  file.path(
    figures_main_dir,
    "Barplot_TumorLike_Epithelial_Subtypes_by_Group.pdf"
  ),
  p2,
  width = 12,
  height = 6
)


############################################################
# Summary log
############################################################

log_file <- file.path(
  logs_dir,
  "Step9_TumorLike_Epithelial.log"
)

log_con <- file(
  log_file,
  open = "wt"
)

writeLines(
  c(
    "TNBC Single-cell RNA-seq Analysis Pipeline",
    "Step 9: Tumor-like epithelial compartment characterization",
    paste("Date:", Sys.time()),
    "",
    "Input object",
    "------------",
    file.path(
      objects_dir,
      "TNBC_Epithelial_CopyKAT.rds"
    ),
    "",
    "Summary",
    "-------",
    paste(
      "Total epithelial cells:",
      ncol(epithelial_obj)
    ),
    paste(
      "Tumor-like epithelial cells:",
      ncol(epithelial_obj_tumor)
    ),
    paste(
      "Number of tumor-like epithelial subtypes:",
      length(
        unique(
          epithelial_obj_tumor$Epithelial_Subtype
        )
      )
    ),
    "",
    "Tumor-like cells per group",
    "-------------------------",
    capture.output(
      print(
        table(
          epithelial_obj_tumor$Class
        )
      )
    ),
    "",
    "Output object",
    "-------------",
    file.path(
      objects_dir,
      "TNBC_Epithelial_TumorLike.rds"
    ),
    "",
    "Key result tables",
    "-----------------",
    file.path(
      tables_main_dir,
      "CNV_by_Group.csv"
    ),
    file.path(
      tables_main_dir,
      "CNV_by_Group_Percent.csv"
    ),
    file.path(
      tables_main_dir,
      "TumorLike_EpithelialSubtype_by_Group.csv"
    ),
    file.path(
      tables_main_dir,
      "TumorLike_EpithelialSubtype_Percent_by_Group.csv"
    ),
    "",
    "Key figures",
    "-----------",
    file.path(
      figures_main_dir,
      "UMAP_TumorLike_Epithelial_Subtypes.pdf"
    ),
    file.path(
      figures_main_dir,
      "Barplot_TumorLike_Epithelial_Subtypes_by_Group.pdf"
    ),
    "",
    "Step 9 completed successfully."
  ),
  con = log_con
)

close(log_con)


############################################################
# Completion message
############################################################

message(
  "Step 9 completed successfully."
)

message(
  "Tumor-like epithelial object saved to: ",
  file.path(
    objects_dir,
    "TNBC_Epithelial_TumorLike.rds"
  )
)

message(
  "Log saved to: ",
  log_file
)