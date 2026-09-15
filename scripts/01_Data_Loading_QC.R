# ============================================================
# Stage 01 — Data Loading and Quality Control
# TNBC / BRCA1 Single-Cell RNA-seq Analysis
# ============================================================

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------

project_dir <- "C:/Users/asus/Desktop/GSE161529"

data_dir <- file.path(
  project_dir,
  "data")

results_dir <- file.path(
  project_dir,
  "results",
  "01_Data_Loading_QC"
)

figures_dir <- file.path(
  results_dir,
  "figures",
  "main"
)

logs_dir <- file.path(
  results_dir,
  "logs"
)

objects_dir <- file.path(
  project_dir,
  "objects"
)

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Logging
# ------------------------------------------------------------

log_file <- file.path(
  logs_dir,
  "Step1_QC.log"
)

sink(log_file, split = TRUE)

cat("============================================================\n")
cat("Stage 01 — Data Loading and Quality Control\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Project directory:", project_dir, "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 4. Sample information
# ------------------------------------------------------------

sample_ids <- c(
  "GSM4909281",
  "GSM4909282",
  "GSM4909283",
  "GSM4909284",
  "GSM4909285",
  "GSM4909286",
  "GSM4909287",
  "GSM4909288"
)

group_labels <- c(
  rep("TotalCell", 4),
  rep("BRCA1_tumour", 4)
)

names(group_labels) <- sample_ids

cat("Sample groups:\n")

print(
  data.frame(
    sample = sample_ids,
    group = group_labels
  )
)

cat("\n")

# ------------------------------------------------------------
# 5. Read each sample independently
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Reading individual samples\n")
cat("------------------------------------------------------------\n")

seurat_objects <- list()

for (sample_id in sample_ids) {
  
  sample_dir <- file.path(
    data_dir,
    sample_id,
    "10X"
  )
  
  cat("\nSample:", sample_id, "\n")
  cat("Directory:", sample_dir, "\n")
  
  if (!dir.exists(sample_dir)) {
    stop(
      "10X directory not found for ",
      sample_id,
      ": ",
      sample_dir
    )
  }
  
  counts <- Read10X(
    data.dir = sample_dir
  )
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample <- sample_id
  obj$group <- group_labels[[sample_id]]
  
  seurat_objects[[sample_id]] <- obj
  
  cat(
    "Cells:", ncol(obj),
    "| Genes:", nrow(obj),
    "| Group:", group_labels[[sample_id]],
    "\n"
  )
}

cat("\nAll samples loaded successfully.\n\n")

# ------------------------------------------------------------
# 6. Merge samples within each biological group
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Merging samples within biological groups\n")
cat("------------------------------------------------------------\n")

# TotalCell

TotalCell <- merge(
  seurat_objects[["GSM4909281"]],
  y = list(
    seurat_objects[["GSM4909282"]],
    seurat_objects[["GSM4909283"]],
    seurat_objects[["GSM4909284"]]
  ),
  add.cell.ids = c(
    "GSM4909281",
    "GSM4909282",
    "GSM4909283",
    "GSM4909284"
  ),
  project = "TotalCell"
)

# BRCA1 tumour

BRCA1tumour <- merge(
  seurat_objects[["GSM4909285"]],
  y = list(
    seurat_objects[["GSM4909286"]],
    seurat_objects[["GSM4909287"]],
    seurat_objects[["GSM4909288"]]
  ),
  add.cell.ids = c(
    "GSM4909285",
    "GSM4909286",
    "GSM4909287",
    "GSM4909288"
  ),
  project = "BRCA1tumour"
)

cat("TotalCell cells:", ncol(TotalCell), "\n")
cat("BRCA1tumour cells:", ncol(BRCA1tumour), "\n\n")

# ------------------------------------------------------------
# 7. Calculate mitochondrial percentage separately
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Calculating mitochondrial percentages\n")
cat("------------------------------------------------------------\n")

TotalCell[["percent.mt"]] <- PercentageFeatureSet(
  TotalCell,
  pattern = "^MT-"
)

BRCA1tumour[["percent.mt"]] <- PercentageFeatureSet(
  BRCA1tumour,
  pattern = "^MT-"
)

cat("Mitochondrial QC metrics calculated.\n\n")

# ------------------------------------------------------------
# 8. Pre-QC violin plots
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Generating pre-QC violin plots\n")
cat("------------------------------------------------------------\n")

p_TotalCell <- VlnPlot(
  TotalCell,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  ncol = 3,
  pt.size = 0
)

p_TotalCell <- p_TotalCell +
  plot_annotation(
    title = "TotalCell — Pre-QC distributions"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "QC_violin_TotalCell_preQC.pdf"
  ),
  plot = p_TotalCell,
  width = 12,
  height = 5
)

p_BRCA1 <- VlnPlot(
  BRCA1tumour,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  ncol = 3,
  pt.size = 0
)

p_BRCA1 <- p_BRCA1 +
  plot_annotation(
    title = "BRCA1 tumour — Pre-QC distributions"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "QC_violin_BRCA1_tumour_preQC.pdf"
  ),
  plot = p_BRCA1,
  width = 12,
  height = 5
)

cat("Pre-QC plots saved.\n\n")

# ------------------------------------------------------------
# 9. QC thresholds
#
# These values are provisional and should be reviewed
# after inspecting the pre-QC distributions.
# ------------------------------------------------------------

# TotalCell

TotalCell_min_features <- 200
TotalCell_max_features <- 5000
TotalCell_max_counts <- 40000
TotalCell_max_mt <- 20

# BRCA1 tumour

BRCA1_min_features <- 200
BRCA1_max_features <- 5000
BRCA1_max_counts <- 40000
BRCA1_max_mt <- 20

# ------------------------------------------------------------
# 10. Cell counts before QC
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Cell counts before QC\n")
cat("------------------------------------------------------------\n")

cat("\nTotalCell:\n")
print(table(TotalCell$sample))

cat("\nBRCA1tumour:\n")
print(table(BRCA1tumour$sample))

# ------------------------------------------------------------
# 11. Apply QC filters
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Applying QC filters\n")
cat("------------------------------------------------------------\n")

TotalCell_QC_filtered <- subset(
  TotalCell,
  subset =
    nFeature_RNA > TotalCell_min_features &
    nFeature_RNA < TotalCell_max_features &
    nCount_RNA < TotalCell_max_counts &
    percent.mt < TotalCell_max_mt
)

BRCA1tumour_QC_filtered <- subset(
  BRCA1tumour,
  subset =
    nFeature_RNA > BRCA1_min_features &
    nFeature_RNA < BRCA1_max_features &
    nCount_RNA < BRCA1_max_counts &
    percent.mt < BRCA1_max_mt
)

# ------------------------------------------------------------
# 12. Join layers after filtering
# ------------------------------------------------------------

TotalCell_QC_filtered <- JoinLayers(
  TotalCell_QC_filtered,
  assay = "RNA"
)

BRCA1tumour_QC_filtered <- JoinLayers(
  BRCA1tumour_QC_filtered,
  assay = "RNA"
)

# ------------------------------------------------------------
# 13. Cell counts after QC
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Cell counts after QC\n")
cat("------------------------------------------------------------\n")

cat("\nTotalCell:\n")
print(table(TotalCell_QC_filtered$sample))

cat("\nBRCA1tumour:\n")
print(table(BRCA1tumour_QC_filtered$sample))

cat("\nTotalCell cells before QC:",
    ncol(TotalCell), "\n")

cat("TotalCell cells after QC:",
    ncol(TotalCell_QC_filtered), "\n")

cat("TotalCell retained:",
    round(
      100 * ncol(TotalCell_QC_filtered) / ncol(TotalCell),
      2
    ),
    "%\n")

cat("\nBRCA1tumour cells before QC:",
    ncol(BRCA1tumour), "\n")

cat("BRCA1tumour cells after QC:",
    ncol(BRCA1tumour_QC_filtered), "\n")

cat("BRCA1tumour retained:",
    round(
      100 * ncol(BRCA1tumour_QC_filtered) /
        ncol(BRCA1tumour),
      2
    ),
    "%\n")

# ------------------------------------------------------------
# 14. Save filtered Seurat objects
# ------------------------------------------------------------
stopifnot(
  all(TotalCell_QC_filtered$group == "TotalCell"),
  all(BRCA1tumour_QC_filtered$group == "BRCA1_tumour")
)

cat("\nCells per sample after QC:\n")

print(table(
  TotalCell_QC_filtered$sample,
  TotalCell_QC_filtered$group
))

print(table(
  BRCA1tumour_QC_filtered$sample,
  BRCA1tumour_QC_filtered$group
))

cat("\n------------------------------------------------------------\n")
cat("Saving filtered Seurat objects\n")
cat("------------------------------------------------------------\n")

saveRDS(
  TotalCell_QC_filtered,
  file = file.path(
    objects_dir,
    "TotalCell_QC_filtered.rds"
  )
)

saveRDS(
  BRCA1tumour_QC_filtered,
  file = file.path(
    objects_dir,
    "BRCA1tumour_QC_filtered.rds"
  )
)

cat("Saved:\n")
cat(
  file.path(
    objects_dir,
    "TotalCell_QC_filtered.rds"
  ),
  "\n"
)

cat(
  file.path(
    objects_dir,
    "BRCA1tumour_QC_filtered.rds"
  ),
  "\n"
)

# ------------------------------------------------------------
# 15. Record final QC thresholds
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("QC thresholds used in this run\n")
cat("------------------------------------------------------------\n")

cat("\nTotalCell:\n")
cat("nFeature_RNA >", TotalCell_min_features, "\n")
cat("nFeature_RNA <", TotalCell_max_features, "\n")
cat("nCount_RNA   <", TotalCell_max_counts, "\n")
cat("percent.mt   <", TotalCell_max_mt, "\n")

cat("\nBRCA1 tumour:\n")
cat("nFeature_RNA >", BRCA1_min_features, "\n")
cat("nFeature_RNA <", BRCA1_max_features, "\n")
cat("nCount_RNA   <", BRCA1_max_counts, "\n")
cat("percent.mt   <", BRCA1_max_mt, "\n")

# ------------------------------------------------------------
# 16. Session information
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("Session information\n")
cat("------------------------------------------------------------\n")

print(sessionInfo())

cat("\n============================================================\n")
cat("Stage 01 completed:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

