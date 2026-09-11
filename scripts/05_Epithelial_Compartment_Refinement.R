############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
# Step 5: Epithelial compartment re-analysis and refinement
#
# This script:
# - Loads the annotated Seurat object
# - Extracts epithelial cells
# - Re-analyzes the epithelial compartment
# - Identifies epithelial cluster marker genes
# - Creates a balanced epithelial reference
# - Generates lineage verification plots
# - Applies refined epithelial subtype annotations
# - Transfers annotations to all epithelial cells
# - Removes the predefined immune-contaminated subtype
# - Saves epithelial objects, figures, tables, and log
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

metadata_dir <- file.path(
  project_dir,
  "metadata"
)

stage_results_dir <- file.path(
  project_dir,
  "results",
  "05_Epithelial_Compartment"
)

figures_main_dir <- file.path(
  stage_results_dir,
  "figures",
  "main"
)

figures_supplementary_dir <- file.path(
  stage_results_dir,
  "figures",
  "supplementary"
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
  figures_main_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_supplementary_dir,
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
# Load annotated Seurat object
############################################################

srobj <- readRDS(
  file.path(
    objects_dir,
    "TNBC_Annotated.rds"
  )
)


############################################################
# Extract epithelial cells
############################################################

epithelial_obj <- subset(
  srobj,
  subset = Manual_CellType == "Epithelial_cells"
)

epithelial_obj <- JoinLayers(
  epithelial_obj,
  assay = "RNA"
)


############################################################
# Re-analysis of epithelial compartment
############################################################

DefaultAssay(epithelial_obj) <- "RNA"

epithelial_obj <- NormalizeData(
  epithelial_obj
)

epithelial_obj <- FindVariableFeatures(
  epithelial_obj,
  selection.method = "vst",
  nfeatures = 3000
)

epithelial_obj <- ScaleData(
  epithelial_obj
)

epithelial_obj <- RunPCA(
  epithelial_obj,
  npcs = 50,
  verbose = FALSE
)


############################################################
# Epithelial PCA elbow plot
############################################################

pdf(
  file.path(
    figures_supplementary_dir,
    "Epithelial_ElbowPlot.pdf"
  ),
  width = 7,
  height = 5
)

ElbowPlot(
  epithelial_obj,
  ndims = 50
)

dev.off()


############################################################
# Epithelial clustering and UMAP
############################################################

epithelial_obj <- FindNeighbors(
  epithelial_obj,
  dims = 1:30,
  verbose = FALSE
)

epithelial_obj <- FindClusters(
  epithelial_obj,
  resolution = 0.5,
  verbose = FALSE
)

epithelial_obj <- RunUMAP(
  epithelial_obj,
  dims = 1:30,
  verbose = FALSE
)


############################################################
# Save epithelial subset object
############################################################

saveRDS(
  epithelial_obj,
  file = file.path(
    objects_dir,
    "TNBC_Epithelial_Subset.rds"
  ),
  compress = TRUE
)


############################################################
# UMAP visualization of epithelial clusters
############################################################

pdf(
  file.path(
    figures_supplementary_dir,
    "Epithelial_UMAP_Clusters.pdf"
  ),
  width = 10,
  height = 7
)

DimPlot(
  epithelial_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  split.by = "Class"
) +
  NoLegend()

dev.off()


############################################################
# Marker genes for epithelial clusters
############################################################

Idents(epithelial_obj) <- "seurat_clusters"

markers <- FindAllMarkers(
  epithelial_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top60 <- markers %>%
  group_by(cluster) %>%
  slice_max(
    avg_log2FC,
    n = 60
  )

write.csv(
  top60,
  file.path(
    tables_supplementary_dir,
    "Epithelial_top60_markers.csv"
  ),
  row.names = FALSE
)


############################################################
# Create balanced epithelial reference
############################################################

set.seed(123)

epithelial_balanced <- subset(
  epithelial_obj,
  cells = unlist(
    lapply(
      split(
        colnames(epithelial_obj),
        epithelial_obj$orig.ident
      ),
      function(x) {
        sample(
          x,
          min(length(x), 1000)
        )
      }
    )
  )
)


############################################################
# Re-process balanced epithelial reference
############################################################

DefaultAssay(epithelial_balanced) <- "RNA"

epithelial_balanced <- NormalizeData(
  epithelial_balanced
)

epithelial_balanced <- FindVariableFeatures(
  epithelial_balanced,
  selection.method = "vst",
  nfeatures = 3000
)

epithelial_balanced <- ScaleData(
  epithelial_balanced
)

epithelial_balanced <- RunPCA(
  epithelial_balanced,
  npcs = 50,
  verbose = FALSE
)


############################################################
# Balanced epithelial clustering
############################################################

epithelial_balanced <- FindNeighbors(
  epithelial_balanced,
  dims = 1:20,
  verbose = FALSE
)

epithelial_balanced <- FindClusters(
  epithelial_balanced,
  resolution = 0.5,
  verbose = FALSE
)

epithelial_balanced <- RunUMAP(
  epithelial_balanced,
  dims = 1:20,
  verbose = FALSE
)


############################################################
# Save balanced epithelial object
############################################################

saveRDS(
  epithelial_balanced,
  file.path(
    objects_dir,
    "TNBC_Epithelial_Balanced.rds"
  ),
  compress = TRUE
)


############################################################
# Balanced epithelial UMAP
############################################################

pdf(
  file.path(
    figures_supplementary_dir,
    "Epithelial_Balanced_UMAP_Clusters.pdf"
  ),
  width = 10,
  height = 7
)

DimPlot(
  epithelial_balanced,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  split.by = "Class"
) +
  NoLegend()

dev.off()


############################################################
# Marker genes for balanced epithelial clusters
############################################################

markers_epi <- FindAllMarkers(
  epithelial_balanced,
  only.pos = TRUE,
  min.pct = 0.25
)

markers_epi_top20 <- markers_epi %>%
  group_by(cluster) %>%
  slice_max(
    avg_log2FC,
    n = 20
  )

write.csv(
  markers_epi_top20,
  file.path(
    tables_supplementary_dir,
    "Epithelial_Balanced_top20_markers.csv"
  ),
  row.names = FALSE
)


############################################################
# Lineage verification markers
############################################################

annotation_markers <- c(
  # Core epithelial
  "EPCAM", "KRT8", "KRT18",
  
  # Basal / myoepithelial
  "KRT14", "KRT5", "COL17A1", "ACTA2",
  
  # EMT / mesenchymal
  "VIM", "VCAN", "CAV1", "DKK3",
  
  # ECM remodeling
  "COL2A1", "COL11A2", "FBLN2", "HSPG2",
  
  # Luminal mature / hormone receptor
  "ESR1", "PGR", "FOXA1", "AR",
  
  # Luminal progenitor / secretory
  "ALDH1A1", "IGFBP5", "RBP4", "CXCL14",
  
  # Apocrine
  "MUCL1", "AGR2", "PIP", "AKR1C2",
  
  # Secretory inflammatory
  "PIGR", "LTF", "CCL28",
  
  # MHC-II epithelial
  "HLA-DQA1", "HLA-DPA1", "HLA-DRB1",
  
  # Cell cycle
  "MKI67", "TOP2A", "PLK1", "CCNB2",
  
  # Immune contamination
  "PTPRC", "CD3D", "CD68"
)


############################################################
# Annotation marker dot plot
############################################################

p_dot <- DotPlot(
  epithelial_balanced,
  features = annotation_markers
) +
  RotatedAxis() +
  ggtitle(
    "Annotation markers across epithelial clusters"
  ) +
  theme(
    axis.text.x = element_text(
      size = 9,
      face = "italic"
    )
  )

ggsave(
  file.path(
    figures_supplementary_dir,
    "Epithelial_DotPlot_AnnotationMarkers.pdf"
  ),
  plot = p_dot,
  width = 15,
  height = 7
)


############################################################
# Load epithelial cluster annotations
############################################################

annotation_file <- file.path(
  metadata_dir,
  "epithelial_cluster_annotations.csv"
)

if (!file.exists(annotation_file)) {
  stop(
    "Annotation file not found: ",
    annotation_file
  )
}

epi_annotation_table <- read.csv(
  annotation_file,
  stringsAsFactors = FALSE
)


############################################################
# Create named annotation vector
# Names = epithelial cluster IDs
############################################################

epi_annotations <- setNames(
  epi_annotation_table$Epithelial_Subtype,
  epi_annotation_table$cluster
)


############################################################
# Add subtype annotations to balanced epithelial reference
############################################################

Idents(epithelial_balanced) <- "seurat_clusters"

epithelial_balanced$Epithelial_Subtype <- as.vector(
  epi_annotations[
    as.character(
      epithelial_balanced$seurat_clusters
    )
  ]
)


############################################################
# Save annotated balanced epithelial reference
############################################################

saveRDS(
  epithelial_balanced,
  file = file.path(
    objects_dir,
    "TNBC_Epithelial_Balanced_Annotated.rds"
  ),
  compress = TRUE
)


############################################################
# UMAP of refined epithelial subtypes
############################################################

pdf(
  file.path(
    figures_main_dir,
    "Epithelial_Balanced_Annotated.pdf"
  ),
  width = 14,
  height = 10
)

DimPlot(
  epithelial_balanced,
  reduction = "umap",
  group.by = "Epithelial_Subtype",
  label = TRUE,
  repel = TRUE,
  split.by = "Class"
) +
  NoLegend()

dev.off()


############################################################
# Transfer refined annotation to all epithelial cells
############################################################

transfer_anchors <- FindTransferAnchors(
  reference = epithelial_balanced,
  query = epithelial_obj,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = epithelial_balanced$Epithelial_Subtype,
  dims = 1:30
)

epithelial_obj <- AddMetaData(
  epithelial_obj,
  metadata = predictions
)

Idents(epithelial_obj) <- epithelial_obj$predicted.id

epithelial_obj$Epithelial_Subtype <- Idents(
  epithelial_obj
)


############################################################
# Save refined epithelial object before purification
############################################################

saveRDS(
  epithelial_obj,
  file = file.path(
    objects_dir,
    "TNBC_Epithelial_Annotated.rds"
  ),
  compress = TRUE
)


############################################################
# Remove likely immune-contaminated epithelial subtype
############################################################

# Remove the predefined subtype representing likely
# immune-contaminated epithelial cells.
#
# The subtype name must match the annotation file exactly.

subtypes_to_remove <- c(
  "Basal epithelial with immune signal"
)

epithelial_obj_pure <- subset(
  epithelial_obj,
  subset = !(
    Epithelial_Subtype %in% subtypes_to_remove
  )
)


############################################################
# Save purified epithelial object
############################################################

saveRDS(
  epithelial_obj_pure,
  file = file.path(
    objects_dir,
    "TNBC_Epithelial_Pure.rds"
  ),
  compress = TRUE
)


############################################################
# Compare cell numbers before and after purification
############################################################

message(
  "Epithelial cells before purification: ",
  ncol(epithelial_obj)
)

message(
  "Epithelial cells after purification: ",
  ncol(epithelial_obj_pure)
)


############################################################
# UMAP of transferred epithelial subtype annotations
############################################################

pdf(
  file.path(
    figures_main_dir,
    "Epithelial_Transferred_Annotations.pdf"
  ),
  width = 12,
  height = 8
)

DimPlot(
  epithelial_obj,
  reduction = "umap",
  group.by = "Epithelial_Subtype",
  label = TRUE,
  repel = TRUE,
  split.by = "Class"
) +
  NoLegend()

dev.off()


############################################################
# Cell-type composition
############################################################

subtype_table <- table(
  epithelial_obj$Epithelial_Subtype
)

write.csv(
  as.data.frame(subtype_table),
  file = file.path(
    tables_main_dir,
    "Epithelial_Subtype_Abundance.csv"
  ),
  row.names = FALSE
)


############################################################
# Summary log
############################################################

log_file <- file.path(
  logs_dir,
  "Step5_Epithelial_Refinement.log"
)

log_con <- file(
  log_file,
  open = "wt"
)

writeLines(
  c(
    "TNBC Single-cell RNA-seq Analysis Pipeline",
    "Step 5: Epithelial compartment re-analysis and refinement",
    paste("Date:", Sys.time()),
    "",
    "Summary",
    "-------",
    paste(
      "Epithelial cells:",
      ncol(epithelial_obj)
    ),
    paste(
      "Balanced epithelial cells:",
      ncol(epithelial_balanced)
    ),
    paste(
      "Number of epithelial clusters:",
      length(
        unique(
          epithelial_balanced$seurat_clusters
        )
      )
    ),
    paste(
      "Purified epithelial cells:",
      ncol(epithelial_obj_pure)
    ),
    paste(
      "Number of epithelial subtypes:",
      length(
        unique(
          epithelial_obj$Epithelial_Subtype
        )
      )
    ),
    "",
    "Cells per epithelial subtype",
    "--------------------------",
    capture.output(
      print(
        table(
          epithelial_obj$Epithelial_Subtype
        )
      )
    ),
    "",
    "Objects generated",
    "-----------------",
    file.path(
      objects_dir,
      "TNBC_Epithelial_Subset.rds"
    ),
    file.path(
      objects_dir,
      "TNBC_Epithelial_Balanced.rds"
    ),
    file.path(
      objects_dir,
      "TNBC_Epithelial_Balanced_Annotated.rds"
    ),
    file.path(
      objects_dir,
      "TNBC_Epithelial_Annotated.rds"
    ),
    file.path(
      objects_dir,
      "TNBC_Epithelial_Pure.rds"
    ),
    "",
    "Key result files",
    "----------------",
    file.path(
      tables_supplementary_dir,
      "Epithelial_top60_markers.csv"
    ),
    file.path(
      tables_supplementary_dir,
      "Epithelial_Balanced_top20_markers.csv"
    ),
    file.path(
      figures_supplementary_dir,
      "Epithelial_DotPlot_AnnotationMarkers.pdf"
    ),
    file.path(
      figures_main_dir,
      "Epithelial_Balanced_Annotated.pdf"
    ),
    file.path(
      figures_main_dir,
      "Epithelial_Transferred_Annotations.pdf"
    ),
    "",
    "Step 5 completed successfully."
  ),
  con = log_con
)

close(log_con)


############################################################
# Completion message
############################################################

message(
  "Step 5 completed successfully."
)

message(
  "Log saved to: ",
  log_file
)