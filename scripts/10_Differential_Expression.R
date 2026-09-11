#!/usr/bin/env Rscript

############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
#
# STEP 10: Differential Expression Analysis
#           Tumor-like Epithelial Cells
#
# Comparison:
# BRCA1_tumour vs TotalCell
#
# IMPORTANT:
# - FindMarkers() is performed in this step
# - Positive avg_log2FC = higher in BRCA1_tumour
# - Negative avg_log2FC = higher in TotalCell
# - FindMarkers() uses min.pct = 0.10
# - FindMarkers() uses logfc.threshold = 0.25
# - Step 10 focuses on differential expression analysis
# - The complete ranked gene list for GSEA is handled separately in Step 11
# - Significant DEGs are defined afterward using FDR < 0.05
#
# NOTE:
# This script is retained for reproducibility and repository
# documentation. Do not rerun during the current refactoring pass.
############################################################


############################################################
# 10.1 Load libraries
############################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)


############################################################
# 10.2 Project directories
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
  "10_Differential_Expression"
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
# 10.3 Load tumor-like epithelial object
############################################################

epithelial_obj_tumor <- readRDS(
  file.path(
    objects_dir,
    "TNBC_Epithelial_TumorLike.rds"
  )
)


############################################################
# 10.4 Verify metadata
############################################################

if (
  !"Class" %in%
  colnames(
    epithelial_obj_tumor@meta.data
  )
) {
  
  stop(
    "Metadata column 'Class' was not found."
  )
}


required_classes <- c(
  "BRCA1_tumour",
  "TotalCell"
)

available_classes <- unique(
  epithelial_obj_tumor$Class
)

if (
  !all(
    required_classes %in%
    available_classes
  )
) {
  
  stop(
    "Required identities were not found. ",
    "Expected: BRCA1_tumour and TotalCell. ",
    "Available identities: ",
    paste(
      available_classes,
      collapse = ", "
    )
  )
}


############################################################
# 10.5 Set identities
############################################################

Idents(
  epithelial_obj_tumor
) <- "Class"


############################################################
# 10.6 Cell counts
############################################################

n_brca1 <- sum(
  epithelial_obj_tumor$Class ==
    "BRCA1_tumour"
)

n_totalcell <- sum(
  epithelial_obj_tumor$Class ==
    "TotalCell"
)


############################################################
# 10.7 Differential expression
############################################################

deg_tumor <- FindMarkers(
  object = epithelial_obj_tumor,
  ident.1 = "BRCA1_tumour",
  ident.2 = "TotalCell",
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  only.pos = FALSE
)


############################################################
# 10.8 Add gene names
############################################################

deg_tumor$gene <- rownames(
  deg_tumor
)

rownames(
  deg_tumor
) <- NULL


############################################################
# 10.9 Verify DEG output
############################################################

required_deg_columns <- c(
  "p_val",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj",
  "gene"
)

missing_columns <- setdiff(
  required_deg_columns,
  colnames(deg_tumor)
)

if (
  length(missing_columns) > 0
) {
  
  stop(
    "Missing expected columns in FindMarkers() output: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


############################################################
# 10.10 Add explicit biological direction
############################################################
#
# Positive avg_log2FC:
# BRCA1_tumour > TotalCell
#
# Negative avg_log2FC:
# TotalCell > BRCA1_tumour
############################################################

deg_tumor <- deg_tumor %>%
  mutate(
    DEG_direction = case_when(
      
      avg_log2FC > 0 ~
        "Higher_in_BRCA1_tumour",
      
      avg_log2FC < 0 ~
        "Higher_in_TotalCell",
      
      TRUE ~
        "No_difference"
    )
  )


############################################################
# 10.11 Detection difference
############################################################
#
# delta_pct = pct.1 - pct.2
############################################################

deg_tumor <- deg_tumor %>%
  mutate(
    delta_pct =
      pct.1 - pct.2
  )


############################################################
# 10.12 Order complete DEG table
############################################################

deg_tumor <- deg_tumor %>%
  arrange(
    desc(avg_log2FC)
  )


############################################################
# 10.13 Significant DEGs
############################################################

deg_sig <- deg_tumor %>%
  filter(
    p_val_adj < 0.05
  )


############################################################
# 10.14 Genes higher in BRCA1_tumour
############################################################

deg_up <- deg_sig %>%
  filter(
    avg_log2FC > 0
  ) %>%
  arrange(
    desc(avg_log2FC)
  )


############################################################
# 10.15 Genes higher in TotalCell
############################################################

deg_down <- deg_sig %>%
  filter(
    avg_log2FC < 0
  ) %>%
  arrange(
    avg_log2FC
  )


############################################################
# 10.16 Save complete DEG table
############################################################

write.csv(
  deg_tumor,
  file.path(
    tables_supp_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_all.csv"
  ),
  row.names = FALSE
)


############################################################
# 10.17 Save significant DEG table
############################################################

write.csv(
  deg_sig,
  file.path(
    tables_main_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_significant.csv"
  ),
  row.names = FALSE
)


############################################################
# 10.18 Save BRCA1_tumour-enriched genes
############################################################

write.csv(
  deg_up,
  file.path(
    tables_main_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_UP.csv"
  ),
  row.names = FALSE
)


############################################################
# 10.19 Save TotalCell-enriched genes
############################################################

write.csv(
  deg_down,
  file.path(
    tables_main_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_DOWN.csv"
  ),
  row.names = FALSE
)


############################################################
# 10.20 Top 20 genes
############################################################

top_up <- deg_up %>%
  slice_head(n = 20)

top_down <- deg_down %>%
  slice_head(n = 20)


############################################################
# 10.21 Save Top 20 genes
############################################################

write.csv(
  top_up,
  file.path(
    tables_supp_dir,
    "DEG_Top20_Upregulated_BRCA1.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_down,
  file.path(
    tables_supp_dir,
    "DEG_Top20_Downregulated_TotalCell.csv"
  ),
  row.names = FALSE
)


############################################################
# 10.22 Volcano plot data
############################################################

volcano_data <- deg_tumor %>%
  mutate(
    
    log10_padj =
      -log10(
        p_val_adj + 1e-300
      ),
    
    Significant = case_when(
      
      p_val_adj < 0.05 &
        avg_log2FC >= 0.5 ~
        "Higher in BRCA1_tumour",
      
      p_val_adj < 0.05 &
        avg_log2FC <= -0.5 ~
        "Higher in TotalCell",
      
      TRUE ~
        "NS"
    )
  )


############################################################
# 10.23 Volcano labels — BRCA1_tumour
############################################################

label_up <- volcano_data %>%
  filter(
    
    Significant ==
      "Higher in BRCA1_tumour",
    
    pct.1 >= 0.20,
    
    delta_pct >= 0.10
    
  ) %>%
  arrange(
    
    desc(
      abs(avg_log2FC)
    ),
    
    desc(delta_pct)
    
  ) %>%
  slice_head(
    n = 8
  )


############################################################
# 10.24 Volcano labels — TotalCell
############################################################

label_down <- volcano_data %>%
  filter(
    
    Significant ==
      "Higher in TotalCell",
    
    pct.2 >= 0.20,
    
    delta_pct <= -0.10
    
  ) %>%
  arrange(
    
    avg_log2FC,
    
    delta_pct
    
  ) %>%
  slice_head(
    n = 8
  )


############################################################
# 10.25 Fallback labels
############################################################

if (
  nrow(label_up) < 4
) {
  
  label_up <- volcano_data %>%
    filter(
      
      Significant ==
        "Higher in BRCA1_tumour",
      
      pct.1 >= 0.10
      
    ) %>%
    arrange(
      
      desc(
        abs(avg_log2FC)
      )
      
    ) %>%
    slice_head(
      n = 8
    )
}


if (
  nrow(label_down) < 4
) {
  
  label_down <- volcano_data %>%
    filter(
      
      Significant ==
        "Higher in TotalCell",
      
      pct.2 >= 0.10
      
    ) %>%
    arrange(
      avg_log2FC
    ) %>%
    slice_head(
      n = 8
    )
}


############################################################
# 10.26 Combine labels
############################################################

volcano_labels <- bind_rows(
  label_up,
  label_down
) %>%
  distinct(
    gene,
    .keep_all = TRUE
  )


############################################################
# 10.27 Volcano plot
############################################################

p1 <- ggplot(
  volcano_data,
  aes(
    x = avg_log2FC,
    y = log10_padj,
    color = Significant
  )
) +
  
  geom_point(
    size = 0.8,
    alpha = 0.75
  ) +
  
  scale_color_manual(
    values = c(
      "Higher in BRCA1_tumour" = "#D73027",
      "Higher in TotalCell" = "#4575B4",
      "NS" = "grey70"
    )
  ) +
  
  geom_vline(
    xintercept = c(
      -0.5,
      0.5
    ),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  
  geom_text_repel(
    data = volcano_labels,
    aes(
      label = gene
    ),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.25,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  
  theme_classic() +
  
  labs(
    title = "BRCA1_tumour vs TotalCell",
    subtitle = "Tumor-like epithelial cells",
    x = "Average log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Expression direction"
  )


############################################################
# 10.28 Save Volcano
############################################################

ggsave(
  file.path(
    figures_main_dir,
    "Volcano_BRCA1_vs_TotalCell_TumorLike.pdf"
  ),
  p1,
  width = 10,
  height = 7
)


############################################################
# 10.29 Heatmap gene selection
############################################################

top_heatmap_genes <- unique(
  c(
    top_up$gene,
    top_down$gene
  )
)

top_heatmap_genes <-
  top_heatmap_genes[
    top_heatmap_genes %in%
      rownames(
        epithelial_obj_tumor
      )
  ]


############################################################
# 10.30 Heatmap
############################################################

if (
  length(top_heatmap_genes) > 0
) {
  
  p2 <- DoHeatmap(
    epithelial_obj_tumor,
    features = top_heatmap_genes,
    group.by = "Class"
  ) +
    ggtitle(
      "Top differentially expressed genes"
    )
  
  ggsave(
    file.path(
      figures_main_dir,
      "Heatmap_Top_DEGs_TumorLike.pdf"
    ),
    p2,
    width = 11,
    height = 10
  )
  
} else {
  
  warning(
    "No genes available for heatmap."
  )
}


############################################################
# 10.31 Summary statistics
############################################################

n_tested <-
  nrow(deg_tumor)

n_sig <-
  sum(
    deg_tumor$p_val_adj < 0.05
  )

n_up <-
  sum(
    deg_tumor$p_val_adj < 0.05 &
      deg_tumor$avg_log2FC > 0
  )

n_down <-
  sum(
    deg_tumor$p_val_adj < 0.05 &
      deg_tumor$avg_log2FC < 0
  )

n_ns <-
  n_tested - n_sig


############################################################
# 10.32 Save summary CSV
############################################################

step10_summary <- data.frame(
  
  Metric = c(
    
    "Comparison",
    
    "BRCA1_tumour cells",
    
    "TotalCell cells",
    
    "Genes tested",
    
    "Significant DEGs",
    
    "Higher in BRCA1_tumour",
    
    "Higher in TotalCell",
    
    "Non-significant genes",
    
    "Volcano labels"
    
  ),
  
  Value = c(
    
    "BRCA1_tumour vs TotalCell",
    
    n_brca1,
    
    n_totalcell,
    
    n_tested,
    
    n_sig,
    
    n_up,
    
    n_down,
    
    n_ns,
    
    nrow(volcano_labels)
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  step10_summary,
  file.path(
    tables_main_dir,
    "Step10_DEG_summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 10.33 Step 10 log
############################################################

log_file <- file.path(
  logs_dir,
  "Step10_DEG_TumorLike.log"
)


log_lines <- c(
  
  "====================================================",
  
  "TNBC Single-cell RNA-seq Analysis Pipeline",
  
  "STEP 10 — Differential Expression Analysis",
  
  "====================================================",
  
  "",
  
  paste(
    "Date/time:",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  "COMPARISON",
  
  "----------------------------------------------------",
  
  "BRCA1_tumour vs TotalCell",
  
  "",
  
  "FINDMARKERS PARAMETERS",
  
  "----------------------------------------------------",
  
  "Test: Wilcoxon rank-sum test",
  
  "min.pct: 0.10",
  
  "logfc.threshold: 0.25",
  
  "only.pos: FALSE",
  
  "",
  
  "DIRECTION",
  
  "----------------------------------------------------",
  
  "Positive avg_log2FC = higher in BRCA1_tumour",
  
  "Negative avg_log2FC = higher in TotalCell",
  
  "Positive delta_pct = greater detection in BRCA1_tumour",
  
  "Negative delta_pct = greater detection in TotalCell",
  
  "",
  
  "DEG SUMMARY",
  
  "----------------------------------------------------",
  
  paste(
    "Genes tested:",
    n_tested
  ),
  
  paste(
    "Significant DEGs (FDR < 0.05):",
    n_sig
  ),
  
  paste(
    "Higher in BRCA1_tumour:",
    n_up
  ),
  
  paste(
    "Higher in TotalCell:",
    n_down
  ),
  
  paste(
    "Non-significant:",
    n_ns
  ),
  
  "",
  
  "VOLCANO",
  
  "----------------------------------------------------",
  
  "Visualization threshold: FDR < 0.05 and |log2FC| >= 0.5",
  
  "Labels incorporate expression prevalence and detection difference.",
  
  "",
  
  "OUTPUTS",
  
  "----------------------------------------------------",
  
  file.path(
    tables_supp_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_all.csv"
  ),
  
  file.path(
    tables_main_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_significant.csv"
  ),
  
  file.path(
    tables_main_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_UP.csv"
  ),
  
  file.path(
    tables_main_dir,
    "DEG_TumorLike_BRCA1_vs_TotalCell_DOWN.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "DEG_Top20_Upregulated_BRCA1.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "DEG_Top20_Downregulated_TotalCell.csv"
  ),
  
  file.path(
    tables_main_dir,
    "Step10_DEG_summary.csv"
  ),
  
  file.path(
    figures_main_dir,
    "Volcano_BRCA1_vs_TotalCell_TumorLike.pdf"
  ),
  
  file.path(
    figures_main_dir,
    "Heatmap_Top_DEGs_TumorLike.pdf"
  ),
  
  "",
  
  "STEP 10 COMPLETED SUCCESSFULLY.",
  
  "===================================================="
  
)


writeLines(
  log_lines,
  con = log_file
)


############################################################
# 10.34 Console summary
############################################################

message(
  "Step 10 completed successfully."
)

message(
  "Genes tested: ",
  n_tested
)

message(
  "Significant DEGs: ",
  n_sig
)

message(
  "Higher in BRCA1_tumour: ",
  n_up
)

message(
  "Higher in TotalCell: ",
  n_down
)

message(
  "Log saved to: ",
  log_file
)


############################################################
# End of Step 10
############################################################