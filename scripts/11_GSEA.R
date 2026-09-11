############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
#
# STEP 11: GSEA
#           Tumor-like Epithelial Cells
#
# Comparison:
# BRCA1_tumour vs TotalCell
#
# Repository-ready canonical version
#
# IMPORTANT:
# - Step 11 is independent from the Step 10 R environment.
# - The DEG table is loaded from the canonical Step 10 output.
# - All eligible genes from Step 10 are used for ranking.
# - No |log2FC| cutoff is applied before GSEA.
# - No p-value filtering is applied before GSEA.
# - Ranking variable: avg_log2FC.
# - Positive ranking = higher in BRCA1_tumour.
# - Negative ranking = higher in TotalCell.
# - Final significance criterion: adjusted p-value < 0.05.
# - Multiple-testing correction: Benjamini-Hochberg (BH).
#
# GENE SET DATABASES:
# - GO Biological Process
# - Reactome
# - KEGG
# - MSigDB Hallmark
#
# VISUALIZATION:
# - Complete GSEA results are saved as CSV files.
# - Significant results are separated by enrichment direction.
# - Dot plots display up to 10 representative significant terms
#   per enrichment direction.
#
# BIOLOGICAL INTERPRETATION:
# - NES > 0 = enrichment toward BRCA1_tumour.
# - NES < 0 = enrichment toward TotalCell.
# - Statistical significance = adjusted p-value < 0.05.
############################################################


############################################################
# 11.1 Load libraries
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(ggplot2)
library(dplyr)
library(stringr)
library(msigdbr)


############################################################
# 11.2 Project directories
############################################################

# ==========================================================
# Repository root
# ==========================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"


# ==========================================================
# Stage 11 directory
# ==========================================================

stage_dir <- file.path(
  project_dir,
  "results",
  "11_GSEA"
)


# ==========================================================
# Stage 11 output directories
# ==========================================================

figures_main_dir <- file.path(
  stage_dir,
  "figures",
  "main"
)

tables_main_dir <- file.path(
  stage_dir,
  "tables",
  "main"
)

tables_supp_dir <- file.path(
  stage_dir,
  "tables",
  "supplementary"
)

tables_validation_dir <- file.path(
  stage_dir,
  "tables",
  "validation"
)

logs_dir <- file.path(
  stage_dir,
  "logs"
)


# ==========================================================
# Create required directories
# ==========================================================

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
  tables_validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 11.3 Load Step 10 DEG table
############################################################
#
# Canonical Stage 10 input:
#
# results/
#   10_Differential_Expression/
#     tables/
#       supplementary/
#         DEG_TumorLike_BRCA1_vs_TotalCell_all.csv
############################################################

deg_file <- file.path(
  project_dir,
  "results",
  "10_Differential_Expression",
  "tables",
  "supplementary",
  "DEG_TumorLike_BRCA1_vs_TotalCell_all.csv"
)


if (!file.exists(deg_file)) {
  
  stop(
    "Stage 10 DEG file was not found:\n",
    deg_file,
    "\nPlease run Stage 10 first."
  )
}


deg_tumor <- read.csv(
  deg_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 11.4 Verify required columns
############################################################

required_columns <- c(
  "gene",
  "avg_log2FC",
  "p_val_adj"
)


missing_columns <- setdiff(
  required_columns,
  colnames(deg_tumor)
)


if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns in Stage 10 DEG table: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


############################################################
# 11.5 Input summary
############################################################

n_input <- nrow(
  deg_tumor
)


message(
  "Genes loaded from Stage 10: ",
  n_input
)


############################################################
# 11.6 Prepare complete ranked gene table
############################################################
#
# ALL eligible genes are retained for GSEA ranking.
#
# No p-value filtering.
# No p.adjust filtering.
# No |log2FC| filtering.
#
# Ranking variable:
# avg_log2FC
#
# Positive = BRCA1_tumour
# Negative = TotalCell
############################################################

gene_rank_df <- deg_tumor %>%
  
  select(
    gene,
    avg_log2FC
  ) %>%
  
  filter(
    !is.na(gene),
    gene != "",
    !is.na(avg_log2FC),
    is.finite(avg_log2FC)
  )


############################################################
# 11.7 Remove duplicated gene symbols
############################################################

gene_rank_df <- gene_rank_df %>%
  
  mutate(
    abs_log2FC = abs(avg_log2FC)
  ) %>%
  
  arrange(
    desc(abs_log2FC)
  ) %>%
  
  distinct(
    gene,
    .keep_all = TRUE
  ) %>%
  
  select(
    -abs_log2FC
  )


############################################################
# 11.8 Map gene symbols to ENTREZ IDs
############################################################

gene_id_map <- bitr(
  
  gene_rank_df$gene,
  
  fromType = "SYMBOL",
  
  toType = "ENTREZID",
  
  OrgDb = org.Hs.eg.db
  
)


############################################################
# 11.9 Merge ENTREZ IDs
############################################################

gene_rank_df <- gene_rank_df %>%
  
  inner_join(
    gene_id_map,
    by = c(
      "gene" = "SYMBOL"
    )
  )


############################################################
# 11.10 Resolve duplicated ENTREZ IDs
############################################################

gene_rank_df <- gene_rank_df %>%
  
  mutate(
    abs_log2FC = abs(avg_log2FC)
  ) %>%
  
  arrange(
    desc(abs_log2FC)
  ) %>%
  
  distinct(
    ENTREZID,
    .keep_all = TRUE
  ) %>%
  
  select(
    -abs_log2FC
  )


############################################################
# 11.11 Sort by log2FC
############################################################

gene_rank_df <- gene_rank_df %>%
  
  arrange(
    desc(avg_log2FC)
  )


############################################################
# 11.12 Create final named numeric vector
############################################################

gene_rank_entrez <- gene_rank_df$avg_log2FC

names(
  gene_rank_entrez
) <- gene_rank_df$ENTREZID


############################################################
# 11.13 Final sorting
############################################################

gene_rank_entrez <- sort(
  gene_rank_entrez,
  decreasing = TRUE
)


############################################################
# 11.14 Validate GSEA input
############################################################

if (!is.numeric(gene_rank_entrez)) {
  
  stop(
    "GSEA ranking is not numeric."
  )
}


if (is.null(names(gene_rank_entrez))) {
  
  stop(
    "GSEA ranking has no ENTREZID names."
  )
}


if (anyDuplicated(names(gene_rank_entrez)) > 0) {
  
  stop(
    "Duplicate ENTREZ IDs remain."
  )
}


if (!all(diff(gene_rank_entrez) <= 0)) {
  
  stop(
    "GSEA ranking is not sorted decreasingly."
  )
}


############################################################
# 11.15 GSEA input summary
############################################################

n_ranked <- length(
  gene_rank_entrez
)

n_positive <- sum(
  gene_rank_entrez > 0
)

n_negative <- sum(
  gene_rank_entrez < 0
)

n_zero <- sum(
  gene_rank_entrez == 0
)

n_mapped <- n_ranked

mapping_rate <- 100 *
  n_mapped /
  n_input


if (n_ranked < 100) {
  
  stop(
    "Too few genes available for GSEA after ENTREZ mapping: ",
    n_ranked
  )
}


message(
  "Genes entering GSEA after ENTREZ mapping: ",
  n_ranked
)

message(
  "Positive-ranked genes: ",
  n_positive
)

message(
  "Negative-ranked genes: ",
  n_negative
)

message(
  "Zero-ranked genes: ",
  n_zero
)

message(
  "ENTREZ mapping rate: ",
  round(
    mapping_rate,
    2
  ),
  "%"
)

message(
  "GSEA ranking validation: PASSED"
)


############################################################
# 11.16 Save ranked gene list
############################################################

ranked_gene_table <- data.frame(
  
  ENTREZID = names(
    gene_rank_entrez
  ),
  
  avg_log2FC = as.numeric(
    gene_rank_entrez
  )
)


write.csv(
  
  ranked_gene_table,
  
  file.path(
    tables_main_dir,
    "GSEA_RankedGenes_BRCA1_vs_TotalCell_ENTREZ.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 11.17 GSEA — GO Biological Process
############################################################
#
# pvalueCutoff = 1 retains the complete GSEA result.
# Final significance is determined independently by:
#
# p.adjust < 0.05
############################################################

gsea_go_bp <- gseGO(
  
  geneList = gene_rank_entrez,
  
  OrgDb = org.Hs.eg.db,
  
  keyType = "ENTREZID",
  
  ont = "BP",
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 1,
  
  pAdjustMethod = "BH",
  
  verbose = FALSE
)


############################################################
# 11.18 GSEA — Reactome
############################################################

gsea_reactome <- gsePathway(
  
  geneList = gene_rank_entrez,
  
  organism = "human",
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 1,
  
  pAdjustMethod = "BH",
  
  verbose = FALSE
)


############################################################
# 11.19 GSEA — KEGG
############################################################
#
# pvalueCutoff = 1 is used to retain the complete result.
# Significance is filtered later using p.adjust < 0.05.
############################################################

gsea_kegg <- gseKEGG(
  
  geneList = gene_rank_entrez,
  
  organism = "hsa",
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 1,
  
  pAdjustMethod = "BH",
  
  verbose = FALSE
)


############################################################
# 11.20 GSEA — Hallmark
############################################################
#
# Complete GSEA results are retained.
############################################################

hallmark_df <- msigdbr(
  
  species = "Homo sapiens",
  
  collection = "H"
  
) %>%
  
  select(
    gs_name,
    gene_symbol
  ) %>%
  
  inner_join(
    
    gene_rank_df %>%
      select(
        gene,
        ENTREZID
      ),
    
    by = c(
      "gene_symbol" = "gene"
    )
  )


gsea_hallmark <- GSEA(
  
  geneList = gene_rank_entrez,
  
  TERM2GENE = hallmark_df %>%
    select(
      gs_name,
      ENTREZID
    ),
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 1,
  
  pAdjustMethod = "BH",
  
  verbose = FALSE
)


############################################################
# 11.21 Extract GSEA results
############################################################

gsea_go_df <- as.data.frame(
  gsea_go_bp
)

gsea_reactome_df <- as.data.frame(
  gsea_reactome
)

gsea_kegg_df <- as.data.frame(
  gsea_kegg
)

gsea_hallmark_df <- as.data.frame(
  gsea_hallmark
)


############################################################
# 11.22 Add biological direction
############################################################

add_direction <- function(df) {
  
  if (nrow(df) > 0) {
    
    df <- df %>%
      
      mutate(
        
        Enrichment_Group = case_when(
          
          NES > 0 ~ "BRCA1_tumour",
          
          NES < 0 ~ "TotalCell",
          
          TRUE ~ "No_direction"
        )
      )
  }
  
  return(df)
}


gsea_go_df <- add_direction(
  gsea_go_df
)

gsea_reactome_df <- add_direction(
  gsea_reactome_df
)

gsea_kegg_df <- add_direction(
  gsea_kegg_df
)

gsea_hallmark_df <- add_direction(
  gsea_hallmark_df
)


############################################################
# 11.23 Significant GSEA results
############################################################

gsea_go_sig <- gsea_go_df %>%
  filter(
    p.adjust < 0.05
  )

gsea_reactome_sig <- gsea_reactome_df %>%
  filter(
    p.adjust < 0.05
  )

gsea_kegg_sig <- gsea_kegg_df %>%
  filter(
    p.adjust < 0.05
  )

gsea_hallmark_sig <- gsea_hallmark_df %>%
  filter(
    p.adjust < 0.05
  )


############################################################
# 11.24 Separate enrichment directions
############################################################

gsea_go_brca1 <- gsea_go_sig %>%
  filter(
    NES > 0
  ) %>%
  arrange(
    p.adjust
  )

gsea_go_totalcell <- gsea_go_sig %>%
  filter(
    NES < 0
  ) %>%
  arrange(
    p.adjust
  )


gsea_reactome_brca1 <- gsea_reactome_sig %>%
  filter(
    NES > 0
  ) %>%
  arrange(
    p.adjust
  )

gsea_reactome_totalcell <- gsea_reactome_sig %>%
  filter(
    NES < 0
  ) %>%
  arrange(
    p.adjust
  )


gsea_kegg_brca1 <- gsea_kegg_sig %>%
  filter(
    NES > 0
  ) %>%
  arrange(
    p.adjust
  )

gsea_kegg_totalcell <- gsea_kegg_sig %>%
  filter(
    NES < 0
  ) %>%
  arrange(
    p.adjust
  )


gsea_hallmark_brca1 <- gsea_hallmark_sig %>%
  filter(
    NES > 0
  ) %>%
  arrange(
    p.adjust
  )

gsea_hallmark_totalcell <- gsea_hallmark_sig %>%
  filter(
    NES < 0
  ) %>%
  arrange(
    p.adjust
  )


############################################################
# 11.25 Save complete GSEA results
############################################################

write.csv(
  gsea_go_df,
  file.path(
    tables_main_dir,
    "GSEA_GO_BP_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_reactome_df,
  file.path(
    tables_main_dir,
    "GSEA_Reactome_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_kegg_df,
  file.path(
    tables_main_dir,
    "GSEA_KEGG_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_hallmark_df,
  file.path(
    tables_main_dir,
    "GSEA_Hallmark_all.csv"
  ),
  row.names = FALSE
)


############################################################
# 11.26 Save significant GSEA results
############################################################

write.csv(
  gsea_go_sig,
  file.path(
    tables_supp_dir,
    "GSEA_GO_BP_significant_padj_0.05.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_reactome_sig,
  file.path(
    tables_supp_dir,
    "GSEA_Reactome_significant_padj_0.05.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_kegg_sig,
  file.path(
    tables_supp_dir,
    "GSEA_KEGG_significant_padj_0.05.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_hallmark_sig,
  file.path(
    tables_supp_dir,
    "GSEA_Hallmark_significant_padj_0.05.csv"
  ),
  row.names = FALSE
)


############################################################
# 11.27 Save directional results
############################################################
#
# All four databases are saved by enrichment direction.
############################################################

write.csv(
  gsea_go_brca1,
  file.path(
    tables_supp_dir,
    "GSEA_GO_BP_BRCA1_tumour.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_go_totalcell,
  file.path(
    tables_supp_dir,
    "GSEA_GO_BP_TotalCell.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_reactome_brca1,
  file.path(
    tables_supp_dir,
    "GSEA_Reactome_BRCA1_tumour.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_reactome_totalcell,
  file.path(
    tables_supp_dir,
    "GSEA_Reactome_TotalCell.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_kegg_brca1,
  file.path(
    tables_supp_dir,
    "GSEA_KEGG_BRCA1_tumour.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_kegg_totalcell,
  file.path(
    tables_supp_dir,
    "GSEA_KEGG_TotalCell.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_hallmark_brca1,
  file.path(
    tables_supp_dir,
    "GSEA_Hallmark_BRCA1_tumour.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_hallmark_totalcell,
  file.path(
    tables_supp_dir,
    "GSEA_Hallmark_TotalCell.csv"
  ),
  row.names = FALSE
)


############################################################
# 11.28 GSEA dot plot function
############################################################
#
# Maximum of 10 representative significant terms
# per enrichment direction.
#
# Selection:
# - Top 10 BRCA1_tumour-enriched terms by adjusted p-value.
# - Top 10 TotalCell-enriched terms by adjusted p-value.
#
# If fewer than 10 significant terms are available in a
# direction, all available terms are displayed.
############################################################

make_gsea_dotplot <- function(
    gsea_df,
    title_text,
    output_file,
    width = 11,
    height = 8,
    wrap_width = 55,
    top_n = 10
) {
  
  if (nrow(gsea_df) == 0) {
    
    message(
      "No significant pathways for: ",
      title_text
    )
    
    return(NULL)
  }
  
  
  plot_df <- gsea_df %>%
    
    filter(
      NES != 0
    ) %>%
    
    mutate(
      Direction = case_when(
        
        NES > 0 ~ "BRCA1_tumour",
        
        NES < 0 ~ "TotalCell",
        
        TRUE ~ "No_direction"
      )
    ) %>%
    
    group_by(
      Direction
    ) %>%
    
    arrange(
      p.adjust,
      .by_group = TRUE
    ) %>%
    
    slice_head(
      n = top_n
    ) %>%
    
    ungroup()
  
  
  if (nrow(plot_df) == 0) {
    
    message(
      "No plottable significant pathways for: ",
      title_text
    )
    
    return(NULL)
  }
  
  
  plot_df <- plot_df %>%
    
    mutate(
      
      neg_log10_padj =
        -log10(
          pmax(
            p.adjust,
            .Machine$double.xmin
          )
        ),
      
      Description_wrapped =
        str_wrap(
          Description,
          width = wrap_width
        )
    )
  
  
  plot_df$Description_wrapped <- factor(
    plot_df$Description_wrapped,
    levels = plot_df$Description_wrapped[
      order(
        plot_df$NES
      )
    ]
  )
  
  
  p <- ggplot(
    
    plot_df,
    
    aes(
      x = NES,
      y = Description_wrapped,
      size = abs(NES),
      color = neg_log10_padj
    )
    
  ) +
    
    geom_point(
      alpha = 0.85
    ) +
    
    scale_color_viridis_c(
      name = "-log10 adjusted p-value"
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      
      axis.text.y =
        element_text(
          size = 9
        ),
      
      plot.title =
        element_text(
          face = "bold",
          size = 14
        ),
      
      plot.subtitle =
        element_text(
          size = 10
        )
    ) +
    
    labs(
      
      title = title_text,
      
      subtitle =
        "Top 10 significant terms per enrichment direction",
      
      x =
        "Normalized Enrichment Score (NES)",
      
      y = NULL,
      
      size = "|NES|"
    )
  
  
  plot_height <- max(
    height,
    min(
      16,
      0.35 *
        nrow(plot_df)
    )
  )
  
  
  ggsave(
    
    filename = output_file,
    
    plot = p,
    
    width = width,
    
    height = plot_height,
    
    units = "in"
  )
  
  
  return(p)
}


############################################################
# 11.29 GO BP dot plot
############################################################

p_go <- make_gsea_dotplot(
  
  gsea_go_sig,
  
  "GSEA — GO Biological Process",
  
  file.path(
    figures_main_dir,
    "GSEA_GO_BP_dotplot.pdf"
  )
)


############################################################
# 11.30 Reactome dot plot
############################################################

p_reactome <- make_gsea_dotplot(
  
  gsea_reactome_sig,
  
  "GSEA — Reactome",
  
  file.path(
    figures_main_dir,
    "GSEA_Reactome_dotplot.pdf"
  )
)


############################################################
# 11.31 KEGG dot plot
############################################################

p_kegg <- make_gsea_dotplot(
  
  gsea_kegg_sig,
  
  "GSEA — KEGG",
  
  file.path(
    figures_main_dir,
    "GSEA_KEGG_dotplot.pdf"
  )
)


############################################################
# 11.32 Hallmark dot plot
############################################################

p_hallmark <- make_gsea_dotplot(
  
  gsea_hallmark_sig,
  
  "GSEA — Hallmark",
  
  file.path(
    figures_main_dir,
    "GSEA_Hallmark_dotplot.pdf"
  )
)


############################################################
# 11.33 Summary statistics
############################################################

# ----------------------------------------------------------
# GO Biological Process
# ----------------------------------------------------------

n_go_total <- nrow(
  gsea_go_df
)

n_go_sig <- nrow(
  gsea_go_sig
)

n_go_brca1 <- nrow(
  gsea_go_brca1
)

n_go_totalcell <- nrow(
  gsea_go_totalcell
)


# ----------------------------------------------------------
# Reactome
# ----------------------------------------------------------

n_reactome_total <- nrow(
  gsea_reactome_df
)

n_reactome_sig <- nrow(
  gsea_reactome_sig
)

n_reactome_brca1 <- nrow(
  gsea_reactome_brca1
)

n_reactome_totalcell <- nrow(
  gsea_reactome_totalcell
)


# ----------------------------------------------------------
# KEGG
# ----------------------------------------------------------

n_kegg_total <- nrow(
  gsea_kegg_df
)

n_kegg_sig <- nrow(
  gsea_kegg_sig
)

n_kegg_brca1 <- nrow(
  gsea_kegg_brca1
)

n_kegg_totalcell <- nrow(
  gsea_kegg_totalcell
)


# ----------------------------------------------------------
# Hallmark
# ----------------------------------------------------------

n_hallmark_total <- nrow(
  gsea_hallmark_df
)

n_hallmark_sig <- nrow(
  gsea_hallmark_sig
)

n_hallmark_brca1 <- nrow(
  gsea_hallmark_brca1
)

n_hallmark_totalcell <- nrow(
  gsea_hallmark_totalcell
)


############################################################
# 11.34 Create Stage 11 summary table
############################################################

summary_table <- data.frame(
  
  Database = c(
    "GO Biological Process",
    "Reactome",
    "KEGG",
    "MSigDB Hallmark"
  ),
  
  Pathways_or_GeneSets_Tested = c(
    n_go_total,
    n_reactome_total,
    n_kegg_total,
    n_hallmark_total
  ),
  
  Significant = c(
    n_go_sig,
    n_reactome_sig,
    n_kegg_sig,
    n_hallmark_sig
  ),
  
  BRCA1_tumour_Enriched = c(
    n_go_brca1,
    n_reactome_brca1,
    n_kegg_brca1,
    n_hallmark_brca1
  ),
  
  TotalCell_Enriched = c(
    n_go_totalcell,
    n_reactome_totalcell,
    n_kegg_totalcell,
    n_hallmark_totalcell
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  summary_table,
  
  file.path(
    tables_validation_dir,
    "Step11_GSEA_summary.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 11.35 Validation summary
############################################################

validation_table <- data.frame(
  
  Metric = c(
    
    "Genes loaded from Stage 10",
    
    "Genes entering GSEA after ENTREZ mapping",
    
    "Positive-ranked genes",
    
    "Negative-ranked genes",
    
    "Zero-ranked genes",
    
    "ENTREZ mapping rate (%)",
    
    "GSEA ranking validation",
    
    "GO BP significant pathways",
    
    "Reactome significant pathways",
    
    "KEGG significant pathways",
    
    "Hallmark significant gene sets"
    
  ),
  
  Value = c(
    
    n_input,
    
    n_ranked,
    
    n_positive,
    
    n_negative,
    
    n_zero,
    
    round(
      mapping_rate,
      2
    ),
    
    "PASSED",
    
    n_go_sig,
    
    n_reactome_sig,
    
    n_kegg_sig,
    
    n_hallmark_sig
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  validation_table,
  
  file.path(
    tables_validation_dir,
    "Step11_GSEA_validation_summary.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 11.36 Verify expected outputs
############################################################

expected_files <- c(
  
  # --------------------------------------------------------
  # Main figures
  # --------------------------------------------------------
  
  file.path(
    figures_main_dir,
    "GSEA_GO_BP_dotplot.pdf"
  ),
  
  file.path(
    figures_main_dir,
    "GSEA_Reactome_dotplot.pdf"
  ),
  
  file.path(
    figures_main_dir,
    "GSEA_KEGG_dotplot.pdf"
  ),
  
  file.path(
    figures_main_dir,
    "GSEA_Hallmark_dotplot.pdf"
  ),
  
  
  # --------------------------------------------------------
  # Main tables
  # --------------------------------------------------------
  
  file.path(
    tables_main_dir,
    "GSEA_RankedGenes_BRCA1_vs_TotalCell_ENTREZ.csv"
  ),
  
  file.path(
    tables_main_dir,
    "GSEA_GO_BP_all.csv"
  ),
  
  file.path(
    tables_main_dir,
    "GSEA_Reactome_all.csv"
  ),
  
  file.path(
    tables_main_dir,
    "GSEA_KEGG_all.csv"
  ),
  
  file.path(
    tables_main_dir,
    "GSEA_Hallmark_all.csv"
  ),
  
  
  # --------------------------------------------------------
  # Significant tables
  # --------------------------------------------------------
  
  file.path(
    tables_supp_dir,
    "GSEA_GO_BP_significant_padj_0.05.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_Reactome_significant_padj_0.05.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_KEGG_significant_padj_0.05.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_Hallmark_significant_padj_0.05.csv"
  ),
  
  
  # --------------------------------------------------------
  # Directional tables
  # --------------------------------------------------------
  
  file.path(
    tables_supp_dir,
    "GSEA_GO_BP_BRCA1_tumour.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_GO_BP_TotalCell.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_Reactome_BRCA1_tumour.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_Reactome_TotalCell.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_KEGG_BRCA1_tumour.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_KEGG_TotalCell.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_Hallmark_BRCA1_tumour.csv"
  ),
  
  file.path(
    tables_supp_dir,
    "GSEA_Hallmark_TotalCell.csv"
  ),
  
  
  # --------------------------------------------------------
  # Validation
  # --------------------------------------------------------
  
  file.path(
    tables_validation_dir,
    "Step11_GSEA_summary.csv"
  ),
  
  file.path(
    tables_validation_dir,
    "Step11_GSEA_validation_summary.csv"
  )
)


missing_outputs <- expected_files[
  !file.exists(
    expected_files
  )
]


if (length(missing_outputs) > 0) {
  
  stop(
    "Stage 11 completed with missing expected outputs:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}


############################################################
# 11.37 Step 11 log
############################################################

log_file <- file.path(
  
  logs_dir,
  
  "Step11_GSEA.log"
)


log_lines <- c(
  
  "====================================================",
  
  "TNBC Single-cell RNA-seq Analysis Pipeline",
  
  "STEP 11 — GSEA",
  
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
  
  "INPUT",
  
  "----------------------------------------------------",
  
  paste(
    "Stage 10 DEG file:",
    deg_file
  ),
  
  paste(
    "Genes loaded from Stage 10:",
    n_input
  ),
  
  paste(
    "Genes entering GSEA after ENTREZ mapping:",
    n_ranked
  ),
  
  paste(
    "Positive-ranked genes:",
    n_positive
  ),
  
  paste(
    "Negative-ranked genes:",
    n_negative
  ),
  
  paste(
    "Zero-ranked genes:",
    n_zero
  ),
  
  paste(
    "ENTREZ mapping rate (%):",
    round(
      mapping_rate,
      2
    )
  ),
  
  "",
  
  "RANKING",
  
  "----------------------------------------------------",
  
  "Ranking variable: avg_log2FC",
  
  "Positive avg_log2FC = higher in BRCA1_tumour",
  
  "Negative avg_log2FC = higher in TotalCell",
  
  "",
  
  "GSEA FILTERING",
  
  "----------------------------------------------------",
  
  "All eligible ranked genes were used.",
  
  "No |log2FC| cutoff was applied.",
  
  "No p-value filtering was applied before GSEA.",
  
  "Complete GSEA results were retained.",
  
  "Final significance criterion: p.adjust < 0.05",
  
  "Multiple-testing correction: Benjamini-Hochberg (BH)",
  
  "",
  
  "GO BIOLOGICAL PROCESS",
  
  "----------------------------------------------------",
  
  paste(
    "Pathways tested:",
    n_go_total
  ),
  
  paste(
    "Significant pathways:",
    n_go_sig
  ),
  
  paste(
    "Enriched in BRCA1_tumour:",
    n_go_brca1
  ),
  
  paste(
    "Enriched in TotalCell:",
    n_go_totalcell
  ),
  
  "",
  
  "REACTOME",
  
  "----------------------------------------------------",
  
  paste(
    "Pathways tested:",
    n_reactome_total
  ),
  
  paste(
    "Significant pathways:",
    n_reactome_sig
  ),
  
  paste(
    "Enriched in BRCA1_tumour:",
    n_reactome_brca1
  ),
  
  paste(
    "Enriched in TotalCell:",
    n_reactome_totalcell
  ),
  
  "",
  
  "KEGG",
  
  "----------------------------------------------------",
  
  paste(
    "Pathways tested:",
    n_kegg_total
  ),
  
  paste(
    "Significant pathways:",
    n_kegg_sig
  ),
  
  paste(
    "Enriched in BRCA1_tumour:",
    n_kegg_brca1
  ),
  
  paste(
    "Enriched in TotalCell:",
    n_kegg_totalcell
  ),
  
  "",
  
  "HALLMARK — MSigDB",
  
  "----------------------------------------------------",
  
  "Gene set collection: MSigDB Hallmark",
  
  paste(
    "Gene sets tested:",
    n_hallmark_total
  ),
  
  paste(
    "Significant gene sets:",
    n_hallmark_sig
  ),
  
  paste(
    "Enriched in BRCA1_tumour:",
    n_hallmark_brca1
  ),
  
  paste(
    "Enriched in TotalCell:",
    n_hallmark_totalcell
  ),
  
  "",
  
  "INTERPRETATION",
  
  "----------------------------------------------------",
  
  "NES > 0 = enrichment toward BRCA1_tumour",
  
  "NES < 0 = enrichment toward TotalCell",
  
  "p.adjust < 0.05 = statistically significant enrichment",
  
  "",
  
  "VISUALIZATION",
  
  "----------------------------------------------------",
  
  "Dot plots display up to 10 significant terms per enrichment direction.",
  
  "Terms are selected by adjusted p-value.",
  
  "",
  
  "OUTPUT",
  
  "----------------------------------------------------",
  
  paste(
    "Stage 11 directory:",
    stage_dir
  ),
  
  paste(
    "Main figures:",
    figures_main_dir
  ),
  
  paste(
    "Main tables:",
    tables_main_dir
  ),
  
  paste(
    "Supplementary tables:",
    tables_supp_dir
  ),
  
  paste(
    "Validation tables:",
    tables_validation_dir
  ),
  
  paste(
    "Log file:",
    log_file
  ),
  
  "",
  
  "OUTPUT VALIDATION",
  
  "----------------------------------------------------",
  
  paste(
    "Expected output files:",
    length(expected_files)
  ),
  
  paste(
    "Missing output files:",
    length(missing_outputs)
  ),
  
  "",
  
  "STAGE 11 COMPLETED SUCCESSFULLY.",
  
  "===================================================="
)


writeLines(
  
  log_lines,
  
  con = log_file
)


############################################################
# 11.38 Step 11 summary log
############################################################
#
# This file is an existing Stage 11 repository output and
# is therefore preserved as a canonical log artifact.
############################################################

summary_log_file <- file.path(
  
  logs_dir,
  
  "Step11_GSEA_summary.log"
)


summary_log_lines <- c(
  
  "====================================================",
  
  "TNBC Single-cell RNA-seq Analysis Pipeline",
  
  "STEP 11 — GSEA SUMMARY",
  
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
  
  "BRCA1_tumour vs TotalCell",
  
  "",
  
  "INPUT",
  
  paste(
    "Genes loaded from Stage 10:",
    n_input
  ),
  
  paste(
    "Genes entering GSEA after ENTREZ mapping:",
    n_ranked
  ),
  
  paste(
    "ENTREZ mapping rate (%):",
    round(
      mapping_rate,
      2
    )
  ),
  
  "",
  
  "GSEA RESULTS",
  
  paste(
    "GO BP significant pathways:",
    n_go_sig
  ),
  
  paste(
    "  BRCA1_tumour enriched:",
    n_go_brca1
  ),
  
  paste(
    "  TotalCell enriched:",
    n_go_totalcell
  ),
  
  "",
  
  paste(
    "Reactome significant pathways:",
    n_reactome_sig
  ),
  
  paste(
    "  BRCA1_tumour enriched:",
    n_reactome_brca1
  ),
  
  paste(
    "  TotalCell enriched:",
    n_reactome_totalcell
  ),
  
  "",
  
  paste(
    "KEGG significant pathways:",
    n_kegg_sig
  ),
  
  paste(
    "  BRCA1_tumour enriched:",
    n_kegg_brca1
  ),
  
  paste(
    "  TotalCell enriched:",
    n_kegg_totalcell
  ),
  
  "",
  
  paste(
    "Hallmark significant gene sets:",
    n_hallmark_sig
  ),
  
  paste(
    "  BRCA1_tumour enriched:",
    n_hallmark_brca1
  ),
  
  paste(
    "  TotalCell enriched:",
    n_hallmark_totalcell
  ),
  
  "",
  
  "INTERPRETATION",
  
  "NES > 0 = enrichment toward BRCA1_tumour",
  
  "NES < 0 = enrichment toward TotalCell",
  
  "p.adjust < 0.05 = statistically significant enrichment",
  
  "",
  
  "VALIDATION",
  
  "GSEA ranking validation: PASSED",
  
  paste(
    "Expected output files:",
    length(expected_files)
  ),
  
  paste(
    "Missing output files:",
    length(missing_outputs)
  ),
  
  "",
  
  "STAGE 11 SUMMARY COMPLETED SUCCESSFULLY.",
  
  "===================================================="
)


writeLines(
  
  summary_log_lines,
  
  con = summary_log_file
)


############################################################
# 11.39 Verify log outputs
############################################################

expected_log_files <- c(
  
  log_file,
  
  summary_log_file
)


missing_log_files <- expected_log_files[
  !file.exists(
    expected_log_files
  )
]


if (length(missing_log_files) > 0) {
  
  stop(
    "Stage 11 completed but one or more expected log files are missing:\n",
    paste(
      missing_log_files,
      collapse = "\n"
    )
  )
}


############################################################
# 11.40 Console summary
############################################################

message(
  "===================================================="
)

message(
  "Stage 11 completed successfully."
)

message(
  "Genes loaded from Stage 10: ",
  n_input
)

message(
  "Genes entering GSEA: ",
  n_ranked
)

message(
  "ENTREZ mapping rate: ",
  round(
    mapping_rate,
    2
  ),
  "%"
)

message(
  "GO BP significant pathways: ",
  n_go_sig
)

message(
  "  BRCA1_tumour: ",
  n_go_brca1
)

message(
  "  TotalCell: ",
  n_go_totalcell
)

message(
  "Reactome significant pathways: ",
  n_reactome_sig
)

message(
  "  BRCA1_tumour: ",
  n_reactome_brca1
)

message(
  "  TotalCell: ",
  n_reactome_totalcell
)

message(
  "KEGG significant pathways: ",
  n_kegg_sig
)

message(
  "  BRCA1_tumour: ",
  n_kegg_brca1
)

message(
  "  TotalCell: ",
  n_kegg_totalcell
)

message(
  "Hallmark significant gene sets: ",
  n_hallmark_sig
)

message(
  "  BRCA1_tumour: ",
  n_hallmark_brca1
)

message(
  "  TotalCell: ",
  n_hallmark_totalcell
)

message(
  "Main figures saved to: ",
  figures_main_dir
)

message(
  "Main tables saved to: ",
  tables_main_dir
)

message(
  "Supplementary tables saved to: ",
  tables_supp_dir
)

message(
  "Validation tables saved to: ",
  tables_validation_dir
)

message(
  "Main log saved to: ",
  log_file
)

message(
  "Summary log saved to: ",
  summary_log_file
)

message(
  "===================================================="
)


############################################################
# End of Step 11
############################################################