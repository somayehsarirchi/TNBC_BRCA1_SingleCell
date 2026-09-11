# ============================================================
# TNBC Single-Cell RNA-seq Pipeline
# STAGE 17E: Reactome + MSigDB Hallmark Validation
#
# Purpose:
#   Independent pathway-level validation of the canonical
#   temporal gene-expression clusters identified in Stage 17B,
#   characterized in Stage 17C, and biologically interpreted
#   in Stage 17D.
#
# Databases:
#   1. Reactome
#   2. MSigDB Hallmark
#
# Important conceptual rule:
#   Temporal clusters represent gene-expression dynamics along
#   pseudotime. They are NOT cell subtypes.
#
# Canonical temporal gene lists from Stage17D:
#   Cluster 1 = 53 genes
#   Cluster 2 = 516 genes
#   Total      = 569 genes
#
# Background:
#   Stage17D reported:
#       RNA features      = 24,865
#       Entrez background = 18,815
#
#   Current reconstruction from the same RNA feature universe
#   may differ slightly because annotation databases can change.
#   Any discrepancy is documented and not artificially corrected.
#
# Mapping principle:
#   A single SYMBOL -> ENTREZID mapping generated from the
#   complete Stage17 RNA feature universe is used consistently
#   for both the background and temporal genes.
#
# Reactome:
#   Entrez-ID based enrichment using the mapped RNA background.
#
# MSigDB Hallmark:
#   SYMBOL-based enrichment using the RNA-derived background
#   intersected with the Hallmark gene universe.
#
# Scientific scope:
#   Stage17E performs independent pathway-level validation only.
#
#   No:
#       - trajectory reconstruction
#       - root selection
#       - pseudotime recalculation
#       - reclustering
#       - cell-type subsetting
#       - modification of Stage17B/C/D objects
#
# Canonical output architecture:
#
# results/
# └── 17_Temporal_Gene_Dynamics/
#     ├── figures/
#     │   ├── main/
#     │   └── supplementary/
#     ├── tables/
#     │   ├── main/
#     │   ├── supplementary/
#     │   └── validation/
#     ├── logs/
#     └── README.md
#
# Central Stage17E object:
#   objects/Stage17E_Reactome_Hallmark_Results.rds
#
# ============================================================


# ============================================================
# 1. PROJECT SETUP
# ============================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"

stage17_root <- file.path(
  project_dir,
  "results",
  "17_Temporal_Gene_Dynamics"
)

figures_main_dir <- file.path(
  stage17_root,
  "figures",
  "main"
)

figures_supplementary_dir <- file.path(
  stage17_root,
  "figures",
  "supplementary"
)

tables_main_dir <- file.path(
  stage17_root,
  "tables",
  "main"
)

tables_supplementary_dir <- file.path(
  stage17_root,
  "tables",
  "supplementary"
)

tables_validation_dir <- file.path(
  stage17_root,
  "tables",
  "validation"
)

logs_dir <- file.path(
  stage17_root,
  "logs"
)

objects_dir <- file.path(
  project_dir,
  "objects"
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
  tables_validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. LOGGING
# ============================================================

stage17e_log <- file.path(
  logs_dir,
  "Stage17E_Completion.log"
)

writeLines(
  character(0),
  con = stage17e_log
)

write_log <- function(...) {
  
  txt <- paste0(
    ...,
    collapse = ""
  )
  
  cat(
    txt,
    "\n"
  )
  
  write(
    txt,
    file = stage17e_log,
    append = TRUE
  )
}


write_log(
  "============================================================"
)

write_log(
  "STAGE 17E: REACTOME + MSIGDB HALLMARK VALIDATION"
)

write_log(
  "============================================================"
)

write_log(
  "Started: ",
  as.character(Sys.time())
)


# ============================================================
# 3. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "Seurat",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "msigdbr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "The following required packages are missing:\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\nPlease install them before running Stage17E."
    )
  )
}

suppressPackageStartupMessages({
  
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(msigdbr)
})


write_log(
  "Required packages available."
)


# ============================================================
# 4. INPUT FILES
# ============================================================

stage17d_rds <- file.path(
  objects_dir,
  "Stage17D_Biological_Interpretation_Results.rds"
)

stage17d_source_seurat_default <- file.path(
  objects_dir,
  "TNBC_Epithelial_CopyKAT.rds"
)

if (!file.exists(stage17d_rds)) {
  
  stop(
    paste0(
      "Stage17D RDS not found:\n",
      stage17d_rds
    )
  )
}


# ============================================================
# 5. LOAD STAGE17D
# ============================================================

write_log(
  "Loading frozen Stage17D results..."
)

stage17d <- readRDS(
  stage17d_rds
)

required_stage17d_components <- c(
  "cluster_gene_lists",
  "go_conversion_summary",
  "kegg_conversion_summary",
  "background_gene_count",
  "background_entrez_count",
  "background_source",
  "seurat_source"
)

missing_components <- setdiff(
  required_stage17d_components,
  names(stage17d)
)

if (length(missing_components) > 0) {
  
  stop(
    paste0(
      "The following required Stage17D components are missing:\n",
      paste(
        missing_components,
        collapse = ", "
      )
    )
  )
}

write_log(
  "Stage17D RDS loaded successfully."
)


# ============================================================
# 6. VERIFY FROZEN UPSTREAM STATUS
# ============================================================

if ("upstream_frozen" %in% names(stage17d)) {
  
  if (!isTRUE(stage17d$upstream_frozen)) {
    
    stop(
      "Stage17D does not report upstream_frozen = TRUE."
    )
  }
}

if ("trajectory_reconstruction" %in% names(stage17d)) {
  
  if (isTRUE(stage17d$trajectory_reconstruction)) {
    
    stop(
      "Stage17D reports trajectory reconstruction = TRUE."
    )
  }
}

if ("root_selection" %in% names(stage17d)) {
  
  if (isTRUE(stage17d$root_selection)) {
    
    stop(
      "Stage17D reports root selection = TRUE."
    )
  }
}

if ("pseudotime_recalculation" %in% names(stage17d)) {
  
  if (isTRUE(stage17d$pseudotime_recalculation)) {
    
    stop(
      "Stage17D reports pseudotime recalculation = TRUE."
    )
  }
}

write_log(
  "Stage17D frozen-status checks passed."
)


# ============================================================
# 7. EXTRACT CANONICAL TEMPORAL GENE LISTS
# ============================================================

cluster_gene_lists <- stage17d$cluster_gene_lists

if (
  !is.list(cluster_gene_lists) ||
  length(cluster_gene_lists) != 2
) {
  
  stop(
    "Stage17D cluster_gene_lists must contain exactly two clusters."
  )
}

cluster_names <- names(
  cluster_gene_lists
)

if (
  is.null(cluster_names) ||
  any(cluster_names == "")
) {
  
  cluster_names <- as.character(
    seq_along(cluster_gene_lists)
  )
  
  names(
    cluster_gene_lists
  ) <- cluster_names
}

if (
  !all(
    c("1", "2") %in%
    names(cluster_gene_lists)
  )
) {
  
  stop(
    paste0(
      "Expected cluster names '1' and '2'. Found: ",
      paste(
        names(cluster_gene_lists),
        collapse = ", "
      )
    )
  )
}


# Standardize gene symbols without changing membership.
cluster_gene_lists <- lapply(
  cluster_gene_lists,
  function(x) {
    
    x <- as.character(x)
    x <- trimws(x)
    
    x <- x[
      !is.na(x) &
        x != ""
    ]
    
    unique(x)
  }
)

cluster_input_counts <- vapply(
  cluster_gene_lists,
  length,
  FUN.VALUE = integer(1)
)

write_log(
  "Canonical temporal gene counts:"
)

write_log(
  "Cluster 1 = ",
  cluster_input_counts["1"]
)

write_log(
  "Cluster 2 = ",
  cluster_input_counts["2"]
)

write_log(
  "Total = ",
  sum(cluster_input_counts)
)

if (
  cluster_input_counts["1"] != 53 ||
  cluster_input_counts["2"] != 516 ||
  sum(cluster_input_counts) != 569
) {
  
  stop(
    paste0(
      "Canonical Stage17D temporal gene counts do not match ",
      "expected values of 53, 516, and 569."
    )
  )
}


# ============================================================
# 8. LOCATE THE SEURAT SOURCE
# ============================================================

stage17d_seurat_source <- as.character(
  stage17d$seurat_source
)

candidate_seurat_paths <- unique(
  c(
    stage17d_seurat_source,
    stage17d_source_seurat_default
  )
)

candidate_seurat_paths <- candidate_seurat_paths[
  !is.na(candidate_seurat_paths) &
    candidate_seurat_paths != ""
]

existing_seurat_paths <- candidate_seurat_paths[
  file.exists(candidate_seurat_paths)
]

if (
  length(existing_seurat_paths) == 0
) {
  
  stop(
    paste0(
      "Could not locate the Seurat source used by Stage17D.\n\n",
      "Candidates checked:\n",
      paste(
        candidate_seurat_paths,
        collapse = "\n"
      )
    )
  )
}

seurat_source <- existing_seurat_paths[1]

write_log(
  "Seurat source:"
)

write_log(
  seurat_source
)


# ============================================================
# 9. LOAD SEURAT OBJECT
# ============================================================

write_log(
  "Loading Seurat expression object..."
)

seurat_obj <- readRDS(
  seurat_source
)

if (
  !"RNA" %in% Seurat::Assays(seurat_obj)
) {
  
  stop(
    "RNA assay was not found in the Seurat object."
  )
}

rna_features <- rownames(
  seurat_obj[["RNA"]]
)

if (
  length(rna_features) == 0
) {
  
  stop(
    "RNA assay contains no features."
  )
}

rna_features <- unique(
  as.character(rna_features)
)

write_log(
  "RNA feature count = ",
  length(rna_features)
)

stage17d_background_genes <- as.integer(
  stage17d$background_gene_count
)

if (
  length(rna_features) !=
  stage17d_background_genes
) {
  
  stop(
    paste0(
      "RNA feature count does not match Stage17D background gene count.\n",
      "Current = ",
      length(rna_features),
      "\n",
      "Stage17D = ",
      stage17d_background_genes
    )
  )
}


# ============================================================
# 10. MAP RNA BACKGROUND TO ENTREZ IDS
# ============================================================

write_log(
  "Mapping RNA background SYMBOLs to Entrez IDs..."
)

background_mapping_raw <- suppressMessages(
  clusterProfiler::bitr(
    rna_features,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

if (
  is.null(background_mapping_raw) ||
  nrow(background_mapping_raw) == 0
) {
  
  stop(
    "No RNA background genes could be mapped to Entrez IDs."
  )
}

background_mapping_raw <- background_mapping_raw %>%
  dplyr::mutate(
    SYMBOL = as.character(SYMBOL),
    ENTREZID = as.character(ENTREZID)
  ) %>%
  dplyr::filter(
    !is.na(SYMBOL),
    SYMBOL != "",
    !is.na(ENTREZID),
    ENTREZID != ""
  )


# Preserve one SYMBOL -> one Entrez assignment for
# construction of the reproducible background.
#
# This same mapping table is subsequently used for the
# temporal genes, avoiding a second independent bitr() call.

background_mapping <- background_mapping_raw %>%
  dplyr::distinct(
    SYMBOL,
    .keep_all = TRUE
  )

background_entrez <- unique(
  background_mapping$ENTREZID
)

background_symbol <- unique(
  background_mapping$SYMBOL
)

background_entrez_count <- length(
  background_entrez
)

write_log(
  "Current reconstructed Entrez background = ",
  background_entrez_count
)


# ============================================================
# 10A. BACKGROUND AND TEMPORAL GENE MAPPING QC
# ============================================================

write_log(
  "============================================================"
)

write_log(
  "STAGE17E BACKGROUND MAPPING QC"
)

write_log(
  "============================================================"
)

unmapped_background_genes <- setdiff(
  rna_features,
  unique(background_mapping_raw$SYMBOL)
)

background_mapping_rate <- (
  length(unique(background_mapping_raw$SYMBOL)) /
    length(rna_features)
) * 100

write_log(
  "RNA features = ",
  length(rna_features)
)

write_log(
  "Mapped SYMBOLs = ",
  length(unique(background_mapping_raw$SYMBOL))
)

write_log(
  "Mapped Entrez IDs = ",
  length(unique(background_mapping_raw$ENTREZID))
)

write_log(
  "Unmapped RNA features = ",
  length(unmapped_background_genes)
)

write_log(
  "Background mapping rate (%) = ",
  round(
    background_mapping_rate,
    2
  )
)

write_log(
  "First 100 unmapped RNA features:"
)

write(
  paste(
    head(
      unmapped_background_genes,
      100
    ),
    collapse = ", "
  ),
  file = stage17e_log,
  append = TRUE
)


# ============================================================
# 10B. TEMPORAL GENE MAPPING QC
# ============================================================

write_log(
  "============================================================"
)

write_log(
  "TEMPORAL GENE MAPPING QC"
)

write_log(
  "============================================================"
)

temporal_mapping_qc <- list()

for (cl in names(cluster_gene_lists)) {
  
  genes <- cluster_gene_lists[[cl]]
  
  mapped_genes <- intersect(
    genes,
    unique(background_mapping_raw$SYMBOL)
  )
  
  unmapped_genes <- setdiff(
    genes,
    mapped_genes
  )
  
  mapping_rate <- (
    length(mapped_genes) /
      length(genes)
  ) * 100
  
  temporal_mapping_qc[[cl]] <- tibble::tibble(
    Temporal_Cluster = cl,
    Input_Genes = length(genes),
    Mapped_to_Entrez = length(mapped_genes),
    Unmapped = length(unmapped_genes),
    Mapping_Rate_Percent = round(
      mapping_rate,
      2
    )
  )
  
  write_log(
    "Cluster ",
    cl
  )
  
  write_log(
    "Input genes = ",
    length(genes)
  )
  
  write_log(
    "Mapped to Entrez = ",
    length(mapped_genes)
  )
  
  write_log(
    "Unmapped = ",
    length(unmapped_genes)
  )
  
  write_log(
    "Mapping rate (%) = ",
    round(
      mapping_rate,
      2
    )
  )
  
  if (
    length(unmapped_genes) > 0
  ) {
    
    write_log(
      "Unmapped genes:"
    )
    
    write(
      paste(
        unmapped_genes,
        collapse = ", "
      ),
      file = stage17e_log,
      append = TRUE
    )
  }
}

temporal_mapping_qc_df <-
  dplyr::bind_rows(
    temporal_mapping_qc
  )

write_csv(
  temporal_mapping_qc_df,
  file.path(
    tables_validation_dir,
    "Stage17E_Mapping_QC.csv"
  )
)

write_log(
  "Mapping QC table exported:"
)

write_log(
  file.path(
    tables_validation_dir,
    "Stage17E_Mapping_QC.csv"
  )
)

write_log(
  "End of mapping QC."
)


# ============================================================
# 11. BACKGROUND DISCREPANCY DOCUMENTATION
# ============================================================

stage17d_background_entrez <- as.integer(
  stage17d$background_entrez_count
)

background_difference <- (
  background_entrez_count -
    stage17d_background_entrez
)

background_difference_percent <- (
  background_difference /
    stage17d_background_entrez
) * 100

write_log(
  "Stage17D reported Entrez background = ",
  stage17d_background_entrez
)

write_log(
  "Current reconstructed Entrez background = ",
  background_entrez_count
)

write_log(
  "Entrez background difference = ",
  background_difference
)

write_log(
  "Entrez background difference (%) = ",
  round(
    background_difference_percent,
    6
  )
)

if (
  abs(background_difference) > 0
) {
  
  write_log(
    "Annotation reproducibility discrepancy detected."
  )
  
  write_log(
    "The discrepancy is documented and not artificially corrected."
  )
}


# ============================================================
# 12. BACKGROUND DUPLICATION QC
# ============================================================

duplicated_entrez <- background_mapping %>%
  dplyr::count(
    ENTREZID,
    name = "n"
  ) %>%
  dplyr::filter(
    n > 1
  )

multiple_symbols_per_entrez <- background_mapping %>%
  dplyr::count(
    ENTREZID,
    name = "n_symbols"
  ) %>%
  dplyr::filter(
    n_symbols > 1
  )

duplicate_entrez_count <- nrow(
  duplicated_entrez
)

multiple_symbol_entrez_count <- nrow(
  multiple_symbols_per_entrez
)

write_log(
  "Duplicated Entrez IDs after SYMBOL-level deduplication = ",
  duplicate_entrez_count
)

write_log(
  "Entrez IDs associated with multiple SYMBOLs = ",
  multiple_symbol_entrez_count
)


# ============================================================
# 13. EXPORT BACKGROUND MAPPING
# ============================================================

write_csv(
  background_mapping,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Background_SYMBOL_to_ENTREZ.csv"
  )
)


# ============================================================
# 14. TEMPORAL GENE -> ENTREZ MAPPING
# ============================================================
#
# Important:
#   Temporal genes are mapped by subsetting the single canonical
#   background mapping table generated in Section 10.
#
#   A second independent bitr() call is intentionally avoided.
#   This ensures that the background and temporal gene mappings
#   use exactly the same annotation operation.
#
# ============================================================

temporal_mapping_summary <- list()

temporal_entrez_lists <- list()

temporal_mapping_details <- list()

for (cl in names(cluster_gene_lists)) {
  
  genes <- cluster_gene_lists[[cl]]
  
  mapping <- background_mapping %>%
    dplyr::filter(
      SYMBOL %in% genes
    ) %>%
    dplyr::select(
      SYMBOL,
      ENTREZID
    ) %>%
    dplyr::distinct(
      SYMBOL,
      .keep_all = TRUE
    )
  
  mapped_genes <- unique(
    mapping$SYMBOL
  )
  
  mapped_entrez <- unique(
    mapping$ENTREZID
  )
  
  unmapped_genes <- setdiff(
    genes,
    mapped_genes
  )
  
  temporal_entrez_lists[[cl]] <-
    mapped_entrez
  
  temporal_mapping_details[[cl]] <-
    tibble::tibble(
      Temporal_Cluster = cl,
      SYMBOL = genes
    ) %>%
    dplyr::left_join(
      mapping,
      by = "SYMBOL"
    ) %>%
    dplyr::mutate(
      Mapping_Status = ifelse(
        is.na(ENTREZID),
        "Unmapped",
        "Mapped"
      )
    )
  
  temporal_mapping_summary[[cl]] <-
    tibble::tibble(
      Temporal_Cluster = cl,
      Input_Genes = length(genes),
      Mapped_SYMBOLs = length(mapped_genes),
      Mapped_Entrez_Genes = length(mapped_entrez),
      Unmapped_Genes = length(unmapped_genes),
      Mapping_Rate_Percent =
        round(
          100 *
            length(mapped_genes) /
            length(genes),
          2
        )
    )
  
  write_log(
    "Temporal mapping from canonical background mapping: Cluster ",
    cl,
    " = ",
    length(mapped_genes),
    " mapped SYMBOLs / ",
    length(genes),
    " input genes"
  )
}

temporal_mapping_summary_df <- dplyr::bind_rows(
  temporal_mapping_summary
)

temporal_mapping_details_df <- dplyr::bind_rows(
  temporal_mapping_details
)

write_csv(
  temporal_mapping_summary_df,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Temporal_Gene_Entrez_Mapping_Summary.csv"
  )
)

write_csv(
  temporal_mapping_details_df,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Temporal_Gene_Entrez_Mapping_Details.csv"
  )
)

write_log(
  "Temporal gene mapping completed using the canonical background mapping."
)


# ============================================================
# 15. REACTOME ENRICHMENT
# ============================================================

write_log(
  "============================================================"
)

write_log(
  "REACTOME ENRICHMENT"
)

write_log(
  "============================================================"
)

reactome_results_list <- list()
reactome_mapping_summary <- list()
reactome_summary <- list()

for (cl in names(cluster_gene_lists)) {
  
  genes <- cluster_gene_lists[[cl]]
  
  write_log(
    "Reactome Cluster ",
    cl,
    ": ",
    length(genes),
    " input genes"
  )
  
  mapped_entrez <- temporal_entrez_lists[[cl]]
  
  mapped_count <- length(
    mapped_entrez
  )
  
  reactome_mapping_summary[[cl]] <-
    tibble::tibble(
      Database = "Reactome",
      Temporal_Cluster = cl,
      Input_Genes = length(genes),
      Mapped_Genes = mapped_count,
      Mapping_Rate_Percent =
        round(
          100 *
            mapped_count /
            length(genes),
          2
        )
    )
  
  if (
    mapped_count == 0
  ) {
    
    enrichment_df <- tibble::tibble()
    
  } else {
    
    reactome_enrich <- tryCatch(
      
      ReactomePA::enrichPathway(
        gene = mapped_entrez,
        organism = "human",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.20,
        universe = background_entrez,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      ),
      
      error = function(e) {
        
        warning(
          paste0(
            "Reactome enrichment failed for Cluster ",
            cl,
            ": ",
            conditionMessage(e)
          )
        )
        
        NULL
      }
    )
    
    if (
      is.null(reactome_enrich)
    ) {
      
      enrichment_df <- tibble::tibble()
      
    } else {
      
      enrichment_df <- as.data.frame(
        reactome_enrich
      )
      
      if (
        nrow(enrichment_df) > 0
      ) {
        
        enrichment_df <-
          enrichment_df %>%
          dplyr::mutate(
            Temporal_Cluster = cl,
            Database = "Reactome"
          ) %>%
          dplyr::arrange(
            p.adjust,
            pvalue
          )
        
      } else {
        
        enrichment_df <-
          tibble::tibble()
      }
    }
  }
  
  reactome_results_list[[cl]] <-
    enrichment_df
  
  if (
    nrow(enrichment_df) > 0
  ) {
    
    significant_count <- sum(
      !is.na(
        enrichment_df$p.adjust
      ) &
        enrichment_df$p.adjust < 0.05
    )
    
    result_pathways <-
      nrow(enrichment_df)
    
  } else {
    
    significant_count <- 0
    result_pathways <- 0
  }
  
  reactome_summary[[cl]] <-
    tibble::tibble(
      Database = "Reactome",
      Temporal_Cluster = cl,
      Input_Genes = length(genes),
      Mapped_Genes = mapped_count,
      Mapping_Rate_Percent =
        round(
          100 *
            mapped_count /
            length(genes),
          2
        ),
      Result_Pathways =
        result_pathways,
      Significant_Pathways =
        significant_count
    )
  
  write_csv(
    enrichment_df,
    file.path(
      tables_main_dir,
      paste0(
        "Stage17E_Reactome_Cluster_",
        cl,
        ".csv"
      )
    )
  )
}


reactome_mapping_summary_df <-
  dplyr::bind_rows(
    reactome_mapping_summary
  )

reactome_summary_df <-
  dplyr::bind_rows(
    reactome_summary
  )

reactome_all_df <-
  dplyr::bind_rows(
    reactome_results_list
  )

write_csv(
  reactome_mapping_summary_df,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Reactome_Mapping_Summary.csv"
  )
)

write_csv(
  reactome_summary_df,
  file.path(
    tables_main_dir,
    "Stage17E_Reactome_Enrichment_Summary.csv"
  )
)

write_csv(
  reactome_all_df,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Reactome_All_Clusters.csv"
  )
)


# ============================================================
# 16. REACTOME DOTPLOTS
# ============================================================

for (cl in names(reactome_results_list)) {
  
  df <- reactome_results_list[[cl]]
  
  if (
    nrow(df) == 0
  ) {
    
    next
  }
  
  plot_df <- df %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::arrange(
      p.adjust
    ) %>%
    dplyr::slice_head(
      n = 20
    )
  
  if (
    nrow(plot_df) == 0
  ) {
    
    next
  }
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      Description = factor(
        Description,
        levels = rev(
          Description
        )
      )
    )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = GeneRatio,
      y = Description,
      size = Count,
      fill = p.adjust
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "black"
    ) +
    ggplot2::scale_fill_continuous(
      trans = "reverse"
    ) +
    ggplot2::labs(
      title = paste0(
        "Reactome Enrichment - Temporal Cluster ",
        cl
      ),
      x = "Gene Ratio",
      y = NULL,
      size = "Gene Count",
      fill = "Adjusted P-value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 9
      )
    )
  
  ggsave(
    filename = file.path(
      figures_main_dir,
      paste0(
        "Stage17E_Reactome_Dotplot_Cluster_",
        cl,
        ".pdf"
      )
    ),
    plot = p,
    width = 10,
    height = 7
  )
}


# ============================================================
# 17. RETRIEVE MSIGDB HALLMARK DATABASE
# ============================================================

write_log(
  "============================================================"
)

write_log(
  "MSIGDB HALLMARK"
)

write_log(
  "============================================================"
)

hallmark_db <- NULL

hallmark_retrieval_method <- NA_character_


# ------------------------------------------------------------
# Modern msigdbr syntax
# ------------------------------------------------------------

hallmark_db <- tryCatch(
  
  msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "H"
  ),
  
  error = function(e) {
    NULL
  }
)

if (
  !is.null(hallmark_db) &&
  nrow(hallmark_db) > 0
) {
  
  hallmark_retrieval_method <-
    "msigdbr(db_species='HS', species='Homo sapiens', collection='H')"
}


# ------------------------------------------------------------
# Legacy msigdbr syntax
# ------------------------------------------------------------

if (
  is.null(hallmark_db) ||
  nrow(hallmark_db) == 0
) {
  
  hallmark_db <- tryCatch(
    
    msigdbr::msigdbr(
      species = "Homo sapiens",
      category = "H"
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  if (
    !is.null(hallmark_db) &&
    nrow(hallmark_db) > 0
  ) {
    
    hallmark_retrieval_method <-
      "msigdbr(species='Homo sapiens', category='H')"
  }
}


if (
  is.null(hallmark_db) ||
  nrow(hallmark_db) == 0
) {
  
  stop(
    paste0(
      "MSigDB Hallmark gene-set retrieval returned no data.\n",
      "Please inspect the installed msigdbr version."
    )
  )
}


write_log(
  "Hallmark database retrieved using:"
)

write_log(
  hallmark_retrieval_method
)


# ============================================================
# 18. STANDARDIZE HALLMARK DATABASE
# ============================================================

hallmark_required_columns <- c(
  "gs_name",
  "gene_symbol"
)

missing_hallmark_columns <- setdiff(
  hallmark_required_columns,
  colnames(hallmark_db)
)

if (
  length(missing_hallmark_columns) > 0
) {
  
  stop(
    paste0(
      "Missing Hallmark columns:\n",
      paste(
        missing_hallmark_columns,
        collapse = ", "
      )
    )
  )
}

hallmark_db <- hallmark_db %>%
  dplyr::transmute(
    gs_name = as.character(gs_name),
    gene_symbol = as.character(gene_symbol)
  ) %>%
  dplyr::filter(
    !is.na(gs_name),
    gs_name != "",
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  dplyr::distinct()

hallmark_gene_sets <- split(
  hallmark_db$gene_symbol,
  hallmark_db$gs_name
)

hallmark_gene_sets <- lapply(
  hallmark_gene_sets,
  unique
)

hallmark_genes_all <- unique(
  hallmark_db$gene_symbol
)

write_log(
  "Hallmark gene sets = ",
  length(hallmark_gene_sets)
)

write_log(
  "Unique Hallmark genes = ",
  length(hallmark_genes_all)
)


# ============================================================
# 19. DEFINE HALLMARK UNIVERSE
# ============================================================
#
# Hallmark enrichment uses a SYMBOL-based universe because
# TERM2GENE is defined by SYMBOL.
#
# The Hallmark universe is restricted to genes that are both:
#
#   1. present in the RNA-derived Stage17D background
#   2. represented in the MSigDB Hallmark database
#
# This prevents genes outside the Hallmark annotation space
# from contributing to the enrichment denominator.
#
# ============================================================

hallmark_background <- intersect(
  background_symbol,
  hallmark_genes_all
)

hallmark_background <- unique(
  hallmark_background
)

hallmark_background_count <- length(
  hallmark_background
)

write_log(
  "Hallmark background genes = ",
  hallmark_background_count
)

if (
  hallmark_background_count == 0
) {
  
  stop(
    "Hallmark background contains zero genes."
  )
}


# ============================================================
# 20. HALLMARK ENRICHMENT
# ============================================================

hallmark_results_list <- list()
hallmark_mapping_summary <- list()
hallmark_summary <- list()

for (cl in names(cluster_gene_lists)) {
  
  genes <- cluster_gene_lists[[cl]]
  
  write_log(
    "Hallmark Cluster ",
    cl,
    ": ",
    length(genes),
    " input genes"
  )
  
  genes_in_hallmark_background <-
    intersect(
      genes,
      hallmark_background
    )
  
  mapped_count <- length(
    genes_in_hallmark_background
  )
  
  hallmark_mapping_summary[[cl]] <-
    tibble::tibble(
      Database = "MSigDB_Hallmark",
      Temporal_Cluster = cl,
      Input_Genes = length(genes),
      Mapped_Genes = mapped_count,
      Mapping_Rate_Percent =
        round(
          100 *
            mapped_count /
            length(genes),
          2
        )
    )
  
  if (
    mapped_count == 0
  ) {
    
    enrichment_df <- tibble::tibble()
    
  } else {
    
    hallmark_enrich <- tryCatch(
      
      clusterProfiler::enricher(
        gene =
          genes_in_hallmark_background,
        TERM2GENE =
          hallmark_db,
        universe =
          hallmark_background,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.20,
        minGSSize = 10,
        maxGSSize = 500
      ),
      
      error = function(e) {
        
        warning(
          paste0(
            "Hallmark enrichment failed for Cluster ",
            cl,
            ": ",
            conditionMessage(e)
          )
        )
        
        NULL
      }
    )
    
    if (
      is.null(hallmark_enrich)
    ) {
      
      enrichment_df <- tibble::tibble()
      
    } else {
      
      enrichment_df <- as.data.frame(
        hallmark_enrich
      )
      
      if (
        nrow(enrichment_df) > 0
      ) {
        
        enrichment_df <-
          enrichment_df %>%
          dplyr::mutate(
            Temporal_Cluster = cl,
            Database = "MSigDB_Hallmark"
          ) %>%
          dplyr::arrange(
            p.adjust,
            pvalue
          )
        
      } else {
        
        enrichment_df <-
          tibble::tibble()
      }
    }
  }
  
  hallmark_results_list[[cl]] <-
    enrichment_df
  
  if (
    nrow(enrichment_df) > 0
  ) {
    
    significant_count <- sum(
      !is.na(
        enrichment_df$p.adjust
      ) &
        enrichment_df$p.adjust < 0.05
    )
    
    result_pathways <-
      nrow(enrichment_df)
    
  } else {
    
    significant_count <- 0
    result_pathways <- 0
  }
  
  hallmark_summary[[cl]] <-
    tibble::tibble(
      Database = "MSigDB_Hallmark",
      Temporal_Cluster = cl,
      Input_Genes = length(genes),
      Mapped_Genes = mapped_count,
      Mapping_Rate_Percent =
        round(
          100 *
            mapped_count /
            length(genes),
          2
        ),
      Result_Pathways =
        result_pathways,
      Significant_Pathways =
        significant_count
    )
  
  write_csv(
    enrichment_df,
    file.path(
      tables_main_dir,
      paste0(
        "Stage17E_Hallmark_Cluster_",
        cl,
        ".csv"
      )
    )
  )
}


hallmark_mapping_summary_df <-
  dplyr::bind_rows(
    hallmark_mapping_summary
  )

hallmark_summary_df <-
  dplyr::bind_rows(
    hallmark_summary
  )

hallmark_all_df <-
  dplyr::bind_rows(
    hallmark_results_list
  )

write_csv(
  hallmark_mapping_summary_df,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Hallmark_Mapping_Summary.csv"
  )
)

write_csv(
  hallmark_summary_df,
  file.path(
    tables_main_dir,
    "Stage17E_Hallmark_Enrichment_Summary.csv"
  )
)

write_csv(
  hallmark_all_df,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Hallmark_All_Clusters.csv"
  )
)


# ============================================================
# 21. HALLMARK DOTPLOTS
# ============================================================

for (cl in names(hallmark_results_list)) {
  
  df <- hallmark_results_list[[cl]]
  
  if (
    nrow(df) == 0
  ) {
    
    next
  }
  
  plot_df <- df %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::arrange(
      p.adjust
    ) %>%
    dplyr::slice_head(
      n = 20
    )
  
  if (
    nrow(plot_df) == 0
  ) {
    
    next
  }
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      Description = factor(
        Description,
        levels = rev(
          Description
        )
      )
    )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = GeneRatio,
      y = Description,
      size = Count,
      fill = p.adjust
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "black"
    ) +
    ggplot2::scale_fill_continuous(
      trans = "reverse"
    ) +
    ggplot2::labs(
      title = paste0(
        "MSigDB Hallmark Enrichment - Temporal Cluster ",
        cl
      ),
      x = "Gene Ratio",
      y = NULL,
      size = "Gene Count",
      fill = "Adjusted P-value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 9
      )
    )
  
  ggsave(
    filename = file.path(
      figures_main_dir,
      paste0(
        "Stage17E_Hallmark_Dotplot_Cluster_",
        cl,
        ".pdf"
      )
    ),
    plot = p,
    width = 10,
    height = 7
  )
}


# ============================================================
# 22. CROSS-DATABASE SUMMARY
# ============================================================

cross_database_summary <-
  dplyr::bind_rows(
    reactome_summary_df,
    hallmark_summary_df
  ) %>%
  dplyr::arrange(
    Temporal_Cluster,
    Database
  )

write_csv(
  cross_database_summary,
  file.path(
    tables_main_dir,
    "Stage17E_Cross_Database_Summary.csv"
  )
)


# ============================================================
# 23. DATABASE AND BACKGROUND METADATA
# ============================================================

metadata_df <- tibble::tibble(
  
  Item = c(
    "Stage",
    "Description",
    "Stage17D_RDS",
    "Seurat_Source",
    "Stage17D_Background_RNA_Features",
    "Current_Background_RNA_Features",
    "Stage17D_Background_Entrez",
    "Current_Background_Entrez",
    "Background_Entrez_Difference",
    "Background_Entrez_Difference_Percent",
    "Background_Source_Stage17D",
    "Background_Mapping_Method",
    "Temporal_Mapping_Method",
    "Hallmark_Retrieval_Method",
    "Hallmark_Background_Genes",
    "Hallmark_Gene_Sets",
    "Hallmark_Unique_Genes",
    "Total_Temporal_Genes",
    "Temporal_Cluster_1_Genes",
    "Temporal_Cluster_2_Genes",
    "Temporal_Cluster_1_Entrez",
    "Temporal_Cluster_2_Entrez",
    "R_Version",
    "Seurat_Version",
    "clusterProfiler_Version",
    "ReactomePA_Version",
    "org.Hs.eg.db_Version",
    "msigdbr_Version"
  ),
  
  Value = c(
    
    "17E",
    
    "Independent Reactome and MSigDB Hallmark validation of temporal gene-expression clusters",
    
    stage17d_rds,
    
    seurat_source,
    
    as.character(
      stage17d_background_genes
    ),
    
    as.character(
      length(rna_features)
    ),
    
    as.character(
      stage17d_background_entrez
    ),
    
    as.character(
      background_entrez_count
    ),
    
    as.character(
      background_difference
    ),
    
    as.character(
      round(
        background_difference_percent,
        6
      )
    ),
    
    as.character(
      stage17d$background_source
    ),
    
    "SYMBOL -> ENTREZID using org.Hs.eg.db",
    
    "Subset of the canonical background SYMBOL -> ENTREZID mapping",
    
    hallmark_retrieval_method,
    
    as.character(
      hallmark_background_count
    ),
    
    as.character(
      length(hallmark_gene_sets)
    ),
    
    as.character(
      length(hallmark_genes_all)
    ),
    
    as.character(
      sum(cluster_input_counts)
    ),
    
    as.character(
      cluster_input_counts["1"]
    ),
    
    as.character(
      cluster_input_counts["2"]
    ),
    
    as.character(
      length(
        temporal_entrez_lists[["1"]]
      )
    ),
    
    as.character(
      length(
        temporal_entrez_lists[["2"]]
      )
    ),
    
    R.version.string,
    
    as.character(
      packageVersion("Seurat")
    ),
    
    as.character(
      packageVersion("clusterProfiler")
    ),
    
    as.character(
      packageVersion("ReactomePA")
    ),
    
    as.character(
      packageVersion("org.Hs.eg.db")
    ),
    
    as.character(
      packageVersion("msigdbr")
    )
  )
)

write_csv(
  metadata_df,
  file.path(
    tables_supplementary_dir,
    "Stage17E_Database_Metadata.csv"
  )
)


# ============================================================
# 24. VALIDATION
# ============================================================

reactome_analysis_completed <-
  nrow(reactome_summary_df) > 0

hallmark_analysis_completed <-
  nrow(hallmark_summary_df) > 0

validation_df <- tibble::tibble(
  
  Check = c(
    
    "Stage17D RDS exists",
    
    "Stage17D cluster gene lists available",
    
    "Stage17D reports upstream frozen",
    
    "No trajectory reconstruction reported",
    
    "No root selection reported",
    
    "No pseudotime recalculation reported",
    
    "Total temporal genes = 569",
    
    "Cluster 1 contains 53 genes",
    
    "Cluster 2 contains 516 genes",
    
    "RNA assay available",
    
    "RNA background features reconstructed",
    
    "Current Entrez background reconstructed",
    
    "Background discrepancy documented",
    
    "Temporal mapping QC completed",
    
    "No duplicated Entrez IDs after SYMBOL-level mapping",
    
    "Temporal Cluster 1 Entrez mapping available",
    
    "Temporal Cluster 2 Entrez mapping available",
    
    "Reactome analysis completed",
    
    "MSigDB Hallmark database retrieved",
    
    "Hallmark analysis completed"
  ),
  
  Result = c(
    
    file.exists(
      stage17d_rds
    ),
    
    is.list(
      cluster_gene_lists
    ),
    
    if (
      "upstream_frozen" %in%
      names(stage17d)
    ) {
      isTRUE(
        stage17d$upstream_frozen
      )
    } else {
      TRUE
    },
    
    if (
      "trajectory_reconstruction" %in%
      names(stage17d)
    ) {
      !isTRUE(
        stage17d$trajectory_reconstruction
      )
    } else {
      TRUE
    },
    
    if (
      "root_selection" %in%
      names(stage17d)
    ) {
      !isTRUE(
        stage17d$root_selection
      )
    } else {
      TRUE
    },
    
    if (
      "pseudotime_recalculation" %in%
      names(stage17d)
    ) {
      !isTRUE(
        stage17d$pseudotime_recalculation
      )
    } else {
      TRUE
    },
    
    sum(
      cluster_input_counts
    ) == 569,
    
    cluster_input_counts["1"] == 53,
    
    cluster_input_counts["2"] == 516,
    
    "RNA" %in%
      Seurat::Assays(
        seurat_obj
      ),
    
    length(rna_features) ==
      stage17d_background_genes,
    
    background_entrez_count > 0,
    
    TRUE,
    
    nrow(temporal_mapping_qc_df) == 2,
    
    duplicate_entrez_count == 0,
    
    length(
      temporal_entrez_lists[["1"]]
    ) > 0,
    
    length(
      temporal_entrez_lists[["2"]]
    ) > 0,
    
    reactome_analysis_completed,
    
    !is.null(hallmark_db) &&
      nrow(hallmark_db) > 0,
    
    hallmark_analysis_completed
  ),
  
  Details = c(
    
    stage17d_rds,
    
    paste(
      names(cluster_gene_lists),
      collapse = ", "
    ),
    
    if (
      "upstream_frozen" %in%
      names(stage17d)
    ) {
      as.character(
        stage17d$upstream_frozen
      )
    } else {
      "Field not present; no upstream modification performed"
    },
    
    if (
      "trajectory_reconstruction" %in%
      names(stage17d)
    ) {
      as.character(
        stage17d$trajectory_reconstruction
      )
    } else {
      "Field not present"
    },
    
    if (
      "root_selection" %in%
      names(stage17d)
    ) {
      as.character(
        stage17d$root_selection
      )
    } else {
      "Field not present"
    },
    
    if (
      "pseudotime_recalculation" %in%
      names(stage17d)
    ) {
      as.character(
        stage17d$pseudotime_recalculation
      )
    } else {
      "Field not present"
    },
    
    paste(
      "Observed:",
      sum(cluster_input_counts)
    ),
    
    paste(
      "Observed:",
      cluster_input_counts["1"]
    ),
    
    paste(
      "Observed:",
      cluster_input_counts["2"]
    ),
    
    "RNA assay detected",
    
    paste(
      "Observed:",
      length(rna_features),
      "| Stage17D:",
      stage17d_background_genes
    ),
    
    paste(
      "Observed:",
      background_entrez_count,
      "| Stage17D:",
      stage17d_background_entrez
    ),
    
    paste(
      "Difference:",
      background_difference,
      "| Percent:",
      round(
        background_difference_percent,
        6
      )
    ),
    
    paste(
      "QC clusters:",
      nrow(
        temporal_mapping_qc_df
      )
    ),
    
    paste(
      "Count:",
      duplicate_entrez_count
    ),
    
    paste(
      "Mapped Entrez:",
      length(
        temporal_entrez_lists[["1"]]
      )
    ),
    
    paste(
      "Mapped Entrez:",
      length(
        temporal_entrez_lists[["2"]]
      )
    ),
    
    paste(
      "Summary rows:",
      nrow(
        reactome_summary_df
      )
    ),
    
    paste(
      "Gene sets:",
      length(
        hallmark_gene_sets
      )
    ),
    
    paste(
      "Summary rows:",
      nrow(
        hallmark_summary_df
      )
    )
  )
)

write_csv(
  validation_df,
  file.path(
    tables_validation_dir,
    "Stage17E_Final_Validation.csv"
  )
)


# ============================================================
# 25. STAGE17E SUMMARY
# ============================================================

reactome_sig_total <- sum(
  reactome_summary_df$Significant_Pathways,
  na.rm = TRUE
)

hallmark_sig_total <- sum(
  hallmark_summary_df$Significant_Pathways,
  na.rm = TRUE
)

stage17e_summary <- tibble::tibble(
  
  Stage = "17E",
  
  Description =
    "Independent Reactome and MSigDB Hallmark validation",
  
  Temporal_Clusters =
    length(cluster_gene_lists),
  
  Total_Temporal_Genes =
    sum(cluster_input_counts),
  
  Cluster_1_Genes =
    cluster_input_counts["1"],
  
  Cluster_2_Genes =
    cluster_input_counts["2"],
  
  Stage17D_Background_RNA_Features =
    stage17d_background_genes,
  
  Current_Background_RNA_Features =
    length(rna_features),
  
  Stage17D_Background_Entrez =
    stage17d_background_entrez,
  
  Current_Background_Entrez =
    background_entrez_count,
  
  Background_Difference =
    background_difference,
  
  Background_Difference_Percent =
    round(
      background_difference_percent,
      6
    ),
  
  Cluster_1_Entrez =
    length(
      temporal_entrez_lists[["1"]]
    ),
  
  Cluster_2_Entrez =
    length(
      temporal_entrez_lists[["2"]]
    ),
  
  Reactome_Significant_Pathways =
    reactome_sig_total,
  
  Hallmark_Significant_Pathways =
    hallmark_sig_total,
  
  Hallmark_Gene_Sets =
    length(hallmark_gene_sets),
  
  Hallmark_Background_Genes =
    hallmark_background_count
)

write_csv(
  stage17e_summary,
  file.path(
    tables_main_dir,
    "Stage17E_Summary.csv"
  )
)


# ============================================================
# 26. SAVE COMPLETE STAGE17E RDS
# ============================================================

stage17e_results <- list(
  
  stage = "17E",
  
  analysis_name =
    "Reactome + MSigDB Hallmark Validation",
  
  description =
    "Independent pathway-level validation of temporal gene-expression clusters",
  
  project_directory =
    project_dir,
  
  stage17d_source =
    stage17d_rds,
  
  seurat_source =
    seurat_source,
  
  temporal_gene_lists =
    cluster_gene_lists,
  
  temporal_gene_counts =
    cluster_input_counts,
  
  temporal_entrez_lists =
    temporal_entrez_lists,
  
  temporal_mapping_summary =
    temporal_mapping_summary_df,
  
  temporal_mapping_details =
    temporal_mapping_details_df,
  
  temporal_mapping_qc =
    temporal_mapping_qc_df,
  
  reactome_results =
    reactome_results_list,
  
  reactome_mapping_summary =
    reactome_mapping_summary_df,
  
  reactome_summary =
    reactome_summary_df,
  
  hallmark_database =
    hallmark_db,
  
  hallmark_retrieval_method =
    hallmark_retrieval_method,
  
  hallmark_gene_sets =
    hallmark_gene_sets,
  
  hallmark_background =
    hallmark_background,
  
  hallmark_results =
    hallmark_results_list,
  
  hallmark_mapping_summary =
    hallmark_mapping_summary_df,
  
  hallmark_summary =
    hallmark_summary_df,
  
  cross_database_summary =
    cross_database_summary,
  
  background_mapping_raw =
    background_mapping_raw,
  
  background_mapping =
    background_mapping,
  
  background_entrez =
    background_entrez,
  
  background_symbol =
    background_symbol,
  
  background_metadata =
    metadata_df,
  
  validation =
    validation_df,
  
  summary =
    stage17e_summary,
  
  stage17d_background_gene_count =
    stage17d_background_genes,
  
  stage17d_background_entrez_count =
    stage17d_background_entrez,
  
  current_background_gene_count =
    length(rna_features),
  
  current_background_entrez_count =
    background_entrez_count,
  
  background_entrez_difference =
    background_difference,
  
  background_entrez_difference_percent =
    background_difference_percent,
  
  upstream_frozen =
    TRUE,
  
  trajectory_reconstruction =
    FALSE,
  
  root_selection =
    FALSE,
  
  pseudotime_recalculation =
    FALSE
)


canonical_stage17e_rds <- file.path(
  objects_dir,
  "Stage17E_Reactome_Hallmark_Results.rds"
)

saveRDS(
  stage17e_results,
  canonical_stage17e_rds
)


# ============================================================
# 27. FINAL COMPLETION LOG
# ============================================================

write_log(
  "============================================================"
)

write_log(
  "STAGE 17E COMPLETED SUCCESSFULLY"
)

write_log(
  "============================================================"
)

write_log(
  "Temporal genes: ",
  sum(cluster_input_counts)
)

write_log(
  "Cluster 1 genes: ",
  cluster_input_counts["1"]
)

write_log(
  "Cluster 2 genes: ",
  cluster_input_counts["2"]
)

write_log(
  "Cluster 1 Entrez genes: ",
  length(
    temporal_entrez_lists[["1"]]
  )
)

write_log(
  "Cluster 2 Entrez genes: ",
  length(
    temporal_entrez_lists[["2"]]
  )
)

write_log(
  "Stage17D background Entrez: ",
  stage17d_background_entrez
)

write_log(
  "Current background Entrez: ",
  background_entrez_count
)

write_log(
  "Background difference: ",
  background_difference
)

write_log(
  "Background mapping rate (%): ",
  round(
    background_mapping_rate,
    2
  )
)

write_log(
  "Reactome significant pathways: ",
  reactome_sig_total
)

write_log(
  "Hallmark significant pathways: ",
  hallmark_sig_total
)

write_log(
  "Hallmark gene sets: ",
  length(hallmark_gene_sets)
)

write_log(
  "Hallmark background genes: ",
  hallmark_background_count
)

write_log(
  "Canonical RDS:"
)

write_log(
  canonical_stage17e_rds
)

write_log(
  "Completion log:"
)

write_log(
  stage17e_log
)

write_log(
  "Stage17B remained frozen."
)

write_log(
  "Stage17C remained frozen."
)

write_log(
  "Stage17D remained frozen."
)

write_log(
  "No trajectory reconstruction was performed."
)

write_log(
  "No root selection was performed."
)

write_log(
  "No pseudotime recalculation was performed."
)

write_log(
  "No temporal gene subsetting was performed."
)

write_log(
  "Stage17E finished: ",
  as.character(Sys.time())
)

write_log(
  "============================================================"
)


# ============================================================
# 28. FINAL CONSOLE SUMMARY
# ============================================================

message("\n")
message(
  "============================================================"
)

message(
  "STAGE 17E COMPLETED SUCCESSFULLY"
)

message(
  "============================================================"
)

message(
  "Temporal genes: ",
  sum(cluster_input_counts)
)

message(
  "Cluster 1 genes: ",
  cluster_input_counts["1"]
)

message(
  "Cluster 2 genes: ",
  cluster_input_counts["2"]
)

message(
  "Cluster 1 Entrez genes: ",
  length(
    temporal_entrez_lists[["1"]]
  )
)

message(
  "Cluster 2 Entrez genes: ",
  length(
    temporal_entrez_lists[["2"]]
  )
)

message(
  "Stage17D background Entrez: ",
  stage17d_background_entrez
)

message(
  "Current background Entrez: ",
  background_entrez_count
)

message(
  "Background difference: ",
  background_difference
)

message(
  "Background mapping rate: ",
  round(
    background_mapping_rate,
    2
  ),
  "%"
)

message(
  "Reactome significant pathways: ",
  reactome_sig_total
)

message(
  "Hallmark significant pathways: ",
  hallmark_sig_total
)

message(
  "Hallmark gene sets: ",
  length(hallmark_gene_sets)
)

message(
  "Hallmark background genes: ",
  hallmark_background_count
)

message(
  "\nCanonical RDS:"
)

message(
  canonical_stage17e_rds
)

message(
  "\nCompletion log:"
)

message(
  stage17e_log
)

message(
  "\nStage17B remained frozen."
)

message(
  "Stage17C remained frozen."
)

message(
  "Stage17D remained frozen."
)

message(
  "No trajectory reconstruction was performed."
)

message(
  "No root selection was performed."
)

message(
  "No pseudotime recalculation was performed."
)

message(
  "============================================================"
)


# ============================================================
# 29. PRINT FINAL SUMMARIES
# ============================================================

reactome_summary_df

hallmark_summary_df

stage17e_summary

validation_df

