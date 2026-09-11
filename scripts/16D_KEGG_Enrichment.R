# ==============================================================================
# Stage 16D: KEGG Enrichment of Strong Pseudotime-associated Genes
#
# Purpose:
#   Functional interpretation of strong trajectory-associated genes identified
#   in Stage 16B using KEGG pathway enrichment.
#
# Frozen inputs:
#   Stage 15  -> Frozen Monocle3 trajectory
#   Stage 16A -> Frozen analysis population
#   Stage 16B -> graph_test trajectory-associated genes
#   Stage 16C -> Gene characterization
#
# Important:
#   - Stage 15 remains completely frozen.
#   - Stage 16A/16B/16C remain unchanged.
#   - KEGG.db is NOT required.
#   - KEGG information is obtained through KEGGREST.
#   - clusterProfiler::enricher() is used with custom TERM2GENE.
# ==============================================================================


# ==============================================================================
# 1. Load Required Libraries
# ==============================================================================

suppressPackageStartupMessages({
  
  library(clusterProfiler)
  library(dplyr)
  library(KEGGREST)
  
})


# ==============================================================================
# 2. Define Project Paths
# ==============================================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"

results_root <- file.path(
  project_dir,
  "results"
)

# ------------------------------------------------------------------------------
# Stage 16B input
# ------------------------------------------------------------------------------

stage16b_dir <- file.path(
  results_root,
  "16B_Pseudotime-associated_Gene_Discovery"
)

# ------------------------------------------------------------------------------
# Stage 16D output
# ------------------------------------------------------------------------------

step_id <- "16D_KEGG_Enrichment"

stage16d_dir <- file.path(
  results_root,
  step_id
)

logs_dir <- file.path(
  stage16d_dir,
  "logs"
)

dir.create(
  stage16d_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. Header
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "STAGE 16D — KEGG ENRICHMENT\n"
)

cat(
  "============================================================\n\n"
)


# ==============================================================================
# 4. Load Stage 16B Results
# ==============================================================================

stage16b_rds <- file.path(
  stage16b_dir,
  "Stage16B_GraphTest_Complete_Results.rds"
)


if (!file.exists(stage16b_rds)) {
  
  stop(
    paste0(
      "Stage 16B result file not found:\n",
      stage16b_rds
    )
  )
  
}


stage16b_results <- readRDS(
  stage16b_rds
)


# ==============================================================================
# 5. Retrieve Frozen Stage 16B Inputs
# ==============================================================================

strong_trajectory_genes <-
  stage16b_results$strong_genes


graph_test_df <-
  stage16b_results$graph_test_dataframe


analysis_cells <-
  stage16b_results$analysis_cells


target_partition <-
  stage16b_results$target_partition


expected_root <-
  stage16b_results$root_node


# ==============================================================================
# 6. Validate Stage 16B Input
# ==============================================================================

if (is.null(strong_trajectory_genes)) {
  
  stop(
    "Strong trajectory gene table is missing from Stage 16B results."
  )
  
}


if (nrow(strong_trajectory_genes) == 0) {
  
  stop(
    "Stage 16B contains zero strong trajectory-associated genes."
  )
  
}


cat(
  "Strong trajectory genes:",
  nrow(strong_trajectory_genes),
  "\n"
)


cat(
  "Analysis cells:",
  length(analysis_cells),
  "\n"
)


cat(
  "Frozen target partition:",
  target_partition,
  "\n"
)


cat(
  "Frozen root:",
  expected_root,
  "\n\n"
)


# ==============================================================================
# 7. Identify Gene Symbols
#
# Stage 16B strong_genes contains gene IDs and gene_short_name.
# We use gene symbols for KEGG ID conversion.
# ==============================================================================

cat(
  "Preparing gene symbols for KEGG conversion...\n"
)


if ("gene_short_name" %in% colnames(strong_trajectory_genes)) {
  
  strong_symbols <- unique(
    as.character(
      strong_trajectory_genes$gene_short_name
    )
  )
  
} else if ("SYMBOL" %in% colnames(strong_trajectory_genes)) {
  
  strong_symbols <- unique(
    as.character(
      strong_trajectory_genes$SYMBOL
    )
  )
  
} else {
  
  stop(
    paste0(
      "Could not identify gene symbols in strong trajectory gene table.\n",
      "Available columns:\n",
      paste(
        colnames(strong_trajectory_genes),
        collapse = ", "
      )
    )
  )
  
}


strong_symbols <- strong_symbols[
  !is.na(strong_symbols) &
    strong_symbols != ""
]


cat(
  "Unique strong trajectory gene symbols:",
  length(strong_symbols),
  "\n\n"
)


# ==============================================================================
# 8. Identify Universe Gene Symbols
#
# The universe is the set of genes actually tested by graph_test.
# This is important because enrichment should be conditioned on the genes
# available to the trajectory analysis rather than the whole genome.
# ==============================================================================

if ("gene_short_name" %in% colnames(graph_test_df)) {
  
  universe_symbols <- unique(
    as.character(
      graph_test_df$gene_short_name
    )
  )
  
} else if ("gene_id" %in% colnames(graph_test_df)) {
  
  universe_symbols <- unique(
    as.character(
      graph_test_df$gene_id
    )
  )
  
} else {
  
  stop(
    paste0(
      "Could not identify universe gene identifiers in graph_test_df.\n",
      "Available columns:\n",
      paste(
        colnames(graph_test_df),
        collapse = ", "
      )
    )
  )
  
}


universe_symbols <- universe_symbols[
  !is.na(universe_symbols) &
    universe_symbols != ""
]


cat(
  "Universe gene identifiers:",
  length(universe_symbols),
  "\n\n"
)


# ==============================================================================
# 9. Convert Strong Genes to Entrez IDs
# ==============================================================================

cat(
  "Converting strong trajectory genes to Entrez IDs...\n"
)


strong_conversion <- tryCatch(
  
  {
    
    bitr(
      strong_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
  },
  
  error = function(e) {
    
    stop(
      paste0(
        "Gene ID conversion failed:\n",
        e$message
      )
    )
    
  }
  
)


strong_conversion <- strong_conversion %>%
  
  dplyr::filter(
    !is.na(ENTREZID)
  ) %>%
  
  dplyr::distinct(
    SYMBOL,
    ENTREZID
  )


cat(
  "Strong genes successfully converted:",
  nrow(strong_conversion),
  "\n"
)


# ==============================================================================
# 10. Convert Universe Genes to Entrez IDs
# ==============================================================================

cat(
  "Converting universe genes to Entrez IDs...\n"
)


universe_conversion <- tryCatch(
  
  {
    
    bitr(
      universe_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
  },
  
  error = function(e) {
    
    stop(
      paste0(
        "Universe gene ID conversion failed:\n",
        e$message
      )
    )
    
  }
  
)


universe_conversion <- universe_conversion %>%
  
  dplyr::filter(
    !is.na(ENTREZID)
  ) %>%
  
  dplyr::distinct(
    SYMBOL,
    ENTREZID
  )


# ==============================================================================
# 11. Create Final Entrez Gene Sets
# ==============================================================================

strong_entrez <- unique(
  as.character(
    strong_conversion$ENTREZID
  )
)


universe_entrez <- unique(
  as.character(
    universe_conversion$ENTREZID
  )
)


# Ensure all tested genes are part of the universe
strong_entrez <- intersect(
  strong_entrez,
  universe_entrez
)


cat(
  "\nFinal KEGG gene sets:\n"
)


cat(
  "Strong trajectory Entrez genes:",
  length(strong_entrez),
  "\n"
)


cat(
  "Universe Entrez genes:",
  length(universe_entrez),
  "\n\n"
)


if (length(strong_entrez) < 10) {
  
  stop(
    "Too few strong trajectory genes available for KEGG enrichment."
  )
}


if (length(universe_entrez) < length(strong_entrez)) {
  
  stop(
    "Universe gene set is smaller than the strong trajectory gene set."
  )
}


# ==============================================================================
# 12. Save Gene Conversion Tables
# ==============================================================================

write.csv(
  
  strong_conversion,
  
  file.path(
    stage16d_dir,
    "00_Strong_Trajectory_Gene_Conversion.csv"
  ),
  
  row.names = FALSE
  
)


write.csv(
  
  universe_conversion,
  
  file.path(
    stage16d_dir,
    "00_Universe_Gene_Conversion.csv"
  ),
  
  row.names = FALSE
  
)


# ==============================================================================
# 13. Retrieve KEGG Gene-Pathway Associations
#
# KEGGREST is used instead of KEGG.db.
# ==============================================================================

cat(
  "Retrieving human KEGG pathway-gene associations...\n"
)


kegg_links_raw <- tryCatch(
  
  {
    
    KEGGREST::keggLink(
      "pathway",
      "hsa"
    )
    
  },
  
  error = function(e) {
    
    stop(
      paste0(
        "KEGG pathway-gene retrieval failed:\n",
        e$message
      )
    )
    
  }
  
)


# ==============================================================================
# 14. Convert KEGG Link Object to Data Frame
# ==============================================================================

kegg_links_df <- data.frame(
  
  ENTREZID = sub(
    "^hsa:",
    "",
    names(kegg_links_raw)
  ),
  
  KEGG_Pathway = sub(
    "^path:",
    "",
    as.character(kegg_links_raw)
  ),
  
  stringsAsFactors = FALSE
  
)


kegg_links_df <- kegg_links_df %>%
  
  dplyr::filter(
    !is.na(ENTREZID),
    !is.na(KEGG_Pathway),
    ENTREZID != "",
    KEGG_Pathway != ""
  ) %>%
  
  dplyr::distinct()


cat(
  "KEGG gene-pathway associations:",
  nrow(kegg_links_df),
  "\n\n"
)


# ==============================================================================
# 15. Retrieve KEGG Pathway Names
# ==============================================================================

cat(
  "Retrieving KEGG pathway names...\n"
)


kegg_names_raw <- tryCatch(
  
  {
    
    KEGGREST::keggList(
      "pathway",
      "hsa"
    )
    
  },
  
  error = function(e) {
    
    stop(
      paste0(
        "KEGG pathway-name retrieval failed:\n",
        e$message
      )
    )
    
  }
  
)


kegg_names <- data.frame(
  
  KEGG_Pathway = sub(
    "^path:",
    "",
    names(kegg_names_raw)
  ),
  
  Description = sub(
    " - Homo sapiens \\(human\\)$",
    "",
    as.character(kegg_names_raw)
  ),
  
  stringsAsFactors = FALSE
  
)


# ==============================================================================
# 16. Clean KEGG Pathway Names
# ==============================================================================

kegg_names <- kegg_names %>%
  
  dplyr::filter(
    !is.na(KEGG_Pathway),
    !is.na(Description)
  ) %>%
  
  dplyr::distinct(
    KEGG_Pathway,
    .keep_all = TRUE
  )


# ==============================================================================
# 17. Restrict KEGG TERM2GENE to Valid Human Pathways
# ==============================================================================

kegg_term2gene <- kegg_links_df %>%
  
  dplyr::inner_join(
    kegg_names,
    by = "KEGG_Pathway"
  ) %>%
  
  dplyr::select(
    KEGG_Pathway,
    ENTREZID
  ) %>%
  
  dplyr::distinct()


cat(
  "KEGG TERM2GENE dimensions:",
  nrow(kegg_term2gene),
  "rows x",
  ncol(kegg_term2gene),
  "columns\n"
)


if (nrow(kegg_term2gene) == 0) {
  
  stop(
    "KEGG TERM2GENE table is empty."
  )
  
}


# ==============================================================================
# 18. Restrict TERM2NAME to Pathways Present in TERM2GENE
# ==============================================================================

kegg_names_enrichment <- kegg_names %>%
  
  dplyr::filter(
    KEGG_Pathway %in%
      unique(kegg_term2gene$KEGG_Pathway)
  )


# ==============================================================================
# 19. Verify Strong Genes are Represented in KEGG
# ==============================================================================

strong_entrez_kegg <- intersect(
  
  strong_entrez,
  
  unique(
    kegg_term2gene$ENTREZID
  )
  
)


universe_entrez_kegg <- intersect(
  
  universe_entrez,
  
  unique(
    kegg_term2gene$ENTREZID
  )
  
)


cat(
  "\nKEGG representation:\n"
)


cat(
  "Strong genes represented in KEGG:",
  length(strong_entrez_kegg),
  "\n"
)


cat(
  "Universe genes represented in KEGG:",
  length(universe_entrez_kegg),
  "\n\n"
)


if (length(strong_entrez_kegg) < 10) {
  
  stop(
    "Too few strong trajectory genes are represented in KEGG."
  )
  
}


# ==============================================================================
# 20. Save KEGG Mapping Tables
# ==============================================================================

write.csv(
  
  kegg_links_df,
  
  file.path(
    stage16d_dir,
    "00_KEGG_Gene_Pathway_Associations.csv"
  ),
  
  row.names = FALSE
  
)


write.csv(
  
  kegg_names,
  
  file.path(
    stage16d_dir,
    "00_KEGG_Pathway_Names.csv"
  ),
  
  row.names = FALSE
  
)


write.csv(
  
  kegg_term2gene,
  
  file.path(
    stage16d_dir,
    "00_KEGG_TERM2GENE.csv"
  ),
  
  row.names = FALSE
  
)


# ==============================================================================
# 21. Run KEGG Enrichment
# ==============================================================================

cat(
  "Running KEGG enrichment...\n"
)


kegg_result <- tryCatch(
  
  {
    
    clusterProfiler::enricher(
      
      gene = strong_entrez_kegg,
      
      universe = universe_entrez_kegg,
      
      TERM2GENE = kegg_term2gene,
      
      TERM2NAME = kegg_names_enrichment,
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      qvalueCutoff = 0.05,
      
      minGSSize = 10,
      
      maxGSSize = 500
      
    )
    
  },
  
  error = function(e) {
    
    stop(
      paste0(
        "KEGG enrichment failed:\n",
        e$message
      )
    )
    
  }
  
)


# ==============================================================================
# 22. Convert Enrichment Result
# ==============================================================================

if (
  is.null(kegg_result) ||
  length(kegg_result) == 0
) {
  
  kegg_df <- data.frame()
  
} else {
  
  kegg_df <- as.data.frame(
    kegg_result
  )
  
}


cat(
  "\nNumber of significant KEGG pathways:",
  nrow(kegg_df),
  "\n\n"
)


# ==============================================================================
# 23. Inspect Results
# ==============================================================================

if (nrow(kegg_df) > 0) {
  
  print(
    head(
      kegg_df,
      20
    )
  )
  
} else {
  
  cat(
    "No significant KEGG pathways detected at the specified cutoff.\n"
  )
  
}


# ==============================================================================
# 24. Save Complete KEGG Result
# ==============================================================================

write.csv(
  
  kegg_df,
  
  file.path(
    stage16d_dir,
    "01_KEGG_Enrichment_Strong_Trajectory_Genes.csv"
  ),
  
  row.names = FALSE
  
)


saveRDS(
  
  kegg_result,
  
  file.path(
    stage16d_dir,
    "01_KEGG_Enrichment_Strong_Trajectory_Genes.rds"
  )
  
)


# ==============================================================================
# 25. Significant KEGG Pathways
# ==============================================================================

if (nrow(kegg_df) > 0) {
  
  kegg_significant <- kegg_df %>%
    
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    
    dplyr::arrange(
      p.adjust,
      pvalue
    )
  
} else {
  
  kegg_significant <- data.frame()
  
}


write.csv(
  
  kegg_significant,
  
  file.path(
    stage16d_dir,
    "02_Significant_KEGG_Pathways.csv"
  ),
  
  row.names = FALSE
  
)


# ==============================================================================
# 26. Top KEGG Pathways
# ==============================================================================

if (nrow(kegg_significant) > 0) {
  
  top_kegg_20 <- kegg_significant %>%
    
    dplyr::slice_head(
      n = 20
    )
  
  top_kegg_50 <- kegg_significant %>%
    
    dplyr::slice_head(
      n = 50
    )
  
} else {
  
  top_kegg_20 <- data.frame()
  
  top_kegg_50 <- data.frame()
  
}


write.csv(
  
  top_kegg_20,
  
  file.path(
    stage16d_dir,
    "03_Top20_KEGG_Pathways.csv"
  ),
  
  row.names = FALSE
  
)


write.csv(
  
  top_kegg_50,
  
  file.path(
    stage16d_dir,
    "04_Top50_KEGG_Pathways.csv"
  ),
  
  row.names = FALSE
  
)


# ==============================================================================
# 27. KEGG Summary
# ==============================================================================

kegg_summary <- data.frame(
  
  Parameter = c(
    
    "Strong_trajectory_genes_Stage16B",
    
    "Strong_Entrez_genes",
    
    "Universe_Entrez_genes",
    
    "Strong_genes_represented_in_KEGG",
    
    "Universe_genes_represented_in_KEGG",
    
    "KEGG_gene_pathway_associations",
    
    "KEGG_pathways_available",
    
    "Significant_KEGG_pathways"
    
  ),
  
  Value = c(
    
    nrow(strong_trajectory_genes),
    
    length(strong_entrez),
    
    length(universe_entrez),
    
    length(strong_entrez_kegg),
    
    length(universe_entrez_kegg),
    
    nrow(kegg_term2gene),
    
    length(
      unique(
        kegg_term2gene$KEGG_Pathway
      )
    ),
    
    nrow(kegg_significant)
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  
  kegg_summary,
  
  file.path(
    stage16d_dir,
    "05_KEGG_Enrichment_Summary.csv"
  ),
  
  row.names = FALSE
  
)


# ==============================================================================
# 28. Save Complete Stage 16D Results
# ==============================================================================

stage16d_results <- list(
  
  strong_symbols = strong_symbols,
  
  universe_symbols = universe_symbols,
  
  strong_conversion = strong_conversion,
  
  universe_conversion = universe_conversion,
  
  strong_entrez = strong_entrez,
  
  universe_entrez = universe_entrez,
  
  strong_entrez_kegg = strong_entrez_kegg,
  
  universe_entrez_kegg = universe_entrez_kegg,
  
  kegg_links_df = kegg_links_df,
  
  kegg_names = kegg_names,
  
  kegg_term2gene = kegg_term2gene,
  
  kegg_result = kegg_result,
  
  kegg_dataframe = kegg_df,
  
  kegg_significant = kegg_significant,
  
  top_kegg_20 = top_kegg_20,
  
  top_kegg_50 = top_kegg_50,
  
  summary = kegg_summary,
  
  target_partition = target_partition,
  
  root_node = expected_root,
  
  analysis_cells = analysis_cells
  
)


saveRDS(
  
  stage16d_results,
  
  file.path(
    stage16d_dir,
    "Stage16D_KEGG_Complete_Results.rds"
  )
  
)


# ==============================================================================
# 29. Completion Log
# ==============================================================================

stage16d_log <- file.path(
  
  logs_dir,
  
  "Stage16D_Completion.log"
  
)


log_lines <- c(
  
  "============================================================",
  
  "TNBC Single-Cell Analysis Pipeline",
  
  "Stage 16D: KEGG Enrichment of Strong Trajectory Genes",
  
  "============================================================",
  
  "",
  
  paste(
    "Completion time:",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  paste(
    "Frozen target partition:",
    target_partition
  ),
  
  paste(
    "Frozen root:",
    expected_root
  ),
  
  paste(
    "Analysis cells:",
    length(analysis_cells)
  ),
  
  paste(
    "Strong Stage 16B genes:",
    nrow(strong_trajectory_genes)
  ),
  
  paste(
    "Strong Entrez genes:",
    length(strong_entrez)
  ),
  
  paste(
    "Strong genes represented in KEGG:",
    length(strong_entrez_kegg)
  ),
  
  paste(
    "Universe Entrez genes:",
    length(universe_entrez)
  ),
  
  paste(
    "KEGG associations:",
    nrow(kegg_term2gene)
  ),
  
  paste(
    "Significant KEGG pathways:",
    nrow(kegg_significant)
  ),
  
  "",
  
  "Stage 15: FROZEN",
  
  "Stage 16A: FROZEN",
  
  "Stage 16B: COMPLETED",
  
  "Stage 16C: COMPLETED",
  
  "Stage 16D: COMPLETED",
  
  "",
  
  "STATUS: STAGE 16D COMPLETED SUCCESSFULLY",
  
  "============================================================"
  
)


writeLines(
  
  log_lines,
  
  con = stage16d_log
  
)


# ==============================================================================
# 30. Final Message
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "STAGE 16D COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n"
)

cat(
  "Strong Stage 16B genes:",
  nrow(strong_trajectory_genes),
  "\n"
)

cat(
  "Strong Entrez genes:",
  length(strong_entrez),
  "\n"
)

cat(
  "Strong genes represented in KEGG:",
  length(strong_entrez_kegg),
  "\n"
)

cat(
  "Universe Entrez genes:",
  length(universe_entrez),
  "\n"
)

cat(
  "KEGG associations:",
  nrow(kegg_term2gene),
  "\n"
)

cat(
  "Significant KEGG pathways:",
  nrow(kegg_significant),
  "\n"
)

cat(
  "\nResults directory:\n",
  stage16d_dir,
  "\n"
)

cat(
  "============================================================\n"
)