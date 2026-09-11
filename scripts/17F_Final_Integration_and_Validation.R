# ============================================================
# STAGE 17F
# FINAL INTEGRATION, VALIDATION, EXPORT AND ARCHIVING
# OF STAGE 17 TEMPORAL GENE DYNAMICS
#
# This is the FINAL Stage17 step.
#
# Integrates:
#   Stage17C = temporal gene-expression characterization
#   Stage17D = GO Biological Process + KEGG
#   Stage17E = Reactome + MSigDB Hallmark
#
# Stage17F performs:
#   1. Final integration
#   2. Biological theme annotation
#   3. Cross-database comparison
#   4. Final biological interpretation
#   5. Final validation
#   6. Canonical export
#   7. Final RDS archiving
#
# IMPORTANT:
#   - No new enrichment analysis is performed.
#   - No new clustering is performed.
#   - No pseudotime recalculation is performed.
#   - No trajectory reconstruction is performed.
#   - No root selection is performed.
#   - No genes are re-selected according to Pattern.
#   - Stage17C, Stage17D and Stage17E remain frozen.
#
# Temporal clusters represent gene-expression dynamics
# along inferred pseudotime and NOT independent cell subtypes.
#
# Biological themes are integration-level annotations only.
# They do not modify the underlying enrichment statistics.
#
# Final canonical object:
#   objects/Stage17F_Temporal_Gene_Dynamics_Final.rds
# ============================================================


# ============================================================
# 1. PROJECT PATHS
# ============================================================

project_dir <- "C:/Users/asus/Desktop/TNBC_BRCA1_SingleCell"

stage17_results_dir <- file.path(
  project_dir,
  "results",
  "17_Temporal_Gene_Dynamics"
)

figures_main_dir <- file.path(
  stage17_results_dir,
  "figures",
  "main"
)

figures_supp_dir <- file.path(
  stage17_results_dir,
  "figures",
  "supplementary"
)

tables_main_dir <- file.path(
  stage17_results_dir,
  "tables",
  "main"
)

tables_supp_dir <- file.path(
  stage17_results_dir,
  "tables",
  "supplementary"
)

tables_validation_dir <- file.path(
  stage17_results_dir,
  "tables",
  "validation"
)

logs_dir <- file.path(
  stage17_results_dir,
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
  figures_supp_dir,
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

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. FROZEN INPUT RDS FILES
# ============================================================

stage17c_rds <- file.path(
  objects_dir,
  "Stage17C_Temporal_Characterization_Results.rds"
)

stage17d_rds <- file.path(
  objects_dir,
  "Stage17D_Biological_Interpretation_Results.rds"
)

stage17e_rds <- file.path(
  objects_dir,
  "Stage17E_Reactome_Hallmark_Results.rds"
)


# ============================================================
# 3. FINAL OUTPUT RDS
# ============================================================

final_rds <- file.path(
  objects_dir,
  "Stage17F_Temporal_Gene_Dynamics_Final.rds"
)


# ============================================================
# 4. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "dplyr",
  "tibble",
  "ggplot2",
  "stringr"
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
    paste(
      "The following required packages are missing:",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 5. START LOG
# ============================================================

start_time <- Sys.time()

log_messages <- character(0)

log_file <- file.path(
  logs_dir,
  "Stage17F_Completion.log"
)

add_log <- function(message) {
  
  timestamp <- format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
  
  log_messages <<- c(
    log_messages,
    paste0(
      "[",
      timestamp,
      "] ",
      message
    )
  )
  
  message
}


add_log(
  "Stage17F final integration started."
)

add_log(
  "Stage17F is the final Stage17 step: integration, validation, export and archiving."
)


# ============================================================
# 6. CHECK FROZEN INPUT FILES
# ============================================================

required_input_files <- c(
  stage17c_rds,
  stage17d_rds,
  stage17e_rds
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  
  stop(
    paste(
      "Required frozen Stage17 input RDS file(s) not found:",
      paste(
        missing_input_files,
        collapse = "\n"
      )
    )
  )
}

add_log(
  "Frozen Stage17C, Stage17D and Stage17E RDS files were found."
)


# ============================================================
# 7. LOAD FROZEN RDS OBJECTS
# ============================================================

stage17c <- readRDS(
  stage17c_rds
)

stage17d <- readRDS(
  stage17d_rds
)

stage17e <- readRDS(
  stage17e_rds
)

add_log(
  "Stage17C RDS loaded."
)

add_log(
  "Stage17D RDS loaded."
)

add_log(
  "Stage17E RDS loaded."
)


# ============================================================
# 8. VALIDATE TOP-LEVEL OBJECT TYPES
# ============================================================

if (!is.list(stage17c)) {
  
  stop(
    "Stage17C RDS does not contain a list object."
  )
}

if (!is.list(stage17d)) {
  
  stop(
    "Stage17D RDS does not contain a list object."
  )
}

if (!is.list(stage17e)) {
  
  stop(
    "Stage17E RDS does not contain a list object."
  )
}

add_log(
  "Stage17C, Stage17D and Stage17E top-level structures validated."
)


# ============================================================
# 9. VALIDATE STAGE IDENTIFIERS
# ============================================================

if (
  "stage" %in% names(stage17c) &&
  !identical(
    as.character(stage17c$stage),
    "17C"
  )
) {
  
  stop(
    paste(
      "Unexpected Stage17C stage identifier:",
      stage17c$stage
    )
  )
}


if (
  "stage" %in% names(stage17d) &&
  !identical(
    as.character(stage17d$stage),
    "17D"
  )
) {
  
  stop(
    paste(
      "Unexpected Stage17D stage identifier:",
      stage17d$stage
    )
  )
}


if (
  "stage" %in% names(stage17e) &&
  !identical(
    as.character(stage17e$stage),
    "17E"
  )
) {
  
  stop(
    paste(
      "Unexpected Stage17E stage identifier:",
      stage17e$stage
    )
  )
}

add_log(
  "Stage identifiers validated."
)


# ============================================================
# 10. VALIDATE REQUIRED STAGE17C STRUCTURE
#
# Actual Stage17C schema contains:
#   cluster_centroids
#   cluster_assignments
#   cluster_summary
#   validation_results
#   validation_checks
#
# Stage17C does NOT use the obsolete temporal_behavior field.
# ============================================================

required_stage17c_fields <- c(
  "cluster_centroids",
  "cluster_assignments",
  "cluster_summary"
)

missing_stage17c_fields <- setdiff(
  required_stage17c_fields,
  names(stage17c)
)

if (length(missing_stage17c_fields) > 0) {
  
  stop(
    paste(
      "Missing required Stage17C field(s):",
      paste(
        missing_stage17c_fields,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 11. VALIDATE REQUIRED STAGE17D STRUCTURE
#
# Actual Stage17D schema contains:
#   cluster_gene_lists
#   cluster_pattern_summary
#   biological_program_description
#   go_results
#   kegg_results
#   go_summary
#   kegg_summary
# ============================================================

required_stage17d_fields <- c(
  "cluster_gene_lists",
  "cluster_pattern_summary",
  "biological_program_description",
  "go_results",
  "kegg_results"
)

missing_stage17d_fields <- setdiff(
  required_stage17d_fields,
  names(stage17d)
)

if (length(missing_stage17d_fields) > 0) {
  
  stop(
    paste(
      "Missing required Stage17D field(s):",
      paste(
        missing_stage17d_fields,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 12. VALIDATE REQUIRED STAGE17E STRUCTURE
#
# Actual Stage17E analytical result objects:
#   reactome_results
#   hallmark_results
# ============================================================

required_stage17e_fields <- c(
  "reactome_results",
  "hallmark_results"
)

missing_stage17e_fields <- setdiff(
  required_stage17e_fields,
  names(stage17e)
)

if (length(missing_stage17e_fields) > 0) {
  
  stop(
    paste(
      "Missing required Stage17E field(s):",
      paste(
        missing_stage17e_fields,
        collapse = ", "
      )
    )
  )
}

add_log(
  "Required Stage17C, Stage17D and Stage17E fields validated."
)


# ============================================================
# 13. EXTRACT CANONICAL TEMPORAL GENE LISTS
# ============================================================

cluster_gene_lists <- stage17d$cluster_gene_lists

if (!is.list(cluster_gene_lists)) {
  
  stop(
    "Stage17D cluster_gene_lists is not a list."
  )
}


cluster_ids <- names(
  cluster_gene_lists
)

if (is.null(cluster_ids)) {
  
  cluster_ids <- as.character(
    seq_along(cluster_gene_lists)
  )
  
  names(
    cluster_gene_lists
  ) <- cluster_ids
}

cluster_ids <- as.character(
  cluster_ids
)

cluster_ids <- cluster_ids[
  nzchar(cluster_ids)
]

if (length(cluster_ids) == 0) {
  
  stop(
    "No temporal clusters were detected."
  )
}


cluster_gene_lists <- lapply(
  cluster_gene_lists,
  function(x) {
    
    x <- as.character(x)
    
    x <- x[
      !is.na(x) &
        nzchar(x)
    ]
    
    unique(x)
  }
)

names(
  cluster_gene_lists
) <- cluster_ids


temporal_gene_counts <- vapply(
  cluster_gene_lists,
  length,
  integer(1)
)


total_temporal_genes <- sum(
  temporal_gene_counts
)


add_log(
  paste(
    "Temporal clusters:",
    paste(
      cluster_ids,
      collapse = ", "
    )
  )
)

add_log(
  paste(
    "Total temporal genes:",
    total_temporal_genes
  )
)


# ============================================================
# 14. VALIDATE EXPECTED TEMPORAL GENE COUNTS
# ============================================================

if (
  !identical(
    as.integer(
      temporal_gene_counts[["1"]]
    ),
    53L
  )
) {
  
  stop(
    paste(
      "Unexpected Cluster 1 gene count:",
      temporal_gene_counts[["1"]],
      "| Expected: 53"
    )
  )
}


if (
  !identical(
    as.integer(
      temporal_gene_counts[["2"]]
    ),
    516L
  )
) {
  
  stop(
    paste(
      "Unexpected Cluster 2 gene count:",
      temporal_gene_counts[["2"]],
      "| Expected: 516"
    )
  )
}


if (
  !identical(
    as.integer(total_temporal_genes),
    569L
  )
) {
  
  stop(
    paste(
      "Unexpected total temporal gene count:",
      total_temporal_genes,
      "| Expected: 569"
    )
  )
}

add_log(
  "Canonical temporal gene counts validated: 53 + 516 = 569."
)


# ============================================================
# 15. CHECK FOR DUPLICATED GENES ACROSS CLUSTERS
# ============================================================

final_temporal_gene_table <- dplyr::bind_rows(
  lapply(
    cluster_ids,
    function(cl) {
      
      tibble::tibble(
        Temporal_Cluster = cl,
        Gene = cluster_gene_lists[[cl]]
      )
    }
  )
)


duplicate_gene_table <- final_temporal_gene_table %>%
  dplyr::count(
    Gene,
    name = "Cluster_Count"
  ) %>%
  dplyr::filter(
    Cluster_Count > 1
  )


if (
  nrow(duplicate_gene_table) > 0
) {
  
  stop(
    paste(
      "Genes are duplicated across temporal clusters:",
      nrow(duplicate_gene_table)
    )
  )
}


add_log(
  "No duplicated genes were detected across temporal clusters."
)


# ============================================================
# 16. EXTRACT TEMPORAL CHARACTERIZATION
#
# Canonical source:
#   Stage17D$cluster_pattern_summary
#
# This preserves the actual Stage17C-derived temporal
# characterization without relying on obsolete field names.
# ============================================================

cluster_pattern_summary <- stage17d$cluster_pattern_summary

if (!is.data.frame(cluster_pattern_summary)) {
  
  stop(
    "Stage17D cluster_pattern_summary is not a data.frame/tibble."
  )
}


cluster_pattern_summary <- cluster_pattern_summary %>%
  dplyr::mutate(
    Temporal_Cluster =
      as.character(
        Temporal_Cluster
      ),
    Pattern =
      as.character(
        Pattern
      )
  ) %>%
  dplyr::arrange(
    match(
      Temporal_Cluster,
      cluster_ids
    )
  )


if (
  !all(
    cluster_ids %in%
    cluster_pattern_summary$Temporal_Cluster
  )
) {
  
  stop(
    "Not all temporal clusters have a Stage17D temporal pattern summary."
  )
}


# ============================================================
# 17. VALIDATE EXPECTED TEMPORAL PATTERNS
# ============================================================

pattern_cluster1 <- cluster_pattern_summary %>%
  dplyr::filter(
    Temporal_Cluster == "1"
  ) %>%
  dplyr::pull(
    Pattern
  )

pattern_cluster2 <- cluster_pattern_summary %>%
  dplyr::filter(
    Temporal_Cluster == "2"
  ) %>%
  dplyr::pull(
    Pattern
  )


if (
  length(pattern_cluster1) != 1 ||
  pattern_cluster1 !=
  "Early_Peaking_Decreasing"
) {
  
  stop(
    paste(
      "Unexpected Cluster 1 temporal pattern:",
      paste(pattern_cluster1, collapse = "; "),
      "| Expected: Early_Peaking_Decreasing"
    )
  )
}


if (
  length(pattern_cluster2) != 1 ||
  pattern_cluster2 !=
  "Late_Increasing"
) {
  
  stop(
    paste(
      "Unexpected Cluster 2 temporal pattern:",
      paste(pattern_cluster2, collapse = "; "),
      "| Expected: Late_Increasing"
    )
  )
}


add_log(
  "Cluster 1 pattern validated as Early_Peaking_Decreasing."
)

add_log(
  "Cluster 2 pattern validated as Late_Increasing."
)


# ============================================================
# 18. CREATE HUMAN-READABLE TEMPORAL DIRECTION
# ============================================================

temporal_summary <- cluster_pattern_summary %>%
  dplyr::mutate(
    Temporal_Direction =
      dplyr::case_when(
        
        Pattern ==
          "Early_Peaking_Decreasing" ~
          "Early-peaking program with subsequent attenuation",
        
        Pattern ==
          "Late_Increasing" ~
          "Late-emerging program with increasing expression",
        
        Pattern ==
          "Increasing" ~
          "Increasing along pseudotime",
        
        Pattern ==
          "Decreasing" ~
          "Decreasing along pseudotime",
        
        TRUE ~
          Pattern
      )
  )


# ============================================================
# 19. EXTRACT STAGE17D BIOLOGICAL PROGRAM DESCRIPTIONS
# ============================================================

biological_program_description <-
  stage17d$biological_program_description


if (
  !is.data.frame(
    biological_program_description
  )
) {
  
  stop(
    "Stage17D biological_program_description is not a data.frame/tibble."
  )
}


biological_program_description <-
  biological_program_description %>%
  dplyr::mutate(
    Temporal_Cluster =
      as.character(
        Temporal_Cluster
      )
  )


# ============================================================
# 20. HELPER FUNCTION:
#     EXTRACT ENRICHMENT RESULTS SAFELY
# ============================================================

extract_enrichment_table <- function(
    enrichment_object,
    database_name,
    cluster_id) {
  
  empty_table <- tibble::tibble(
    Database = character(),
    Temporal_Cluster = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = numeric()
  )
  
  
  if (
    is.null(
      enrichment_object
    )
  ) {
    
    return(
      empty_table
    )
  }
  
  
  result_table <- tryCatch(
    {
      
      if (
        is.data.frame(
          enrichment_object
        )
      ) {
        
        enrichment_object
        
      } else if (
        !is.null(
          enrichment_object@result
        )
      ) {
        
        enrichment_object@result
        
      } else {
        
        NULL
      }
      
    },
    error = function(e) {
      
      NULL
    }
  )
  
  
  if (
    is.null(result_table) ||
    nrow(result_table) == 0
  ) {
    
    return(
      empty_table
    )
  }
  
  
  result_table <- as.data.frame(
    result_table,
    stringsAsFactors = FALSE
  )
  
  
  required_columns <- c(
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "qvalue",
    "geneID",
    "Count"
  )
  
  
  for (
    column_name in required_columns
  ) {
    
    if (
      !column_name %in%
      colnames(result_table)
    ) {
      
      result_table[[column_name]] <-
        NA
    }
  }
  
  
  tibble::tibble(
    
    Database =
      database_name,
    
    Temporal_Cluster =
      as.character(
        cluster_id
      ),
    
    ID =
      as.character(
        result_table$ID
      ),
    
    Description =
      as.character(
        result_table$Description
      ),
    
    GeneRatio =
      as.character(
        result_table$GeneRatio
      ),
    
    BgRatio =
      as.character(
        result_table$BgRatio
      ),
    
    pvalue =
      suppressWarnings(
        as.numeric(
          result_table$pvalue
        )
      ),
    
    p.adjust =
      suppressWarnings(
        as.numeric(
          result_table$p.adjust
        )
      ),
    
    qvalue =
      suppressWarnings(
        as.numeric(
          result_table$qvalue
        )
      ),
    
    geneID =
      as.character(
        result_table$geneID
      ),
    
    Count =
      suppressWarnings(
        as.numeric(
          result_table$Count
        )
      )
  )
}


# ============================================================
# 21. EXTRACT GO RESULTS
# ============================================================

go_results_list <- stage17d$go_results

go_all <- tibble::tibble()


for (
  cl in cluster_ids
) {
  
  go_object <- NULL
  
  
  if (
    !is.null(
      go_results_list[[cl]]
    )
  ) {
    
    go_object <-
      go_results_list[[cl]]
    
  } else {
    
    index <- match(
      cl,
      names(
        go_results_list
      )
    )
    
    if (
      !is.na(index)
    ) {
      
      go_object <-
        go_results_list[[index]]
    }
  }
  
  
  current_go <-
    extract_enrichment_table(
      enrichment_object =
        go_object,
      database_name =
        "GO_BP",
      cluster_id =
        cl
    )
  
  
  if (
    nrow(current_go) > 0
  ) {
    
    go_all <-
      dplyr::bind_rows(
        go_all,
        current_go
      )
  }
}


if (
  nrow(go_all) > 0
) {
  
  go_all <-
    go_all %>%
    dplyr::arrange(
      Temporal_Cluster,
      p.adjust
    )
}


add_log(
  paste(
    "GO BP pathways preserved:",
    nrow(go_all)
  )
)


# ============================================================
# 22. EXTRACT KEGG RESULTS
# ============================================================

kegg_results_list <- stage17d$kegg_results

kegg_all <- tibble::tibble()


for (
  cl in cluster_ids
) {
  
  kegg_object <- NULL
  
  
  if (
    !is.null(
      kegg_results_list[[cl]]
    )
  ) {
    
    kegg_object <-
      kegg_results_list[[cl]]
    
  } else {
    
    index <- match(
      cl,
      names(
        kegg_results_list
      )
    )
    
    if (
      !is.na(index)
    ) {
      
      kegg_object <-
        kegg_results_list[[index]]
    }
  }
  
  
  current_kegg <-
    extract_enrichment_table(
      enrichment_object =
        kegg_object,
      database_name =
        "KEGG",
      cluster_id =
        cl
    )
  
  
  if (
    nrow(current_kegg) > 0
  ) {
    
    kegg_all <-
      dplyr::bind_rows(
        kegg_all,
        current_kegg
      )
  }
}


if (
  nrow(kegg_all) > 0
) {
  
  kegg_all <-
    kegg_all %>%
    dplyr::arrange(
      Temporal_Cluster,
      p.adjust
    )
}


add_log(
  paste(
    "KEGG pathways preserved:",
    nrow(kegg_all)
  )
)


# ============================================================
# 23. EXTRACT REACTOME RESULTS
# ============================================================

reactome_results_list <-
  stage17e$reactome_results

reactome_all <-
  tibble::tibble()


for (
  cl in cluster_ids
) {
  
  reactome_object <- NULL
  
  
  if (
    !is.null(
      reactome_results_list[[cl]]
    )
  ) {
    
    reactome_object <-
      reactome_results_list[[cl]]
    
  } else {
    
    index <- match(
      cl,
      names(
        reactome_results_list
      )
    )
    
    if (
      !is.na(index)
    ) {
      
      reactome_object <-
        reactome_results_list[[index]]
    }
  }
  
  
  current_reactome <-
    extract_enrichment_table(
      enrichment_object =
        reactome_object,
      database_name =
        "Reactome",
      cluster_id =
        cl
    )
  
  
  if (
    nrow(current_reactome) > 0
  ) {
    
    reactome_all <-
      dplyr::bind_rows(
        reactome_all,
        current_reactome
      )
  }
}


if (
  nrow(reactome_all) > 0
) {
  
  reactome_all <-
    reactome_all %>%
    dplyr::arrange(
      Temporal_Cluster,
      p.adjust
    )
}


add_log(
  paste(
    "Reactome pathways preserved:",
    nrow(reactome_all)
  )
)


# ============================================================
# 24. EXTRACT HALLMARK RESULTS
# ============================================================

hallmark_results_list <-
  stage17e$hallmark_results

hallmark_all <-
  tibble::tibble()


for (
  cl in cluster_ids
) {
  
  hallmark_object <- NULL
  
  
  if (
    !is.null(
      hallmark_results_list[[cl]]
    )
  ) {
    
    hallmark_object <-
      hallmark_results_list[[cl]]
    
  } else {
    
    index <- match(
      cl,
      names(
        hallmark_results_list
      )
    )
    
    if (
      !is.na(index)
    ) {
      
      hallmark_object <-
        hallmark_results_list[[index]]
    }
  }
  
  
  current_hallmark <-
    extract_enrichment_table(
      enrichment_object =
        hallmark_object,
      database_name =
        "Hallmark",
      cluster_id =
        cl
    )
  
  
  if (
    nrow(current_hallmark) > 0
  ) {
    
    hallmark_all <-
      dplyr::bind_rows(
        hallmark_all,
        current_hallmark
      )
  }
}


if (
  nrow(hallmark_all) > 0
) {
  
  hallmark_all <-
    hallmark_all %>%
    dplyr::arrange(
      Temporal_Cluster,
      p.adjust
    )
}


add_log(
  paste(
    "Hallmark pathways preserved:",
    nrow(hallmark_all)
  )
)


# ============================================================
# 25. COMBINE ALL ENRICHMENT DATABASES
#
# No enrichment is performed here.
# Results are copied from frozen Stage17D/17E objects.
# ============================================================

all_enrichment <-
  dplyr::bind_rows(
    go_all,
    kegg_all,
    reactome_all,
    hallmark_all
  )


if (
  nrow(all_enrichment) > 0
) {
  
  all_enrichment <-
    all_enrichment %>%
    dplyr::arrange(
      Temporal_Cluster,
      Database,
      p.adjust
    )
}


add_log(
  paste(
    "Combined enrichment rows:",
    nrow(all_enrichment)
  )
)


# ============================================================
# 26. DEFINE BIOLOGICAL THEMES
#
# Theme assignment is an integration-level annotation.
# It does not alter enrichment statistics.
#
# The classifier is intentionally conservative.
# ============================================================

assign_biological_theme <- function(
    pathway_description) {
  
  x <- tolower(
    as.character(
      pathway_description
    )
  )
  
  
  dplyr::case_when(
    
    # --------------------------------------------------------
    # Immune / antigen presentation / inflammatory signaling
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "interferon",
          "antigen processing",
          "antigen presentation",
          "mhc",
          "class ii",
          "immunological synapse",
          "t cell",
          "t-cell",
          "tcr",
          "cd3",
          "cd28",
          "zap-70",
          "zap70",
          "pd-1",
          "pd1",
          "pd-l1",
          "pd1",
          "complement",
          "tnfr",
          "nf-kb",
          "nfkb",
          "host defense",
          "lymphocyte",
          "immune response"
        ),
        collapse = "|"
      )
    ) ~
      "Immune / antigen-presentation / inflammatory processes",
    
    
    # --------------------------------------------------------
    # Protein synthesis / translation
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "translation",
          "ribosome",
          "ribosomal",
          "ribonucleoprotein",
          "translational",
          "ribosome quality",
          "nonsense-mediated",
          "nmd",
          "srp-dependent",
          "40s",
          "60s",
          "peptide chain elongation",
          "translation initiation",
          "translation termination",
          "translation elongation",
          "mrna translation"
        ),
        collapse = "|"
      )
    ) ~
      "Protein synthesis / translational processes",
    
    
    # --------------------------------------------------------
    # Mitochondrial / oxidative metabolism
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "oxidative phosphorylation",
          "mitochondrial translation",
          "mitochondrial respiration",
          "respiratory electron transport",
          "aerobic respiration",
          "electron transport chain",
          "respiratory chain",
          "mitochondrial",
          "cristae",
          "proton motive force",
          "atp synthesis coupled electron transport",
          "atp biosynthetic"
        ),
        collapse = "|"
      )
    ) ~
      "Mitochondrial / oxidative energy metabolism",
    
    
    # --------------------------------------------------------
    # Cell cycle / proliferation
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "cell cycle",
          "g1/s",
          "g1 s",
          "s phase",
          "g2/m",
          "g2 m",
          "mitotic",
          "cyclin",
          "cdk",
          "apc/c",
          "dna replication",
          "replication fork"
        ),
        collapse = "|"
      )
    ) ~
      "Cell cycle / proliferation",
    
    
    # --------------------------------------------------------
    # MYC / mTOR / growth
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "myc",
          "mtorc",
          "mtor",
          "growth factor",
          "growth signaling"
        ),
        collapse = "|"
      )
    ) ~
      "Growth / MYC / mTOR signaling",
    
    
    # --------------------------------------------------------
    # Stress / hypoxia / redox
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "reactive oxygen",
          "oxidative stress",
          "hypoxia",
          "hif",
          "nfe2l2",
          "nrf2",
          "keap1",
          "starvation",
          "stress response"
        ),
        collapse = "|"
      )
    ) ~
      "Cellular stress / hypoxia / redox",
    
    
    # --------------------------------------------------------
    # DNA damage / repair
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "p53",
          "dna damage",
          "dna repair",
          "dna damage checkpoint",
          "nucleotide excision",
          "checkpoint",
          "repair"
        ),
        collapse = "|"
      )
    ) ~
      "DNA damage / repair / checkpoint",
    
    
    # --------------------------------------------------------
    # Adhesion / extracellular interactions / EMT
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "epithelial-mesenchymal",
          "epithelial mesenchymal",
          "emt",
          "cell adhesion",
          "cell-cell adhesion",
          "cell matrix",
          "extracellular matrix",
          "integrin",
          "cadherin",
          "cell migration",
          "cell junction",
          "junction organization"
        ),
        collapse = "|"
      )
    ) ~
      "Cell adhesion / extracellular interactions",
    
    
    # --------------------------------------------------------
    # General metabolism
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "metabolism",
          "metabolic",
          "amino acid",
          "nucleotide metabolism",
          "lipid metabolism",
          "carbon metabolism"
        ),
        collapse = "|"
      )
    ) ~
      "Metabolic processes",
    
    
    # --------------------------------------------------------
    # Cell death
    # --------------------------------------------------------
    
    stringr::str_detect(
      x,
      paste(
        c(
          "apoptosis",
          "programmed cell death",
          "necroptosis",
          "cell death"
        ),
        collapse = "|"
      )
    ) ~
      "Cell death / apoptosis",
    
    
    # --------------------------------------------------------
    # Default
    # --------------------------------------------------------
    
    TRUE ~
      "Other biological processes"
  )
}


# ============================================================
# 27. ADD BIOLOGICAL THEMES
# ============================================================

if (
  nrow(all_enrichment) > 0
) {
  
  all_enrichment <-
    all_enrichment %>%
    dplyr::mutate(
      Biological_Theme =
        assign_biological_theme(
          Description
        )
    )
}


# ============================================================
# 28. CREATE DATABASE-LEVEL SUMMARY
# ============================================================

database_names <- c(
  "GO_BP",
  "KEGG",
  "Reactome",
  "Hallmark"
)

database_summary <- tibble::tibble()


for (
  cl in cluster_ids
) {
  
  for (
    db in database_names
  ) {
    
    subset_data <-
      all_enrichment %>%
      dplyr::filter(
        Temporal_Cluster == cl,
        Database == db
      )
    
    
    significant_count <-
      sum(
        !is.na(
          subset_data$p.adjust
        ) &
          subset_data$p.adjust < 0.05
      )
    
    
    database_summary <-
      dplyr::bind_rows(
        database_summary,
        tibble::tibble(
          
          Temporal_Cluster =
            cl,
          
          Database =
            db,
          
          Total_Pathways =
            nrow(
              subset_data
            ),
          
          Significant_Pathways =
            significant_count
        )
      )
  }
}


# ============================================================
# 29. CREATE THEME-LEVEL SUMMARY
# ============================================================

theme_summary <-
  all_enrichment %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  ) %>%
  dplyr::group_by(
    Temporal_Cluster,
    Biological_Theme
  ) %>%
  dplyr::summarise(
    
    Database_Count =
      dplyr::n_distinct(
        Database
      ),
    
    Pathway_Count =
      dplyr::n(),
    
    Best_Adjusted_P =
      min(
        p.adjust,
        na.rm = TRUE
      ),
    
    Representative_Pathway =
      Description[
        which.min(
          p.adjust
        )
      ],
    
    .groups =
      "drop"
  ) %>%
  dplyr::arrange(
    Temporal_Cluster,
    Best_Adjusted_P
  )


# ============================================================
# 30. CREATE CROSS-DATABASE THEME SUMMARY
#
# "Database support" describes recurrence across databases.
# It does NOT imply statistical independence.
# ============================================================

cross_database_theme <-
  all_enrichment %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  ) %>%
  dplyr::group_by(
    Temporal_Cluster,
    Biological_Theme
  ) %>%
  dplyr::summarise(
    
    Supporting_Databases =
      paste(
        sort(
          unique(
            Database
          )
        ),
        collapse = "; "
      ),
    
    Number_of_Databases =
      dplyr::n_distinct(
        Database
      ),
    
    Number_of_Pathways =
      dplyr::n(),
    
    Best_Adjusted_P =
      min(
        p.adjust,
        na.rm = TRUE
      ),
    
    Representative_Pathway =
      Description[
        which.min(
          p.adjust
        )
      ],
    
    .groups =
      "drop"
  ) %>%
  dplyr::arrange(
    Temporal_Cluster,
    dplyr::desc(
      Number_of_Databases
    ),
    Best_Adjusted_P
  )


# ============================================================
# 31. CREATE REPRESENTATIVE PATHWAY TABLE
#
# This is a compact publication-oriented view.
# Complete enrichment results remain preserved separately.
# ============================================================

representative_pathways <-
  all_enrichment %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  ) %>%
  dplyr::group_by(
    Temporal_Cluster,
    Biological_Theme
  ) %>%
  dplyr::arrange(
    p.adjust,
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = 5
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    Temporal_Cluster,
    Biological_Theme,
    p.adjust
  )


# ============================================================
# 32. CREATE INTEGRATED CLUSTER INTERPRETATION
#
# IMPORTANT:
# Biological interpretation remains conservative.
#
# Cluster 1:
#   Early-peaking immune/antigen-presentation-associated
#   program with mitochondrial respiratory activity.
#
# Cluster 2:
#   Broad late-emerging cellular-state program with
#   increased translational and mitochondrial energy
#   metabolism and an epithelial-associated component
#   supported by temporal gene composition rather than
#   direct epithelial GO enrichment.
# ============================================================

get_pattern_value <- function(
    cluster_id
) {
  
  result <-
    temporal_summary %>%
    dplyr::filter(
      Temporal_Cluster ==
        cluster_id
    )
  
  if (
    nrow(result) == 0
  ) {
    
    return(
      NA_character_
    )
  }
  
  as.character(
    result$Pattern[[1]]
  )
}


get_direction_value <- function(
    cluster_id
) {
  
  result <-
    temporal_summary %>%
    dplyr::filter(
      Temporal_Cluster ==
        cluster_id
    )
  
  if (
    nrow(result) == 0
  ) {
    
    return(
      NA_character_
    )
  }
  
  as.character(
    result$Temporal_Direction[[1]]
  )
}


get_theme_names <- function(
    cluster_id
) {
  
  result <-
    cross_database_theme %>%
    dplyr::filter(
      Temporal_Cluster ==
        cluster_id
    ) %>%
    dplyr::arrange(
      dplyr::desc(
        Number_of_Databases
      ),
      Best_Adjusted_P
    )
  
  if (
    nrow(result) == 0
  ) {
    
    return(
      NA_character_
    )
  }
  
  paste(
    result$Biological_Theme,
    collapse = "; "
  )
}


integrated_interpretation <-
  tibble::tibble()


for (
  cl in cluster_ids
) {
  
  pattern <-
    get_pattern_value(
      cl
    )
  
  direction <-
    get_direction_value(
      cl
    )
  
  gene_count <-
    temporal_gene_counts[[cl]]
  
  
  themes <-
    get_theme_names(
      cl
    )
  
  
  if (
    identical(
      pattern,
      "Early_Peaking_Decreasing"
    )
  ) {
    
    interpretation <-
      paste(
        "An early-peaking temporal gene-expression",
        "program is followed by attenuation along",
        "the inferred pseudotime. Multi-database",
        "enrichment supports antigen processing and",
        "presentation, MHC-associated immune functions",
        "and immune-related processes, together with",
        "mitochondrial respiratory activity. This pattern",
        "describes a trajectory-associated transcriptional",
        "program and does not establish a transition",
        "from an immune cell state to another cell type."
      )
    
  } else if (
    identical(
      pattern,
      "Late_Increasing"
    )
  ) {
    
    interpretation <-
      paste(
        "A broad late-emerging temporal gene-expression",
        "program increases toward later pseudotime.",
        "Enrichment is dominated by protein synthesis,",
        "translation-related processes and mitochondrial",
        "oxidative energy metabolism, with additional",
        "growth, stress or cell-cycle-associated signals",
        "where supported by the enrichment results.",
        "The temporal gene composition also contains",
        "epithelial-associated genes, but the enrichment",
        "results do not by themselves establish a distinct",
        "epithelial cell-state transition."
      )
    
  } else {
    
    interpretation <-
      paste(
        "A trajectory-associated transcriptional program",
        "with the observed temporal pattern:",
        direction,
        "."
      )
  }
  
  
  integrated_interpretation <-
    dplyr::bind_rows(
      integrated_interpretation,
      tibble::tibble(
        
        Temporal_Cluster =
          cl,
        
        Gene_Count =
          gene_count,
        
        Pattern =
          pattern,
        
        Temporal_Direction =
          direction,
        
        Major_Biological_Themes =
          themes,
        
        Integrated_Interpretation =
          interpretation
      )
    )
}


# ============================================================
# 33. CREATE FINAL INTEGRATED TABLE
# ============================================================

get_database_themes <- function(
    cluster_id,
    database_name
) {
  
  result <-
    all_enrichment %>%
    dplyr::filter(
      Temporal_Cluster ==
        cluster_id,
      Database ==
        database_name,
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::group_by(
      Biological_Theme
    ) %>%
    dplyr::summarise(
      Best_P =
        min(
          p.adjust,
          na.rm = TRUE
        ),
      .groups =
        "drop"
    ) %>%
    dplyr::arrange(
      Best_P
    )
  
  
  if (
    nrow(result) == 0
  ) {
    
    return(
      NA_character_
    )
  }
  
  
  paste(
    result$Biological_Theme,
    collapse = "; "
  )
}


get_representative_pathways <- function(
    cluster_id,
    database_name,
    n = 5
) {
  
  result <-
    all_enrichment %>%
    dplyr::filter(
      Temporal_Cluster ==
        cluster_id,
      Database ==
        database_name,
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::arrange(
      p.adjust
    ) %>%
    dplyr::slice_head(
      n = n
    )
  
  
  if (
    nrow(result) == 0
  ) {
    
    return(
      NA_character_
    )
  }
  
  
  paste(
    result$Description,
    collapse = "; "
  )
}


final_integration_table <-
  tibble::tibble()


for (
  cl in cluster_ids
) {
  
  current_interpretation <-
    integrated_interpretation %>%
    dplyr::filter(
      Temporal_Cluster ==
        cl
    )
  
  
  current_row <-
    tibble::tibble(
      
      Temporal_Cluster =
        cl,
      
      Gene_Count =
        temporal_gene_counts[[cl]],
      
      Temporal_Pattern =
        get_pattern_value(
          cl
        ),
      
      Temporal_Direction =
        get_direction_value(
          cl
        ),
      
      GO_BP_Themes =
        get_database_themes(
          cl,
          "GO_BP"
        ),
      
      KEGG_Themes =
        get_database_themes(
          cl,
          "KEGG"
        ),
      
      Reactome_Themes =
        get_database_themes(
          cl,
          "Reactome"
        ),
      
      Hallmark_Themes =
        get_database_themes(
          cl,
          "Hallmark"
        ),
      
      GO_BP_Representative_Pathways =
        get_representative_pathways(
          cl,
          "GO_BP"
        ),
      
      KEGG_Representative_Pathways =
        get_representative_pathways(
          cl,
          "KEGG"
        ),
      
      Reactome_Representative_Pathways =
        get_representative_pathways(
          cl,
          "Reactome"
        ),
      
      Hallmark_Representative_Pathways =
        get_representative_pathways(
          cl,
          "Hallmark"
        ),
      
      Integrated_Interpretation =
        current_interpretation$
        Integrated_Interpretation[[1]]
    )
  
  
  final_integration_table <-
    dplyr::bind_rows(
      final_integration_table,
      current_row
    )
}


# ============================================================
# 34. CREATE COMPACT CLUSTER COMPARISON
# ============================================================

cluster_comparison <-
  final_integration_table %>%
  dplyr::select(
    Temporal_Cluster,
    Gene_Count,
    Temporal_Pattern,
    Temporal_Direction,
    Integrated_Interpretation
  )


# ============================================================
# 35. CREATE FINAL DATABASE SUMMARY
# ============================================================

integrated_significant_pathways <-
  sum(
    !is.na(
      all_enrichment$p.adjust
    ) &
      all_enrichment$p.adjust < 0.05
  )


number_of_themes <-
  dplyr::n_distinct(
    all_enrichment$Biological_Theme
  )


number_of_cross_database_themes <-
  if (
    nrow(cross_database_theme) > 0
  ) {
    
    sum(
      cross_database_theme$
        Number_of_Databases >= 2,
      na.rm = TRUE
    )
    
  } else {
    
    0L
  }


final_summary <-
  tibble::tibble(
    
    Stage =
      "17F",
    
    Description =
      "Final integration, validation, export and archiving of Stage17 temporal gene dynamics",
    
    Temporal_Genes =
      total_temporal_genes,
    
    Number_of_Temporal_Clusters =
      length(
        cluster_ids
      ),
    
    Cluster_1_Genes =
      as.integer(
        temporal_gene_counts[["1"]]
      ),
    
    Cluster_2_Genes =
      as.integer(
        temporal_gene_counts[["2"]]
      ),
    
    GO_Pathways =
      nrow(
        go_all
      ),
    
    KEGG_Pathways =
      nrow(
        kegg_all
      ),
    
    Reactome_Pathways =
      nrow(
        reactome_all
      ),
    
    Hallmark_Pathways =
      nrow(
        hallmark_all
      ),
    
    Integrated_Significant_Pathways =
      integrated_significant_pathways,
    
    Biological_Themes =
      number_of_themes,
    
    Cross_Database_Supported_Themes =
      number_of_cross_database_themes,
    
    New_Enrichment_Performed =
      FALSE,
    
    New_Clustering_Performed =
      FALSE,
    
    Pseudotime_Recalculated =
      FALSE,
    
    Trajectory_Reconstructed =
      FALSE,
    
    Root_Selected =
      FALSE
  )


# ============================================================
# 36. CREATE HUMAN-READABLE FINAL BIOLOGICAL SUMMARY
# ============================================================

summary_lines <- c(
  
  "STAGE 17 FINAL BIOLOGICAL SUMMARY",
  "=================================",
  "",
  
  "Stage17 integrates trajectory-associated temporal",
  "gene-expression dynamics with GO Biological Process,",
  "KEGG, Reactome and MSigDB Hallmark enrichment.",
  "",
  
  "Temporal clusters represent gene-expression dynamics",
  "along inferred pseudotime and should NOT be interpreted",
  "as independent cell subtypes or direct evidence of",
  "cell-state conversion.",
  "",
  
  paste0(
    "Total temporal genes: ",
    total_temporal_genes
  ),
  
  paste0(
    "Temporal clusters: ",
    length(cluster_ids)
  ),
  
  paste0(
    "Cluster 1 genes: ",
    temporal_gene_counts[["1"]]
  ),
  
  paste0(
    "Cluster 2 genes: ",
    temporal_gene_counts[["2"]]
  ),
  
  paste0(
    "GO pathways: ",
    nrow(go_all)
  ),
  
  paste0(
    "KEGG pathways: ",
    nrow(kegg_all)
  ),
  
  paste0(
    "Reactome pathways: ",
    nrow(reactome_all)
  ),
  
  paste0(
    "Hallmark pathways: ",
    nrow(hallmark_all)
  ),
  
  paste0(
    "Integrated significant pathways: ",
    integrated_significant_pathways
  ),
  
  paste0(
    "Biological themes: ",
    number_of_themes
  ),
  
  paste0(
    "Themes supported by at least two databases: ",
    number_of_cross_database_themes
  ),
  
  ""
)


for (
  cl in cluster_ids
) {
  
  current <-
    final_integration_table %>%
    dplyr::filter(
      Temporal_Cluster ==
        cl
    )
  
  
  if (
    nrow(current) == 0
  ) {
    
    next
  }
  
  
  summary_lines <-
    c(
      summary_lines,
      
      paste0(
        "Cluster ",
        cl
      ),
      
      "-------------------------",
      
      paste0(
        "Gene count: ",
        current$Gene_Count[[1]]
      ),
      
      paste0(
        "Temporal pattern: ",
        current$Temporal_Pattern[[1]]
      ),
      
      paste0(
        "Temporal direction: ",
        current$Temporal_Direction[[1]]
      ),
      
      "",
      
      "Integrated interpretation:",
      
      current$
        Integrated_Interpretation[[1]],
      
      "",
      
      paste0(
        "GO themes: ",
        current$GO_BP_Themes[[1]]
      ),
      
      paste0(
        "KEGG themes: ",
        current$KEGG_Themes[[1]]
      ),
      
      paste0(
        "Reactome themes: ",
        current$Reactome_Themes[[1]]
      ),
      
      paste0(
        "Hallmark themes: ",
        current$Hallmark_Themes[[1]]
      ),
      
      ""
    )
}


summary_lines <-
  c(
    summary_lines,
    
    "Scientific caution:",
    "These results identify transcriptional programs",
    "associated with inferred pseudotime.",
    "They do not establish causality or direct biological",
    "conversion between cell states or cell types.",
    "",
    
    "Cross-database recurrence indicates concordance of",
    "annotations across databases but does not imply that",
    "the databases are statistically independent.",
    "",
    
    "Reactome pathway counts should not be interpreted as",
    "a direct measure of biological activity because",
    "pathway sizes and hierarchical redundancy differ.",
    "",
    
    "Viral-related pathway annotations should not be",
    "interpreted as evidence of viral infection without",
    "independent experimental evidence."
  )


# ============================================================
# 37. EXPORT FINAL TABLES
# ============================================================

write.csv(
  temporal_summary,
  file.path(
    tables_main_dir,
    "Stage17F_Temporal_Dynamics_Summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  final_temporal_gene_table,
  file.path(
    tables_main_dir,
    "Stage17F_Temporal_Gene_List.csv"
  ),
  row.names = FALSE
)


write.csv(
  final_integration_table,
  file.path(
    tables_main_dir,
    "Stage17F_Final_Integrated_Temporal_Pathway_Table.csv"
  ),
  row.names = FALSE
)


write.csv(
  cluster_comparison,
  file.path(
    tables_main_dir,
    "Stage17F_Final_Cluster_Comparison.csv"
  ),
  row.names = FALSE
)


write.csv(
  integrated_interpretation,
  file.path(
    tables_main_dir,
    "Stage17F_Integrated_Cluster_Interpretation.csv"
  ),
  row.names = FALSE
)


write.csv(
  final_summary,
  file.path(
    tables_main_dir,
    "Stage17F_Final_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 38. EXPORT COMPLETE ENRICHMENT TABLES
# ============================================================

write.csv(
  go_all,
  file.path(
    tables_supp_dir,
    "Stage17F_GO_BP_All_Clusters.csv"
  ),
  row.names = FALSE
)


write.csv(
  kegg_all,
  file.path(
    tables_supp_dir,
    "Stage17F_KEGG_All_Clusters.csv"
  ),
  row.names = FALSE
)


write.csv(
  reactome_all,
  file.path(
    tables_supp_dir,
    "Stage17F_Reactome_All_Clusters.csv"
  ),
  row.names = FALSE
)


write.csv(
  hallmark_all,
  file.path(
    tables_supp_dir,
    "Stage17F_Hallmark_All_Clusters.csv"
  ),
  row.names = FALSE
)


write.csv(
  all_enrichment,
  file.path(
    tables_supp_dir,
    "Stage17F_All_Database_Enrichment_With_Themes.csv"
  ),
  row.names = FALSE
)


write.csv(
  database_summary,
  file.path(
    tables_supp_dir,
    "Stage17F_Database_Level_Summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  theme_summary,
  file.path(
    tables_supp_dir,
    "Stage17F_Biological_Theme_Summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  cross_database_theme,
  file.path(
    tables_supp_dir,
    "Stage17F_Cross_Database_Biological_Themes.csv"
  ),
  row.names = FALSE
)


write.csv(
  representative_pathways,
  file.path(
    tables_main_dir,
    "Stage17F_Representative_Pathways_By_Theme.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 39. EXPORT FINAL BIOLOGICAL SUMMARY TEXT
# ============================================================

writeLines(
  summary_lines,
  con = file.path(
    tables_main_dir,
    "Stage17F_Final_Biological_Summary.txt"
  )
)


# ============================================================
# 40. CREATE FINAL TEMPORAL CLUSTER SIZE FIGURE
# ============================================================

cluster_plot_data <-
  tibble::tibble(
    
    Temporal_Cluster =
      names(
        temporal_gene_counts
      ),
    
    Gene_Count =
      as.integer(
        temporal_gene_counts
      )
  )


p_cluster_size <-
  ggplot2::ggplot(
    cluster_plot_data,
    ggplot2::aes(
      x =
        Temporal_Cluster,
      y =
        Gene_Count
    )
  ) +
  
  ggplot2::geom_col() +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label =
        Gene_Count
    ),
    vjust =
      -0.4,
    size =
      4
  ) +
  
  ggplot2::labs(
    title =
      "Temporal Gene Counts by Cluster",
    
    x =
      "Temporal Cluster",
    
    y =
      "Number of Temporal Genes"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face =
          "bold",
        size =
          13
      )
  )


ggplot2::ggsave(
  filename =
    file.path(
      figures_main_dir,
      "Stage17F_Temporal_Cluster_Gene_Counts.pdf"
    ),
  plot =
    p_cluster_size,
  width =
    7,
  height =
    5
)


# ============================================================
# 41. CREATE CROSS-DATABASE THEME FIGURE
# ============================================================

if (
  nrow(cross_database_theme) > 0
) {
  
  theme_plot_data <-
    cross_database_theme
  
  
  theme_plot_data$Biological_Theme <-
    factor(
      theme_plot_data$Biological_Theme,
      levels =
        rev(
          unique(
            theme_plot_data$Biological_Theme
          )
        )
    )
  
  
  p_theme <-
    ggplot2::ggplot(
      theme_plot_data,
      ggplot2::aes(
        x =
          Number_of_Databases,
        y =
          Biological_Theme,
        size =
          Number_of_Pathways
      )
    ) +
    
    ggplot2::geom_point() +
    
    ggplot2::facet_wrap(
      ~ Temporal_Cluster,
      scales =
        "free_y"
    ) +
    
    ggplot2::scale_x_continuous(
      breaks =
        1:4,
      limits =
        c(
          0.5,
          4.5
        )
    ) +
    
    ggplot2::labs(
      title =
        "Cross-Database Recurrence of Biological Themes",
      
      subtitle =
        "Database recurrence is descriptive and does not imply statistical independence",
      
      x =
        "Number of Enrichment Databases",
      
      y =
        "Biological Theme",
      
      size =
        "Number of Significant Pathways"
    ) +
    
    ggplot2::theme_bw() +
    
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          face =
            "bold",
          size =
            13
        ),
      
      plot.subtitle =
        ggplot2::element_text(
          size =
            9
        ),
      
      axis.text.y =
        ggplot2::element_text(
          size =
            9
        ),
      
      strip.text =
        ggplot2::element_text(
          face =
            "bold"
        )
    )
  
  
  ggplot2::ggsave(
    filename =
      file.path(
        figures_main_dir,
        "Stage17F_Cross_Database_Biological_Themes.pdf"
      ),
    plot =
      p_theme,
    width =
      11,
    height =
      7
  )
}


# ============================================================
# 42. CREATE STAGE17C-D-E PROVENANCE TABLE
# ============================================================

provenance_table <-
  tibble::tibble(
    
    Final_Stage =
      "17F",
    
    Stage17C_Source =
      normalizePath(
        stage17c_rds,
        winslash =
          "/",
        mustWork =
          FALSE
      ),
    
    Stage17D_Source =
      normalizePath(
        stage17d_rds,
        winslash =
          "/",
        mustWork =
          FALSE
      ),
    
    Stage17E_Source =
      normalizePath(
        stage17e_rds,
        winslash =
          "/",
        mustWork =
          FALSE
      ),
    
    Analysis_Status =
      "Final integration, validation, export and archiving",
    
    New_Enrichment_Performed =
      FALSE,
    
    New_Clustering_Performed =
      FALSE,
    
    Pseudotime_Recalculated =
      FALSE,
    
    Trajectory_Reconstructed =
      FALSE,
    
    Root_Selected =
      FALSE,
    
    Pattern_Based_Gene_Subsetting =
      FALSE,
    
    Final_RDS =
      normalizePath(
        final_rds,
        winslash =
          "/",
        mustWork =
          FALSE
      )
  )


write.csv(
  provenance_table,
  file.path(
    tables_validation_dir,
    "Stage17F_Provenance.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 43. FINAL VALIDATION
# ============================================================

validation_table <-
  tibble::tibble(
    Check =
      character(),
    
    Value =
      character(),
    
    Status =
      character()
  )


add_validation <- function(
    check,
    value,
    status
) {
  
  validation_table <<-
    dplyr::bind_rows(
      validation_table,
      
      tibble::tibble(
        
        Check =
          check,
        
        Value =
          as.character(
            value
          ),
        
        Status =
          status
      )
    )
}


# ------------------------------------------------------------
# Check 1
# ------------------------------------------------------------

add_validation(
  "Stage17C RDS exists",
  file.exists(
    stage17c_rds
  ),
  ifelse(
    file.exists(
      stage17c_rds
    ),
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 2
# ------------------------------------------------------------

add_validation(
  "Stage17D RDS exists",
  file.exists(
    stage17d_rds
  ),
  ifelse(
    file.exists(
      stage17d_rds
    ),
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 3
# ------------------------------------------------------------

add_validation(
  "Stage17E RDS exists",
  file.exists(
    stage17e_rds
  ),
  ifelse(
    file.exists(
      stage17e_rds
    ),
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 4
# ------------------------------------------------------------

add_validation(
  "Stage17C structure validated",
  length(
    missing_stage17c_fields
  ) == 0,
  ifelse(
    length(
      missing_stage17c_fields
    ) == 0,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 5
# ------------------------------------------------------------

add_validation(
  "Stage17D structure validated",
  length(
    missing_stage17d_fields
  ) == 0,
  ifelse(
    length(
      missing_stage17d_fields
    ) == 0,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 6
# ------------------------------------------------------------

add_validation(
  "Stage17E structure validated",
  length(
    missing_stage17e_fields
  ) == 0,
  ifelse(
    length(
      missing_stage17e_fields
    ) == 0,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 7
# ------------------------------------------------------------

add_validation(
  "Temporal clusters",
  length(
    cluster_ids
  ),
  ifelse(
    length(
      cluster_ids
    ) == 2,
    "PASS",
    "CHECK"
  )
)


# ------------------------------------------------------------
# Check 8
# ------------------------------------------------------------

add_validation(
  "Total temporal genes",
  total_temporal_genes,
  ifelse(
    total_temporal_genes == 569,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 9
# ------------------------------------------------------------

add_validation(
  "Cluster 1 gene count",
  temporal_gene_counts[["1"]],
  ifelse(
    temporal_gene_counts[["1"]] == 53,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 10
# ------------------------------------------------------------

add_validation(
  "Cluster 2 gene count",
  temporal_gene_counts[["2"]],
  ifelse(
    temporal_gene_counts[["2"]] == 516,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 11
# ------------------------------------------------------------

add_validation(
  "Cluster 1 temporal pattern",
  pattern_cluster1,
  ifelse(
    pattern_cluster1 ==
      "Early_Peaking_Decreasing",
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 12
# ------------------------------------------------------------

add_validation(
  "Cluster 2 temporal pattern",
  pattern_cluster2,
  ifelse(
    pattern_cluster2 ==
      "Late_Increasing",
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 13
# ------------------------------------------------------------

add_validation(
  "Temporal gene table row count",
  nrow(
    final_temporal_gene_table
  ),
  ifelse(
    nrow(
      final_temporal_gene_table
    ) ==
      total_temporal_genes,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 14
# ------------------------------------------------------------

add_validation(
  "Duplicated genes across clusters",
  nrow(
    duplicate_gene_table
  ),
  ifelse(
    nrow(
      duplicate_gene_table
    ) == 0,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 15
# ------------------------------------------------------------

add_validation(
  "GO enrichment table available",
  nrow(
    go_all
  ),
  ifelse(
    is.data.frame(
      go_all
    ),
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 16
# ------------------------------------------------------------

add_validation(
  "KEGG enrichment table available",
  nrow(
    kegg_all
  ),
  ifelse(
    is.data.frame(
      kegg_all
    ),
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 17
# ------------------------------------------------------------

add_validation(
  "Reactome enrichment table available",
  nrow(
    reactome_all
  ),
  ifelse(
    nrow(
      reactome_all
    ) > 0,
    "PASS",
    "CHECK"
  )
)


# ------------------------------------------------------------
# Check 18
# ------------------------------------------------------------

add_validation(
  "Hallmark enrichment table available",
  nrow(
    hallmark_all
  ),
  ifelse(
    nrow(
      hallmark_all
    ) > 0,
    "PASS",
    "CHECK"
  )
)


# ------------------------------------------------------------
# Check 19
# ------------------------------------------------------------

add_validation(
  "Combined enrichment table available",
  nrow(
    all_enrichment
  ),
  ifelse(
    nrow(
      all_enrichment
    ) > 0,
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 20
# ------------------------------------------------------------

add_validation(
  "Final integration table has one row per cluster",
  nrow(
    final_integration_table
  ),
  ifelse(
    nrow(
      final_integration_table
    ) ==
      length(
        cluster_ids
      ),
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 21
# ------------------------------------------------------------

add_validation(
  "Integrated interpretation has one row per cluster",
  nrow(
    integrated_interpretation
  ),
  ifelse(
    nrow(
      integrated_interpretation
    ) ==
      length(
        cluster_ids
      ),
    "PASS",
    "FAIL"
  )
)


# ------------------------------------------------------------
# Check 22
# ------------------------------------------------------------

add_validation(
  "No new enrichment performed",
  TRUE,
  "PASS"
)


# ------------------------------------------------------------
# Check 23
# ------------------------------------------------------------

add_validation(
  "No new clustering performed",
  TRUE,
  "PASS"
)


# ------------------------------------------------------------
# Check 24
# ------------------------------------------------------------

add_validation(
  "No pseudotime recalculation performed",
  TRUE,
  "PASS"
)


# ------------------------------------------------------------
# Check 25
# ------------------------------------------------------------

add_validation(
  "No trajectory reconstruction performed",
  TRUE,
  "PASS"
)


# ------------------------------------------------------------
# Check 26
# ------------------------------------------------------------

add_validation(
  "No root selection performed",
  TRUE,
  "PASS"
)


# ------------------------------------------------------------
# Check 27
# ------------------------------------------------------------

add_validation(
  "No Pattern-based gene subsetting performed",
  TRUE,
  "PASS"
)


# ------------------------------------------------------------
# Check 28
# ------------------------------------------------------------

add_validation(
  "Final biological summary created",
  length(
    summary_lines
  ),
  ifelse(
    length(
      summary_lines
    ) > 0,
    "PASS",
    "FAIL"
  )
)


# ============================================================
# 44. SAVE VALIDATION TABLE
# ============================================================

write.csv(
  validation_table,
  file.path(
    tables_validation_dir,
    "Stage17F_Final_Validation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 45. CREATE FINAL STAGE17F RDS
# ============================================================

stage17f_results <- list(
  
  stage =
    "17F",
  
  description =
    "Final integration, validation, export and archiving of Stage17 temporal gene dynamics",
  
  temporal_gene_lists =
    cluster_gene_lists,
  
  temporal_gene_counts =
    temporal_gene_counts,
  
  temporal_gene_table =
    final_temporal_gene_table,
  
  temporal_behavior =
    temporal_summary,
  
  cluster_pattern_summary =
    cluster_pattern_summary,
  
  biological_program_description =
    biological_program_description,
  
  go_all =
    go_all,
  
  kegg_all =
    kegg_all,
  
  reactome_all =
    reactome_all,
  
  hallmark_all =
    hallmark_all,
  
  all_enrichment =
    all_enrichment,
  
  database_summary =
    database_summary,
  
  theme_summary =
    theme_summary,
  
  cross_database_theme =
    cross_database_theme,
  
  representative_pathways =
    representative_pathways,
  
  integrated_interpretation =
    integrated_interpretation,
  
  final_integration_table =
    final_integration_table,
  
  cluster_comparison =
    cluster_comparison,
  
  final_summary =
    final_summary,
  
  validation =
    validation_table,
  
  provenance =
    provenance_table,
  
  stage17c_source =
    normalizePath(
      stage17c_rds,
      winslash =
        "/",
      mustWork =
        FALSE
    ),
  
  stage17d_source =
    normalizePath(
      stage17d_rds,
      winslash =
        "/",
      mustWork =
        FALSE
    ),
  
  stage17e_source =
    normalizePath(
      stage17e_rds,
      winslash =
        "/",
      mustWork =
        FALSE
    ),
  
  final_output_directory =
    normalizePath(
      stage17_results_dir,
      winslash =
        "/",
      mustWork =
        FALSE
    ),
  
  trajectory_reconstruction =
    FALSE,
  
  root_selection =
    FALSE,
  
  pseudotime_recalculation =
    FALSE,
  
  new_enrichment =
    FALSE,
  
  new_clustering =
    FALSE,
  
  pattern_based_gene_subsetting =
    FALSE,
  
  start_time =
    start_time
)


# ============================================================
# 46. SAVE FINAL RDS
# ============================================================

saveRDS(
  stage17f_results,
  final_rds
)


add_log(
  paste(
    "Final canonical Stage17F RDS saved:",
    final_rds
  )
)


# ============================================================
# 47. RELOAD FINAL RDS FOR INTEGRITY CHECK
# ============================================================

reloaded_final <-
  readRDS(
    final_rds
  )


required_final_fields <- c(
  
  "stage",
  
  "temporal_gene_lists",
  
  "temporal_gene_counts",
  
  "temporal_gene_table",
  
  "temporal_behavior",
  
  "cluster_pattern_summary",
  
  "go_all",
  
  "kegg_all",
  
  "reactome_all",
  
  "hallmark_all",
  
  "all_enrichment",
  
  "database_summary",
  
  "theme_summary",
  
  "cross_database_theme",
  
  "representative_pathways",
  
  "integrated_interpretation",
  
  "final_integration_table",
  
  "cluster_comparison",
  
  "final_summary",
  
  "validation",
  
  "provenance"
)


missing_final_fields <-
  setdiff(
    required_final_fields,
    names(
      reloaded_final
    )
  )


final_rds_integrity_status <-
  if (
    length(
      missing_final_fields
    ) == 0
  ) {
    
    "PASS"
    
  } else {
    
    "FAIL"
  }


if (
  final_rds_integrity_status ==
  "FAIL"
) {
  
  stop(
    paste(
      "Final RDS integrity check failed. Missing field(s):",
      paste(
        missing_final_fields,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 48. ADD FINAL RDS INTEGRITY VALIDATION
# ============================================================

validation_table <-
  dplyr::bind_rows(
    validation_table,
    
    tibble::tibble(
      
      Check =
        "Final RDS integrity",
      
      Value =
        "All required fields present after reload",
      
      Status =
        final_rds_integrity_status
    )
  )


# ============================================================
# 49. SAVE UPDATED VALIDATION
# ============================================================

write.csv(
  validation_table,
  file.path(
    tables_validation_dir,
    "Stage17F_Final_Validation.csv"
  ),
  row.names = FALSE
)


# Update validation inside the final RDS.

stage17f_results$validation <-
  validation_table


saveRDS(
  stage17f_results,
  final_rds
)


# ============================================================
# 50. FINAL RUNTIME
# ============================================================

end_time <- Sys.time()

duration_minutes <-
  as.numeric(
    difftime(
      end_time,
      start_time,
      units =
        "mins"
    )
  )


stage17f_results$end_time <-
  end_time

stage17f_results$runtime_minutes <-
  duration_minutes


saveRDS(
  stage17f_results,
  final_rds
)


# ============================================================
# 51. FINAL LOG
# ============================================================

add_log(
  paste(
    "Temporal genes:",
    total_temporal_genes
  )
)

add_log(
  paste(
    "Cluster 1 genes:",
    temporal_gene_counts[["1"]]
  )
)

add_log(
  paste(
    "Cluster 2 genes:",
    temporal_gene_counts[["2"]]
  )
)

add_log(
  paste(
    "GO pathways:",
    nrow(go_all)
  )
)

add_log(
  paste(
    "KEGG pathways:",
    nrow(kegg_all)
  )
)

add_log(
  paste(
    "Reactome pathways:",
    nrow(reactome_all)
  )
)

add_log(
  paste(
    "Hallmark pathways:",
    nrow(hallmark_all)
  )
)

add_log(
  paste(
    "Integrated significant pathways:",
    integrated_significant_pathways
  )
)

add_log(
  paste(
    "Biological themes:",
    number_of_themes
  )
)

add_log(
  paste(
    "Themes supported by at least two databases:",
    number_of_cross_database_themes
  )
)

add_log(
  paste(
    "Final RDS integrity:",
    final_rds_integrity_status
  )
)

add_log(
  paste(
    "Runtime:",
    round(
      duration_minutes,
      3
    ),
    "minutes"
  )
)

add_log(
  "No upstream Stage17 analysis was rerun."
)

add_log(
  "Stage17F completed successfully."
)


writeLines(
  log_messages,
  con =
    log_file
)


# ============================================================
# 52. FINAL CONSOLE REPORT
# ============================================================

cat(
  
  "\n",
  
  "============================================================\n",
  
  "STAGE 17F COMPLETED SUCCESSFULLY\n",
  
  "FINAL INTEGRATION + VALIDATION + EXPORT + ARCHIVING\n",
  
  "============================================================\n",
  
  "\n",
  
  "Temporal genes:",
  total_temporal_genes,
  "\n",
  
  "Temporal clusters:",
  length(cluster_ids),
  "\n",
  
  "Cluster 1 genes:",
  temporal_gene_counts[["1"]],
  "\n",
  
  "Cluster 2 genes:",
  temporal_gene_counts[["2"]],
  "\n",
  
  "\n",
  
  "GO pathways:",
  nrow(go_all),
  "\n",
  
  "KEGG pathways:",
  nrow(kegg_all),
  "\n",
  
  "Reactome pathways:",
  nrow(reactome_all),
  "\n",
  
  "Hallmark pathways:",
  nrow(hallmark_all),
  "\n",
  
  "\n",
  
  "Integrated significant pathways:",
  integrated_significant_pathways,
  "\n",
  
  "Biological themes:",
  number_of_themes,
  "\n",
  
  "Themes supported by >=2 databases:",
  number_of_cross_database_themes,
  "\n",
  
  "\n",
  
  "Final RDS integrity:",
  final_rds_integrity_status,
  "\n",
  
  "\n",
  
  "Final canonical RDS:\n",
  final_rds,
  "\n",
  
  "\n",
  
  "Results directory:\n",
  stage17_results_dir,
  "\n",
  
  "\n",
  
  "No new enrichment performed.\n",
  
  "No new clustering performed.\n",
  
  "No pseudotime recalculation performed.\n",
  
  "No trajectory reconstruction performed.\n",
  
  "No root selection performed.\n",
  
  "No Pattern-based gene subsetting performed.\n",
  
  "\n",
  
  "============================================================\n",
  
  "\n"
)
