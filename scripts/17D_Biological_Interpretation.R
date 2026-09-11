# ==============================================================================
# STAGE 17D
# BIOLOGICAL INTERPRETATION OF TEMPORAL GENE CLUSTERS
# FINAL OPTIMIZED VERSION
# ==============================================================================
#
# Purpose:
#   Biological interpretation of the two temporal gene-expression programs
#   identified in Stage 17B and characterized in Stage 17C.
#
# Important:
#   Temporal clusters represent gene-expression dynamics along pseudotime.
#   They are NOT interpreted as independent cell subtypes.
#
# Inputs:
#   Stage17B_Temporal_Clustering_Results.rds
#   Stage17C_Temporal_Characterization_Results.rds
#
# Optional source:
#   TNBC_Epithelial_CopyKAT.rds
#
# Main analyses:
#   1. Recover validated temporal cluster assignments from Stage 17B
#   2. Recover temporal profiles from Stage 17B centroids
#   3. Summarize cluster composition and representative genes
#   4. Define cautious temporal program descriptions
#   5. GO Biological Process enrichment
#   6. KEGG enrichment
#   7. Generate temporal profile and enrichment figures
#   8. Save integrated Stage 17D object
#   9. Validate output and provenance
#
# No:
#   - trajectory reconstruction
#   - root selection
#   - pseudotime recalculation
#   - reclustering
#   - modification of Stage 17B or Stage 17C
#
# ==============================================================================


# ==============================================================================
# 1. PROJECT CONFIGURATION
# ==============================================================================

project_dir <- "YOUR_PROJECT_DIRECTORY"


results_dir <-
  file.path(
    project_dir,
    "results"
  )


objects_dir <-
  file.path(
    project_dir,
    "objects"
  )


stage17_dir <-
  file.path(
    results_dir,
    "17_Temporal_Gene_Dynamics"
  )


figures_main_dir <-
  file.path(
    stage17_dir,
    "figures",
    "main"
  )


figures_supp_dir <-
  file.path(
    stage17_dir,
    "figures",
    "supplementary"
  )


tables_main_dir <-
  file.path(
    stage17_dir,
    "tables",
    "main"
  )


tables_supp_dir <-
  file.path(
    stage17_dir,
    "tables",
    "supplementary"
  )


validation_dir <-
  file.path(
    stage17_dir,
    "tables",
    "validation"
  )


logs_dir <-
  file.path(
    stage17_dir,
    "logs"
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
  validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. LOGGING
# ==============================================================================

stage17d_log <-
  file.path(
    logs_dir,
    "Stage17D_Completion.log"
  )


writeLines(
  character(0),
  con = stage17d_log
)


write_log <- function(...) {
  
  txt <-
    paste0(
      ...,
      collapse = ""
    )
  
  cat(
    txt,
    "\n"
  )
  
  write(
    txt,
    file = stage17d_log,
    append = TRUE
  )
  
}


write_log(
  "============================================================"
)

write_log(
  "STAGE 17D: BIOLOGICAL INTERPRETATION"
)

write_log(
  "============================================================"
)

write_log(
  "Project directory: ",
  project_dir
)


# ==============================================================================
# 3. REQUIRED PACKAGES
# ==============================================================================

required_packages <-
  c(
    "dplyr",
    "ggplot2",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot"
  )


for (
  pkg in required_packages
) {
  
  if (
    !requireNamespace(
      pkg,
      quietly = TRUE
    )
  ) {
    
    stop(
      paste(
        "Required package is not installed:",
        pkg
      )
    )
    
  }
  
}


suppressPackageStartupMessages({
  
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
})


# ==============================================================================
# 4. INPUT FILES
# ==============================================================================

stage17b_object_rds <-
  file.path(
    objects_dir,
    "Stage17B_Temporal_Clustering_Results.rds"
  )


stage17c_object_rds <-
  file.path(
    objects_dir,
    "Stage17C_Temporal_Characterization_Results.rds"
  )


seurat_object_rds <-
  file.path(
    objects_dir,
    "TNBC_Epithelial_CopyKAT.rds"
  )


# ==============================================================================
# 5. INPUT FILE VALIDATION
# ==============================================================================

if (
  !file.exists(
    stage17b_object_rds
  )
) {
  
  stop(
    paste(
      "Stage 17B object not found:",
      stage17b_object_rds
    )
  )
  
}


if (
  !file.exists(
    stage17c_object_rds
  )
) {
  
  stop(
    paste(
      "Stage 17C object not found:",
      stage17c_object_rds
    )
  )
  
}


write_log(
  "Stage 17B object found."
)


write_log(
  "Stage 17C object found."
)


# ==============================================================================
# 6. LOAD STAGE 17B AND STAGE 17C
# ==============================================================================

stage17b_results <-
  readRDS(
    stage17b_object_rds
  )


stage17c_results <-
  readRDS(
    stage17c_object_rds
  )


if (
  !is.list(
    stage17b_results
  )
) {
  
  stop(
    "Stage 17B object is not a list."
  )
  
}


if (
  !is.list(
    stage17c_results
  )
) {
  
  stop(
    "Stage 17C object is not a list."
  )
  
}


write_log(
  "Stage 17B loaded successfully."
)


write_log(
  "Stage 17C loaded successfully."
)


# ==============================================================================
# 7. VALIDATE STAGE 17B FROZEN STRUCTURE
# ==============================================================================

required_17b_fields <-
  c(
    "stage",
    "gene_clusters",
    "cluster_assignments",
    "cluster_centroids",
    "cluster_summary",
    "gene_cluster_correlation_table",
    "strong_temporal_genes",
    "strong_temporal_genes_in_expression",
    "n_bins",
    "bin_breaks",
    "pseudotime",
    "synchronized_cells",
    "upstream_frozen",
    "trajectory_reconstruction",
    "root_selection",
    "pseudotime_recalculation"
  )


missing_17b_fields <-
  setdiff(
    required_17b_fields,
    names(
      stage17b_results
    )
  )


if (
  length(
    missing_17b_fields
  ) > 0
) {
  
  stop(
    paste(
      "Stage 17B is missing required fields:",
      paste(
        missing_17b_fields,
        collapse = ", "
      )
    )
  )
  
}


if (
  !isTRUE(
    stage17b_results$upstream_frozen
  )
) {
  
  stop(
    "Stage 17B upstream_frozen flag is not TRUE."
  )
  
}


if (
  isTRUE(
    stage17b_results$trajectory_reconstruction
  )
) {
  
  stop(
    "Stage 17B indicates trajectory reconstruction was performed."
  )
  
}


if (
  isTRUE(
    stage17b_results$root_selection
  )
) {
  
  stop(
    "Stage 17B indicates root selection was performed."
  )
  
}


if (
  isTRUE(
    stage17b_results$pseudotime_recalculation
  )
) {
  
  stop(
    "Stage 17B indicates pseudotime recalculation was performed."
  )
  
}


# ==============================================================================
# 8. VALIDATE STAGE 17C
# ==============================================================================

if (
  !(
    "stage" %in%
    names(
      stage17c_results
    )
  )
) {
  
  stop(
    "Stage 17C object does not contain a stage field."
  )
  
}


if (
  !(
    "stage17c" %in%
    tolower(
      as.character(
        stage17c_results$stage
      )
    )
  )
) {
  
  if (
    !(
      "17c" %in%
      tolower(
        as.character(
          stage17c_results$stage
        )
      )
    )
  ) {
    
    write_log(
      "Warning: Stage 17C stage label is ",
      stage17c_results$stage
    )
    
  }
  
}


# ==============================================================================
# 9. RECOVER THE ACTUAL STAGE 17B GENE CLUSTER TABLE
# ==============================================================================

temporal_gene_table <-
  as.data.frame(
    stage17b_results$gene_cluster_correlation_table,
    stringsAsFactors = FALSE
  )


required_gene_cluster_columns <-
  c(
    "Gene",
    "Cluster",
    "Centroid_Pearson_Correlation",
    "Silhouette_Width",
    "Rank_Within_Cluster"
  )


missing_gene_cluster_columns <-
  setdiff(
    required_gene_cluster_columns,
    colnames(
      temporal_gene_table
    )
  )


if (
  length(
    missing_gene_cluster_columns
  ) > 0
) {
  
  stop(
    paste(
      "Stage 17B gene cluster table is missing:",
      paste(
        missing_gene_cluster_columns,
        collapse = ", "
      )
    )
  )
  
}


# ==============================================================================
# 10. STANDARDIZE GENE CLUSTER INFORMATION
# ==============================================================================

temporal_gene_table$Gene <-
  as.character(
    temporal_gene_table$Gene
  )


temporal_gene_table$Cluster <-
  as.character(
    temporal_gene_table$Cluster
  )


temporal_gene_table$Centroid_Pearson_Correlation <-
  as.numeric(
    temporal_gene_table$Centroid_Pearson_Correlation
  )


temporal_gene_table$Silhouette_Width <-
  as.numeric(
    temporal_gene_table$Silhouette_Width
  )


temporal_gene_table$Rank_Within_Cluster <-
  as.numeric(
    temporal_gene_table$Rank_Within_Cluster
  )


temporal_gene_table <-
  temporal_gene_table[
    !is.na(
      temporal_gene_table$Gene
    ) &
      temporal_gene_table$Gene != "" &
      !is.na(
        temporal_gene_table$Cluster
      ),
    ,
    drop = FALSE
  ]


temporal_gene_table <-
  temporal_gene_table[
    !duplicated(
      temporal_gene_table$Gene
    ),
    ,
    drop = FALSE
  ]


temporal_gene_table$Temporal_Cluster <-
  temporal_gene_table$Cluster


# ==============================================================================
# 11. BASIC TEMPORAL CLUSTER VALIDATION
# ==============================================================================

if (
  nrow(
    temporal_gene_table
  ) != 569
) {
  
  write_log(
    "Warning: expected 569 temporal genes; observed ",
    nrow(
      temporal_gene_table
    )
  )
  
}


cluster_ids <-
  sort(
    unique(
      temporal_gene_table$Temporal_Cluster
    )
  )


if (
  length(
    cluster_ids
  ) < 2
) {
  
  stop(
    "Fewer than two temporal clusters detected."
  )
  
}


cluster_gene_counts <-
  table(
    temporal_gene_table$Temporal_Cluster
  )


# ==============================================================================
# 12. RECOVER VALIDATED CLUSTER CENTROIDS
# ==============================================================================

cluster_centroids <-
  as.matrix(
    stage17b_results$cluster_centroids
  )


if (
  nrow(
    cluster_centroids
  ) != length(
    cluster_ids
  )
) {
  
  stop(
    "Number of centroid rows does not match number of temporal clusters."
  )
  
}


if (
  ncol(
    cluster_centroids
  ) != stage17b_results$n_bins
) {
  
  stop(
    "Number of centroid columns does not match Stage 17B bin count."
  )
  
}


# ==============================================================================
# 13. STANDARDIZE CENTROID ROW NAMES
# ==============================================================================

if (
  is.null(
    rownames(
      cluster_centroids
    )
  )
) {
  
  rownames(
    cluster_centroids
  ) <-
    cluster_ids
  
}


if (
  !all(
    cluster_ids %in%
    rownames(
      cluster_centroids
    )
  )
) {
  
  rownames(
    cluster_centroids
  ) <-
    as.character(
      seq_len(
        nrow(
          cluster_centroids
        )
      )
    )
  
}


# Reorder centroid rows according to cluster IDs

cluster_centroids <-
  cluster_centroids[
    match(
      cluster_ids,
      rownames(
        cluster_centroids
      )
    ),
    ,
    drop = FALSE
  ]


# ==============================================================================
# 14. RECOVER PSEUDOTIME BIN INFORMATION
# ==============================================================================

bin_breaks <-
  as.numeric(
    stage17b_results$bin_breaks
  )


n_bins <-
  as.integer(
    stage17b_results$n_bins
  )


if (
  length(
    bin_breaks
  ) !=
  n_bins + 1
) {
  
  stop(
    "Stage 17B bin_breaks length is inconsistent with n_bins."
  )
  
}


bin_centers <-
  (
    bin_breaks[
      -length(
        bin_breaks
      )
    ] +
      bin_breaks[
        -1
      ]
  ) / 2


bin_labels <-
  paste0(
    "Bin_",
    seq_len(
      n_bins
    )
  )


colnames(
  cluster_centroids
) <-
  bin_labels


# ==============================================================================
# 15. CREATE CENTROID TABLE
# ==============================================================================

cluster_centroid_table <-
  data.frame(
    Cluster =
      cluster_ids,
    cluster_centroids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


write.csv(
  cluster_centroid_table,
  file.path(
    tables_main_dir,
    "01_Temporal_Cluster_Centroids.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. CREATE LONG-FORM CENTROID TABLE
# ==============================================================================

cluster_centroid_long <-
  cluster_centroid_table %>%
  
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        bin_labels
      ),
    names_to =
      "Bin",
    values_to =
      "Centroid_Z"
  ) %>%
  
  dplyr::mutate(
    Bin_Number =
      as.integer(
        sub(
          "Bin_",
          "",
          Bin
        )
      ),
    Pseudotime_Bin_Center =
      bin_centers[
        Bin_Number
      ]
  ) %>%
  
  dplyr::arrange(
    Cluster,
    Bin_Number
  )


write.csv(
  cluster_centroid_long,
  file.path(
    tables_supp_dir,
    "01_Temporal_Cluster_Centroids_Long.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. DEFINE TEMPORAL BEHAVIOR FROM VALIDATED CENTROIDS
# ==============================================================================

classify_temporal_pattern <- function(
    centroid,
    early_n = 4,
    late_n = 4
) {
  
  n <-
    length(
      centroid
    )
  
  
  if (
    n < 8
  ) {
    
    return(
      "Unclassified"
    )
    
  }
  
  
  early_mean <-
    mean(
      centroid[
        seq_len(
          early_n
        )
      ],
      na.rm = TRUE
    )
  
  
  late_mean <-
    mean(
      centroid[
        (n - late_n + 1):n
      ],
      na.rm = TRUE
    )
  
  
  peak_index <-
    which.max(
      centroid
    )
  
  
  trough_index <-
    which.min(
      centroid
    )
  
  
  range_value <-
    max(
      centroid,
      na.rm = TRUE
    ) -
    min(
      centroid,
      na.rm = TRUE
    )
  
  
  if (
    range_value == 0
  ) {
    
    return(
      "Stable"
    )
    
  }
  
  
  if (
    peak_index <= early_n &&
    late_mean < early_mean
  ) {
    
    return(
      "Early_Peaking_Decreasing"
    )
    
  }
  
  
  if (
    trough_index <= early_n &&
    late_mean > early_mean
  ) {
    
    return(
      "Late_Increasing"
    )
    
  }
  
  
  if (
    late_mean > early_mean
  ) {
    
    return(
      "Late_Increasing"
    )
    
  }
  
  
  if (
    early_mean > late_mean
  ) {
    
    return(
      "Early_Decreasing"
    )
    
  }
  
  
  return(
    "Complex_or_Intermediate"
  )
  
}


cluster_pattern_summary <-
  data.frame(
    
    Temporal_Cluster =
      cluster_ids,
    
    Pattern =
      vapply(
        seq_along(
          cluster_ids
        ),
        function(i) {
          
          classify_temporal_pattern(
            cluster_centroids[
              i,
            ]
          )
          
        },
        character(
          1
        )
      ),
    
    Early_Mean_Z =
      apply(
        cluster_centroids[
          ,
          seq_len(
            min(
              4,
              n_bins
            )
          ),
          drop = FALSE
        ],
        1,
        mean,
        na.rm = TRUE
      ),
    
    Late_Mean_Z =
      apply(
        cluster_centroids[
          ,
          max(
            1,
            n_bins - 3
          ):n_bins,
          drop = FALSE
        ],
        1,
        mean,
        na.rm = TRUE
      ),
    
    Peak_Bin =
      apply(
        cluster_centroids,
        1,
        which.max
      ),
    
    Trough_Bin =
      apply(
        cluster_centroids,
        1,
        which.min
      ),
    
    stringsAsFactors =
      FALSE
    
  )


cluster_pattern_summary$Early_Mean_Z <-
  as.numeric(
    cluster_pattern_summary$Early_Mean_Z
  )


cluster_pattern_summary$Late_Mean_Z <-
  as.numeric(
    cluster_pattern_summary$Late_Mean_Z
  )


write.csv(
  cluster_pattern_summary,
  file.path(
    tables_main_dir,
    "02_Temporal_Pattern_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. CREATE CLUSTER GENE LISTS
# ==============================================================================

cluster_gene_lists <-
  split(
    temporal_gene_table$Gene,
    temporal_gene_table$Temporal_Cluster
  )


cluster_gene_lists <-
  lapply(
    cluster_gene_lists,
    unique
  )


# ==============================================================================
# 19. CREATE TOP REPRESENTATIVE GENE TABLE
# ==============================================================================

top_genes_per_cluster <-
  temporal_gene_table %>%
  
  dplyr::group_by(
    Temporal_Cluster
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Centroid_Pearson_Correlation
    ),
    dplyr::desc(
      Silhouette_Width
    ),
    .by_group = TRUE
  ) %>%
  
  dplyr::slice_head(
    n = 30
  ) %>%
  
  dplyr::ungroup()


write.csv(
  top_genes_per_cluster,
  file.path(
    tables_main_dir,
    "03_Top30_Representative_Genes_Per_Temporal_Cluster.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 20. BIOLOGICAL PROGRAM DESCRIPTIONS
# ==============================================================================

biological_program_description <-
  data.frame(
    
    Temporal_Cluster =
      cluster_ids,
    
    Pattern =
      cluster_pattern_summary$Pattern[
        match(
          cluster_ids,
          cluster_pattern_summary$Temporal_Cluster
        )
      ],
    
    Gene_Count =
      as.integer(
        cluster_gene_counts[
          cluster_ids
        ]
      ),
    
    Biological_Description =
      NA_character_,
    
    stringsAsFactors =
      FALSE
    
  )


if (
  "1" %in%
  cluster_ids
) {
  
  biological_program_description$Biological_Description[
    biological_program_description$Temporal_Cluster == "1"
  ] <-
    "Early-peaking immune/antigen-presentation-associated temporal program"
  
}


if (
  "2" %in%
  cluster_ids
) {
  
  biological_program_description$Biological_Description[
    biological_program_description$Temporal_Cluster == "2"
  ] <-
    "Broad late-emerging epithelial/cellular-state-associated temporal program"
  
}


# Generic fallback for unexpected cluster IDs

missing_descriptions <-
  is.na(
    biological_program_description$Biological_Description
  )


if (
  any(
    missing_descriptions
  )
) {
  
  biological_program_description$Biological_Description[
    missing_descriptions
  ] <-
    "Temporally distinct gene-expression program; biological interpretation based on representative genes and pathway enrichment"
  
}


write.csv(
  biological_program_description,
  file.path(
    tables_main_dir,
    "04_Biological_Program_Descriptions.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 21. DETERMINE STAGE 17 TRANSCRIPTOME-WIDE BACKGROUND
# ==============================================================================

stage17_association_table <-
  NULL


if (
  "association_table" %in%
  names(
    stage17b_results
  )
) {
  
  stage17_association_table <-
    as.data.frame(
      stage17b_results$association_table,
      stringsAsFactors = FALSE
    )
  
}


background_genes <-
  NULL


background_source <-
  NULL


if (
  !is.null(
    stage17_association_table
  ) &&
  "Gene" %in%
  colnames(
    stage17_association_table
  )
) {
  
  background_genes <-
    unique(
      as.character(
        stage17_association_table$Gene
      )
    )
  
  
  background_genes <-
    background_genes[
      !is.na(
        background_genes
      ) &
        background_genes != ""
    ]
  
  
  if (
    length(
      background_genes
    ) > 0
  ) {
    
    background_source <-
      "Stage17_transcriptome_wide_association_table"
    
  }
  
}


# If Stage17 association table is unavailable, use Stage17B expression-tested genes
# stored in the Stage17B object as the next-best source.

if (
  is.null(
    background_genes
  ) ||
  length(
    background_genes
  ) == 0
) {
  
  if (
    "strong_temporal_genes_in_expression" %in%
    names(
      stage17b_results
    ) &&
    "source_seurat_expression_object" %in%
    names(
      stage17b_results
    )
  ) {
    
    source_seurat_path <-
      stage17b_results$source_seurat_expression_object
    
  } else {
    
    source_seurat_path <-
      seurat_object_rds
    
  }
  
  
  if (
    file.exists(
      source_seurat_path
    )
  ) {
    
    seurat_obj <-
      tryCatch(
        
        readRDS(
          source_seurat_path
        ),
        
        error = function(e) {
          NULL
        }
        
      )
    
    
    if (
      !is.null(
        seurat_obj
      ) &&
      "RNA" %in%
      names(
        seurat_obj@assays
      )
    ) {
      
      background_genes <-
        unique(
          rownames(
            seurat_obj[["RNA"]]
          )
        )
      
      
      background_genes <-
        background_genes[
          !is.na(
            background_genes
          ) &
            background_genes != ""
        ]
      
      
      if (
        length(
          background_genes
        ) > 0
      ) {
        
        background_source <-
          "RNA_assay_features_from_Stage17_source_Seurat_object"
        
      }
      
    }
    
  }
  
}


if (
  is.null(
    background_genes
  ) ||
  length(
    background_genes
  ) == 0
) {
  
  stop(
    "Unable to establish an appropriate transcriptome-wide background."
  )
  
}


# ==============================================================================
# 22. MAP TEMPORAL GENES TO ENTREZ IDS
# ==============================================================================

all_cluster_genes <-
  unique(
    temporal_gene_table$Gene
  )


all_conversion <-
  suppressMessages(
    
    tryCatch(
      
      clusterProfiler::bitr(
        all_cluster_genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      ),
      
      error = function(e) {
        NULL
      }
      
    )
    
  )


if (
  is.null(
    all_conversion
  ) ||
  nrow(
    all_conversion
  ) == 0
) {
  
  stop(
    "No temporal genes could be mapped to Entrez IDs."
  )
  
}


all_conversion <-
  all_conversion[
    !duplicated(
      all_conversion$SYMBOL
    ),
    ,
    drop = FALSE
  ]


# ==============================================================================
# 23. MAP BACKGROUND GENES TO ENTREZ IDS
# ==============================================================================

background_conversion <-
  suppressMessages(
    
    tryCatch(
      
      clusterProfiler::bitr(
        background_genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      ),
      
      error = function(e) {
        NULL
      }
      
    )
    
  )


if (
  is.null(
    background_conversion
  ) ||
  nrow(
    background_conversion
  ) == 0
) {
  
  stop(
    "No background genes could be mapped to Entrez IDs."
  )
  
}


background_conversion <-
  background_conversion[
    !duplicated(
      background_conversion$ENTREZID
    ),
    ,
    drop = FALSE
  ]


background_entrez <-
  unique(
    background_conversion$ENTREZID
  )


write_log(
  "Temporal genes: ",
  length(
    all_cluster_genes
  )
)


write_log(
  "Temporal genes mapped to Entrez: ",
  nrow(
    all_conversion
  )
)


write_log(
  "Background genes: ",
  length(
    background_genes
  )
)


write_log(
  "Background genes mapped to Entrez: ",
  length(
    background_entrez
  )
)


write_log(
  "Background source: ",
  background_source
)


# ==============================================================================
# 24. GO BIOLOGICAL PROCESS ENRICHMENT
# ==============================================================================

go_results <-
  list()


go_conversion_summary <-
  data.frame()


for (
  cl in names(
    cluster_gene_lists
  )
) {
  
  genes <-
    unique(
      cluster_gene_lists[[cl]]
    )
  
  
  conversion <-
    suppressMessages(
      
      tryCatch(
        
        clusterProfiler::bitr(
          genes,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Hs.eg.db
        ),
        
        error = function(e) {
          NULL
        }
        
      )
      
    )
  
  
  if (
    is.null(
      conversion
    ) ||
    nrow(
      conversion
    ) == 0
  ) {
    
    go_conversion_summary <-
      dplyr::bind_rows(
        
        go_conversion_summary,
        
        data.frame(
          Temporal_Cluster =
            cl,
          Input_Genes =
            length(
              genes
            ),
          Mapped_Genes =
            0,
          stringsAsFactors =
            FALSE
        )
        
      )
    
    next
    
  }
  
  
  conversion <-
    conversion[
      !duplicated(
        conversion$ENTREZID
      ),
      ,
      drop = FALSE
    ]
  
  
  gene_ids <-
    unique(
      conversion$ENTREZID
    )
  
  
  go_conversion_summary <-
    dplyr::bind_rows(
      
      go_conversion_summary,
      
      data.frame(
        Temporal_Cluster =
          cl,
        Input_Genes =
          length(
            genes
          ),
        Mapped_Genes =
          length(
            gene_ids
          ),
        stringsAsFactors =
          FALSE
      )
      
    )
  
  
  ego <-
    tryCatch(
      
      clusterProfiler::enrichGO(
        gene =
          gene_ids,
        universe =
          background_entrez,
        OrgDb =
          org.Hs.eg.db,
        keyType =
          "ENTREZID",
        ont =
          "BP",
        pAdjustMethod =
          "BH",
        pvalueCutoff =
          0.05,
        qvalueCutoff =
          0.20,
        readable =
          TRUE
      ),
      
      error = function(e) {
        NULL
      }
      
    )
  
  
  if (
    !is.null(
      ego
    )
  ) {
    
    go_df <-
      as.data.frame(
        ego
      )
    
    
    if (
      nrow(
        go_df
      ) > 0
    ) {
      
      go_df$Temporal_Cluster <-
        cl
      
      
      go_results[[cl]] <-
        ego
      
      
      write.csv(
        go_df,
        file.path(
          tables_supp_dir,
          paste0(
            "GO_BP_Cluster_",
            cl,
            ".csv"
          )
        ),
        row.names = FALSE
      )
      
    }
    
  }
  
}


write.csv(
  go_conversion_summary,
  file.path(
    tables_supp_dir,
    "02_GO_Gene_Mapping_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 25. KEGG ENRICHMENT
# ==============================================================================

kegg_results <-
  list()


kegg_conversion_summary <-
  data.frame()


for (
  cl in names(
    cluster_gene_lists
  )
) {
  
  genes <-
    unique(
      cluster_gene_lists[[cl]]
    )
  
  
  conversion <-
    suppressMessages(
      
      tryCatch(
        
        clusterProfiler::bitr(
          genes,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Hs.eg.db
        ),
        
        error = function(e) {
          NULL
        }
        
      )
      
    )
  
  
  if (
    is.null(
      conversion
    ) ||
    nrow(
      conversion
    ) == 0
  ) {
    
    kegg_conversion_summary <-
      dplyr::bind_rows(
        
        kegg_conversion_summary,
        
        data.frame(
          Temporal_Cluster =
            cl,
          Input_Genes =
            length(
              genes
            ),
          Mapped_Genes =
            0,
          stringsAsFactors =
            FALSE
        )
        
      )
    
    next
    
  }
  
  
  conversion <-
    conversion[
      !duplicated(
        conversion$ENTREZID
      ),
      ,
      drop = FALSE
    ]
  
  
  gene_ids <-
    unique(
      conversion$ENTREZID
    )
  
  
  kegg_conversion_summary <-
    dplyr::bind_rows(
      
      kegg_conversion_summary,
      
      data.frame(
        Temporal_Cluster =
          cl,
        Input_Genes =
          length(
            genes
          ),
        Mapped_Genes =
          length(
            gene_ids
          ),
        stringsAsFactors =
          FALSE
      )
      
    )
  
  
  ekegg <-
    tryCatch(
      
      clusterProfiler::enrichKEGG(
        gene =
          gene_ids,
        universe =
          background_entrez,
        organism =
          "hsa",
        pAdjustMethod =
          "BH",
        pvalueCutoff =
          0.05
      ),
      
      error = function(e) {
        NULL
      }
      
    )
  
  
  if (
    !is.null(
      ekegg
    )
  ) {
    
    kegg_df <-
      as.data.frame(
        ekegg
      )
    
    
    if (
      nrow(
        kegg_df
      ) > 0
    ) {
      
      kegg_df$Temporal_Cluster <-
        cl
      
      
      kegg_results[[cl]] <-
        ekegg
      
      
      write.csv(
        kegg_df,
        file.path(
          tables_supp_dir,
          paste0(
            "KEGG_Cluster_",
            cl,
            ".csv"
          )
        ),
        row.names = FALSE
      )
      
    }
    
  }
  
}


write.csv(
  kegg_conversion_summary,
  file.path(
    tables_supp_dir,
    "03_KEGG_Gene_Mapping_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 26. COMBINE GO RESULTS
# ==============================================================================

go_summary <-
  data.frame()


if (
  length(
    go_results
  ) > 0
) {
  
  go_summary <-
    dplyr::bind_rows(
      
      lapply(
        
        names(
          go_results
        ),
        
        function(cl) {
          
          df <-
            as.data.frame(
              go_results[[cl]]
            )
          
          
          if (
            nrow(
              df
            ) == 0
          ) {
            
            return(
              NULL
            )
            
          }
          
          
          df$Temporal_Cluster <-
            cl
          
          
          df
          
        }
        
      )
      
    )
  
}


if (
  nrow(
    go_summary
  ) > 0
) {
  
  write.csv(
    go_summary,
    file.path(
      tables_main_dir,
      "05_GO_BP_All_Temporal_Clusters.csv"
    ),
    row.names = FALSE
  )
  
}


# ==============================================================================
# 27. COMBINE KEGG RESULTS
# ==============================================================================

kegg_summary <-
  data.frame()


if (
  length(
    kegg_results
  ) > 0
) {
  
  kegg_summary <-
    dplyr::bind_rows(
      
      lapply(
        
        names(
          kegg_results
        ),
        
        function(cl) {
          
          df <-
            as.data.frame(
              kegg_results[[cl]]
            )
          
          
          if (
            nrow(
              df
            ) == 0
          ) {
            
            return(
              NULL
            )
            
          }
          
          
          df$Temporal_Cluster <-
            cl
          
          
          df
          
        }
        
      )
      
    )
  
}


if (
  nrow(
    kegg_summary
  ) > 0
) {
  
  write.csv(
    kegg_summary,
    file.path(
      tables_main_dir,
      "06_KEGG_All_Temporal_Clusters.csv"
    ),
    row.names = FALSE
  )
  
}


# ==============================================================================
# 28. INTEGRATED CLUSTER SUMMARY
# ==============================================================================

cluster_summary <-
  temporal_gene_table %>%
  
  dplyr::group_by(
    Temporal_Cluster
  ) %>%
  
  dplyr::summarise(
    
    Gene_Count =
      dplyr::n(),
    
    Mean_Centroid_Correlation =
      mean(
        Centroid_Pearson_Correlation,
        na.rm = TRUE
      ),
    
    Median_Centroid_Correlation =
      median(
        Centroid_Pearson_Correlation,
        na.rm = TRUE
      ),
    
    Mean_Silhouette =
      mean(
        Silhouette_Width,
        na.rm = TRUE
      ),
    
    Median_Silhouette =
      median(
        Silhouette_Width,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  )


integrated_cluster_summary <-
  dplyr::left_join(
    cluster_summary,
    cluster_pattern_summary,
    by =
      "Temporal_Cluster"
  )


integrated_cluster_summary <-
  dplyr::left_join(
    integrated_cluster_summary,
    biological_program_description[
      ,
      c(
        "Temporal_Cluster",
        "Biological_Description"
      ),
      drop = FALSE
    ],
    by =
      "Temporal_Cluster"
  )


write.csv(
  integrated_cluster_summary,
  file.path(
    tables_main_dir,
    "07_Integrated_Temporal_Cluster_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 29. TEMPORAL CLUSTER PROFILE PLOT
# ==============================================================================

cluster_profile_plot_data <-
  cluster_centroid_long


p_profiles <-
  ggplot(
    cluster_profile_plot_data,
    aes(
      x =
        Pseudotime_Bin_Center,
      y =
        Centroid_Z,
      group =
        Cluster,
      color =
        Cluster
    )
  ) +
  
  geom_line(
    linewidth =
      1.2
  ) +
  
  geom_point(
    size =
      2
  ) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    alpha =
      0.5
  ) +
  
  theme_classic() +
  
  labs(
    title =
      "Temporal Dynamics of Gene-Expression Programs",
    subtitle =
      "Stage 17B temporal cluster centroids",
    x =
      "Pseudotime",
    y =
      "Mean gene-wise Z-score",
    color =
      "Temporal Cluster"
  )


ggsave(
  file.path(
    figures_main_dir,
    "Stage17D_Temporal_Cluster_Profiles.pdf"
  ),
  p_profiles,
  width =
    10,
  height =
    7
)


# ==============================================================================
# 30. TEMPORAL CLUSTER CENTROID HEATMAP
# ==============================================================================

centroid_heatmap_df <-
  cluster_centroid_long %>%
  
  dplyr::mutate(
    Cluster =
      factor(
        Cluster,
        levels =
          cluster_ids
      ),
    Bin_Number =
      factor(
        Bin_Number,
        levels =
          seq_len(
            n_bins
          )
      )
  )


p_centroid_heatmap <-
  ggplot(
    centroid_heatmap_df,
    aes(
      x =
        Bin_Number,
      y =
        Cluster,
      fill =
        Centroid_Z
    )
  ) +
  
  geom_tile() +
  
  theme_classic() +
  
  labs(
    title =
      "Temporal Cluster Centroid Heatmap",
    x =
      "Pseudotime bin",
    y =
      "Temporal cluster",
    fill =
      "Mean Z-score"
  )


ggsave(
  file.path(
    figures_supp_dir,
    "Stage17D_Temporal_Cluster_Centroid_Heatmap.pdf"
  ),
  p_centroid_heatmap,
  width =
    10,
  height =
    4
)


# ==============================================================================
# 31. TEMPORAL CLUSTER SIZE PLOT
# ==============================================================================

p_cluster_size <-
  ggplot(
    integrated_cluster_summary,
    aes(
      x =
        Temporal_Cluster,
      y =
        Gene_Count
    )
  ) +
  
  geom_col() +
  
  theme_classic() +
  
  labs(
    title =
      "Temporal Cluster Gene Counts",
    x =
      "Temporal cluster",
    y =
      "Number of genes"
  )


ggsave(
  file.path(
    figures_main_dir,
    "Stage17D_Temporal_Cluster_Gene_Counts.pdf"
  ),
  p_cluster_size,
  width =
    8,
  height =
    6
)


# ==============================================================================
# 32. GO DOTPLOTS
# ==============================================================================

if (
  length(
    go_results
  ) > 0
) {
  
  for (
    cl in names(
      go_results
    )
  ) {
    
    ego <-
      go_results[[cl]]
    
    
    if (
      nrow(
        as.data.frame(
          ego
        )
      ) == 0
    ) {
      
      next
      
    }
    
    
    p_go <-
      tryCatch(
        
        enrichplot::dotplot(
          ego,
          showCategory =
            15
        ) +
          
          ggtitle(
            paste(
              "GO Biological Process - Temporal Cluster",
              cl
            )
          ) +
          
          theme_classic() +
          
          theme(
            plot.title =
              element_text(
                hjust =
                  0.5
              )
          ),
        
        error = function(e) {
          NULL
        }
        
      )
    
    
    if (
      !is.null(
        p_go
      )
    ) {
      
      ggsave(
        file.path(
          figures_main_dir,
          paste0(
            "Stage17D_GO_BP_Dotplot_Cluster_",
            cl,
            ".pdf"
          )
        ),
        p_go,
        width =
          10,
        height =
          7
      )
      
    }
    
  }
  
}


# ==============================================================================
# 33. KEGG DOTPLOTS
# ==============================================================================

if (
  length(
    kegg_results
  ) > 0
) {
  
  for (
    cl in names(
      kegg_results
    )
  ) {
    
    ekegg <-
      kegg_results[[cl]]
    
    
    if (
      nrow(
        as.data.frame(
          ekegg
        )
      ) == 0
    ) {
      
      next
      
    }
    
    
    p_kegg <-
      tryCatch(
        
        enrichplot::dotplot(
          ekegg,
          showCategory =
            15
        ) +
          
          ggtitle(
            paste(
              "KEGG Pathways - Temporal Cluster",
              cl
            )
          ) +
          
          theme_classic() +
          
          theme(
            plot.title =
              element_text(
                hjust =
                  0.5
              )
          ),
        
        error = function(e) {
          NULL
        }
        
      )
    
    
    if (
      !is.null(
        p_kegg
      )
    ) {
      
      ggsave(
        file.path(
          figures_main_dir,
          paste0(
            "Stage17D_KEGG_Dotplot_Cluster_",
            cl,
            ".pdf"
          )
        ),
        p_kegg,
        width =
          10,
        height =
          7
      )
      
    }
    
  }
  
}


# ==============================================================================
# 34. FINAL BIOLOGICAL INTERPRETATION TABLE
# ==============================================================================

biological_interpretation <-
  integrated_cluster_summary


biological_interpretation$Interpretation_Caution <-
  paste(
    "Temporal gene-expression program along inferred pseudotime;",
    "not an independent cell subtype or evidence of cell-type conversion."
  )


write.csv(
  biological_interpretation,
  file.path(
    tables_main_dir,
    "08_Biological_Interpretation_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 35. ENRICHMENT SUMMARY
# ==============================================================================

go_result_counts <-
  data.frame(
    Temporal_Cluster =
      cluster_ids,
    Significant_GO_Terms =
      vapply(
        cluster_ids,
        function(cl) {
          
          if (
            cl %in%
            names(
              go_results
            )
          ) {
            
            nrow(
              as.data.frame(
                go_results[[cl]]
              )
            )
            
          } else {
            
            0
            
          }
          
        },
        numeric(
          1
        )
      ),
    stringsAsFactors =
      FALSE
  )


kegg_result_counts <-
  data.frame(
    Temporal_Cluster =
      cluster_ids,
    Significant_KEGG_Pathways =
      vapply(
        cluster_ids,
        function(cl) {
          
          if (
            cl %in%
            names(
              kegg_results
            )
          ) {
            
            nrow(
              as.data.frame(
                kegg_results[[cl]]
              )
            )
            
          } else {
            
            0
            
          }
          
        },
        numeric(
          1
        )
      ),
    stringsAsFactors =
      FALSE
  )


enrichment_summary <-
  dplyr::left_join(
    go_result_counts,
    kegg_result_counts,
    by =
      "Temporal_Cluster"
  )


write.csv(
  enrichment_summary,
  file.path(
    tables_main_dir,
    "09_Enrichment_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 36. STAGE 17D VALIDATION
# ==============================================================================

validation_checks <-
  data.frame(
    
    Check = c(
      
      "Stage 17B loaded",
      
      "Stage 17C loaded",
      
      "Stage 17B gene cluster table available",
      
      "Gene column available",
      
      "Cluster column available",
      
      "Silhouette values available",
      
      "Temporal genes retained",
      
      "Expected 569 strong temporal genes",
      
      "At least two temporal clusters",
      
      "Cluster centroids available",
      
      "20 pseudotime bins available",
      
      "Cluster centroids have 20 bins",
      
      "Transcriptome-wide background available",
      
      "Temporal genes mapped to Entrez",
      
      "Background genes mapped to Entrez",
      
      "GO analysis executed",
      
      "KEGG analysis executed",
      
      "Main biological interpretation table saved",
      
      "Temporal profile figure saved"
      
    ),
    
    Passed = c(
      
      is.list(
        stage17b_results
      ),
      
      is.list(
        stage17c_results
      ),
      
      "gene_cluster_correlation_table" %in%
        names(
          stage17b_results
        ),
      
      "Gene" %in%
        colnames(
          temporal_gene_table
        ),
      
      "Cluster" %in%
        colnames(
          temporal_gene_table
        ),
      
      "Silhouette_Width" %in%
        colnames(
          temporal_gene_table
        ),
      
      nrow(
        temporal_gene_table
      ) > 0,
      
      nrow(
        temporal_gene_table
      ) == 569,
      
      length(
        cluster_ids
      ) >= 2,
      
      !is.null(
        cluster_centroids
      ),
      
      n_bins == 20,
      
      ncol(
        cluster_centroids
      ) == 20,
      
      length(
        background_genes
      ) > 0,
      
      nrow(
        all_conversion
      ) > 0,
      
      length(
        background_entrez
      ) > 0,
      
      TRUE,
      
      TRUE,
      
      file.exists(
        file.path(
          tables_main_dir,
          "08_Biological_Interpretation_Summary.csv"
        )
      ),
      
      file.exists(
        file.path(
          figures_main_dir,
          "Stage17D_Temporal_Cluster_Profiles.pdf"
        )
      )
      
    ),
    
    stringsAsFactors =
      FALSE
    
  )


# ==============================================================================
# 37. VALIDATION FAILURE CHECK
# ==============================================================================

if (
  !all(
    validation_checks$Passed
  )
) {
  
  failed_checks <-
    validation_checks$Check[
      !validation_checks$Passed
    ]
  
  
  write_log(
    "Validation failures detected:"
  )
  
  
  write_log(
    paste(
      failed_checks,
      collapse = " | "
    )
  )
  
}


write.csv(
  validation_checks,
  file.path(
    validation_dir,
    "Stage17D_Final_Validation.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 38. FINAL STAGE 17D OBJECT
# ==============================================================================

stage17d_complete <-
  list(
    
    stage =
      "17D",
    
    analysis_name =
      "Biological Interpretation of Temporal Gene Programs",
    
    project_directory =
      project_dir,
    
    stage17b_source =
      stage17b_object_rds,
    
    stage17c_source =
      stage17c_object_rds,
    
    seurat_source =
      if (
        "source_seurat_expression_object" %in%
        names(
          stage17b_results
        )
      )
        stage17b_results$source_seurat_expression_object
    else
      seurat_object_rds,
    
    temporal_gene_table =
      temporal_gene_table,
    
    cluster_gene_lists =
      cluster_gene_lists,
    
    cluster_summary =
      cluster_summary,
    
    cluster_pattern_summary =
      cluster_pattern_summary,
    
    biological_program_description =
      biological_program_description,
    
    biological_interpretation =
      biological_interpretation,
    
    top_genes_per_cluster =
      top_genes_per_cluster,
    
    cluster_centroids =
      cluster_centroids,
    
    cluster_centroid_long =
      cluster_centroid_long,
    
    go_results =
      go_results,
    
    go_summary =
      go_summary,
    
    kegg_results =
      kegg_results,
    
    kegg_summary =
      kegg_summary,
    
    go_conversion_summary =
      go_conversion_summary,
    
    kegg_conversion_summary =
      kegg_conversion_summary,
    
    background_gene_count =
      length(
        background_genes
      ),
    
    background_entrez_count =
      length(
        background_entrez
      ),
    
    background_source =
      background_source,
    
    mapped_temporal_gene_count =
      nrow(
        all_conversion
      ),
    
    enrichment_summary =
      enrichment_summary,
    
    validation =
      validation_checks,
    
    upstream_frozen =
      TRUE,
    
    trajectory_reconstruction =
      FALSE,
    
    root_selection =
      FALSE,
    
    pseudotime_recalculation =
      FALSE,
    
    description =
      paste(
        "Biological interpretation of temporally distinct",
        "gene-expression programs identified along frozen Stage 15",
        "pseudotime."
      )
    
  )


# ==============================================================================
# 39. SAVE CENTRAL STAGE 17D RDS
# ==============================================================================

stage17d_object_rds <-
  file.path(
    objects_dir,
    "Stage17D_Biological_Interpretation_Results.rds"
  )


saveRDS(
  stage17d_complete,
  stage17d_object_rds
)


# ==============================================================================
# 40. RELOAD VERIFICATION
# ==============================================================================

stage17d_test <-
  readRDS(
    stage17d_object_rds
  )


reload_success <-
  is.list(
    stage17d_test
  )


if (
  !reload_success
) {
  
  stop(
    "Stage 17D reload verification failed."
  )
  
}


# ==============================================================================
# 41. FINAL VALIDATION UPDATE
# ==============================================================================

validation_checks$Passed[
  validation_checks$Check ==
    "Main biological interpretation table saved"
] <-
  file.exists(
    file.path(
      tables_main_dir,
      "08_Biological_Interpretation_Summary.csv"
    )
  )


validation_checks$Passed[
  validation_checks$Check ==
    "Temporal profile figure saved"
] <-
  file.exists(
    file.path(
      figures_main_dir,
      "Stage17D_Temporal_Cluster_Profiles.pdf"
    )
  )


write.csv(
  validation_checks,
  file.path(
    validation_dir,
    "Stage17D_Final_Validation.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 42. COMPLETION LOG
# ==============================================================================

write_log(
  "============================================================"
)

write_log(
  "STAGE 17D COMPLETED SUCCESSFULLY"
)

write_log(
  "============================================================"
)

write_log(
  "Temporal genes: ",
  nrow(
    temporal_gene_table
  )
)

write_log(
  "Temporal clusters: ",
  length(
    cluster_ids
  )
)

write_log(
  "Cluster 1 genes: ",
  sum(
    temporal_gene_table$Temporal_Cluster == "1"
  )
)

write_log(
  "Cluster 2 genes: ",
  sum(
    temporal_gene_table$Temporal_Cluster == "2"
  )
)

write_log(
  "Pseudotime bins: ",
  n_bins
)

write_log(
  "Temporal genes mapped to Entrez: ",
  nrow(
    all_conversion
  )
)

write_log(
  "Background genes: ",
  length(
    background_genes
  )
)

write_log(
  "Background genes mapped to Entrez: ",
  length(
    background_entrez
  )
)

write_log(
  "GO clusters with significant results: ",
  length(
    go_results
  )
)

write_log(
  "KEGG clusters with significant results: ",
  length(
    kegg_results
  )
)

write_log(
  "Validation passed: ",
  sum(
    validation_checks$Passed
  ),
  "/",
  nrow(
    validation_checks
  )
)

write_log(
  "Central RDS: ",
  stage17d_object_rds
)

write_log(
  "Completion log: ",
  stage17d_log
)

write_log(
  "Stage 17B remained frozen."
)

write_log(
  "Stage 17C remained frozen."
)

write_log(
  "No trajectory reconstruction, root selection, or pseudotime recalculation was performed."
)

write_log(
  "============================================================"
)


# ==============================================================================
# 43. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "STAGE 17D COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n"
)

cat(
  "Temporal genes: ",
  nrow(
    temporal_gene_table
  ),
  "\n"
)

cat(
  "Temporal clusters: ",
  length(
    cluster_ids
  ),
  "\n"
)

cat(
  "Cluster 1 genes: ",
  sum(
    temporal_gene_table$Temporal_Cluster == "1"
  ),
  "\n"
)

cat(
  "Cluster 2 genes: ",
  sum(
    temporal_gene_table$Temporal_Cluster == "2"
  ),
  "\n"
)

cat(
  "Pseudotime bins: ",
  n_bins,
  "\n"
)

cat(
  "Temporal genes mapped to Entrez: ",
  nrow(
    all_conversion
  ),
  "\n"
)

cat(
  "Background genes: ",
  length(
    background_genes
  ),
  "\n"
)

cat(
  "GO clusters with results: ",
  length(
    go_results
  ),
  "\n"
)

cat(
  "KEGG clusters with results: ",
  length(
    kegg_results
  ),
  "\n"
)

cat(
  "Validation: ",
  sum(
    validation_checks$Passed
  ),
  "/",
  nrow(
    validation_checks
  ),
  " passed\n"
)

cat(
  "\nCentral RDS:\n",
  stage17d_object_rds,
  "\n"
)

cat(
  "\nCompletion log:\n",
  stage17d_log,
  "\n"
)

cat(
  "\nStage 17B remained frozen.\n"
)

cat(
  "Stage 17C remained frozen.\n"
)

cat(
  "No trajectory reconstruction, root selection, or pseudotime recalculation was performed.\n"
)

cat(
  "============================================================\n"
)


# ==============================================================================
# 44. RETURN FINAL OBJECT
# ==============================================================================

invisible(
  stage17d_complete
)

# ============================================================
# STAGE 17D: EXPORT AND ARCHIVING OF BIOLOGICAL INTERPRETATION
# ============================================================
#
# Purpose:
#   Export the already-computed Stage 17D results stored in the
#   central RDS object into the canonical repository structure.
#
# Important:
#   - No enrichment analysis is rerun.
#   - No upstream stage is rerun.
#   - Stage 17B and Stage 17C remain frozen.
#   - Stage 17D results are read-only during this export.
#   - The central Stage 17D RDS is not modified.
#
# ============================================================



objects_dir <- file.path(
  project_dir,
  "objects"
)

results_dir <- file.path(
  project_dir,
  "results",
  "17_Temporal_Gene_Dynamics"
)

tables_main_dir <- file.path(
  results_dir,
  "tables",
  "main"
)

tables_supp_dir <- file.path(
  results_dir,
  "tables",
  "supplementary"
)

tables_validation_dir <- file.path(
  results_dir,
  "tables",
  "validation"
)

logs_dir <- file.path(
  results_dir,
  "logs"
)

stage17d_rds <- file.path(
  objects_dir,
  "Stage17D_Biological_Interpretation_Results.rds"
)

stage17d_log <- file.path(
  logs_dir,
  "Stage17D_Export_Completion.log"
)

# ------------------------------------------------------------
# Create canonical directories
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Logging
# ------------------------------------------------------------

writeLines(
  character(0),
  con = stage17d_log
)

write_log <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(txt, "\n")
  write(
    txt,
    file = stage17d_log,
    append = TRUE
  )
}

write_log("============================================================")
write_log("STAGE 17D: EXPORT AND ARCHIVING")
write_log("============================================================")
write_log("Project directory: ", project_dir)

# ------------------------------------------------------------
# Verify source RDS
# ------------------------------------------------------------

if (!file.exists(stage17d_rds)) {
  stop(
    paste0(
      "Stage17D RDS not found: ",
      stage17d_rds
    )
  )
}

write_log("Stage17D RDS found.")

# ------------------------------------------------------------
# Load existing Stage17D object
# ------------------------------------------------------------

stage17d <- readRDS(stage17d_rds)

write_log("Stage17D RDS loaded successfully.")

# ------------------------------------------------------------
# Verify required fields
# ------------------------------------------------------------

required_fields <- c(
  "temporal_gene_table",
  "cluster_gene_lists",
  "cluster_summary",
  "cluster_pattern_summary",
  "biological_program_description",
  "biological_interpretation",
  "top_genes_per_cluster",
  "go_results",
  "go_summary",
  "kegg_results",
  "kegg_summary",
  "go_conversion_summary",
  "kegg_conversion_summary",
  "enrichment_summary",
  "validation"
)

missing_fields <- setdiff(
  required_fields,
  names(stage17d)
)

if (length(missing_fields) > 0) {
  stop(
    paste0(
      "Required Stage17D fields are missing: ",
      paste(missing_fields, collapse = ", ")
    )
  )
}

write_log("Required Stage17D fields verified.")

# ------------------------------------------------------------
# Basic structural validation before export
# ------------------------------------------------------------

if (nrow(stage17d$temporal_gene_table) != 569) {
  stop("Unexpected number of temporal genes.")
}

if (length(stage17d$cluster_gene_lists) != 2) {
  stop("Unexpected number of temporal clusters.")
}

if (nrow(stage17d$go_summary) == 0) {
  stop("GO summary is unexpectedly empty.")
}

if (nrow(stage17d$validation) == 0) {
  stop("Validation table is unexpectedly empty.")
}

write_log(
  "Temporal genes: ",
  nrow(stage17d$temporal_gene_table)
)

write_log(
  "Temporal clusters: ",
  length(stage17d$cluster_gene_lists)
)

write_log(
  "GO enrichment rows: ",
  nrow(stage17d$go_summary)
)

write_log(
  "KEGG enrichment rows: ",
  nrow(stage17d$kegg_summary)
)

# ------------------------------------------------------------
# MAIN TABLES
# ------------------------------------------------------------

write.csv(
  as.data.frame(stage17d$biological_interpretation),
  file.path(
    tables_main_dir,
    "Stage17D_Biological_Interpretation.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$biological_program_description),
  file.path(
    tables_main_dir,
    "Stage17D_Biological_Program_Description.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$cluster_summary),
  file.path(
    tables_main_dir,
    "Stage17D_Cluster_Summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# GO RESULTS BY CLUSTER
# ------------------------------------------------------------

go_cluster_1 <- stage17d$go_results[["1"]]
go_cluster_2 <- stage17d$go_results[["2"]]

if (!inherits(go_cluster_1, "enrichResult")) {
  stop("GO result for Cluster 1 is not an enrichResult object.")
}

if (!inherits(go_cluster_2, "enrichResult")) {
  stop("GO result for Cluster 2 is not an enrichResult object.")
}

go_cluster_1_df <- as.data.frame(go_cluster_1)
go_cluster_2_df <- as.data.frame(go_cluster_2)

if (nrow(go_cluster_1_df) > 0) {
  go_cluster_1_df$Temporal_Cluster <- "1"
  
  write.csv(
    go_cluster_1_df,
    file.path(
      tables_main_dir,
      "Stage17D_GO_BP_Cluster_1.csv"
    ),
    row.names = FALSE
  )
}

if (nrow(go_cluster_2_df) > 0) {
  go_cluster_2_df$Temporal_Cluster <- "2"
  
  write.csv(
    go_cluster_2_df,
    file.path(
      tables_main_dir,
      "Stage17D_GO_BP_Cluster_2.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# SUPPLEMENTARY TABLES
# ------------------------------------------------------------

write.csv(
  as.data.frame(stage17d$temporal_gene_table),
  file.path(
    tables_supp_dir,
    "Stage17D_Temporal_Gene_Table.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$top_genes_per_cluster),
  file.path(
    tables_supp_dir,
    "Stage17D_Top_Genes_Per_Cluster.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$cluster_pattern_summary),
  file.path(
    tables_supp_dir,
    "Stage17D_Cluster_Pattern_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$go_summary),
  file.path(
    tables_supp_dir,
    "Stage17D_GO_BP_All_Clusters.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$go_conversion_summary),
  file.path(
    tables_supp_dir,
    "Stage17D_GO_Conversion_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$kegg_conversion_summary),
  file.path(
    tables_supp_dir,
    "Stage17D_KEGG_Conversion_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(stage17d$enrichment_summary),
  file.path(
    tables_supp_dir,
    "Stage17D_Enrichment_Summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# VALIDATION TABLE
# ------------------------------------------------------------

write.csv(
  as.data.frame(stage17d$validation),
  file.path(
    tables_validation_dir,
    "Stage17D_Final_Validation.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# KEGG STATUS
# ------------------------------------------------------------

kegg_rows <- nrow(stage17d$kegg_summary)

if (kegg_rows == 0) {
  
  write_log(
    "KEGG enrichment results: no significant pathways."
  )
  
} else {
  
  write.csv(
    as.data.frame(stage17d$kegg_summary),
    file.path(
      tables_supp_dir,
      "Stage17D_KEGG_All_Clusters.csv"
    ),
    row.names = FALSE
  )
  
  write_log(
    "KEGG enrichment rows exported: ",
    kegg_rows
  )
}

# ------------------------------------------------------------
# Verify exported files
# ------------------------------------------------------------

expected_files <- c(
  file.path(
    tables_main_dir,
    "Stage17D_Biological_Interpretation.csv"
  ),
  file.path(
    tables_main_dir,
    "Stage17D_Biological_Program_Description.csv"
  ),
  file.path(
    tables_main_dir,
    "Stage17D_Cluster_Summary.csv"
  ),
  file.path(
    tables_main_dir,
    "Stage17D_GO_BP_Cluster_1.csv"
  ),
  file.path(
    tables_main_dir,
    "Stage17D_GO_BP_Cluster_2.csv"
  ),
  file.path(
    tables_supp_dir,
    "Stage17D_Temporal_Gene_Table.csv"
  ),
  file.path(
    tables_supp_dir,
    "Stage17D_Top_Genes_Per_Cluster.csv"
  ),
  file.path(
    tables_supp_dir,
    "Stage17D_Cluster_Pattern_Summary.csv"
  ),
  file.path(
    tables_supp_dir,
    "Stage17D_GO_BP_All_Clusters.csv"
  ),
  file.path(
    tables_supp_dir,
    "Stage17D_GO_Conversion_Summary.csv"
  ),
  file.path(
    tables_supp_dir,
    "Stage17D_KEGG_Conversion_Summary.csv"
  ),
  file.path(
    tables_supp_dir,
    "Stage17D_Enrichment_Summary.csv"
  ),
  file.path(
    tables_validation_dir,
    "Stage17D_Final_Validation.csv"
  )
)

export_check <- file.exists(expected_files)

if (!all(export_check)) {
  
  missing_exports <- expected_files[!export_check]
  
  stop(
    paste0(
      "One or more expected export files are missing:\n",
      paste(missing_exports, collapse = "\n")
    )
  )
}

# ------------------------------------------------------------
# Final verification
# ------------------------------------------------------------

exported_main <- list.files(
  tables_main_dir,
  pattern = "^Stage17D_.*\\.csv$",
  full.names = FALSE
)

exported_supplementary <- list.files(
  tables_supp_dir,
  pattern = "^Stage17D_.*\\.csv$",
  full.names = FALSE
)

exported_validation <- list.files(
  tables_validation_dir,
  pattern = "^Stage17D_.*\\.csv$",
  full.names = FALSE
)

validation_passed <- all(
  stage17d$validation$Passed
)

write_log("============================================================")
write_log("STAGE 17D EXPORT COMPLETED SUCCESSFULLY")
write_log("============================================================")
write_log(
  "Temporal genes: ",
  nrow(stage17d$temporal_gene_table)
)
write_log(
  "GO rows exported: ",
  nrow(stage17d$go_summary)
)
write_log(
  "GO Cluster 1 rows: ",
  nrow(go_cluster_1_df)
)
write_log(
  "GO Cluster 2 rows: ",
  nrow(go_cluster_2_df)
)
write_log(
  "KEGG rows: ",
  nrow(stage17d$kegg_summary)
)
write_log(
  "Validation: ",
  sum(stage17d$validation$Passed),
  "/",
  nrow(stage17d$validation)
)
write_log(
  "Main tables exported: ",
  length(exported_main)
)
write_log(
  "Supplementary tables exported: ",
  length(exported_supplementary)
)
write_log(
  "Validation tables exported: ",
  length(exported_validation)
)
write_log(
  "Central RDS preserved: ",
  stage17d_rds
)
write_log(
  "No enrichment analysis was rerun."
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
  "Stage 17B remained frozen."
)
write_log(
  "Stage 17C remained frozen."
)
write_log("============================================================")

cat(
  "\nSTAGE 17D EXPORT COMPLETED SUCCESSFULLY\n",
  "GO rows exported: ", nrow(stage17d$go_summary), "\n",
  "GO Cluster 1: ", nrow(go_cluster_1_df), "\n",
  "GO Cluster 2: ", nrow(go_cluster_2_df), "\n",
  "KEGG rows: ", nrow(stage17d$kegg_summary), "\n",
  "Validation: ",
  sum(stage17d$validation$Passed),
  "/",
  nrow(stage17d$validation), "\n",
  "\nMain tables: ",
  tables_main_dir,
  "\nSupplementary tables: ",
  tables_supp_dir,
  "\nValidation tables: ",
  tables_validation_dir,
  "\nExport log: ",
  stage17d_log,
  "\n"
)

