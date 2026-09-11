############################################################
# TNBC Single-cell RNA-seq Analysis Pipeline
#
# STEP 12: Integrated GSEA Pathway Prioritization
#          Tumor-like Epithelial Cells
#
# Comparison:
# BRCA1_tumour vs TotalCell
#
# IMPORTANT:
# - Step 12 is independent from the Step 11 environment
# - Significant GSEA results are loaded from Step 11 CSV files
# - GO BP, Reactome, KEGG, and MSigDB Hallmark are integrated
# - No new differential expression analysis is performed
# - No new GSEA is performed
# - No additional p-value filtering is applied
# - Step 12 focuses on integration, prioritization, and
#   cross-database recurrence
# - Existing adjusted p-values from Step 11 are retained
#
# SIGNIFICANCE:
# - Step 11 significance criterion:
#   p.adjust < 0.05
#
# BIOLOGICAL DIRECTION:
# - NES > 0 = enrichment toward BRCA1_tumour
# - NES < 0 = enrichment toward TotalCell
#
# DATABASES:
# - GO Biological Process
# - Reactome
# - KEGG
# - MSigDB Hallmark
#
# OUTPUT:
# - Integrated significant GSEA table
# - Direction-specific integrated tables
# - Ranked representative terms
# - Exact recurring terms across databases
# - Integrated summary plots
# - Step 12 summary and log
############################################################


############################################################
# 12.1 Load libraries
############################################################

library(dplyr)
library(ggplot2)
library(stringr)
library(readr)


############################################################
# 12.2 Project and Stage 12 directories
############################################################

project_dir <- "YOUR_PROJECT_DIRECTORY"

stage12_dir <- file.path(
  project_dir,
  "results",
  "12_GSEA_Integration"
)

figures_main_dir <- file.path(
  stage12_dir,
  "figures",
  "main"
)

figures_supp_dir <- file.path(
  stage12_dir,
  "figures",
  "supplementary"
)

tables_main_dir <- file.path(
  stage12_dir,
  "tables",
  "main"
)

tables_supp_dir <- file.path(
  stage12_dir,
  "tables",
  "supplementary"
)

tables_validation_dir <- file.path(
  stage12_dir,
  "tables",
  "validation"
)

logs_dir <- file.path(
  stage12_dir,
  "logs"
)


############################################################
# Create Stage 12 directories
############################################################

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


############################################################
# 12.3 Define Step 11 input files
############################################################

stage11_tables_dir <- file.path(
  project_dir,
  "results",
  "11_GSEA",
  "tables",
  "supplementary"
)

input_files <- list(
  
  GO_BP = file.path(
    stage11_tables_dir,
    "GSEA_GO_BP_significant_padj_0.05.csv"
  ),
  
  Reactome = file.path(
    stage11_tables_dir,
    "GSEA_Reactome_significant_padj_0.05.csv"
  ),
  
  KEGG = file.path(
    stage11_tables_dir,
    "GSEA_KEGG_significant_padj_0.05.csv"
  ),
  
  Hallmark = file.path(
    stage11_tables_dir,
    "GSEA_Hallmark_significant_padj_0.05.csv"
  )
)


############################################################
# 12.4 Verify input files
############################################################

input_files_vec <- unlist(
  input_files,
  use.names = TRUE
)

missing_files <- input_files_vec[
  !file.exists(input_files_vec)
]

if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "The following required Step 11 input files were not found:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
  
}

message(
  "All required Step 11 input files were found."
)


############################################################
# 12.5 Read GSEA results
############################################################

gsea_go <- read.csv(
  input_files$GO_BP,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gsea_reactome <- read.csv(
  input_files$Reactome,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gsea_kegg <- read.csv(
  input_files$KEGG,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gsea_hallmark <- read.csv(
  input_files$Hallmark,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 12.6 Add database labels
############################################################

gsea_go <- gsea_go %>%
  mutate(
    Database = "GO_BP"
  )

gsea_reactome <- gsea_reactome %>%
  mutate(
    Database = "Reactome"
  )

gsea_kegg <- gsea_kegg %>%
  mutate(
    Database = "KEGG"
  )

gsea_hallmark <- gsea_hallmark %>%
  mutate(
    Database = "Hallmark"
  )


############################################################
# 12.7 Verify required columns
############################################################

required_columns <- c(
  "Description",
  "NES",
  "p.adjust",
  "Enrichment_Group"
)


check_columns <- function(
    df,
    database_name
) {
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        "Missing required columns in ",
        database_name,
        ": ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
    
  }
  
}


check_columns(
  gsea_go,
  "GO Biological Process"
)

check_columns(
  gsea_reactome,
  "Reactome"
)

check_columns(
  gsea_kegg,
  "KEGG"
)

check_columns(
  gsea_hallmark,
  "Hallmark"
)


############################################################
# 12.8 Standardize GSEA tables
############################################################

standardize_gsea <- function(
    df
) {
  
  df %>%
    
    mutate(
      
      Description = as.character(
        Description
      ),
      
      NES = as.numeric(
        NES
      ),
      
      p.adjust = as.numeric(
        p.adjust
      ),
      
      Enrichment_Group = as.character(
        Enrichment_Group
      )
      
    ) %>%
    
    filter(
      
      !is.na(Description),
      
      Description != "",
      
      !is.na(NES),
      
      !is.na(p.adjust)
      
    )
  
}


gsea_go <- standardize_gsea(
  gsea_go
)

gsea_reactome <- standardize_gsea(
  gsea_reactome
)

gsea_kegg <- standardize_gsea(
  gsea_kegg
)

gsea_hallmark <- standardize_gsea(
  gsea_hallmark
)


############################################################
# 12.9 Combine all GSEA results
############################################################

gsea_integrated <- bind_rows(
  gsea_go,
  gsea_reactome,
  gsea_kegg,
  gsea_hallmark
)


############################################################
# 12.10 Validate integrated dataset
############################################################

if (
  nrow(gsea_integrated) == 0
) {
  
  stop(
    "No valid GSEA records remained after standardization."
  )
  
}

expected_databases <- c(
  "GO_BP",
  "Reactome",
  "KEGG",
  "Hallmark"
)

observed_databases <- sort(
  unique(
    gsea_integrated$Database
  )
)

missing_databases <- setdiff(
  expected_databases,
  observed_databases
)

if (
  length(missing_databases) > 0
) {
  
  stop(
    paste0(
      "The following expected databases are missing from the integrated GSEA data: ",
      paste(
        missing_databases,
        collapse = ", "
      )
    )
  )
  
}


############################################################
# 12.11 Add significance and ranking metrics
############################################################

gsea_integrated <- gsea_integrated %>%
  
  mutate(
    
    NegLog10_Padj =
      -log10(
        p.adjust
      ),
    
    Abs_NES =
      abs(
        NES
      ),
    
    Direction = case_when(
      
      NES > 0 ~
        "BRCA1_tumour",
      
      NES < 0 ~
        "TotalCell",
      
      TRUE ~
        "No_direction"
      
    )
    
  )


############################################################
# 12.12 Verify biological direction
############################################################

direction_check <- gsea_integrated %>%
  
  filter(
    Direction != "No_direction"
  ) %>%
  
  mutate(
    
    Direction_Match =
      Direction ==
      Enrichment_Group
    
  )


n_direction_mismatches <- sum(
  !direction_check$Direction_Match,
  na.rm = TRUE
)

if (
  n_direction_mismatches > 0
) {
  
  warning(
    paste0(
      "Some GSEA entries have inconsistent direction labels. ",
      "Number of mismatches: ",
      n_direction_mismatches
    )
  )
  
} else {
  
  message(
    "Directional NES validation passed."
  )
  
}


############################################################
# 12.13 Define output file paths
############################################################

integrated_file <- file.path(
  tables_main_dir,
  "GSEA_Integrated_significant.csv"
)

brca1_file <- file.path(
  tables_supp_dir,
  "GSEA_Integrated_BRCA1_tumour.csv"
)

totalcell_file <- file.path(
  tables_supp_dir,
  "GSEA_Integrated_TotalCell.csv"
)

top10_file <- file.path(
  tables_supp_dir,
  "GSEA_Integrated_top10_per_database_direction.csv"
)

recurring_file <- file.path(
  tables_supp_dir,
  "GSEA_Integrated_recurring_terms.csv"
)

database_summary_file <- file.path(
  tables_supp_dir,
  "GSEA_Integrated_database_summary.csv"
)

integration_summary_file <- file.path(
  tables_validation_dir,
  "Step12_GSEA_integration_summary.csv"
)

database_validation_file <- file.path(
  tables_validation_dir,
  "Step12_Database_Validation.csv"
)

top10_plot_file <- file.path(
  figures_main_dir,
  "GSEA_Integrated_top10_dotplot.pdf"
)

recurrence_plot_file <- file.path(
  figures_main_dir,
  "GSEA_Integrated_recurrence_plot.pdf"
)

log_file <- file.path(
  logs_dir,
  "Step12_GSEA_integration.log"
)


############################################################
# 12.14 Save integrated significant results
############################################################

write.csv(
  gsea_integrated,
  integrated_file,
  row.names = FALSE
)


############################################################
# 12.15 Separate enrichment directions
############################################################

gsea_integrated_brca1 <- gsea_integrated %>%
  
  filter(
    NES > 0
  ) %>%
  
  arrange(
    p.adjust,
    desc(Abs_NES)
  )


gsea_integrated_totalcell <- gsea_integrated %>%
  
  filter(
    NES < 0
  ) %>%
  
  arrange(
    p.adjust,
    desc(Abs_NES)
  )


############################################################
# 12.16 Save directional integrated results
############################################################

write.csv(
  gsea_integrated_brca1,
  brca1_file,
  row.names = FALSE
)

write.csv(
  gsea_integrated_totalcell,
  totalcell_file,
  row.names = FALSE
)


############################################################
# 12.17 Rank representative GSEA terms
############################################################

gsea_ranked <- gsea_integrated %>%
  
  arrange(
    p.adjust,
    desc(Abs_NES)
  )


gsea_ranked_brca1 <- gsea_integrated_brca1 %>%
  
  arrange(
    p.adjust,
    desc(Abs_NES)
  )


gsea_ranked_totalcell <- gsea_integrated_totalcell %>%
  
  arrange(
    p.adjust,
    desc(Abs_NES)
  )


############################################################
# 12.18 Select top representative terms
############################################################
#
# Maximum 10 terms per database and direction.
#
# This is for visualization/summary only.
# It does NOT alter the complete integrated results.
############################################################

top_representative <- gsea_integrated %>%
  
  group_by(
    Database,
    Direction
  ) %>%
  
  arrange(
    p.adjust,
    desc(Abs_NES)
  ) %>%
  
  slice_head(
    n = 10
  ) %>%
  
  ungroup()


############################################################
# 12.19 Save representative terms
############################################################

write.csv(
  top_representative,
  top10_file,
  row.names = FALSE
)


############################################################
# 12.20 Identify exact recurring descriptions
############################################################
#
# Exact recurrence is assessed using pathway/gene-set
# descriptions.
#
# This is intentionally conservative:
# biologically related but differently named terms are NOT
# automatically considered identical.
############################################################

recurring_terms <- gsea_integrated %>%
  
  group_by(
    Description,
    Direction
  ) %>%
  
  summarise(
    
    Database_Count =
      n_distinct(
        Database
      ),
    
    Databases =
      paste(
        sort(
          unique(
            Database
          )
        ),
        collapse = "; "
      ),
    
    Best_Padj =
      min(
        p.adjust,
        na.rm = TRUE
      ),
    
    Max_Abs_NES =
      max(
        Abs_NES,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  ) %>%
  
  filter(
    Database_Count >= 2
  ) %>%
  
  arrange(
    Direction,
    desc(Database_Count),
    Best_Padj
  )


############################################################
# 12.21 Save recurring terms
############################################################

write.csv(
  recurring_terms,
  recurring_file,
  row.names = FALSE
)


############################################################
# 12.22 Create pathway count summary
############################################################

database_summary <- gsea_integrated %>%
  
  group_by(
    Database,
    Direction
  ) %>%
  
  summarise(
    
    Significant_Terms =
      n(),
    
    Best_Padj =
      min(
        p.adjust,
        na.rm = TRUE
      ),
    
    Maximum_Abs_NES =
      max(
        Abs_NES,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )


############################################################
# 12.23 Save database summary
############################################################

write.csv(
  database_summary,
  database_summary_file,
  row.names = FALSE
)


############################################################
# 12.24 Integrated dot plot
############################################################
#
# Top 10 representative terms per database and direction.
############################################################

if (
  nrow(top_representative) > 0
) {
  
  plot_data <- top_representative %>%
    
    mutate(
      
      Description_wrapped =
        str_wrap(
          Description,
          width = 55
        )
      
    ) %>%
    
    mutate(
      
      Description_wrapped =
        reorder(
          Description_wrapped,
          NES
        )
      
    )
  
  
  p_integrated <- ggplot(
    plot_data,
    aes(
      x = NES,
      y = Description_wrapped,
      size = Abs_NES,
      color = NegLog10_Padj
    )
  ) +
    
    geom_point(
      alpha = 0.85
    ) +
    
    facet_grid(
      Direction ~ Database,
      scales = "free_y",
      space = "free_y"
    ) +
    
    scale_color_viridis_c(
      name = "-log10 adjusted p-value"
    ) +
    
    theme_classic() +
    
    theme(
      
      axis.text.y =
        element_text(
          size = 8
        ),
      
      strip.text =
        element_text(
          face = "bold"
        ),
      
      panel.spacing =
        grid::unit(
          1,
          "lines"
        )
      
    ) +
    
    labs(
      
      title =
        "Integrated GSEA — Representative Significant Terms",
      
      subtitle =
        "Top 10 terms per database and enrichment direction",
      
      x =
        "Normalized Enrichment Score (NES)",
      
      y =
        NULL,
      
      size =
        "|NES|"
      
    )
  
  
  ggsave(
    top10_plot_file,
    p_integrated,
    width = 16,
    height = 14
  )
  
} else {
  
  warning(
    "No representative GSEA terms were available for plotting."
  )
  
}


############################################################
# 12.25 Cross-database recurrence plot
############################################################
#
# IMPORTANT:
# This plot is conditional.
#
# If no exact recurring terms exist across >=2 databases,
# the plot is intentionally NOT generated.
############################################################

n_recurring_terms <- nrow(
  recurring_terms
)


if (
  n_recurring_terms > 0
) {
  
  recurrence_plot_data <-
    recurring_terms %>%
    
    mutate(
      
      Description_wrapped =
        str_wrap(
          Description,
          width = 60
        )
      
    )
  
  
  p_recurrence <- ggplot(
    recurrence_plot_data,
    aes(
      x = Database_Count,
      y = reorder(
        Description_wrapped,
        Database_Count
      ),
      size = Max_Abs_NES
    )
  ) +
    
    geom_point(
      alpha = 0.85
    ) +
    
    facet_wrap(
      ~ Direction,
      scales = "free_y"
    ) +
    
    scale_x_continuous(
      breaks = 2:4
    ) +
    
    theme_classic() +
    
    labs(
      
      title =
        "Cross-Database GSEA Recurrence",
      
      subtitle =
        "Terms recurring across multiple enrichment databases",
      
      x =
        "Number of databases",
      
      y =
        NULL,
      
      size =
        "Maximum |NES|"
      
    )
  
  
  ggsave(
    recurrence_plot_file,
    p_recurrence,
    width = 12,
    height = 10
  )
  
  message(
    "Cross-database recurrence plot generated."
  )
  
} else {
  
  message(
    "No exact recurring GSEA terms were identified. ",
    "Cross-database recurrence plot was intentionally skipped."
  )
  
}


############################################################
# 12.26 Summary statistics
############################################################

n_integrated_total <-
  nrow(
    gsea_integrated
  )

n_integrated_brca1 <-
  nrow(
    gsea_integrated_brca1
  )

n_integrated_totalcell <-
  nrow(
    gsea_integrated_totalcell
  )

n_recurrent_brca1 <-
  sum(
    recurring_terms$Direction ==
      "BRCA1_tumour",
    na.rm = TRUE
  )

n_recurrent_totalcell <-
  sum(
    recurring_terms$Direction ==
      "TotalCell",
    na.rm = TRUE
  )


############################################################
# 12.27 Save Step 12 summary
############################################################

step12_summary <- data.frame(
  
  Metric = c(
    
    "Integrated significant GSEA terms",
    
    "Significant terms enriched in BRCA1_tumour",
    
    "Significant terms enriched in TotalCell",
    
    "Exact recurring terms across >=2 databases",
    
    "Recurring terms in BRCA1_tumour",
    
    "Recurring terms in TotalCell"
    
  ),
  
  Value = c(
    
    n_integrated_total,
    
    n_integrated_brca1,
    
    n_integrated_totalcell,
    
    n_recurring_terms,
    
    n_recurrent_brca1,
    
    n_recurrent_totalcell
    
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  step12_summary,
  integration_summary_file,
  row.names = FALSE
)


############################################################
# 12.28 Database-level validation table
############################################################

database_validation <- gsea_integrated %>%
  
  group_by(
    Database
  ) %>%
  
  summarise(
    
    Total_Terms =
      n(),
    
    BRCA1_tumour_Terms =
      sum(
        NES > 0,
        na.rm = TRUE
      ),
    
    TotalCell_Terms =
      sum(
        NES < 0,
        na.rm = TRUE
      ),
    
    No_Direction_Terms =
      sum(
        NES == 0,
        na.rm = TRUE
      ),
    
    Minimum_Padj =
      min(
        p.adjust,
        na.rm = TRUE
      ),
    
    Maximum_Abs_NES =
      max(
        Abs_NES,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )


write.csv(
  database_validation,
  database_validation_file,
  row.names = FALSE
)


############################################################
# 12.29 Final output validation
############################################################
#
# Mandatory outputs:
# - integrated significant table
# - directional tables
# - representative terms
# - recurring terms table
# - database summary
# - validation tables
# - integrated dot plot
#
# Conditional output:
# - recurrence plot is required ONLY when recurring terms exist.
############################################################

mandatory_outputs <- c(
  
  integrated_file,
  
  brca1_file,
  
  totalcell_file,
  
  top10_file,
  
  recurring_file,
  
  database_summary_file,
  
  integration_summary_file,
  
  database_validation_file,
  
  top10_plot_file
  
)


if (
  n_recurring_terms > 0
) {
  
  mandatory_outputs <- c(
    mandatory_outputs,
    recurrence_plot_file
  )
  
}


missing_outputs <- mandatory_outputs[
  !file.exists(
    mandatory_outputs
  )
]


if (
  length(missing_outputs) > 0
) {
  
  stop(
    paste0(
      "Stage 12 finished with missing expected outputs:\n",
      paste(
        missing_outputs,
        collapse = "\n"
      )
    )
  )
  
}


############################################################
# 12.30 Additional validation
############################################################

if (
  length(
    unique(
      gsea_integrated$Database
    )
  ) != 4
) {
  
  stop(
    "Stage 12 validation failed: expected four GSEA databases."
  )
  
}


if (
  any(
    !is.finite(
      gsea_integrated$NES
    )
  )
) {
  
  stop(
    "Stage 12 validation failed: non-finite NES values detected."
  )
  
}


if (
  any(
    gsea_integrated$p.adjust < 0
  )
) {
  
  stop(
    "Stage 12 validation failed: invalid adjusted p-values detected."
  )
  
}


############################################################
# 12.31 Write final Step 12 log
############################################################

log_lines <- c(
  
  "====================================================",
  
  "TNBC Single-cell RNA-seq Analysis Pipeline",
  
  "STEP 12 — Integrated GSEA Pathway Prioritization",
  
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
  
  "INPUT DATABASES",
  
  "----------------------------------------------------",
  
  "GO Biological Process",
  
  "Reactome",
  
  "KEGG",
  
  "MSigDB Hallmark",
  
  "",
  
  "SIGNIFICANCE",
  
  "----------------------------------------------------",
  
  "Step 11 adjusted p-value criterion: p.adjust < 0.05",
  
  "No additional significance filtering was applied.",
  
  "No new differential expression analysis was performed.",
  
  "No new GSEA was performed.",
  
  "",
  
  "INTEGRATED RESULTS",
  
  "----------------------------------------------------",
  
  paste(
    "Total significant terms:",
    n_integrated_total
  ),
  
  paste(
    "BRCA1_tumour-associated terms:",
    n_integrated_brca1
  ),
  
  paste(
    "TotalCell-associated terms:",
    n_integrated_totalcell
  ),
  
  "",
  
  "CROSS-DATABASE RECURRENCE",
  
  "----------------------------------------------------",
  
  paste(
    "Exact recurring terms across >=2 databases:",
    n_recurring_terms
  ),
  
  paste(
    "Recurring terms in BRCA1_tumour:",
    n_recurrent_brca1
  ),
  
  paste(
    "Recurring terms in TotalCell:",
    n_recurrent_totalcell
  ),
  
  "",
  
  "INTERPRETATION",
  
  "----------------------------------------------------",
  
  "NES > 0 = pathway enriched toward BRCA1_tumour",
  
  "NES < 0 = pathway enriched toward TotalCell",
  
  "p.adjust < 0.05 = statistically significant term",
  
  "",
  
  "QUALITY CONTROL",
  
  "----------------------------------------------------",
  
  paste(
    "Expected databases:",
    paste(
      expected_databases,
      collapse = ", "
    )
  ),
  
  paste(
    "Observed databases:",
    paste(
      observed_databases,
      collapse = ", "
    )
  ),
  
  paste(
    "Directional NES mismatches:",
    n_direction_mismatches
  ),
  
  paste(
    "Integrated dataset rows:",
    n_integrated_total
  ),
  
  paste(
    "Output validation:",
    "PASSED"
  ),
  
  paste(
    "Additional validation:",
    "PASSED"
  ),
  
  "",
  
  if (
    n_recurring_terms > 0
  ) {
    
    "Recurrence plot: GENERATED"
    
  } else {
    
    "Recurrence plot: NOT GENERATED — no recurring terms identified"
    
  },
  
  "",
  
  "STEP 12 COMPLETED SUCCESSFULLY.",
  
  "===================================================="
  
)


writeLines(
  log_lines,
  con = log_file
)


############################################################
# 12.32 Console summary
############################################################

message(
  "===================================================="
)

message(
  "STEP 12 COMPLETED SUCCESSFULLY"
)

message(
  "===================================================="
)

message(
  "Integrated significant GSEA terms: ",
  n_integrated_total
)

message(
  "  BRCA1_tumour: ",
  n_integrated_brca1
)

message(
  "  TotalCell: ",
  n_integrated_totalcell
)

message(
  "Exact recurring terms across >=2 databases: ",
  n_recurring_terms
)

message(
  "  BRCA1_tumour: ",
  n_recurrent_brca1
)

message(
  "  TotalCell: ",
  n_recurrent_totalcell
)

if (
  n_recurring_terms > 0
) {
  
  message(
    "Recurrence plot: generated."
  )
  
} else {
  
  message(
    "Recurrence plot: intentionally skipped because ",
    "no recurring terms were identified."
  )
  
}

message(
  "Stage 12 output directory: ",
  stage12_dir
)

message(
  "Log saved to: ",
  log_file
)

message(
  "===================================================="
)


############################################################
# 12.33 Stage Check
############################################################

list.files(
  path = stage12_dir,
  recursive = TRUE,
  full.names = FALSE
)

table(
  gsea_integrated$Database,
  gsea_integrated$Enrichment_Group
)

table(
  gsea_integrated$Enrichment_Group
)

gsea_integrated %>%
  
  count(
    Description,
    Enrichment_Group
  ) %>%
  
  filter(
    n > 1
  ) %>%
  
  arrange(
    desc(n)
  )

gsea_integrated %>%
  
  group_by(
    Enrichment_Group
  ) %>%
  
  summarise(
    
    n = n(),
    
    min_NES =
      min(
        NES,
        na.rm = TRUE
      ),
    
    max_NES =
      max(
        NES,
        na.rm = TRUE
      )
    
  )

