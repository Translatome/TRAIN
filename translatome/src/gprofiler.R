#!/bin/env Rscript
#'
#'=============================================================================
#'                         Analysis with gProfiler2
#'=============================================================================
#'
#' Description:
#' ------------
#' This script allows you to call the API of gprofileR for gene ontology
#' and pathway enrichments.
#'
#' Requirements:
#' -------------
#' - r-base: 4.1.3
#' - data.table: 1.14.8
#' - ggplot2: 3.3.6
#' - gprofiler2: 0.2.1
#' - src/some_functions.R: 0.1.0
#' All these packages are managed by the conda environment file
#'
#' Executions:
#' -----------
#' Execution in snakemake rule:
#'    shell:
#'        "src/gprofiler.R"
#'
#' Execution without snakemake:
#'    Replace all paths provided by snakemake in the snakemake value section
#'
#' Credits:
#' --------
#' - Author: Céline Mandier and Julie Ripoll
#' - Date: 2022-07-22
#' - License: CeCILL
#' - Version: 0.1.0
#'
#' Inputs:
#' -------
#' - Count table
#'    TSV file, output from explore_results.py
#'    contains all counts for transcriptome and translatome
#'
#' Outputs:
#' --------
#' - results tables
#'    TSV files
#'    one file by sources, i.e. pathway, gene ontology and protein complex
#'    and by category of transcription / translation / both / divergent
#' - Manhattan enrichment plots
#'    PNG files
#'    top terms selected by adjusted p.value and user threshold in a manhattan
#'    plot. One plot by sources and category.
#' - dotplots
#'    PNG files
#'    top terms selected by term size and user threshold in a dotplot plot.
#'    One plot by sources and category.
#' - session information
#'    RData archive
#'
#' Log:
#' ----
#' - log file
#'    LOG catch terminal stdout
#'

#==============================================================================
# Snakemake variables and libraries
#==============================================================================
# Snakemake values
## Log file
LOG_FILE        <- snakemake@log

## Input file
INPUT_FILE      <- snakemake@input[["data"]]
BACKGROUND      <- snakemake@input[["genome"]]
ORGANISM_NAMES  <- snakemake@input[["gpro_names"]]

## Params
TOOLS           <- snakemake@params[["stat_tools"]]
USER_ORGANISM   <- snakemake@params[["organism"]]
CORRECTION      <- snakemake@params[["correction"]]
PVAL_THRESHOLD  <- snakemake@params[["pval_threshold"]]
PATHWAY_SOURCES <- snakemake@params[["pathways"]]
GO_SOURCES      <- snakemake@params[["go_sources"]]
COMPLEX         <- snakemake@params[["complex_sources"]]
TERM_NUMBER     <- snakemake@params[["term_number"]]

## Output path
PATHWAY_OUTPUT  <- snakemake@output[["pathway"]]
GO_OUTPUT       <- snakemake@output[["go_output"]]
CORUM_OUTPUT    <- snakemake@output[["corum_output"]]
SESSION         <- snakemake@output[["session_info"]]

# Import libraries
library("data.table")
library("gprofiler2")
library("ggplot2")
library("dplyr")

# Import functions
source("translatome/src/some_functions.R")

#==============================================================================

# Save session information
log <- file(as.character(LOG_FILE), open = "wt")
sink(log)

# Start analyse
print("Start_time")
Sys.time()


#==============================================================================
# Variables and data information
#==============================================================================

# Variables selected by end-user
print(paste0("Multiple testing method: ", CORRECTION))
print(paste0("p-value threshold for the multiple testing: ", PVAL_THRESHOLD))
print("Sources: ")
print(PATHWAY_SOURCES)
print(GO_SOURCES)
print(COMPLEX)
print(paste0("Top X terms to plot: ", TERM_NUMBER))

# Import the file that you want to analyse
data <- fread(
  INPUT_FILE,
  header = TRUE,
  data.table = FALSE
)
print(str(data))

# Import file for background
genome <- fread(
  BACKGROUND,
  header = TRUE,
  data.table = FALSE
)
print(str(genome))

#==============================================================================
# Convert scientific name to id for gprofiler
#==============================================================================

gpro_names <- read.delim(
  as.character(ORGANISM_NAMES),
  sep = "\t",
  header = TRUE
)

gpro_names$ensembl_name <- sub(" ", "_", tolower(gpro_names$scientific.name))

gpro_orgn_var <- gpro_names[gpro_names$ensembl_name == USER_ORGANISM, ]
print("Organism information available on gprofiler: ")
print(gpro_orgn_var)

# keep first item if multiple
if (length(gpro_orgn_var) > 1) {
  gpro_orgn_var <- gpro_orgn_var[1, ]
}
print("Keep first item if multiple")
print(gpro_orgn_var)

#==============================================================================
# Call API
#==============================================================================

# Pol-seq or kinetic data (subfolder by time allows to don't change code)
# Create output directories
dir.create(PATHWAY_OUTPUT, showWarnings = FALSE)
dir.create(GO_OUTPUT, showWarnings = FALSE)
dir.create(CORUM_OUTPUT, showWarnings = FALSE)

# Subset by categories and pass the identifier to the call of the API
for (cat in c(
  "Both_mRNA_UP",
  "Both_mRNA_DOWN",
  "Transcription_UP",
  "Transcription_DOWN",
  "Translation_UP",
  "Translation_DOWN",
  "Divergent_UPDN",
  "Divergent_DNUP"
)) {
  if (cat %in% unique(data$Categories)) {
    # Split data by category
    sig_genes_df <- data[data$Categories == cat, ]
    print(paste(
      "Number of genes in", cat, ":",
      length(sig_genes_df$identifier)
    ))
    # Vector of identifiers
    genes <- as.vector(sig_genes_df$identifier)
    # Functions for gprofiler2
    ## Pathways
    if (!is.null(PATHWAY_SOURCES)) {
      pathway_res <- callAPI(
        genes,
        genome$identifier,
        gpro_orgn_var$id,
        PVAL_THRESHOLD,
        CORRECTION,
        PATHWAY_SOURCES,
        paste0(PATHWAY_OUTPUT, "/Pathways_table_", cat, ".tsv", sep = "")
      )
      print(table(pathway_res$result$source))
      gproGraph(
        pathway_res,
        paste0(PATHWAY_OUTPUT, "/Pathways_graph_", cat, ".png", sep = ""),
        paste0(PATHWAY_OUTPUT, "/Pathways_dotplot_", cat, ".png", sep = ""),
        term_number = TERM_NUMBER
      )
    }
  } else {
    tab_zero <- print("No enrichment")
    write.table(
      tab_zero, 
      paste0(PATHWAY_OUTPUT, "/Pathways_table_", cat, ".tsv", sep = "")
    )
  }
}


for (cat in c(
  "Both_mRNA_UP",
  "Both_mRNA_DOWN",
  "Transcription_UP",
  "Transcription_DOWN",
  "Translation_UP",
  "Translation_DOWN",
  "Divergent_UPDN",
  "Divergent_DNUP"
)) {
  if (cat %in% unique(data$Categories)) {
    # Split data by category
    sig_genes_df <- data[data$Categories == cat, ]
    print(paste(
      "Number of genes in", cat, ":",
      length(sig_genes_df$identifier)
    ))
    # Vector of identifiers
    genes <- as.vector(sig_genes_df$identifier)
    # Functions for gprofiler2
    ## Gene ontology
    if (!is.null(GO_SOURCES)) {
      go_res <- callAPI(
        genes,
        genome$identifier,
        gpro_orgn_var$id,
        PVAL_THRESHOLD,
        CORRECTION,
        GO_SOURCES,
        paste0(
          GO_OUTPUT, "/GO_table_",
          cat, ".tsv",
          sep = ""
        )
      )
      print(table(go_res$result$source))
      gproGraph(
        go_res,
        paste0(
          GO_OUTPUT, "/GO_graph_",
          cat, ".png",
          sep = ""
        ),
        paste0(
          GO_OUTPUT, "/GO_dotplot_",
          cat, ".png",
          sep = ""
        ),
        term_number = TERM_NUMBER
      )
    }
  } else {
    tab_zero <- print("No enrichment")
    write.table(
      tab_zero, 
      paste0(
        GO_OUTPUT, "/GO_table_",
        cat, ".tsv",
        sep = ""
      )
    )
  }
}


for (cat in c(
  "Both_mRNA_UP",
  "Both_mRNA_DOWN",
  "Transcription_UP",
  "Transcription_DOWN",
  "Translation_UP",
  "Translation_DOWN",
  "Divergent_UPDN",
  "Divergent_DNUP"
)) {
  if (cat %in% unique(data$Categories)) {
    # Split data by category
    sig_genes_df <- data[data$Categories == cat, ]
    print(paste(
      "Number of genes in", cat, ":",
      length(sig_genes_df$identifier)
    ))
    # Vector of identifiers
    genes <- as.vector(sig_genes_df$identifier)
    # Functions for gprofiler2
    ## Protein complex
    if (!is.null(COMPLEX)) {
      complex_res <- callAPI(
        genes,
        genome$identifier,
        gpro_orgn_var$id,
        PVAL_THRESHOLD,
        CORRECTION,
        COMPLEX,
        paste0(
          CORUM_OUTPUT, "/Complex_table_",
          cat, ".tsv",
          sep = ""
        )
      )
      print(table(complex_res$result$source))
      gproGraph(
        complex_res,
        paste0(
          CORUM_OUTPUT, "/Complex_graph_",
          cat, ".png",
          sep = ""
        ),
        paste0(
          CORUM_OUTPUT, "/Complex_dotplot_",
          cat, ".png",
          sep = ""
        ),
        term_number = TERM_NUMBER
      )
    }
  } else {
    tab_zero <- print("No enrichment")
    write.table(
      tab_zero,
      paste0(CORUM_OUTPUT, "/Complex_table_",
             cat, ".tsv",
             sep = ""
      )
    )
  }
}


#==============================================================================

# Export session
save.image(SESSION)

# End time
print("End_time")
print(Sys.time())

# End save info session
sink()
