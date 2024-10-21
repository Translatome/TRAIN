#!/bin/env Rscript

#'=============================================================================
#'                            Summary of gprofiler
#'=============================================================================
#'
#' Description:
#' ------------
#' The purpose of this script is to make a summary of the number of terms for
#' each gprofiler source (KEGG, REAC, WP, GO:BP, GO:MF, GO:CC and CORUM) for
#' each category file (Both_mRNA_DOWN, Both_mRNA_UP, Divergent_DNUP, 
#' Divergent_UPDN, Transcription_DOWN, Transcription_UP, Translation_DOWN, 
#' Translation_UP).
#'
#' Requirements:
#' -------------
#' All these packages are managed by the conda environment file
#'
#' Executions:
#' -----------
#' Execution in snakemake rule:
#'    script:
#'        "src/summary_gprofiler.R"
#'
#' Execution without snakemake:
#'    Replace all paths provided by snakemake in the snakemake value section
#'
#' Credits:
#' --------
#' - Authors: Céline Mandier and Julie Ripoll
#' - Created: 2023-04-27
#' - License: CeCILL
#' - Version: 0.1.0
#'
#' Inputs:
#' -------
#' - Pathways enrichment folder
#'    Includes all result files for each source and category
#'
#' Outputs:
#' --------
#'  - gProfiler summary
#'    CSV file
#'    Number of terms for each category by source
#' - Session info
#'    RData file
#'    Collect Information About the Current R Session
#'
#' Log:
#' ----
#' - log file
#'    LOG catch terminal stdout
#'

#==============================================================================
# Snakemake variables
#==============================================================================

# Inputs
PATH_INPUT_GPRO <- snakemake@input[["input_gpro"]]
ORGANISM_NAMES  <- snakemake@input[["gpro_names"]]

# params
VERSION         <- snakemake@params["gpro_version"]
USER_ORGANISM   <- snakemake@params[["organism"]]

# Output
LOG_FILE        <- snakemake@log
OUTPUT_GPRO     <- snakemake@output[["output_gpro"]]
SESSION         <- snakemake@output[["session_info"]]


#------------------------------------------------------------------------------

# save session information
log <- file(as.character(LOG_FILE), open = "wt")
sink(log)

# Start analyse
print("Start_time")
Sys.time()


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

# Export databases API version of grpofiler
version_info <- paste0(
  "curl 'https://biit.cs.ut.ee/gprofiler/api/util/data_versions?organism=",
  gpro_orgn_var$id,
  "' > ",
  VERSION
)
system(version_info)


#==============================================================================
# Prepare summary table for gprofiler
#==============================================================================

folder <- c("/Gene_ontology/", "/Corum/", "/Pathway/")

# create output file according to catgeories for pol-seq and kinetic analyses
catego <- c(
  "Both_mRNA_DOWN", "Both_mRNA_UP", "Divergent_DNUP",
  "Divergent_UPDN", "Transcription_DOWN", "Transcription_UP",
  "Translation_DOWN", "Translation_UP"
)
dt <- data.frame(
  "Categorie" = catego,
  "CORUM" = "",
  "GOBP"  = "",
  "GOCC"  = "",
  "GOMF"  = "",
  "KEGG"  = "",
  "REAC"  = "",
  "WP"    = ""
)
# Length for loop
loop_length <- nrow(dt)


#==============================================================================
# Filling of summary table
#==============================================================================

root_path <- getwd()

for (i in folder) {
  # Move to the folder of interest
  print(paste0(root_path, "/", PATH_INPUT_GPRO, i))
  setwd(paste0(root_path, "/", PATH_INPUT_GPRO, i))
  for (j in 1:loop_length) {
    # For each category file and each source count the number of terms present
    if (i == "/Corum/") {
      # File name construction
      file <- paste0("Complex_table_", as.character(dt[j, 1]), ".tsv")
      print(as.character(file))
      dt_gpro <- read.delim(file, header = TRUE)
      if (colnames(dt_gpro[1]) == "x") {
        # If there are no results for this file then add 0 in the summary table
        dt$CORUM[j] <- 0
      } else {
        # For the CORUM source there is only one source listed in the category
        # files, so we can directly list the number of lines
        nb_corum    <- nrow(dt_gpro)
        dt$CORUM[j] <- nb_corum
      }
    } else if (i == "/Gene_ontology/") {
      # File name construction
      file <- paste0("GO_table_", as.character(dt[j, 1]), ".tsv")
      print(as.character(file))
      dt_gpro <- read.delim(file, header = TRUE)
      if (colnames(dt_gpro[1]) == "x") {
        # If there are no results for this file then add 0 in the summary table
        dt$GOBP[j] <- 0
        dt$GOCC[j] <- 0
        dt$GOMF[j] <- 0
      } else {
        # For the Gene_ontology folder there are several sources present in the
        # category files, so we list the lines corresponding to each source
        nb_gobp <- nrow(subset(dt_gpro, source == "GO:BP"))
        nb_gocc <- nrow(subset(dt_gpro, source == "GO:CC"))
        nb_gomf <- nrow(subset(dt_gpro, source == "GO:MF"))
        dt$GOBP[j] <- nb_gobp
        dt$GOCC[j] <- nb_gocc
        dt$GOMF[j] <- nb_gomf
      }
    } else if (i == "/Pathway/") {
      # File name construction
      file <- paste0("Pathways_table_", as.character(dt[j, 1]), ".tsv")
      print(as.character(file))
      dt_gpro <- read.delim(file, header = TRUE)
      if (colnames(dt_gpro[1]) == "x") {
        # If there are no results for this file then add 0 in the summary table
        dt$KEGG[j] <- 0
        dt$REAC[j] <- 0
        dt$WP[j]   <- 0
      } else {
        # For the Pathway folder there are several sources present in the
        # category files, so we list the lines corresponding to each source
        nb_kegg <- nrow(subset(dt_gpro, source == "KEGG"))
        nb_reac <- nrow(subset(dt_gpro, source == "REAC"))
        nb_wp   <- nrow(subset(dt_gpro, source == "WP"))
        dt$KEGG[j] <- nb_kegg
        dt$REAC[j] <- nb_reac
        dt$WP[j]   <- nb_wp
      }
    } else {
      print("No enrichment provided, please check your gprofiler results. 
            They should be saved in at least 3 folders: Pathway, 
            Corum and Gene_ontology.")
    }
  }
}

#==============================================================================
# Writing of summary table
#==============================================================================

write.csv(
  dt,
  paste0(root_path, "/", OUTPUT_GPRO),
  row.names = FALSE
)

#==============================================================================

# Export session
save.image(as.character(paste0(root_path, "/", SESSION)))

# End time
print("End_time")
print(Sys.time())

# End save info session
sink()
