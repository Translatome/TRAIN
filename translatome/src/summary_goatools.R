#!/bin/env Rscript
#'=============================================================================
#'                             Summary of goatools
#'=============================================================================
#'
#' Description:
#' ------------
#' Summary of the number of terms according to the depth or level of hierarchy
#' information for the Gene Ontologies and the categories and sources
#' GO: BP, MF & CC.
#'
#' Requirements:
#' -------------
#' All these packages are managed by the conda environment file
#'
#' Executions:
#' -----------
#' Execution in snakemake rule:
#'    script:
#'        "src/summary_goatools.R"
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
#'    Includes all result files for GO: BP, MF & CC and category
#'
#' Outputs:
#' --------
#' - goatool summary depth
#'    CSV file
#'    Number of terms according to the depth of the hierarchy,
#'    categories and source GO: BP, MF and CC
#' - goatool summary level
#'    CSV file
#'    Number of terms according to the level of the hierarchy,
#'    categories and source GO: BP, MF and CC
#' - goatool summary depth graph
#'    PNG file
#'    Number of terms according to the depth of the hierarchy,
#'    categories and source GO: BP, MF and CC
#' - goatool summary level graph
#'    PNG file
#'    Number of terms according to the level of the hierarchy,
#'    categories and source GO: BP, MF and CC
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
# Snakemake variables and libraries
#==============================================================================
# Snakemake variables
## Input
PATH_INPUT_GOA          <- snakemake@params[["hierarchy_table"]]

## Outputs
OUTPUT_GOA_DEPTH        <- snakemake@output[["output_goa_depth"]]
OUTPUT_GOA_LEVEL        <- snakemake@output[["output_goa_level"]]
OUTPUT_GOA_DEPTH_GRAPH  <- snakemake@output[["output_goa_depth_graph"]]
OUTPUT_GOA_LEVEL_GRAPH  <- snakemake@output[["output_goa_level_graph"]]
OUTPUT_GOA_DEPTH_GRAPH2 <- snakemake@output[["output_goa_depth_graph2"]]
OUTPUT_GOA_LEVEL_GRAPH2 <- snakemake@output[["output_goa_level_graph2"]]
SESSION                 <- snakemake@output[["session_info"]]

## Log file
LOG_FILE                <- snakemake@log

# Import libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

#==============================================================================

# Save session information
log <- file(as.character(LOG_FILE), open = "wt")
sink(log)

# Start analyse
print("Start_time")
Sys.time()


#==============================================================================
# Summary depth table for goatools
#==============================================================================

# Depth
dt_BP <- data.frame("Depth" = c(
  "D00", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10",
  "D11", "D12", "D13", "D14", "D15"
))

dt_MF <- data.frame("Depth" = c(
  "D00", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10",
  "D11", "D12", "D13", "D14", "D15"
))

dt_CC <- data.frame("Depth" = c(
  "D00", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10",
  "D11", "D12", "D13", "D14", "D15"
))

# Pol-seq or Kinetic cases with categories
catego <- c(
  "Both_mRNA_DOWN", "Both_mRNA_UP", "Divergent_DNUP", "Divergent_UPDN",
  "Transcription_DOWN", "Transcription_UP", "Translation_DOWN", "Translation_UP"
)

dt_ob <- data.frame(
  "Depth"              = "obsolete_terms",
  "Both_mRNA_DOWN"     = "",
  "Both_mRNA_UP"       = "",
  "Divergent_DNUP"     = "",
  "Divergent_UPDN"     = "",
  "Transcription_DOWN" = "",
  "Transcription_UP"   = "",
  "Translation_DOWN"   = "",
  "Translation_UP"     = ""
)

# Filling of summary depth table
for (cat in catego) {
  # File name construction
  inputfile_depth <- paste0(
    PATH_INPUT_GOA,
    "summary_depth_",
    as.character(cat),
    "_gene_ontology.tsv"
  )
  if (file.size(inputfile_depth) <= 3) {
    # If there are no results for this file then add 0 in the summary table
    print("File is empty")
    dt_MF[[cat]] <- rep(0, 16)
    dt_BP[[cat]] <- rep(0, 16)
    dt_CC[[cat]] <- rep(0, 16)
    dt_ob[[cat]] <- 0
  } else {
    file <- read.delim(inputfile_depth, header = TRUE)
    colnames(file) <- c("Source", "Depth", cat)
    print(head(file))
    MF <- subset(file, Source == "MF")
    MF <- MF[, 2:3]
    BP <- subset(file, Source == "BP")
    BP <- BP[, 2:3]
    CC <- subset(file, Source == "CC")
    CC <- CC[, 2:3]
    dt_MF <- merge(dt_MF, MF, by = "Depth", all = TRUE)
    dt_BP <- merge(dt_BP, BP, by = "Depth", all = TRUE)
    dt_CC <- merge(dt_CC, CC, by = "Depth", all = TRUE)
    obsolete <- subset(file, Source == "obsolete_terms")
    if (nrow(obsolete) == 0) {
      obsolete <- 0
      dt_ob[[cat]] <- obsolete
    } else {
      obsolete     <- obsolete[, 3]
      dt_ob[[cat]] <- obsolete
    }
  }
}
dt_tot <- rbind(dt_BP, dt_CC, dt_MF, dt_ob)
dt_tot$Source <- c(
  rep("BP", 16),
  rep("CC", 16),
  rep("MF", 16),
  "obsolete_terms"
)

# Writing of summary depth table

dt_tot[is.na(dt_tot)] <- 0
write.csv(dt_tot, OUTPUT_GOA_DEPTH, row.names = FALSE)
print(str(dt_tot))

# Graph of summary depth
## Processing table for graph
dt_tot <- as.data.frame(dt_tot)
table_gather <- gather(
  dt_tot,
  key="Categories",
  value="Value",
  2:9
)

table_gather$Source <- str_replace(
  table_gather$Source,
  "obsolete_terms",
  "OT"
)
print(table_gather)
table_gather_without_zero <- subset(table_gather, Value != 0)

dict <- list(
  "Not_significant"    = "#808080",
  "Both_mRNA_UP"       = "#0000cc",
  "Both_mRNA_DOWN"     = "#99ccff",
  "Translation_UP"     = "#990000",
  "Translation_DOWN"   = "#ff9999",
  "Transcription_UP"   = "#cccc00",
  "Transcription_DOWN" = "#ffff66",
  "Divergent_UPDN"     = "#4B0082",
  "Divergent_DNUP"     = "#8B4513"
)

plt <- ggplot(
  table_gather_without_zero, 
  aes(x = Depth, y = Value, fill = Categories)
  ) +
  geom_bar(stat = "identity") +
  facet_grid(
    .~Source,
    scale = "free",
    space= "free"
  ) +
  theme_classic() +
  theme(
    strip.placement = "outside", 
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    legend.title = element_text(color = "black", size = 6),
    legend.text = element_text(color = "black", size = 5),
    strip.text.x = element_text(size = 8)
  ) +
  labs(
    x = "Depth",
    y = "Number of terms"
  ) +
  scale_fill_manual(values = dict)

ggsave(
  plot     = plt,
  width    = 8,
  height   = 4,
  dpi      = 300,
  filename = OUTPUT_GOA_DEPTH_GRAPH
)

ggsave(
  plot     = plt,
  filename = OUTPUT_GOA_DEPTH_GRAPH2
)


#==============================================================================
# Summary level table for goatools
#==============================================================================

# Level
dt_BP_l <- data.frame("Level" = c(
  "L00", "L01", "L02", "L03", "L04", "L05", "L06", "L07", "L08", "L09", "L10",
  "L11","L12", "L13", "L14", "L15"
))

dt_MF_l <- data.frame("Level" = c(
  "L00", "L01", "L02", "L03", "L04", "L05", "L06", "L07", "L08", "L09", "L10",
  "L11","L12", "L13", "L14", "L15"
))

dt_CC_l <- data.frame("Level" = c(
  "L00", "L01", "L02", "L03", "L04", "L05", "L06", "L07", "L08", "L09", "L10",
  "L11","L12", "L13", "L14", "L15"
))

# Pol-seq or Kinetic cases with categories
catego <- c(
  "Both_mRNA_DOWN", "Both_mRNA_UP", "Divergent_DNUP", "Divergent_UPDN",
  "Transcription_DOWN", "Transcription_UP", "Translation_DOWN", "Translation_UP"
)
dt_ob_l <- data.frame(
  "Level"              = "obsolete_terms",
  "Both_mRNA_DOWN"     = "",
  "Both_mRNA_UP"       = "",
  "Divergent_DNUP"     = "",
  "Divergent_UPDN"     = "",
  "Transcription_DOWN" = "",
  "Transcription_UP"   = "",
  "Translation_DOWN"   = "",
  "Translation_UP"     = ""
)
# Filling of summary level table
for (cat in catego) {
  inputfile_level <- paste0(
    PATH_INPUT_GOA,
    "summary_level_",
    as.character(cat),
    "_gene_ontology.tsv"
  ) # File name construction
  if (file.size(inputfile_level) <= 3) {
    # If there are no results for this file then add 0 in the summary table
    print("File is empty")
    dt_MF_l[[cat]] <- rep(0, 16)
    dt_BP_l[[cat]] <- rep(0, 16)
    dt_CC_l[[cat]] <- rep(0, 16)
    dt_ob_l[[cat]] <- 0
  } else {
    file_l <- read.delim(inputfile_level, header = TRUE)
    colnames(file_l) <- c("Source", "Level", cat)
    print(file_l)
    MF_l <- subset(file_l, Source == "MF")
    MF_l <- MF_l[, 2:3]
    BP_l <- subset(file_l, Source == "BP")
    BP_l <- BP_l[, 2:3]
    CC_l <- subset(file_l, Source == "CC")
    CC_l <- CC_l[, 2:3]
    dt_MF_l    <- merge(dt_MF_l, MF_l, by = "Level", all = TRUE)
    dt_BP_l    <- merge(dt_BP_l, BP_l, by = "Level", all = TRUE)
    dt_CC_l    <- merge(dt_CC_l, CC_l, by = "Level", all = TRUE)
    obsolete_l <- subset(file_l, Source == "obsolete_terms")
    if (nrow(obsolete_l) == 0) {
      obsolete_l     <- 0
      dt_ob_l[[cat]] <- obsolete_l
      print(dt_ob_l)
    } else {
      obsolete_l     <- obsolete_l[, 3]
      dt_ob_l[[cat]] <- obsolete_l
    }
  }
}
dt_tot_l <- rbind(dt_BP_l, dt_CC_l, dt_MF_l, dt_ob_l)
dt_tot_l$Source <- c(
  rep("BP", 16),
  rep("CC", 16),
  rep("MF", 16),
  "obsolete_terms"
)

# Writing of summary level table
dt_tot_l[is.na(dt_tot_l)] <- 0

write.csv(
  dt_tot_l,
  OUTPUT_GOA_LEVEL,
  row.names = FALSE
)

# Graph of summary level
## Processing table for graph
dt_tot_l <- as.data.frame(dt_tot_l)
table_gather_l <- gather(
  dt_tot_l,
  key = "Categories",
  value = "Value",
  2:9
)
table_gather_l$Source <- str_replace(
  table_gather_l$Source,
  "obsolete_terms",
  "OT"
)
print(table_gather_l)
table_gather_without_zero_l <- subset(table_gather_l, Value != 0)

plt_l <- ggplot(
  table_gather_without_zero_l,
  aes(x = Level, y = Value, fill = Categories)
  ) +
  geom_bar(stat = "identity") +
  facet_grid(
    .~Source,
    scale = "free",
    space= "free"
  ) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    legend.title = element_text(color = "black", size = 6),
    legend.text = element_text(color = "black", size = 5),
    strip.text.x = element_text(size = 8)
  ) +
  labs(
    x = "Level",
    y = "Number of terms"
  ) +
  scale_fill_manual(values = dict)

ggsave(
  plot = plt_l,
  width = 8,
  height = 4,
  dpi = 300,
  filename = OUTPUT_GOA_LEVEL_GRAPH
)

ggsave(
  plot = plt_l,
  filename = OUTPUT_GOA_LEVEL_GRAPH2
)


#==============================================================================
# Export session
save.image(as.character(SESSION))

# End time
print("End_time")
print(Sys.time())

# End save info session
sink()
