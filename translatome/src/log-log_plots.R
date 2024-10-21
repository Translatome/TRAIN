#!/bin/env Rscript
#'=============================================================================
#'                                Log-Log plots R 
#'=============================================================================
#'
#' Description:
#' ------------
#' This script allows you to create a log-log plot.
#'
#' Requirements:
#' -------------
#' - dplyr: 1.1.1
#' - gplots: 3.1.3
#' - ggplot2: 3.4.1
#' - gghighlight: 0.4.0
#' - tidyr: 1.3.0
#' - repr: 1.1.6
#' All these packages are managed by the conda environment file
#'
#' Executions:
#' -----------
#' Execution in snakemake rule:
#'    shell:
#'        "src/log-log_plots.R"
#'
#' Execution without snakemake:
#'    Replace all paths provided by snakemake in the snakemake value section
#'
#' Credits:
#' --------
#' - Author: Julie Ripoll
#' - Date: 2022-09-23
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
#' - density plots
#'    PNG & PDF files
#'    density curves according to density of genes
#' - log-log plots
#'    PNG & PDF files
#'    scatter-plot of logFC of transcriptome in x and translatome in y
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
LOG_FILE    <- snakemake@log

## Input file
INPUT_PATH  <- snakemake@input["input_path"]

## Params
LOG_FC      <- snakemake@params[["logFC"]]
STAT_TOOL   <- snakemake@params[["stat_tool"]]

## Output path
OUTPUT_PATH <- snakemake@output[["output_path"]]
SESSION     <- snakemake@output[["sessioninfo"]]

# Import libraries
library("dplyr")
library("tidyr")
library(gplots)
library(ggplot2)
library(gghighlight)
library(repr) # plot size

#==============================================================================

# Save session information
log <- file(as.character(LOG_FILE), open = "wt")
sink(log)

# Start analyse
print("Start_time")
Sys.time()

#==============================================================================
# Import working files
#==============================================================================
temp <- list.files(
  path = as.character(INPUT_PATH),
  pattern = "Categorization*",
  full.names = TRUE
)

list2env(
    lapply(
        setNames(temp, make.names(gsub("*.tsv$", "", temp))),
        read.delim
    ),
    envir = .GlobalEnv
)

Pattern_list <- list(mget(ls(pattern = "Categorization")))[[1]]


#==============================================================================
# Iterate on working files
#==============================================================================

for (file in 1:length(Pattern_list)) {
    print("Info text: Import of transcriptome and translatome data from 
          Explore_results")
    data      <- as.data.frame(Pattern_list[file])
    file_name <- names(Pattern_list[file][1])
    print(str(data))
    print(file_name)
    # Split filename
    file_name_end <- strsplit(file_name, split = "Categorization_")[[1]][2]
    print(file_name_end)

    # Rename columns of dataframe
    colnames(data) <- sub(file_name, "", colnames(data))
    colnames(data) <- sub(".", "", colnames(data))

    # Keep good column names for plots
    if (STAT_TOOL == "deseq2") {
        aes_de <- aes(
            x = log2FoldChange_transcription,
            y = log2FoldChange_translation
        )
        aes_cat <- aes(
            x = log2FoldChange_transcription,
            y = log2FoldChange_translation,
            fill = Categories
        )
    } else if (STAT_TOOL == "limma") {
        aes_de <- aes(
          x = logFC_transcription,
          y = logFC_translation
        )
        aes_cat <- aes(
            x = logFC_transcription,
            y = logFC_translation,
            fill = Categories
        )
    } else {
        print("Please add information on the tool used for statistics")
    }

    # Delete NA for log-log plot
    data_wo_na <- na.omit(data)
    print("Info text: Omit NA if exists")
    print(str(data_wo_na))

    #==========================================================================
    # Plot logFC density
    #==========================================================================

    RawPlot <- ggplot(data, aes_de) +
        geom_point(size = 2, shape = 23) +
        theme_classic() +
        geom_rug() + # Marginal density
        geom_density_2d() + # 2D density
        labs(
            title = "Log2 fold Change plot",
            x = "Transcribed mRNA log2FC",
            y = "Translated mRNA log2FC"
        )
    RawPlot

    # Save plot
    dir.create(as.character(OUTPUT_PATH), showWarnings = FALSE)

    ggsave(
        file = paste0("Density_plot_log2FC_", file_name_end, ".pdf"),
        RawPlot,
        path = OUTPUT_PATH,
        height = 6, width = 5 * 1.5
    ) # PDF
    ggsave(
        file = paste0("Density_plot_log2FC_", file_name_end, ".png"),
        RawPlot,
        path = OUTPUT_PATH,
        height = 6, width = 5 * 1.5
    ) # PNG


    #==========================================================================
    # Color dict
    #==========================================================================

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

    #==========================================================================
    # Log-log plot alias "graphique en fuseau"
    #==========================================================================

    scat <- ggplot(data_wo_na, aes_cat) +
        geom_point(
          size = 1.7,
          shape = 21
        ) +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            plot.background = element_rect(
                fill = "white",
                color = "white"
            ),
            panel.background = element_rect(
                fill = "white",
                linetype = "solid"
            )
        ) +
        labs(
            title = "Log2 fold Change plot",
            x = "Transcribed mRNA log2FC",
            y = "Translated mRNA log2FC"
        ) +
        geom_vline(
            xintercept = c(-as.numeric(LOG_FC), as.numeric(LOG_FC)),
            linewidth = 0.5, linetype = "dotted"
        ) +
        geom_hline(
            yintercept = c(-as.numeric(LOG_FC), as.numeric(LOG_FC)),
            linewidth = 0.5, linetype = "dotted"
        ) +
        geom_vline(
          xintercept = 0,
          linewidth = 0.6,
          linetype = "longdash"
        ) +
        geom_hline(
          yintercept = 0,
          linewidth = 0.6,
          linetype = "longdash"
        ) +
        theme(legend.text = element_text(size = 10)) +
        scale_fill_manual(values = dict)
    scat

    # Save plot
    ggsave(
        file = paste0("Scatter_plot_log2FC_", file_name_end, ".pdf"),
        scat,
        path = OUTPUT_PATH,
        height = 6,
        width = 5 * 1.5
    ) # PDF
    ggsave(
        file = paste0("/Scatter_plot_log2FC_", file_name_end, ".png"),
        scat,
        path = OUTPUT_PATH,
        height = 6,
        width = 5 * 1.5
    ) # PNG
}


#==========================================================================

# Export session
save.image(as.character(SESSION))

# End time
print("End_time")
print(Sys.time())

# End save info session
sink()
