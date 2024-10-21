#!/bin/env Rscript
#'=============================================================================
#'               Multidimensional scaling plot for raw data
#'=============================================================================
#'
#' Description:
#' ------------
#' This script allows you to plot a multidimensional scaling graph on raw data.
#' Depending on the results, users can choose between the deseq2 and limma R
#' packages for their statistical analysis.
#' We use limma when replicates present an important bias.
#'
#' Requirements:
#' -------------
#' - data.table: 1.14.8
#' - edgeR: 3.38.1
#' - ggplot2: 3.4.1
#' - ggpubr: 0.6.0
#' - src/some_functions: 0.1.0
#' All these packages are managed by the conda environment file
#'
#' Executions:
#' -----------
#' Execution in snakemake rule:
#'    shell:
#'        "src/explore_counts.R"
#'
#' Execution without snakemake:
#'    Replace all paths provided by snakemake in the snakemake value section
#'
#' Credits:
#' --------
#' - Authors: Céline Mandier and Julie Ripoll
#' - Created: 2023-04-06
#' - License: CeCILL
#' - Version: 0.1.0
#'
#' Inputs:
#' -------
#' - Count table
#'    TSV file & CSV file are supported
#'    output from HTSeq-counts or featureCounts
#' - metadata
#'    TSV file & CSV file are supported
#'    please adapt and complete the template provided in ressources
#'    contains: sample names, fraction, treatment, replicate
#'
#' Outputs:
#' --------
#' - MDS plot
#'    PNG & PDF files, multidimensional scaling plots
#' - session information
#'    RData archive
#'

#==============================================================================
# Snakemake variables and libraries
#==============================================================================
# Snakemake values

## Inputs
FILE_DATA            <- snakemake@input[["data"]]
FILE_METADATA        <- snakemake@input[["metadata"]]

## Params
ANALYSIS             <- snakemake@params["analysis"]
COUNT_TOOLS          <- snakemake@params[["count_tools"]]
NB_REP               <- snakemake@params[["replicat"]]

## Outputs
OUTPUT_FILE_PNG      <- snakemake@output[["mdsplot"]]
OUTPUT_FILE_PDF      <- snakemake@output[["mdsplotpdf"]]
OUTPUT_SCREE_PNG     <- snakemake@output[["screeplot"]]
OUTPUT_SCREE_PDF     <- snakemake@output[["screeplotpdf"]]
SESSION              <- snakemake@output[["sessioninfo"]]

## Special outputs
OUTPUT_FRACTION_PNG  <- snakemake@output[["fractionplot"]]
OUTPUT_FRACTION_PDF  <- snakemake@output[["fractionplotpdf"]]
OUTPUT_TREATMENT_PNG <- snakemake@output[["treatmentplot"]]
OUTPUT_TREATMENT_PDF <- snakemake@output[["treatmentplotpdf"]]


# Import libraries
library("data.table") # version 1.14.8
library("edgeR") # version 3.38.1
library("ggplot2") # version 3.4.1
library("ggpubr") # version 0.6.0 for plots

# Import functions
source("translatome/src/some_functions.R")


#==============================================================================
# Prepare Data
#==============================================================================

# Import count file
data <- ImportFileAsDF(FILE_DATA, COUNT_TOOLS)

# Order by names
orderData <- data[, order(names(data))]
print("Colnames of count data:")
print(colnames(orderData))

#==============================================================================
# Metadata
#==============================================================================

# Import metadata file
metadata <- fread(FILE_METADATA, header = TRUE, data.table = FALSE)

# Order by names
ordermetadata <- metadata[order(metadata$sample_name), ]
print("Metadata file:")
print(ordermetadata)


#==============================================================================
# Pre-processing
#==============================================================================

# Create a new variable combine that combine Fraction and Treatment / Time
if (ANALYSIS == "kinetic") {
  combine <- interaction(
    ordermetadata$fraction,
    ordermetadata$treatment,
    ordermetadata$time
  )
} else {
  combine <- interaction(
    ordermetadata$fraction,
    ordermetadata$treatment
  )
}

#==============================================================================
# Multidimensional scaling plot
#==============================================================================

# Creation of the DGE object and normalize on library size
DGE_objt <- DGEList(
  counts = orderData,
  norm.factors = rep(1,ncol(orderData)),
  remove.zeros = TRUE
)

# MDS plot
mds <- plotMDS(DGE_objt, col = as.numeric(combine))
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y)

fig <- ggplot(
  toplot,
  aes(
    Dim1,
    Dim2,
    color = as.factor(combine),
    shape = as.factor(ordermetadata$replicate)
  )
) +
  labs(
    color = "Factors",
    shape = "Replicates"
  ) +
  geom_point(size = 3) +
  xlab(paste0(
    "PC1: ", 
    round(mds$var.explained[1] * 100, 2),
    "% variance"
  )) +
  ylab(paste0(
    "PC2: ",
    round(mds$var.explained[2] * 100, 2), 
    "% variance"
  )) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(aspect.ratio = 1)

png(
  file = OUTPUT_FILE_PNG,
  width = 900,
  height = 650
)
plot(fig)
dev.off()

pdf(file = OUTPUT_FILE_PDF)
plot(fig)
dev.off()


#==============================================================================
# ScreeHist plot
#==============================================================================

# Select variance values
mds_axis <- as.data.frame(round(mds$var.explained * 100, 2))
mds_axis$components <- seq.int(nrow(mds_axis))
colnames(mds_axis) <- c("variance", "component_number")

# Perform cumulative variance
mds_axis$cumulative_variance <- cumsum(mds_axis$variance)

# Scree plot
histplot <- ScreeHist(mds_axis)

png(
  file = OUTPUT_SCREE_PNG,
  width = 900,
  height = 650
)
plot(histplot)
dev.off()

pdf(file = OUTPUT_SCREE_PDF)
plot(histplot)
dev.off()


#==============================================================================
# One factor MDS plot
#==============================================================================

# MDS plot for fraction
mds <- plotMDS(
  DGE_objt, 
  col = as.numeric(ordermetadata$fraction)
)
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y)

fig <- ggplot(
  toplot,
  aes(
    Dim1,
    Dim2,
    color = as.factor(ordermetadata$fraction),
    shape = as.factor(ordermetadata$replicate)
  )
) +
  labs(
    color = "Factors",
    shape = "Replicates"
  ) +
  geom_point(size = 3) +
  xlab(paste0(
    "PC1: ",
    round(mds$var.explained[1] * 100, 2),
    "% variance"
  )) +
  ylab(paste0(
    "PC2: ",
    round(mds$var.explained[2] * 100, 2),
    "% variance"
  )) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(aspect.ratio = 1)

png(
  file = OUTPUT_FRACTION_PNG,
  width = 900,
  height = 650
)
plot(fig)
dev.off()

pdf(file = OUTPUT_FRACTION_PDF)
plot(fig)
dev.off()


# MDS plot for treatment
mds <- plotMDS(
  DGE_objt, 
  col = as.numeric(ordermetadata$treatment)
)
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y)

fig <- ggplot(
  toplot,
  aes(
    Dim1,
    Dim2,
    color = as.factor(ordermetadata$treatment),
    shape = as.factor(ordermetadata$replicate)
  )
) +
  labs(
    color = "Factors",
    shape = "Replicates"
  ) +
  geom_point(size = 3) +
  xlab(paste0(
    "PC1: ",
    round(mds$var.explained[1] * 100, 2),
    "% variance"
  )) +
  ylab(paste0(
    "PC2: ",
    round(mds$var.explained[2] * 100, 2),
    "% variance"
  )) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(aspect.ratio = 1)

png(
  file = OUTPUT_TREATMENT_PNG,
  width = 900,
  height = 650
)
plot(fig)
dev.off()

pdf(file = OUTPUT_TREATMENT_PDF)
plot(fig)
dev.off()


if (ANALYSIS == "kinetic") {

  OUTPUT_TIME_PNG <- snakemake@output[["timeplot"]]
  OUTPUT_TIME_PDF <- snakemake@output[["timeplotpdf"]]

  # MDS plot for time
  mds <- plotMDS(
    DGE_objt,
    col = as.numeric(ordermetadata$time)
  )
  toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y)

  fig <- ggplot(
    toplot,
    aes(
      Dim1,
      Dim2,
      color = as.factor(ordermetadata$time),
      shape = as.factor(ordermetadata$replicate)
    )
  ) +
    labs(
      color = "Factors",
      shape = "Replicates"
    ) +
    geom_point(size = 3) +
    xlab(paste0(
      "PC1: ",
      round(mds$var.explained[1] * 100, 2),
      "% variance"
    )) +
    ylab(paste0(
      "PC2: ",
      round(mds$var.explained[2] * 100, 2),
      "% variance"
    )) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    theme(aspect.ratio = 1)

  png(
    file = OUTPUT_TIME_PNG,
    width = 900,
    height = 650
  )
  plot(fig)
  dev.off()

  pdf(file = OUTPUT_TIME_PDF)
  plot(fig)
  dev.off()
}


#==============================================================================

# Export session
save.image(SESSION)
