#!/bin/env Rscript
#'=============================================================================
#'                  Statistical kinetic analyses with DESeq2
#'=============================================================================
#'
#' Description:
#' ------------
#' This script allows you to perform DESeq2 differential expression analysis
#' on RNA-seq data.
#'
#' Requirements:
#' -------------
#' - deseq2: 1.30.1
#' - dplyr: 1.1.1
#' - gplots: 3.1.3
#' - ggplot2: 3.4.1
#' - ggpubr: 0.6.0
#' - matrixStats: 0.63.0
#' - stringr: 1.5.0
#' - RColorBrewer: 1.1_3
#' - edgeR: 3.40.0
#' - data.table: 1.14.8
#' - src/some_functions: 0.1.0
#' All these packages are managed by the conda environment file
#'
#' Executions:
#' -----------
#' Execution in snakemake rule:
#'    shell:
#'        "src/deseq2.R"
#'
#' Execution without snakemake:
#'    Replace all paths provided by snakemake in the snakemake value section
#'
#' Credits:
#' --------
#' - Author: Julie Ripoll & Céline Mandier
#' - Date: 2021-12-20
#' - License: CeCILL
#' - Version: 0.1.0
#'
#' Inputs:
#' -------
#' - Count table
#'    TSV file, output from HTSeq-counts or featureCounts
#'    CSV file is also supported
#' - metadata
#'    TSV file & CSV file are supported
#'    Please adapt and complete the template provided in ressources
#'    Contains: sample name, fraction, treatment, time, replicate
#' - metacontrast
#'    TSV file & CSV file are supported
#'    Please adapt and complete the template provided in ressources
#'    Contains: contrasts to use for comparison
#'
#' Outputs:
#' --------
#' - normalized counts
#'    TSV file with normalized counts
#' - dispersion plots
#'    PNG files dispersion and trend of the normalization model
#' - dendrogram plots
#'    PNG files euclidean distance between samples
#' - PCA and scree histograms
#'    - PNG & PDF files
#'      Scree histograms present a red line at 80% (Pareto's threshold)
#'    - TSV files, PCA axes values and cumulative variance for scree
#' - pval plots
#'    PNG files, histogram of the p-value repartition for the comparison
#' - summary of comparisons
#'    TSV file, number of DE genes in the different comparisons
#' - deseq2 table results by comparisons
#'    TSV file, named with the comparison
#'    Contains baseMean mean value of all samples
#'             log2FoldChange value
#'             lfcSE standard error of log2FC
#'             stat value of the Wald test
#'             pvalue p-value of the Wald test
#'             padj adjusted p-value of the Wald test
#'             identifier of the gene, Ensembl identifier
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
LOG_FILE <- snakemake@log

## Input file
INPUT_FILE <- snakemake@input[["data"]]
METADATA <- snakemake@input[["metadata"]]
METACONTRAST <- snakemake@input[["metacontrast"]]

## Params
ADJ_PVAL <- snakemake@params[["pval"]]
LOG_FC <- snakemake@params[["logFC"]]
COUNT_TOOL <- snakemake@params[["count_tool"]]
FILTERING_METHOD <- snakemake@params[["filtering_method"]]

## Output path
RESULTS_PATH <- snakemake@output

# Import libraries
library("edgeR") # for filtering, call this before DESeq2
library("DESeq2") # for RNA-seq analysis
library("data.table") # for fread function
library("dplyr") # for data transformation
library(gplots) # for plots
library(ggplot2) # for plots
library(ggpubr) # for plots
library(matrixStats) # for scree plot
library("stringr")
library(RColorBrewer)

# Import functions
source("translatome/src/some_functions.R")


#==============================================================================

# Session
## Save session information
log <- file(as.character(LOG_FILE), open = "wt")
sink(log)

## Start time
print("Start_time")
print(Sys.time())


#==============================================================================
# Import part
#==============================================================================

# Import metadata
metadata <- fread(
  METADATA, 
  header = TRUE, 
  data.table = FALSE
)
myconditions <- metadata[order(metadata$sample_name), ]
print("Info text: Contrats for test:")
print(myconditions)

## Combine fraction and treatment in one label
myconditions$combine <- interaction(myconditions$fraction,
  myconditions$treatment,
  myconditions$time,
  sep = "_"
)

# Import count file
## Function(data, count_tool_used)
## Count_tool can be "htseq" or "featurecounts" or "None"
## "None" only if colnames are identic to sample names
data <- ImportFileAsDF(INPUT_FILE, COUNT_TOOL)

## Order by names
OrderData <- data[, order(names(data))]


#==============================================================================
# Prepare count
#==============================================================================

# Create DESeq S4 dataset
dataDESeq <- data.frame(lapply(OrderData, as.integer))
rownames(dataDESeq) <- rownames(OrderData)

# Delete row with NA value or replace by zero
print("Info text: Total number of rows")
print(nrow(dataDESeq))

# Is NA
isna <- is.na(dataDESeq)
test <- table(isna)
## Replace NA with 0 if NA
if (isTRUE(test[2] == TRUE)) {
  print("Info text: We choose to replace NA by zero")
  dataDESeq <- mutate_all(dataDESeq, ~ replace(., is.na(.), 0))
  print("Number of rows after replacing NA by zero:")
  print(nrow(dataDESeq))
  print("Note: should be the same than previous total number of rows")
} else {
  test
}

## Create deseq2 Matrix
dsqMat <- DESeqDataSetFromMatrix(
  countData = dataDESeq,
  colData = myconditions,
  design = ~combine
)
print(summary(dsqMat))


#==============================================================================
# Filtering and normalization of counts
#==============================================================================
# Delete null counts
ddsout <- dsqMat[rowSums(counts(dsqMat)) > 1, ]
print("Info text: Number of rows without empty rows")
print(nrow(ddsout))

# Normalization by factor size
ddsout <- estimateSizeFactors(ddsout, type = "ratio")
png(paste(
  RESULTS_PATH, 
  "sparsityPlot_sizefactor.png",
  sep = "/")
)
plotSparsity(ddsout)
dev.off()

# Filter rows:
## Function(lib_size_norm_data, metadata, filter_method, verbose_mode)
## can be edger, deseq or default
## deseq and edger filter low counts
## i.e. counts > at 10 according to replicates or matrix design respectively
## default keep rows with at least 1 counts per million per replicate groups
ddsout2 <- FilterMethods(
  ddsout, 
  myconditions, 
  FILTERING_METHOD, 
  verbose = TRUE
)

# Normalized dispersion with a negative binomial distribution
ddsout2 <- estimateDispersions(ddsout2, fitType = "parametric")
residual <- mcols(ddsout2)$dispGeneEst - mcols(ddsout2)$dispFit
print(paste0("Min of residuals: ", min(residual)))
print(paste0("Max of residuals: ", max(residual)))
print(paste0("Mean of residuals: ", mean(residual)))
print("Info text: Mean of residuals, better if close to zero")

# Mean of normalized counts (x axis) and dispersion estimate for each gene
png(paste(RESULTS_PATH, "dispersionPlot.png", sep = "/"))
plotDispEsts(ddsout2)
dev.off()

# Extract normalized count table
counts <- counts(ddsout2, normalized = TRUE)
## Export table of normalized counts
write.table(counts,
  file = paste(RESULTS_PATH, "Normalized_counts.tsv", sep = "/"),
  quote = FALSE, row.names = TRUE, sep = "\t"
)


#==============================================================================
# PCA
#==============================================================================
# PCA plot for treatment and replicates
## function(matrix, metadata, design, factor1,
## factor2, filterMethods, fitmodel, output)
tryCatch(
  expr = {
    res_pca <- ComputePCA2cond(
      dataDESeq,
      myconditions,
      ~ treatment + replicate,
      "treatment",
      "replicate",
      "default",
      "parametric",
      RESULTS_PATH
    )
    ### function(res_pca$pcaDATA, res_pca$eig_val, factor1, factor2, output)
    PCAplots2cond(
      res_pca$pcaDATA,
      res_pca$eig_val,
      "treatment",
      "replicate",
      RESULTS_PATH
    )
  },
  error = function(e) {
    message("An unexpected error occurred with constrats for treatment 
            and replicate, design should be equilibrate.")
    print(e)
  }
)

## PCA plot for fraction and replicates
tryCatch(
  expr = {
    res_pca2 <- ComputePCA2cond(
      dataDESeq,
      myconditions,
      ~ fraction + replicate,
      "fraction",
      "replicate",
      "default",
      "parametric",
      RESULTS_PATH
    )
    PCAplots2cond(
      res_pca2$pcaDATA,
      res_pca2$eig_val,
      "fraction",
      "replicate",
      RESULTS_PATH
    )
  },
  error = function(e) {
    message("An unexpected error occurred with constrats for fraction 
            and replicate, design should be equilibrate.")
    print(e)
  }
)

# PCA plot for time and replicates
tryCatch(
  expr = {
    res_pca3 <- ComputePCA2cond(
      dataDESeq,
      myconditions,
      ~ time + replicate,
      "time",
      "replicate",
      "default",
      "parametric",
      RESULTS_PATH
    )
    PCAplots2cond(
      res_pca3$pcaDATA,
      res_pca3$eig_val,
      "time",
      "replicate",
      RESULTS_PATH
    )
  },
  error = function(e) {
    message("An unexpected error occurred with constrats for time 
            and replicate, design should be equilibrate.")
    print(e)
  }
)

# PCA plot for treatment and time
tryCatch(
  expr = {
    res_pca4 <- ComputePCA2cond(
      dataDESeq,
      myconditions,
      ~ treatment + time,
      "treatment",
      "time",
      "default",
      "parametric",
      RESULTS_PATH
    )
    PCAplots2cond(
      res_pca4$pcaDATA,
      res_pca4$eig_val,
      "treatment",
      "time",
      RESULTS_PATH
    )
  },
  error = function(e) {
    message("An unexpected error occurred with constrats for treatment 
            and time, design should be equilibrate.")
    print(e)
  }
)

# PCA plot for fraction and time
tryCatch(
  expr = {
    res_pca5 <- ComputePCA2cond(
      dataDESeq,
      myconditions,
      ~ fraction + time,
      "fraction",
      "time",
      "default",
      "parametric",
      RESULTS_PATH
    )
    PCAplots2cond(
      res_pca5$pcaDATA,
      res_pca5$eig_val,
      "fraction",
      "time",
      RESULTS_PATH
    )
  },
  error = function(e) {
    message("An unexpected error occurred with constrats for fraction 
            and time, design should be equilibrate.")
    print(e)
  }
)

# PCA plot for combine and replicates
tryCatch(
  expr = {
    res_pca6 <- ComputePCA2cond(
      dataDESeq,
      myconditions,
      ~ combine + replicate,
      "combine",
      "replicate",
      "default",
      "parametric",
      RESULTS_PATH
    )

    PCAplots2cond(
      res_pca6$pcaDATA,
      res_pca6$eig_val,
      "combine",
      "replicate",
      RESULTS_PATH
    )
    # Export PCA data
    write.table(
      res_pca6$pcaDATA,
      file = paste0(RESULTS_PATH, "/PCA_values.tsv"),
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )

    # Export all dimension value
    write.table(
      res_pca6$eig_val,
      file = paste0(RESULTS_PATH, "/Scree_dimensions.tsv"),
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )

    # Scree plot as hist
    histplot <- ScreeHist(res_pca6$eig_val)

    png(paste0(RESULTS_PATH, "/ScreeHist_plot.png"))
    plot(histplot)
    dev.off()

    pdf(paste0(RESULTS_PATH, "/ScreeHist_plot.pdf"))
    plot(histplot)
    dev.off()
  },
  error = function(e) {
    message("An unexpected error occurred with constrats for combine 
            and replicate, design should be equilibrate.")
    print(e)
  }
)


#==============================================================================
# Differential expression computation
#==============================================================================

# Differential analysis
contrasts <- DESeq(
  ddsout2, 
  test = c("Wald"),
  fitType = "parametric"
)

# Import metacontrast
metacontrast <- fread(
  METACONTRAST,
  header = TRUE,
  data.table = FALSE
)

# Initialization of the vectors for the summary
acc_tot  <- c()
acc_up   <- c()
acc_down <- c()

# Perform differential expression
for (l in 1:length(metacontrast$fraction)) {
  ## Constructing the contrast matrix
  metacontrast$union1[l] <- paste(
    metacontrast$fraction[l],
    metacontrast$treatment1[l],
    metacontrast$time[l],
    sep = "_"
  )
  print(metacontrast$union1[l])
  metacontrast$union2[l] <- paste(
    metacontrast$fraction_bis[l],
    metacontrast$treatment2[l],
    metacontrast$time_bis[l],
    sep = "_"
  )
  print(metacontrast$union2[l])
  ## Test
  comp_res <- ComparisonToTest(
    contrasts,
    "combine",
    metacontrast$union1[l],
    metacontrast$union2[l],
    LOG_FC,
    ADJ_PVAL,
    RESULTS_PATH
  )
  ### na omit
  comp_res_woNA <- na.omit(comp_res)
  ### filter on adjusted p-value
  filt_comp_res <- comp_res_woNA[comp_res_woNA$padj <= as.numeric(ADJ_PVAL), ]
  ### filter on logFC value
  filtlogup <- filt_comp_res[filt_comp_res$log2FoldChange >= as.numeric(LOG_FC), ]
  filtlogdn <- filt_comp_res[filt_comp_res$log2FoldChange <= -as.numeric(LOG_FC), ]
  ### accumulation of results
  acc_up   <- c(acc_up, nrow(filtlogup))
  acc_down <- c(acc_down, nrow(filtlogdn))
  acc_tot  <- c(acc_tot, nrow(filtlogup) + nrow(filtlogdn))
}

## Summary of results
tab <- data.frame(comparison = interaction(
  metacontrast$union1,
  metacontrast$union2,
  sep = "_vs_"
))
tab$total_present_gene  <- rep(nrow(ddsout), each = length(tab$comparison))
tab$total_filtered_gene <- rep(nrow(ddsout2), each = length(tab$comparison))
tab$significant_genes_total <- acc_tot
tab$significant_genes_UP    <- acc_up
tab$significant_genes_DOWN  <- acc_down
write.table(
  tab,
  file = paste(RESULTS_PATH, "Summary_significant_genes.tsv", sep = "/"),
  row.names = FALSE
)


#==============================================================================

# Export session
save.image(paste0(RESULTS_PATH, "/SessionDESeq2.RData"))

# End time
print("End_time")
print(Sys.time())

# End save info session
sink()
