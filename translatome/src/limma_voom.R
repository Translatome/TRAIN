#!/bin/env Rscript
#'=============================================================================
#'                   Statistical analysis with limma-voom
#'=============================================================================
#'
#' Description:
#' ------------
#' In the limma approach to RNA-seq, read counts are converted to
#' log2-counts-per-million (logCPM) and the mean-variance relationship is
#' modeled using either with precision weights or with an empirical Bayesian
#' prior trend. The precision weights approach is called “voom” and the prior
#' trend approach is called “limma-trend”.
#' In both cases, the RNA-seq data can be analyzed as if it was microarray data.
#' This means that any of the linear modeling or gene set testing methods in the
#' limma package can be applied to RNA-seq data.
#'
#' Documentation:
#' --------------
#' - limma:
#'  https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
#' - edgeR:
#'  https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf
#'
#' Requirements:
#' -------------
#' - edgeR: 3.40.0
#' - ggplot2: 3.4.1
#' - data.table: 1.14.8
#' - ggpubr: 0.6.0
#' - src/some_functions: 0.1.0
#' All these packages are managed by the conda environment file
#'
#' Executions:
#' -----------
#' Execution in snakemake rule:
#'    shell:
#'        "src/limma_voom.R"
#'
#' Execution without snakemake:
#'    Replace all paths provided by snakemake in the snakemake value section
#'
#' Credits:
#' --------
#' - Authors: Céline Mandier and Julie Ripoll
#' - Created: 2023-03-27
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
#' - metacontrast
#'    TSV file & CSV file are supported
#'    please adapt and complete the template provided in ressources
#'    The metacontrast file will allow to inform about the desired contrasts.
#'    The goal is to indicate in the first two columns the fraction and the
#'    treatment that we want as a reference for the comparison.
#'    This one will be compared to the same fraction but with a different
#'    treatment. The fractions may differ from one row to another, but
#'    the treatments themselves differ from column to column.
#'
#' Outputs:
#' --------
#' - Result table:
#'    CSV file, contains:
#'    - logFC: log2 fold change
#'    - AveExpr: The average log2-expression level for that gene across
#'      all the arrays and channels in the experiment
#'    - t: is the empirical Bayes moderated t-statistic
#'    - P.Value: p-value associated with the static test
#'    - adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
#'    - B: log-odds that gene is DE
#' - Summary_significant_gene table
#'    TSV file, number of DE genes in the different comparisons
#' - Mean-variance_trend
#'    PNG file
#' - MDS plot
#'    PNG & PDF files, multidimensional scaling plots
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

## LOG
LOG_FILE          <- snakemake@log

## INPUT
FILE_DATA         <- snakemake@input[["data"]]
FILE_METADATA     <- snakemake@input[["metadata"]]
FILE_METACONTRAST <- snakemake@input[["metacontrast"]]

## PARAMS
THRESH_LOGFC      <- snakemake@params[["logFC"]]
THRESH_PVAL       <- snakemake@params[["pval"]]
COUNT_TOOL        <- snakemake@params[["count_tool"]]

## OUTPUT
PATH_OUTPUT       <- snakemake@output

# Import libraries
library("edgeR") # version 3.40.0
library("ggplot2") # version 3.4.1
library("data.table") # version 1.14.8
library("ggpubr") # version 0.6.0

# Import functions
source("translatome/src/some_functions.R")

#==============================================================================

# Save session information
log <- file(as.character(LOG_FILE), open = "wt")
sink(log)

# Start time
print("Start_time")
print(Sys.time())

#==============================================================================
# Prepare Data
#==============================================================================
# Import count file
data <- ImportFileAsDF(FILE_DATA, COUNT_TOOL)

# Order by names
OrderData <- data[, order(names(data))]


#==============================================================================
# Metadata
#==============================================================================

# Import metadata file
metadata <- fread(
  FILE_METADATA,
  header = TRUE,
  data.table = FALSE
)

# Order by names
Ordermetadata <- metadata[order(metadata$sample_name), ]

#==============================================================================
# Design
#==============================================================================

# Create a new variable combine that combine Fraction and Treatment
combine <- interaction(Ordermetadata$fraction, Ordermetadata$treatment)
cat("Levels : ", levels(combine))

# Create design
replicate <- Ordermetadata$replicate
mm <- model.matrix(~ 0 + combine + replicate)


#==============================================================================
# Filtering and normalization factors
#==============================================================================

# Create DGElist object
DGE_objt <- DGEList(OrderData)

# Filtering of genes with lower expression
## No null expression
DGE_objt_no_null <- DGE_objt[rowSums(DGE_objt$counts) > 1, ]
nb_no_null <- nrow(DGE_objt_no_null)

## Lower expression
keep_gene <- filterByExpr(DGE_objt, design = mm, group = NULL)
DGE_objt_keep <- DGE_objt[keep_gene, , keep.lib.sizes = FALSE]
print(str(DGE_objt_keep))

# Calcul of normalization factors
## CalcNormFactors doesn’t normalize the data,
## it just calculates normalization factors for use downstream.
DGE_objt_factor <- calcNormFactors(DGE_objt_keep, method = "TMM")
print(DGE_objt_factor$samples)

cpms <- cpm(DGE_objt_factor, log=FALSE)
write.table(
  cpms,
  file = paste(PATH_OUTPUT, "Normalized_counts.csv", sep = "/"),
  row.names = TRUE, sep = ";"
)

#==============================================================================
# Multidimensional scaling plot
#==============================================================================

mds <- plotMDS(DGE_objt_factor, col = as.numeric(combine))
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y)

fig <- ggplot(
  toplot,
  aes(
    Dim1,
    Dim2,
    color = as.factor(combine),
    shape = as.factor(replicate)
  )
) +
  labs(
    color = "Factors", 
    shape = "Replicates"
  ) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(mds$var.explained[1] * 100, 2), "% variance")) +
  ylab(paste0("PC2: ", round(mds$var.explained[2] * 100, 2), "% variance")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(aspect.ratio = 1)

png(
  file = paste(PATH_OUTPUT, "/plotMDS_norm_data.png", sep = "/"),
  width = 900,
  height = 650
)
plot(fig)
dev.off()

pdf(file = paste(PATH_OUTPUT, "/plotMDS_norm_data.pdf", sep = "/"))
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
  file = paste(PATH_OUTPUT, "/ScreeHist_plot.png", sep = "/"),
  width = 900,
  height = 650
)
plot(histplot)
dev.off()

pdf(file = paste(PATH_OUTPUT, "/ScreeHist_plot.pdf", sep = "/"))
plot(histplot)
dev.off()


#==============================================================================
# Voom transformation and calculation of variance weights
#==============================================================================

# Voom
png(
  file = paste(PATH_OUTPUT, "/mean-variance_trend.png", sep = "/"),
  width = 900,
  height = 650
)
y <- voom(DGE_objt_factor, mm, plot = TRUE)
dev.off()


#==============================================================================
# Fitting linear models in limma
#==============================================================================

# Linear model established using weighted least squares for each gene
fit <- lmFit(y, mm)

# Renaming of lmfit() table columns to match the metacontrast nomenclature
colnames(fit$coefficients) <- sub(
  "^replicate",
  "",
  colnames(fit$coefficients)
)
colnames(fit$coefficients) <- sub(
  "^combine",
  "",
  colnames(fit$coefficients)
)

colnames(fit$stdev.unscaled) <- sub(
  "^replicate",
  "",
  colnames(fit$stdev.unscaled)
)
colnames(fit$stdev.unscaled) <- sub(
  "^combine",
  "",
  colnames(fit$stdev.unscaled)
)

colnames(fit$design) <- sub(
  "^replicate",
  "",
  colnames(fit$design)
)
colnames(fit$design) <- sub(
  "^combine",
  "",
  colnames(fit$design)
)

colnames(fit$cov.coefficients) <- sub(
  "^replicate",
  "",
  colnames(fit$cov.coefficients)
)
colnames(fit$cov.coefficients) <- sub(
  "^combine",
  "",
  colnames(fit$cov.coefficients)
)
rownames(fit$cov.coefficients) <- sub(
  "^replicate",
  "",
  rownames(fit$cov.coefficients)
)
rownames(fit$cov.coefficients) <- sub(
  "^combine",
  "",
  rownames(fit$cov.coefficients)
)

print(head(coef(fit)))


#==============================================================================
# Contrast
#==============================================================================

# Import metacontrast
metacontrast <- fread(
  FILE_METACONTRAST,
  header = TRUE,
  data.table = FALSE
)

# Initialization of the vectors for the summary
lin_tot  <- c()
lin_up   <- c()
lin_down <- c()

# Construction of the required contrast
lvl <- colnames(coef(fit))
for (l in 1:length(metacontrast$fraction)) {
  ## Construction of the contrast matrix for each contrast provided 
  ## in the metacontrast file
  metacontrast$union1[l] <- paste(
    metacontrast$fraction[l], ".",
    metacontrast$treatment1[l],
    sep = ""
  )
  ### Creation of combination fraction1 treatment
  metacontrast$union2[l] <- paste(
    metacontrast$fraction_bis[l], ".",
    metacontrast$treatment2[l],
    sep = ""
  )
  ### Creation of combination fraction2 treatment
  vec <- c()
  for (i in lvl) {
    if (i == metacontrast$union1[l]) {
      vec <- c(vec, 1)
      ## Assign the value 1 for the first combination
    } else if (i == metacontrast$union2[l]) {
      vec <- c(vec, -1)
      ## Assign the value -1 for the second combination
    } else {
      vec <- c(vec, 0)
      ## Assign a value of 0 for levels not present in the requested contrast
    }
  }
  ## Matrix produced
  ### Equivalent to that produced by the function 
  ### makeContrasts(combination1 - combination2, levels = lvl))
  mat <- matrix(
    vec,
    dimnames = list(
      colnames(coef((fit))),
      paste0(
        metacontrast$union1[l], " - ",
        metacontrast$union2[l], sep = ""
      )
    )
  )

  # Estimate contrast for each gene
  tmp <- contrasts.fit(fit, mat)

  # Empirical Bayes smoothing of standard errors
  tmp_ebayes <- eBayes(tmp)

  # Formatting the table
  top_table <- topTable(tmp_ebayes, sort.by = "P", n = Inf)
  identifier <- row.names(top_table)
  top_table <- cbind(identifier, top_table)

  # Writing the result for the requested contrast
  write.table(
    top_table, 
    file = paste(PATH_OUTPUT, "/",
                 metacontrast$fraction[l], "_",
                 metacontrast$treatment1[l], "_vs_",
                 metacontrast$fraction_bis[l], "_",
                 metacontrast$treatment2[l], ".tsv",
                 sep = ""
    ), 
    sep = "\t",
    dec = ".",
    quote = FALSE,
    row.names = FALSE
  )

  # Number of significant genes
  lgFC_data <- subset(
    top_table,
    logFC <= -as.numeric(THRESH_LOGFC) | logFC >= as.numeric(THRESH_LOGFC)
  )
  signif_data <- subset(
    lgFC_data,
    adj.P.Val <= as.numeric(THRESH_PVAL)
  )
  lgFC_UP <- subset(
    signif_data,
    logFC >= as.numeric(THRESH_LOGFC)
  )
  lgFC_DOWN <- subset(
    signif_data,
    logFC <= -as.numeric(THRESH_LOGFC)
  )

  lin_up   <- c(lin_up, nrow(lgFC_UP))
  lin_down <- c(lin_down, nrow(lgFC_DOWN))
  lin_tot  <- c(lin_tot, nrow(signif_data))

  print(metacontrast$union1[l])
  print(metacontrast$union2[l])
  print(nrow(signif_data))
}

#==============================================================================
# Summary of the number of significant genes
#==============================================================================

tab <- data.frame(comparison = interaction(
  metacontrast$union1,
  metacontrast$union2,
  sep = "_vs_")
)
tab$total_present_gene <- rep(
  nb_no_null,
  each = length(tab$comparison)
)
tab$total_filtered_gene <- rep(
  nrow(DGE_objt_keep),
  each = length(tab$comparison)
)
tab$significant_genes_total <- lin_tot
tab$significant_genes_UP    <- lin_up
tab$significant_genes_DOWN  <- lin_down
print(tab)

# Writing the summary of the number of significant genes
write.table(
  tab,
  file = paste(PATH_OUTPUT, "/Summary_significant_genes.tsv", sep = "/"),
  row.names = FALSE
)

#==============================================================================

# Export session
save.image(paste0(PATH_OUTPUT, "/SessionLimma.RData"))

# End time
print("End_time")
print(Sys.time())

# End save info session
sink()
