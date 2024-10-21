..  _Rdocs:

R scripts 
=========

deseq2_kinetic.R
----------------

Statistical analyses with DESeq2

Description
~~~~~~~~~~~
This script allows you to perform DESeq2 differential expression analysis on RNA-seq data.

Requirements
~~~~~~~~~~~~

    * deseq2: 1.30.1
    * dplyr: 1.1.1
    * gplots: 3.1.3
    * ggplot2: 3.4.1
    * ggpubr: 0.6.0
    * matrixStats: 0.63.0
    * stringr: 1.5.0
    * RColorBrewer: 1.1_3
    * edgeR: 3.40.0
    * data.table: 1.14.8
    * src/some_functions: 0.1.0

All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    shell:
        "src/deseq2.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


    Input
        * Count table
            TSV file, output from HTSeq-counts or featureCounts
            CSV file is also supported
        * metadata
            TSV file & CSV file are supported
            Please adapt and complete the template provided in ressources
            Contains: sample name, fraction, treatment, time, replicate
        * metacontrast
            TSV file & CSV file are supported
            Please adapt and complete the template provided in ressources
            Contains: contrasts to use for comparison

    Output
        * normalized counts
            TSV file with normalized counts
        * dispersion plots
            PNG files dispersion and trend of the normalization model
        * dendrogram plots
            PNG files euclidean distance between samples
        * PCA and scree histograms
            * PNG & PDF files
                Scree histograms present a red line at 80% (Pareto's threshold)
            * TSV files, PCA axes values and cumulative variance for scree
        * pval plots
            PNG files, histogram of the p-value repartition for the comparison
        * summary of comparisons
            TSV file, number of DE genes in the different comparisons
        * deseq2 table results by comparisons
            TSV file, named with the comparison
                Contains: 
                    * baseMean mean value of all samples
                    * log2FoldChange value
                    * lfcSE standard error of log2FC
                    * stat value of the Wald test
                    * pvalue p-value of the Wald test
                    * padj adjusted p-value of the Wald test
                    * identifier of the gene, Ensembl identifier
        * session information
            RData archive

    Log
        * log file
            LOG catch terminal stdout



deseq2.R
--------

Statistical analyses with DESeq2

Description
~~~~~~~~~~~
This script allows you to perform DESeq2 differential expression analysis on RNA-seq data.

Requirements
~~~~~~~~~~~~

    * deseq2: 1.30.1
    * dplyr: 1.1.1
    * gplots: 3.1.3
    * ggplot2: 3.4.1
    * ggpubr: 0.6.0
    * matrixStats: 0.63.0
    * stringr: 1.5.0
    * RColorBrewer: 1.1_3
    * edgeR: 3.40.0
    * data.table: 1.14.8
    * src/some_functions: 0.1.0

All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    shell:
        "src/deseq2.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Count table
        TSV file, output from HTSeq-counts or featureCounts
        CSV file is also supported
    * metadata
        TSV file & CSV file are supported
        Please adapt and complete the template provided in ressources
        Contains: sample name, fraction, treatment, replicate
    * metacontrast
        TSV file & CSV file are supported
        Please adapt and complete the template provided in ressources
        Contains: contrasts to use for comparison

Output
~~~~~~

    * normalized counts
        TSV file with normalized counts
    * dispersion plots
        PNG files dispersion and trend of the normalization model
    * dendrogram plots
        PNG files euclidean distance between samples
    * PCA and scree histograms
        * PNG & PDF files
            Scree histograms present a red line at 80% (Pareto's threshold)
        * TSV files, PCA axes values and cumulative variance for scree
    * pval plots
        PNG files, histogram of the p-value repartition for the comparison
    * summary of comparisons
        TSV file, number of DE genes in the different comparisons
    * deseq2 table results by comparisons
        TSV file, named with the comparison
            Contains:
                * baseMean mean value of all samples
                * log2FoldChange value
                * lfcSE standard error of log2FC
                * stat value of the Wald test
                * pvalue p-value of the Wald test
                * padj adjusted p-value of the Wald test
                * identifier of the gene, Ensembl identifier
    * session information
        RData archive

Log
~~~

    * log file
        LOG catch terminal stdout


explore_counts.R
----------------

Multidimensional scaling plot for raw data

Description
~~~~~~~~~~~
This script allows you to plot a multidimensional scaling graph on raw data.
Depending on the results, users can choose between the deseq2 and limma R packages for their statistical analysis.
We use limma when replicates present an important bias.

Requirements
~~~~~~~~~~~~

    * data.table: 1.14.8
    * edgeR: 3.38.1
    * ggplot2: 3.4.1
    * ggpubr: 0.6.0
    * src/some_functions: 0.1.0

All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    shell:
        "src/explore_counts.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Count table
        TSV file & CSV file are supported
        Output from HTSeq-counts or featureCounts
    * metadata
        TSV file & CSV file are supported
        Please adapt and complete the template provided in ressources
        Contains: sample names, fraction, treatment, replicate

Output
~~~~~~

    * MDS plot
        PNG & PDF files, multidimensional scaling plots
    * session information
        RData archive


gprofiler.R
-----------

Analysis with gProfiler2

Description
~~~~~~~~~~~
This script allows you to call the API of gprofileR for gene ontology and pathway enrichments.

Requirements
~~~~~~~~~~~~

    * r-base: 4.1.3
    * data.table: 1.14.8
    * ggplot2: 3.3.6
    * gprofiler2: 0.2.1
    * src/some_functions.R: 0.1.0

All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    shell:
        "src/gprofiler.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Count table
        TSV file, output from explore_results.py
        Contains all counts for transcriptome and translatome

Output
~~~~~~

    * results tables
        TSV files
        One file by sources, i.e. pathway, gene ontology and protein complex, by category of transcription / translation / both / divergent
    * Manhattan enrichment plots
        PNG files
        Top terms selected by adjusted p.value and user threshold in a manhattan plot. One plot by sources and category.
    * dotplots
        PNG files
        Top terms selected by term size and user threshold in a dotplot plot.
        One plot by sources and category.
    * session information
        RData archive

Log
~~~

    * log file
        LOG catch terminal stdout


limma_voom_kinetic.R
--------------------

Kinetic limma-voom statistical analysis

Description
~~~~~~~~~~~
In the limma approach to RNA-seq, read counts are converted to log2-counts-per-million (logCPM) and the mean-variance relationship is modeled using either with precision weights or with an empirical Bayesian prior trend. The precision weights approach is called “voom” and the prior trend approach is called “limma-trend”. In both cases, the RNA-seq data can be analyzed as if it was microarray data. This means that any of the linear modeling or gene set testing methods in the limma package can be applied to RNA-seq data.

Documentation
~~~~~~~~~~~~~

    * limma: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
    * edgeR: https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf

Requirements
~~~~~~~~~~~~

    * edgeR: 3.40.0
    * ggplot2: 3.4.1
    * data.table: 1.14.8
    * ggpubr: 0.6.0
    * src/some_functions: 0.1.0

All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    shell:
        "src/limma_voom_kinetic.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Count table
        TSV file & CSV file are supported
        Output from HTSeq-counts or featureCounts
    * metadata
        TSV file & CSV file are supported
        Please adapt and complete the template provided in ressources
        Contains: sample names, fraction, treatment, time, replicate
    * metacontrast
        TSV file & CSV file are supported
        Please adapt and complete the template provided in ressources
        Contains: contrasts to use for comparison

Output
~~~~~~

    * Results table
        CSV files, contains:
            * logFC: log2 fold change
            * AveExpr: The average log2-expression level for that gene across all the arrays and channels in the experiment
            * t: is the empirical Bayes moderated t-statistic
            * P.Value: p-value associated with the static test
            * adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
            * B: log-odds that gene is DE

    * Summary_significant_gene table
        TSV file, number of DE genes in the different comparisons
    * Mean-variance_trend
        PNG file
    * MDS plot
        PNG & PDF files, multidimensional scaling plots
    * session information
        RData archive

Log
~~~

    * log file
        LOG catch terminal stdout


limma_voom.R
------------

Statistical analysis with limma-voom

Description
~~~~~~~~~~~
In the limma approach to RNA-seq, read counts are converted to log2-counts-per-million (logCPM) and the mean-variance relationship is modeled using either with precision weights or with an empirical Bayesian prior trend. The precision weights approach is called “voom” and the prior trend approach is called “limma-trend”. In both cases, the RNA-seq data can be analyzed as if it was microarray data. This means that any of the linear modeling or gene set testing methods in the limma package can be applied to RNA-seq data.

Documentation
~~~~~~~~~~~~~~

    * limma: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf
    * edgeR: https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf

Requirements
~~~~~~~~~~~~

    * edgeR: 3.40.0
    * ggplot2: 3.4.1
    * data.table: 1.14.8
    * ggpubr: 0.6.0
    * src/some_functions: 0.1.0

All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    shell:
        "src/limma_voom.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Count table
        TSV file & CSV file are supported
        Output from HTSeq-counts or featureCounts
    * metadata
        TSV file & CSV file are supported
        Please adapt and complete the template provided in ressources
        Contains: sample names, fraction, treatment, replicate
    * metacontrast
        TSV file & CSV file are supported
        Please adapt and complete the template provided in ressources
        Contains: contrasts to use for comparison

Output
~~~~~~

    * Results table
        CSV file, contains:
            * logFC: log2 fold change
            * AveExpr: The average log2-expression level for that gene across all the arrays and channels in the experiment
            * t: is the empirical Bayes moderated t-statistic
            * P.Value: p-value associated with the static test
            * adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
            * B: log-odds that gene is DE

    * Summary_significant_gene table
        TSV file, number of DE genes in the different comparisons
    * Mean-variance_trend
        PNG file
    * MDS plot
        PNG & PDF files, multidimensional scaling plots
    * session information
        RData archive

Log
~~~

    * log file
        LOG catch terminal stdout


log-log_plots.R
---------------
Log-Log plots R 

Description
~~~~~~~~~~~
This script allows you to create a log-log plot.

Requirements
~~~~~~~~~~~~

    * dplyr: 1.1.1
    * gplots: 3.1.3
    * ggplot2: 3.4.1
    * gghighlight: 0.4.0
    * tidyr: 1.3.0
    * repr: 1.1.6

All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    shell:
        "src/log-log_plots.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Count table
        TSV file, output from explore_results.py
        Contains all counts for transcriptome and translatome

Output
~~~~~~

    * density plots
        PNG & PDF files
        Density curves according to density of genes
    * log-log plots
        PNG & PDF files
        Scatter-plot of logFC of transcriptome in x and translatome in y
    * session information
        RData archive

Log
~~~

    * log file
        LOG catch terminal stdout


summary_goatools.R
------------------

Summary of goatools

Description
~~~~~~~~~~~
Summary of the number of terms according to the depth or level of hierarchy information for the Gene Ontologies and the categories and sources GO: BP, MF & CC.

Requirements
~~~~~~~~~~~~
All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    script:
        "src/summary_goatools.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Pathways enrichment folder
        Includes all result files for GO: BP, MF & CC and category

Output
~~~~~~

    * goatool summary depth
        CSV file
        Number of terms according to the depth of the hierarchy, categories and source GO: BP, MF and CC
    * goatool summary level
        CSV file
        Number of terms according to the level of the hierarchy, categories and source GO: BP, MF and CC
    * goatool summary depth graph
        PNG file
        Number of terms according to the depth of the hierarchy, categories and source GO: BP, MF and CC
    * goatool summary level graph
        PNG file
        Number of terms according to the level of the hierarchy, categories and source GO: BP, MF and CC
    * Session info
        RData file
        Collect Information About the Current R Session

Log
~~~
    * log file
        LOG catch terminal stdout


summary_gprofiler.R
-------------------

Summary of gprofileR

Description
~~~~~~~~~~~
The purpose of this script is to make a summary of the number of terms for each gprofiler source (KEGG, REAC, WP, 'GO:BP', 'GO:MF', 'GO:CC' and CORUM) for
each category file (Both_mRNA_DOWN, Both_mRNA_UP, Divergent_DNUP, Divergent_UPDN, Transcription_DOWN, Transcription_UP, Translation_DOWN, Translation_UP).

Requirements
~~~~~~~~~~~~
All these packages are managed by the conda environment file

Execution
~~~~~~~~~
Execution in snakemake rule:

::

    script:
        "src/summary_gprofiler.R"

Execution without snakemake:
Replace all paths provided by snakemake in the snakemake value section


Input
~~~~~

    * Pathways enrichment folder
        Includes all result files for each source and category

Output
~~~~~~

    * gProfiler summary
        CSV file
        Number of terms for each category by source
    * Session info
        RData file
        Collect Information About the Current R Session

Log
~~~

    * log file
        LOG catch terminal stdout

