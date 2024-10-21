Configuration
=============

This workflow is managed by a global configuration defined in a yaml configuration file. This file contains several informations required to run the workflow.

Several example configuration files are available in the *configs/* folder:
    * config_PE.yml
    * config_PE_kinetic.yml
    * config_SE.yml
    * config_SE_kinetic.yml

The configuration options will change depending on whether the data is kinetic or not, and whether the data is paired end or single end.

These configuration files are divided into several sections, where each option available to the user is specified. Before each run of the workflow, the user must check that the options are correct according to their project.

The first option in this workflow is the choice of analysis type: pol-seq (classical POL-seq analysis combined with RNA-seq) or kinetic (kinetic POL-seq analysis combined with kinetic RNA-seq).


Data information
----------------

**Path and and sequencing information**

The user should input the following primary paths in the analysis:
    * "inpfol": the pathway to the raw sequencing data
    * "resfol": the pathway to the directory where the analysis outcomes will be saved
    * "ref_fol": the pathway for the reference files that will be obtained through the *00_download_references.smk* program.

We suggest users avoid storing their data and results in the same directory as the pipeline to enable easier re-use.

The filepaths for metadata and metacontrast must be specified in the "metadata" and "metacontrast" parameters and customized to the user's project and data. For information on adapting the files, refer to the `readme <https://github.com/Translatome/TRAIN/blob/main/README.md>`_.

The user should provide information about their samples and assign "YES" to the "paired" parameter if the data is paired end, otherwise "NO".

The number of replicates needs to be entered through the "replicate" parameter.


**Raw filenames**

Next steps: complete the "extension_files" parameter to determine whether the data is compressed or not.

In cases of paired end data, identify the read 1 and read 2 files through the details provided in the "paired_info_1" and "paired_info_2" parameters. The annotation should follow this format: "_1" for the read 1 file, and "_2" for the read 2 file.

If the sample name contains any unique characteristics between the paired-end information and file extension, use the "end_filename" parameter to include this information (e.g. _001.fastq.gz). Otherwise, simply append the correct file extension ('.fastq.gz' or '.fq', etc).


References files
----------------

The reference files section includes all necessary information for downloading genome or transcriptome, and annotation versions for performing analysis. These references are obtained from the ENSEMBL database.

Fill in the following details:
    * "organism": enter the ENSEMBL ID of the species (e.g. homo_sapiens, mus_musculus, ...) in lowercase and without spaces.
    * "organismcode": enter the associated taxonomic code (e.g. 9606, 10090, ...)
    * "release": enter the desired reference version (e.g. 107)
    * "assembly_reference": for downloading reference genome or transcriptome fasta file according to arguments **"dna"** or **"cdna"** respectively.
    * "extensionAnnotation": extension of annotation, available options are **"gtf"** or **"gff"**.


After downloading the reference files using the *00_download_references.smk* snakefile, the user must specify the names of the downloaded genome and annotation files.

Example for human:
    * "assembly_ref": **"Homo_sapiens.GRCh38.dna.primary_assembly.fa"**
    * "annotation": **"Homo_sapiens.GRCh38.107.gtf"**


Available tools
---------------

This pipeline enables users to select specific tools for each step, while other tools remain unchangeable.

Users have the option to specify their tool preference for trimming, counting, and statistical analysis. Changes to tool selections can be made by re-running the snakefile after modifying the configuration file. Output results are saved in separate directories.

It is important to note that if these selections are altered, only the tool defined in the configuration file will be employed for subsequent pipeline steps.

Enrichment can be included or excluded based on the user's request. Selecting "None" will generate a report after the statistical section instead of after the enrichment section.

Available option for tools:
    * Trimming: "cutadapt" / "fastp" 
    * Counting: "featurecounts" / "htseq"
    * Statistics: "deseq2" / "limma"
    * Enrichment: "gprofiler" / "None"


Pipeline parameters
-------------------

The final section of the configuration file pertains to the parameters required for the analysis.

**Threads**

To indicate the number of threads needed, use the "threads" parameter.
Please be aware that the available number of threads may vary depending on the computer or server.

**Trimming**

For the trimming section, there are several available parameters based on the tool used.

Both counting tools require the "quality" and "minilength" parameters to indicate the minimum acceptable Phred score and the minimum length of a read after trimming, respectively.

For Cutadapt, you can use the following parameters:
    * "fref": to specify the fasta file that contains the sequencing adapters for their detection and removal. If there are unknown adapters, we offer a file containing all known adapters in the "ressources" folder.
    * "cutphred": to set the minimum quality value. If one value is given, only the 3' end will be trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff and the 3' end with the second.
    * "cutlength1": to cut the beginning (if value is positive) or tail (if value is negative) of first read 
    * "cutlength2": to cut the beginning (if value is positive) or tail (if value is negative) of second read 
    * "morelength1": if you want to cut at both ends of read 1 then use this option in addition to the cutlength1 option to cut either at the beginning or the end depending on what has been chosen in the cutlength1 option.
    * "morelength2": if you want to cut at both ends of read 2 then use this option in addition to the cutlength2 option to cut either at the beginning or the end depending on what has been chosen in the cutlength2 option.

Please be aware that for "cutlength1" and "morelength1", if the length is positive it removes bases from the beginning. If length is negative, it removes bases from the end.

For fastp, the available parameters are:
    * "trimfront1": to cut the beginning of first read 
    * "trimfront2": to cut the beginning of second read  if paired end data
    * "trimtail1": to cut the tail of first read 
    * "trimtail2": to cut the tail of second read  if paired end data


**Mapping**

The "rlength" parameter is required when indexing the reference *genome/transcriptome/* with STAR.

This parameter must be equivalent to the length of read cuts minus one.

**Counts**

There is a parameter called "stranded" that can be used for counting with featurecounts or HTSeq-count.

RNA-Seq libraries can be either stranded or unstranded. The way the libraries are prepared affects the data generated from next generation sequencing (NGS) and how it is interpreted. Stranded RNA-Seq (also known as strand-specific or directional RNA-Seq) allows you to determine the orientation of the transcript, while this information is not available in unstranded or standard RNA-Seq.
If sequences from read 1 align with the RNA strand, the library is considered "stranded". If sequences from read 2 align with the RNA strand, the library is considered "reverse stranded".

This is important as Feature Counts uses a scale of "0" (unstranded), "1" (stranded), and "2" (reverse stranded).

For HTSeq Count, the options are "yes", "no", or "reverse", with "reverse" indicating a "yes" with a reversed strand interpretation.


**Statistics**

The filtering method is a parameter for DESeq2 analysis, allowing for "deseq", "edger" or "default" options:
    * "edger": uses the filterByExpr() function from the edgeR R package, as detailed in the `documentation <https://rdrr.io/bioc/edgeR/man/filterByExpr.html>`_. 
    * "deseq":  the raw counts are filtered for all pairs of replicates that are less than 10, as recommended by Love `documentation <https://support.bioconductor.org/p/110833/>`_. 
    * "default": the standard approach involves less strict filtering based on raw counts, where each pair of replicates should be less than 1.

The parameters "pval" and "logFC" establish the selection thresholds for p-value and log2FoldChange in the statistical analysis section.

**Enrichment**

Enrichment parameters are specific to gProfiler2. Gene Ontology sources (GO:BP, GO:CC or GO:MF), pathways (KEGG, REAC, WP) and CORUM complexes can be entered.

To determine significant terms, a p-value correction is applied. Therefore, the "correction" parameter requires the user to specify the desired correction type from "g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS" or "analytical", and to specify the p-value threshold in the "pval_threshold" parameter.

The "term-nb" parameter indicates the number of terms needed to display enrichment dotplots.

**Final report**

The required parameters for generating an automated report through Rmarkdown to communicate the results of the statistical and enrichment sections are as follows:
    * "context": users can provide the context of their analysis
    * "comments": users can include any desired comments
    * "acknowledgements": users can acknowledge their partners
    * "author": add author names
    * "filename": report filename, modify as per the project
