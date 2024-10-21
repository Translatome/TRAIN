import os
import sys
import subprocess
from snakemake.utils import validate

from wrappers.wrapper_system import inject, val_mappy
from wrappers.multiqc import create_multiqc
from wrappers.featureCounts import create_featurecounts
from wrappers.htseq_count import create_htseqcount


# Input folder
INPFOL = config["inpfol"] # specify path of inputs
RESFOL = config["resfol"] # specify pathway for outputs
REFFOL = config["ref_fol"] # specify path to references files

# params
ANALYSIS = config["analysis"]
COUNT_TOOL = config["tools"]["count"]
RELEASE = config['release'] # release of the organism specified in config yaml file
PAIRED = config["paired"] # paired-end or single-end data
THREADS = config["threads"] # low: 4 CPUs, medium: 8, high: 16, Ultra: 32
ORGANISM = config["organism"] # organism specified in config yaml file
ANNOT = config["references"]["annotation"]
FORMATREF = config['extensionAnnotation'] # extension of annotation file specified in config yaml file
STRANDED = config["stranded"]
METADATA = config["metadata"] # specify path of metadata file
REPLICATE = config["replicat"]

# Wildcards
PAIRED_INFO_1 = config["paired_info_1"]
PAIRED_INFO_2 = config["paired_info_2"]
END = config["end_filename"]
SAMPLE, = glob_wildcards(INPFOL+"{sample}"+(PAIRED_INFO_1 if PAIRED == "YES" else "")+END)


#######################################################################

rule all:
    """
    Count reads after mapping
    =========================

    Requirements
    ------------
    - snakemake with python 3.10 min. & conda / miniconda or mamba
    - featureCount (subread) / HTSeq-count
    - pandas
    - multiQC
    - r ggplot2
    Except snakemake & conda, programs can be installed using environment yaml files

    Execution
    ---------
    ``
    snakemake 
        -s [name_snakefile]
        --use-conda [use conda env specified in script&config]
        --conda-frontend mamba [use mamba for faster deployment]
        --configfile [path to config file]
        -j [Number of cores]
        -p [Print rules script]
        -r [Print the reason for each executed rule]
        -k [Go on with independent jobs if a job fails]
        -n [Check integrity: dry mode before execution]
    ``

    Examples
    --------
    ``
    >snakemake -s translatome/workflows/04_counts.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 1 -p -r -k -n
    
    # print doc
    >snakemake -s translatome/workflows/04_counts.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
    ``

    """
    input:
        counts = RESFOL+"Counts_"+COUNT_TOOL+"/All_counts.tsv",
        mdsplot = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_raw_data.png"


#######################################################################

print("------------------------------------------------------")
print("Pipeline will be launch for:", ANALYSIS, " analysis")
print("Tool choosen for counting:", COUNT_TOOL)
print("Available tools: featurecounts or htseq")
print("------------------------------------------------------")

#######################################################################

match COUNT_TOOL:
    case "featurecounts":
        if STRANDED not in ["0", "1", "2"]:
            print("""
            WARNING: change the declaration of 'stranded' value in the configuration file.
            Feature Counts takes: '0', '1' or '2'.
            """)
        else:
            rule FeatureCounts:
                """
                Aim of the rule: Count exon in read
                    Why exon ? See: https://www.biostars.org/p/61042/
                    A typical RNS-seq experiment measures counts over exons and exon junctions (if there are any).
                    From that we may expand to transcript/gene models but that process will involve approximations and potential inconsistencies with respect of the real phenomena.
                    The transcript is fragmented into small pieces and these align to exons or junctions.
                    So that is what we can quantify and not the entire transcript.
                """
                message: 
                    """--- featureCounts from subread counts of {wildcards.sample} ---"""
                conda:
                    "../envs/post-alignment/subread.yml"
                log:
                    RESFOL+"Counts_"+COUNT_TOOL+"/{sample}.log"
                input:
                    RESFOL+"Mapping/STAR/Samtools/{sample}.bam"
                params:
                    paired = PAIRED,
                    annotation = REFFOL+RELEASE+"/"+ORGANISM+"/annotation/"+ANNOT,
                    threads = THREADS,
                    formatREF = str(FORMATREF.upper()), # uppercase required GTF or GFF or SAF
                    featuretype = "exon",
                    strand = STRANDED,
                output:
                    RESFOL+"Counts_"+COUNT_TOOL+"/{sample}_counts.tsv"
                shell:
                    inject(create_featurecounts, locals())


            rule MultiQC_FeatureCounts:
                """
                Aim of the rule: Agregate results from featureCounts into a single report
                See: MultiQC wrapper documentation for others parameters
                """
                message:
                    """--- MultiQC on FeatureCounts log files ---"""
                conda:
                    "../envs/pre-alignment/multiqc.yml"
                input:
                    expand(RESFOL+"Counts_"+COUNT_TOOL+"/{sample}_counts.tsv", sample = SAMPLE)
                params:
                    module = "'featureCounts'",
                    title = "'Report on featureCounts outputs'",
                    filename = "'report_MQC_Counts.html'",
                    dirs = RESFOL+"Counts_"+COUNT_TOOL+"/",
                    fullnames = True,
                    viewtags = False,
                    zipdatadir = True,
                    export = True,
                    verbose = True,
                    interactive = True
                output:
                    outQC = directory(RESFOL+"MQC_reports/MQC_FCounts"),
                    done = touch(RESFOL+"MQC_reports/MQC_FCounts.done")
                shell:
                    inject(create_multiqc, locals())


            rule MergeFeatureCounts:
                """
                Aim of the rule: merge all individual counts in one file
                """
                message: 
                    """---merge counts from featureCounts---"""
                conda:
                    "../envs/stats/pandas.yml"
                input:
                    check = RESFOL+"MQC_reports/MQC_FCounts.done"
                params:
                    inputdir = RESFOL+"Counts_"+COUNT_TOOL+"/",
                    annotation = REFFOL+RELEASE+"/"+ORGANISM+"/annotation/"+ANNOT
                output:
                    counts = RESFOL+"Counts_"+COUNT_TOOL+"/All_counts.tsv",
                    heatmap = RESFOL+"explore_counts_"+COUNT_TOOL+"/Raw_correlation_matrix.png"
                script:
                    "../src/merge_counts_feature.py"


    case "htseq":
        if STRANDED not in ["yes", "no", "reverse"]:
            print("""
            WARNING: change the declaration of 'stranded' value in the configuration file.
            HTSeq Counts takes: 'yes', 'no' or 'reverse'.
            """)
        else:
            rule HTSEQ:
                """
                Aim of the rule: Count exon in read
                    Why exon ? See: https://www.biostars.org/p/61042/
                    A typical RNS-seq experiment measures counts over exons and exon junctions (if there are any).
                    From that we may expand to transcript/gene models but that process will involve approximations and 
                    potential inconsistencies with respect of the real phenomena.
                    The transcript is fragmented into small pieces and these align to exons or junctions.
                    So that is what we can quantify and not the entire transcript.
                """
                message: 
                    """---HTSEQ Count on {wildcards.sample}---"""
                conda:
                    "../envs/post-alignment/htseq.yml"
                log:
                    RESFOL+"Counts_"+COUNT_TOOL+"/{sample}.log"
                input:
                    RESFOL+"Mapping/STAR/Samtools/{sample}_sorted.bam"
                params:
                    threads = THREADS,
                    annotation = REFFOL+RELEASE+"/"+ORGANISM+"/annotation/"+ANNOT,
                    stranded = STRANDED,
                    featuretype = "exon",
                    mode = "intersection-nonempty",
                    supplementaryalignments = "ignore"
                output:
                    RESFOL+"Counts_"+COUNT_TOOL+"/{sample}_counts.tsv"
                shell:
                    inject(create_htseqcount, locals())


            rule MultiQC_HTSeq:
                """
                Aim of the rule: Agregate results from HTSeqcounts into a single report
                See: MultiQC wrapper documentation for others parameters
                    WARNINGS: REQUIRED OPTION -a e.g. *__too_low_aQual
                    CAN ALSO BE PERFORMED ON _BasedOnGene reports
                """
                message:
                    """--- MultiQC on HTSeq log files ---"""
                conda:
                    "../envs/pre-alignment/multiqc.yml"
                log:
                    RESFOL+'MQC_reports/MQC_HTSeq/MQC_HTSeq.done'
                input:
                    expand(RESFOL+"Counts_"+COUNT_TOOL+"/{sample}_counts.tsv", sample = SAMPLE)
                params:
                    module = "'htseq'",
                    title = "'Report on htseq outputs'",
                    filename = "'report_MQC_Counts.html'",
                    dirs = RESFOL+"Counts_"+COUNT_TOOL+"/",
                    fullnames = True,
                    viewtags = False,
                    zipdatadir = True,
                    export = True,
                    verbose = True,
                    interactive = True
                output:
                    outQC = directory(RESFOL+'MQC_reports/MQC_HTSeq'),
                    done = touch(RESFOL+'MQC_reports/MQC_HTSeq.done')
                shell:
                    inject(create_multiqc, locals())


            rule MergeHTSeqcounts:
                """
                Aim of the rule: merge all individual counts in one file
                """
                message: 
                    """---merge counts from HTSeq counts---"""
                conda:
                    "../envs/stats/pandas.yml"
                input:
                    check = RESFOL+'MQC_reports/MQC_HTSeq.done'
                params:
                    inputdir = RESFOL+"Counts_"+COUNT_TOOL+"/",
                    annotation = REFFOL+RELEASE+"/"+ORGANISM+"/annotation/"+ANNOT
                output:
                    counts = RESFOL+"Counts_"+COUNT_TOOL+"/All_counts.tsv",
                    heatmap = RESFOL+"explore_counts_"+COUNT_TOOL+"/Raw_correlation_matrix.png"
                script:
                    "../src/merge_counts_htseq.py"


if int(REPLICATE) < 3:
    print("""
    WARNING: replicates should be at least >= 3, you cannot continue with this pipeline.
    """)
else:
    rule MDS_on_raw_data:
        """
        Aim of the rule: Plot a Multi-dimensional scaling plot MDS on raw data before statistical analysis
        """
        message:
            """--- Perform MDS on raw data ---"""
        conda: 
            "../envs/stats/limma.yml"
        input:
            data = RESFOL+"Counts_"+COUNT_TOOL+"/All_counts.tsv",
            metadata = METADATA
        params:
            analysis = ANALYSIS,
            count_tools = COUNT_TOOL,
            replicat = REPLICATE
        output:
            mdsplot = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_raw_data.png",
            mdsplotpdf = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_raw_data.pdf",
            fractionplot = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_fraction_raw_data.png",
            fractionplotpdf = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_fraction_raw_data.pdf",
            treatmentplot = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_treatment_raw_data.png",
            treatmentplotpdf = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_treatment_raw_data.pdf",
            timeplot = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_time_raw_data.png" if ANALYSIS == "kinetic" else [],
            timeplotpdf = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_time_raw_data.pdf" if ANALYSIS == "kinetic" else [],
            screeplot = RESFOL+"explore_counts_"+COUNT_TOOL+"/ScreeHist_raw_data.png",
            screeplotpdf = RESFOL+"explore_counts_"+COUNT_TOOL+"/ScreeHist_raw_data.pdf",
            sessioninfo = RESFOL+"explore_counts_"+COUNT_TOOL+"/SessionMDSrawData.RData"
        script:
            "../src/explore_counts.R"


#######################################################################
# Author: Julie Ripoll
# Contact: julie.ripoll87@gmail.com
# Created: 2021-02-26
# License: CeCILL
# Update: 2023-04-04
