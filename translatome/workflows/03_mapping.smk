import os
import sys
import re
import subprocess

from wrappers.wrapper_system import inject, val_mappy
from wrappers.STAR import create_STARIndex, create_STAR
from wrappers.samtools import create_samtools_index, create_samtools_sort, create_samtools_flagstat, create_samtools_stats, create_samtools_view
from wrappers.multiqc import create_multiqc


# Input folder
INPFOL = config["inpfol"] # specify path of inputs
RESFOL = config["resfol"] # specify pathway for outputs
REFFOL = config["ref_fol"] # specify path to references files

# Params
ANALYSIS = config["analysis"]
TRIM_TOOL = config["tools"]["trimming"]
RLENGTH = config["rlength"] # basis of read length -1
PAIRED = config["paired"] # paired data or not
RELEASE = config["release"] # release version of annotations files
THREADS = config["threads"]  # low, medium, high or Ultra
ORGANISM = config["organism"] # organism specified in config yaml file
GENREF = config["references"]["assembly_ref"]
ANNOT = config["references"]["annotation"]

# Wildcards
PAIRED_INFO_1 = config["paired_info_1"]
PAIRED_INFO_2 = config["paired_info_2"]
END = config["end_filename"]
SAMPLE, = glob_wildcards(INPFOL+"{sample}"+(PAIRED_INFO_1 if PAIRED == "YES" else "")+END)


#######################################################################

rule all:
    """
    Alignment on reference
    ======================

    Requirements
    ------------
    - snakemake with python 3.10 min. & conda / miniconda or mamba
    - star
    - samtools
    - multiQC
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
    # WARNING: use -j 2 due to pipe of rules STAR_alignment and Samtools_view
    >snakemake -s translatome/workflows/03_mapping.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 2 -p -r -k -n

    # print doc
    >snakemake -s translatome/workflows/03_mapping.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
    ``
    
    """
    input:
        mqc_end = RESFOL+'MQC_reports/MQC_STAR'


#######################################################################

print("------------------------------------------------------")
print("Pipeline will be launch for:", ANALYSIS, " analysis")
print("------------------------------------------------------")

#######################################################################

# input files
fareadss = val_mappy(PAIRED, 
                            {"YES":
                            {"fareads1": RESFOL+"trimming_"+TRIM_TOOL+"/{sample}"+PAIRED_INFO_1+END,
                            "fareads2": RESFOL+"trimming_"+TRIM_TOOL+"/{sample}"+PAIRED_INFO_2+END},
                            "NO": {"fareads": RESFOL+"trimming_"+TRIM_TOOL+"/{sample}"+END}})


rule IndexSTAR:
    """
    Aim of the rule: prepare reference for alignment (indexing)
    """
    message:
        """--- Indexing reference genome ---"""
    conda:
        "../envs/aligner/star.yml"
    input:
        faref = REFFOL+RELEASE+"/"+ORGANISM+"/assembly/"+GENREF,
        annot = REFFOL+RELEASE+"/"+ORGANISM+"/annotation/"+ANNOT
    params:
        readlength = RLENGTH,
        threads = THREADS
    output:
        refg = directory(REFFOL+RELEASE+"/"+ORGANISM+"/star_index/")
    shell:
        inject(create_STARIndex, locals())


rule STAR_alignment:
    """
    Aim of the rule: alignment of reads on a reference using STAR
    """
    message:
        """--- Alignment of isoforms on reference fasta (genome or transcriptome) {wildcards.sample} ---"""
    conda:
        "../envs/aligner/star.yml"
    log:
        RESFOL+"Mapping/STAR/{sample}.done"
    input:
        **fareadss,
        refg  = expand(REFFOL+RELEASE+"/"+ORGANISM+"/star_index/", release = RELEASE)
    params:
        outfilename = RESFOL+"Mapping/STAR/{sample}",
        annot = REFFOL+RELEASE+"/"+ORGANISM+"/annotation/"+ANNOT,
        paired = PAIRED,
        readFilesCommand = "zcat", # to uncompress
        threads = THREADS,
        create_folder = True
    output:
        pipe(RESFOL+"Mapping/STAR/{sample}Aligned.out.sam")
    shell:
        inject(create_STAR, locals())


rule Samtools_view:
    """
    Aim of the rule: convert sam in bam
    """
    message:
        """--- Convert in bam files {wildcards.sample} ---"""
    conda:
        "../envs/post-alignment/samtools.yml"
    log:
        RESFOL+"Samtools/samtools_view_{sample}.log"
    input:
        RESFOL+"Mapping/STAR/{sample}Aligned.out.sam"
    params:
        threads = THREADS
    output:
        out = RESFOL+"Mapping/STAR/Samtools/{sample}.bam"
    shell:
        inject(create_samtools_view, locals())


rule Sort_bam:
    """
    Aim of the rule: Sort bam output files
    """
    message:
        """--- Sort bam files {wildcards.sample} ---"""
    conda:
        "../envs/post-alignment/samtools.yml"
    log:
        RESFOL+"Mapping/STAR/Samtools/samtools_sort_{sample}.log"
    input:
        RESFOL+"Mapping/STAR/Samtools/{sample}.bam"
    params:
        threads = THREADS
    output:
        out = RESFOL+"Mapping/STAR/Samtools/{sample}_sorted.bam",
        console = RESFOL+"Mapping/STAR/Samtools/samtools_sort_{sample}.done"
    shell:
        inject(create_samtools_sort, locals())


rule Samtools_index:
    """
    Aim of the rule: indexes SAM/BAM/CRAM files
    """
    message:
        """--- Index bam files with Samtools {wildcards.sample} ---"""
    conda:
        "../envs/post-alignment/samtools.yml"
    log:
        RESFOL+"Mapping/STAR/Samtools/samtools_index_{sample}.log"
    input:
        RESFOL+"Mapping/STAR/Samtools/{sample}_sorted.bam"
    params:
        threads = THREADS
    output:
        RESFOL+"Mapping/STAR/Samtools/{sample}_sorted.bam.bai"
    shell:
        inject(create_samtools_index, locals())


rule MultiQC_STAR:
    """
    Aim of the rule: Generate a uniq report for STAR log files.
    """
    message:
        """--- MultiQC on STAR log files ---"""
    conda:
        "../envs/pre-alignment/multiqc.yml"
    input:
        expand(RESFOL+"Mapping/STAR/Samtools/{sample}_sorted.bam.bai", sample = SAMPLE)
    params:
        module = "star",
        title = "'Report for STAR alignments'",
        filename = "report_star.html",
        dirs = RESFOL+"/Mapping/STAR/", # directory to check
        fullnames = True, # Do not clean the sample names (leave as full file name)
        viewtags = False, # View the available tags and which modules they load
        zipdatadir = True, # Compress the data directory
        export = True, # Export plots as static images in addition to the report
        verbose = True, # Increase output verbosity
        interactive = True # interactive js plots
    output:
        outQC = directory(RESFOL+'MQC_reports/MQC_STAR/'),
        done = touch(RESFOL+'MQC_reports/MQC_star.done')
    shell:
        inject(create_multiqc, locals()) 


#######################################################################
# Author: RipollJ
# Contact: julie.ripoll87@gmail.com
# Created: 2020-08-26
# License: CeCILL
# Update: 2023-03-27
