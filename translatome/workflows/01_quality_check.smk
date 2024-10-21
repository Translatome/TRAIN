import os
import sys
import re
import subprocess

from wrappers.wrapper_system import inject
from wrappers.fastqc import create_fastqc
from wrappers.multiqc import create_multiqc
from wrappers.fastqscreen import create_FastqScreen


# Input folder
INPFOL = config["inpfol"] # specify path of inputs
RESFOL = config["resfol"] # specify pathway for outputs
REFFOL = config["ref_fol"] # specify path to references files

# Params
ANALYSIS = config["analysis"]
THREADS = config["threads"] # low, medium, high or Ultra
EXTFiles = config["extension_files"] # extension of raw files can be fastq or fastq.gz


#######################################################################

rule all:
    """
    Check quality of raw data
    =========================

    Requirements
    ------------
    - snakemake with python 3.10 min. & conda / miniconda or mamba
    - FastQC
    - FastqScreen 
    - MultiQC
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
    >snakemake -s translatome/workflows/01_quality_check.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 1 -p -r -k -n

    # print doc
    >snakemake -s translatome/workflows/01_quality_check.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
    ``

    """
    input:
        MQC_FQC = RESFOL+'MQC_reports/MQC_fastqc.done',
        MQC_FQS = RESFOL+'MQC_reports/MQC_fastqscreen.done'


#######################################################################

print("------------------------------------------------------")
print("Pipeline will be launch for:", ANALYSIS, " analysis")
print("------------------------------------------------------")

#######################################################################

rule FastQC:
    """
    Aim: Checks sequence quality using FastQC
    See: FastQC documentation for others parameters
    Link: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    """
    message:
        """--- FastQC on samples ---"""
    conda:
        "../envs/pre-alignment/fastqc.yml"
    input:
        FReads = INPFOL
    params:
        threads = THREADS,
        create_folder = True,
        realformat = EXTFiles, # specify exact format of your file if fq.gz or fq
        noextract = True # don't uncompress output
    output:
        outQC = directory(RESFOL+'Quality/FQC'),
        done = RESFOL+'Quality/FQC.done'
    shell:
        inject(create_fastqc, locals())


#######################################################################

rule MultiQC:
    """
    Aim: Agrgegate results from FastQC into a single report
    See: MultiQC documentation for others parameters
    Link: https://multiqc.info/docs/
    """
    message:
        """--- MultiQC report ---"""
    conda:
        "../envs/pre-alignment/multiqc.yml"
    input:
        FReads = INPFOL,
        check = RESFOL+'Quality/FQC.done'
    params:
        create_folder = True,
        config = "translatome/src/wrappers/config_tools/multiqc_config.yml", # MultiQC config YAML for fastqc
        module = "'fastqc'", # [module name] Use only this module. Can specify multiple times.
        title = "'Report of FQC on raw data'", # TEXT Report title.
        filename = "'report_fastqc.html'", # report name
        dirs = RESFOL+"Quality/FQC/", # Prepend directory to sample names
        fullnames = True, # Do not clean the sample names (leave as full file name)
        zipdatadir = True, # Compress the data directory.
        export = True, # Export plots as static images in addition to the report
        verbose = True, # Increase output verbosity.
        interactive = True # interactive js plots
    output:
        outQC = directory(RESFOL+'MQC_reports/MQC_fastqc/'),
        done = RESFOL+'MQC_reports/MQC_fastqc.done'
    shell:
        inject(create_multiqc, locals())


######################################################################

rule FastqScreen_mapping:
    """
    Aim: Check for contamination
    Requirements: step 000_DL_references.smk rule FastqScreen_download
    See: FastqScreen documentation for others parameters
    Link: https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html#
    """
    message:
        """--- Check for contamination using Fastq-screen - 2 - mapping ---"""
    conda:
        "../envs/pre-alignment/fastqscreen.yml"
    log:
        RESFOL+'Quality/fastq_screen.done'
    input:
        genomes = REFFOL+'fastq_screen/FastQ_Screen_Genomes/',
        FReads = INPFOL
    params:
        threads = THREADS,
        realformat = EXTFiles,
        aligner = "bowtie2", # default for databases download
        bowtie2 = "'-k 2 --very-fast-local'", # here change name to bowtie or bwa according to wanted aligner and adapt params
        conf = REFFOL+'fastq_screen/FastQ_Screen_Genomes/fastq_screen.conf'
    output:
        outQC = directory(RESFOL+'Quality/fastq_screen'),
        done = RESFOL+'Quality/fastq_screen.done'
    shell:
        inject(create_FastqScreen, locals())


#######################################################################

rule MultiQC_fqscreen:
    """
    Aim: Agrgegate results from FastqScreen into a single report
    See: MultiQC documentation for others parameters
    Link: https://multiqc.info/docs/
    """
    message:
        """--- MultiQC reports for fastq_screen output ---"""
    conda:
        "../envs/pre-alignment/multiqc.yml"
    input:
        RESFOL+'Quality/fastq_screen'
    params:
        module = "'fastq_screen'", # [module name] Use only this module. Can specify multiple times.
        title = "'Report of FastqScreen on raw data'", # TEXT Report title.
        filename = "'report_fastq_screen.html'", # report name
        dirs = RESFOL+"Quality/fastq_screen/*_screen.txt", # Prepend directory to sample names
        fullnames = True, # Do not clean the sample names (leave as full file name)
        zipdatadir = True, # Compress the data directory.
        export = False, # Export plots as static images in addition to the report
        verbose = True, # Increase output verbosity.
        interactive = True # interactive js plots
    output:
        outQC = directory(RESFOL+'MQC_reports/MQC_fastqscreen/'),
        done = RESFOL+'MQC_reports/MQC_fastqscreen.done'
    shell:
        inject(create_multiqc, locals())


#######################################################################
# Author: RipollJ
# Contact: julie.ripoll87@gmail.com
# Created: 2018-08-01
# License: CeCILL
# Update: 2023-03-21
