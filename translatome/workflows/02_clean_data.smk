import os
import sys
import re
import subprocess

from wrappers.wrapper_system import inject, val_mappy, mappy
from wrappers.cutadapt import create_cutadapt
from wrappers.fastp import create_fastp
from wrappers.bbtools import create_bbmerge
from wrappers.multiqc import create_multiqc
from wrappers.fastqc import create_fastqc


# Input folder
INPFOL = config["inpfol"] # specify path of inputs
RESFOL = config["resfol"] # specify pathway for outputs
REFFOL = config["ref_fol"] # specify path to references files

# Params
ANALYSIS = config["analysis"]
TRIM_TOOL = config["tools"]["trimming"]
THREADS = config["threads"] # low, medium, high or Ultra
EXTFiles = config["extension_files"] # extension of raw files can be fastq or fastq.gz
PAIRED = config["paired"] # paired data or not

QUALITY = config["trim"]["quality"] # quality thtreshold
FREF = config["trim"]["fref"] # Comma-delimited list of fasta reference for filtering
MINILENGTH = config["trim"]["minilength"]

CUTLENGTH1 = config['trim']["cutlength1"]
CUTLENGTH2 = config['trim']["cutlength2"]
MORELENGTH1 = config['trim']["morelength1"]
MORELENGTH2 = config['trim']["morelength2"]
CUTPHRED = config['trim']["cutphred"]

TRIMFRONT1 = config['trim']["trimfront1"]
TRIMFRONT2 = config['trim']["trimfront2"]
TRIMTAIL1 = config['trim']["trimtail1"]
TRIMTAIL2 = config['trim']["trimtail2"]


# Wildcards
PAIRED_INFO_1 = config["paired_info_1"]
PAIRED_INFO_2 = config["paired_info_2"]
END = config["end_filename"]
SAMPLE, = glob_wildcards(INPFOL+"{sample}"+(PAIRED_INFO_1 if PAIRED == "YES" else "")+END)


#######################################################################

def get_all_inputs():
    inputs = {}

    if TRIM_TOOL == "cutadapt" and  PAIRED == "YES":
        inputs["bbmerge"] = expand(RESFOL+'Quality/BBmerge/{sample}.done', sample = SAMPLE)

    if TRIM_TOOL == "fastp":
        inputs["fastp_end"] = RESFOL+"MQC_reports/MQC_fastqc2_trim_"+TRIM_TOOL+".done"

    inputs["mqc_end"] = RESFOL+"MQC_reports/MQC_trim_"+TRIM_TOOL+".done"

    return inputs


rule all:
    """
    Cleaning of raw data
    ====================

    Requirements
    ------------
    - snakemake with python 3.10 min. & conda / miniconda or mamba
    - cutadapt / fastp
    - bbmap
    - fastqc
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
    >snakemake -s translatome/workflows/02_clean_data.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 1 -p -r -k -n

    # print doc
    >snakemake -s translatome/workflows/02_clean_data.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
    ``
    
    """
    input:
        **get_all_inputs()


#######################################################################

print("------------------------------------------------------")
print("Pipeline will be launch for:", ANALYSIS, " analysis")
print("Tool choosen for trimming:", TRIM_TOOL)
print("Available tools: cutadapt or fastp")
print("------------------------------------------------------")

#######################################################################

fareadss = val_mappy(PAIRED,
                    {"YES": {"fareads1": INPFOL+"{sample}"+PAIRED_INFO_1+END, 
                            "fareads2": INPFOL+"{sample}"+PAIRED_INFO_2+END},
                    "NO": {"fareadse": INPFOL+"{sample}"+END}})

cutreadss = val_mappy(PAIRED,
                    {"YES": {"cutreads1": RESFOL+"trimming_"+TRIM_TOOL+"/{sample}"+PAIRED_INFO_1+END, 
                            "cutreads2": RESFOL+"trimming_"+TRIM_TOOL+"/{sample}"+PAIRED_INFO_2+END},
                    "NO": {"cutreadse": RESFOL+"trimming_"+TRIM_TOOL+"/{sample}"+END}})


match TRIM_TOOL:
    case "cutadapt":
        rule Cutadapt:
            """
            Aim of the rule: Cleanning of adapter sequences, primers, poly-A tails and other types of unwanted sequence in high-throughput sequencing reads
            See: Cutadapt wrapper documentation for others parameters
            """
            message:
                """--- Cutadapt on {wildcards.sample} / support compressed files --- """
            conda:
                "../envs/pre-alignment/cutadapt.yml"
            input:
                **fareadss
            params:
                paired = PAIRED,
                strandtrim = True,
                strand = "35",
                # list adapters to remove for paired-end data
                adapter1 = FREF,
                adapter2 = FREF,
                # number of cores to use
                cores = THREADS,
                # Remove bases from each read (first read only if paired)
                cut = True,
                # length to cut, positive value start at the beginning of the read whereas negative value at the end
                cutlength1 = CUTLENGTH1,
                cutlength2 = CUTLENGTH2,
                morelength1 = MORELENGTH1,
                morelength2 = MORELENGTH2,
                # debug
                debug = False,
                # Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal.
                qualitycutoff = CUTPHRED,
                # Trim N's and poly-a on ends of reads
                trimn = True,
                polya = True,
                # discard trimmed reads that are shorter than length
                minimumlength = MINILENGTH
            output:
                **cutreadss,
                # write information about each read and its adapter matches to file
                infos = RESFOL+"trimming_"+TRIM_TOOL+"/{sample}_trim.done"
            shell:
                inject(create_cutadapt, locals())


        rule MultiQC_cutadapt:
            """
            Aim of the rule: Agregate results from Cutadapt into a single report
            See: MultiQC wrapper documentation for others parameters
            """
            message:
                """--- MultiQC reports for cutadapt output ---"""
            conda:
                "../envs/pre-alignment/multiqc.yml"
            input:
                expand(RESFOL+"trimming_"+TRIM_TOOL+"/{sample}_trim.done", sample = SAMPLE),
                RESFOL+"MQC_reports/MQC_fastqc2_trim_"+TRIM_TOOL+".done"
            params:
                module = "'cutadapt'",
                title = "'Report for Cutadapt'",
                filename = "'report_cutadapt.html'",
                dirs = RESFOL+"trimming_"+TRIM_TOOL+"/*.done",
                fullnames = True,
                zipdatadir = True,
                export = True,
                verbose = True,
                interactive = True
            output:
                outQC = directory(RESFOL+"MQC_reports/MQC_trim_"+TRIM_TOOL+"/"),
                done = RESFOL+"MQC_reports/MQC_trim_"+TRIM_TOOL+".done"
            shell:
                inject(create_multiqc, locals())


        if PAIRED == "YES":
            rule AverageInsertSize:
                """
                Aim of the rule: Check average insert size of paired-end sample
                See: bbtools wrapper documentation for others parameters
                """
                message:
                    """--- BBmerge check average insert size ---"""
                conda:
                    "../envs/pre-alignment/bbmap.yml"
                input:
                    inputR1 = INPFOL+'{sample}'+PAIRED_INFO_1+END,
                    inputR2 = INPFOL+'{sample}'+PAIRED_INFO_2+END,
                    check = RESFOL+"trimming_"+TRIM_TOOL+"/{sample}_trim.done"
                params:
                    paired = PAIRED,
                    reads = '2m'
                output:
                    merged = temp(RESFOL+'Quality/BBmerge/{sample}.fq'),
                    hist = RESFOL+'Quality/BBmerge/ihist_{sample}.txt',
                    done = RESFOL+'Quality/BBmerge/{sample}.done'
                shell:
                    inject(create_bbmerge, locals())


    case "fastp":
        rule Fastp:
            """
            Aim of the rule: Cleanning of adapter sequences, primers, poly-A tails and other types of unwanted sequence in high-throughput sequencing reads
            See: FASTP wrapper documentation for more informations on parameters
            """
            message:
                """--- fastp on {wildcards.sample} / support compressed files ---"""
            conda:
                "../envs/pre-alignment/fastp.yml"
            log:
                RESFOL+"trimming_"+TRIM_TOOL+"/{sample}_trim.done"
            input:
                **fareadss
            params:
                paired = PAIRED,
                # Adapters options 
                disableadaptertrimming = False, 
                adaptersequence = False,
                adaptersequencer2 = False,
                detectadapterforpe = True,
                adapterfasta = False,
                # Number of cores to use
                thread = THREADS,
                # Quality 
                qualityphred = QUALITY,
                # Discard trimmed reads that are shorter than length
                lengthrequired = MINILENGTH,
                # Global trimming
                trimfront1 = TRIMFRONT1 if TRIMFRONT1 != "" else [],
                trimfront2 = TRIMFRONT2 if TRIMFRONT2 != "" else [],
                trimtail1 = TRIMTAIL1 if TRIMTAIL1 != "" else [],
                trimtail2 = TRIMTAIL2 if TRIMTAIL2 != "" else [],
                # PolyX/G tail trimming
                trimpolyg = True,
                polygminlen = False,
                disabletrimpolyg = False,
                trimpolyx = True,
                polyxminlen = False
            output:
                **cutreadss,
                json =  RESFOL+"trimming_"+TRIM_TOOL+"/Json/{sample}_fastp.json",
                html =  RESFOL+"trimming_"+TRIM_TOOL+"/Html/{sample}_fastp.html"
            shell:
                inject(create_fastp, locals())


        rule MultiQC_Fastp:
            """
            Aim of the rule: Agrgegate results from Fastp (trimmed samples) into a single report
            See: MultiQC wrapper documentation for others parameters
            """
            message:
                """--- MultiQC reports for fastp output ---"""
            conda:
                "../envs/pre-alignment/multiqc.yml"
            input:
                check = RESFOL+"MQC_reports/MQC_fastqc2_trim_"+TRIM_TOOL+".done"
            params:
                module = "'fastp'",
                title = "'Report for Fastp'",
                filename = "'report_fastp.html'",
                dirs = RESFOL+"trimming_"+TRIM_TOOL+"/Json/*_fastp.json",
                export = True,
                verbose = True
            output:
                outQC = directory(RESFOL+"MQC_reports/MQC_trim_"+TRIM_TOOL+"/"),
                done = RESFOL+"MQC_reports/MQC_trim_"+TRIM_TOOL+".done"
            shell:
                inject(create_multiqc, locals())


#######################################################################

def check_trimmed_input():
    if TRIM_TOOL == "cutadapt":
        return expand(RESFOL+"trimming_"+TRIM_TOOL+"/{sample}_trim.done", sample = SAMPLE)
    else:
        return expand(RESFOL+"trimming_"+TRIM_TOOL+"/Json/{sample}_fastp.json", sample = SAMPLE)


rule FastQC_after_trim:
    """
    Aim of the rule: Check sequence quality using FastQC
    See: FastQC wrapper documentation for others parameters
    """
    message:
        """--- FastQC on trimmed samples ---"""
    conda:
        "../envs/pre-alignment/fastqc.yml"
    input:
        check_trimmed_input()
    params:
        FReads = RESFOL+"trimming_"+TRIM_TOOL,
        threads = THREADS,
        create_folder = True,
        realformat = EXTFiles, # specify exact format of your file if fq.gz or fq
        noextract = True # don't uncompress output
    output:
        outQC = directory(RESFOL+"Quality/FQC2_trim_"+TRIM_TOOL+"/"),
        done = RESFOL+"Quality/FQC2_trim_"+TRIM_TOOL+".done"
    shell:
        inject(create_fastqc, locals())


rule MultiQC_after_trim:
    """
    Aim of the rule: Agregate results from FastQC of trimmed samples after trimming into a single report
    See: MultiQC wrapper documentation for others parameters
    """
    message:
        """--- MultiQC reports on FastQC outputs ---"""
    conda:
        "../envs/pre-alignment/multiqc.yml"
    input:
        RESFOL+"Quality/FQC2_trim_"+TRIM_TOOL+".done"
    params:
        config = "translatome/src/wrappers/config_tools/multiqc_config.yml", # MultiQC config YAML for fastqc
        module = "'fastqc'", # [module name] Use only this module. Can specify multiple times.
        title = "'Report FastQC on trimmed outputs'", # TEXT Report title.
        filename = "'report_fastqc2_trim.html'", # report name
        dirs = RESFOL+"Quality/FQC2_trim_"+TRIM_TOOL+"/", # Prepend directory to sample names
        fullnames = True, # Do not clean the sample names (leave as full file name)
        zipdatadir = True, # Compress the data directory.
        export = True, # Export plots as static images in addition to the report
        verbose = True, # Increase output verbosity.
        interactive = True # interactive js plots
    output:
        outQC = directory(RESFOL+"MQC_reports/MQC_FQC2_"+TRIM_TOOL+"/"),
        done = RESFOL+"MQC_reports/MQC_fastqc2_trim_"+TRIM_TOOL+".done"
    shell:
        inject(create_multiqc, locals())


#######################################################################
# Author: RipollJ
# Contact: julie.ripoll87@gmail.com
# Created: 2022-09-08
# License: CeCILL
# Update: 2023-07-17
