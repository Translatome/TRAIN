import os
import sys
import re
import subprocess

from wrappers.wrapper_system import inject
from wrappers.ensembl_download import injectable_ensembl_download, injectable_ensembl_download_simplify
from wrappers.fastqscreen import create_FastqScreenDownload


# Input folder
ANALYSIS = config["analysis"]
INPFOL = config["inpfol"] # specify path of inputs
RESFOL = config["resfol"] # specify pathway for outputs
REFFOL = config["ref_fol"] # specify path to references files
METADATA = config["metadata"] # specify path of metadata file
METACONTRAST = config["metacontrast"] # specify path of metacontrast file

# Params
ORGANISM = config["organism"] # organism specified in config yaml file
ORGANISMCODE = config["organismcode"]
RELEASE = config["release"] # release of the organism specified in config yaml file
ASSEMBLY = config["assembly_reference"] # genome or transcriptome cf config yaml file
EXTANNOT = config["extensionAnnotation"] # extension of annotation file specified in config yaml file
THREADS = config["threads"] # low, medium, high or Ultra


#######################################################################

rule all:
    """
    Download references files
    =========================

    Download references files required for the pipeline
    Sources: fastq_screen and Ensembl APIs

    Requirements
    ------------
    - snakemake with python 3.10 min. & conda / miniconda or mamba
    - fastqscreen
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
    >snakemake -s translatome/workflows/00_download_references.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 1 -p -r -k -n
    
    # print doc
    >snakemake -s translatome/workflows/00_download_references.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
    ``
    
    """
    input:
        RESFOL+"biotype_of_genome_"+RELEASE+"/"+ORGANISM+".done"


#######################################################################

print("------------------------------------------------------")
print("Pipeline will be launch for:", ANALYSIS, " analysis")
print("Download of reference files for:", ORGANISM)
print("Release choosen:", RELEASE)
print("------------------------------------------------------")

#######################################################################

rule FastqScreen_download:
    """
    Download reference genomes used by fastq_screen for check of contamination
        Note: default use bowtie2 indexation for reference genomes
    """
    message:
        """--- Check for contamination using Fastq-screen - 1 - download of genomes ---"""
    conda:
        "../envs/pre-alignment/fastqscreen.yml"
    input:
        REFFOL
    params:
        get_genomes = True,
        threads = THREADS
    output:
        outQC = directory(REFFOL+'fastq_screen/'),
        done = REFFOL+'download_genomes_fastqscreen.done'
    shell:
        inject(create_FastqScreenDownload, locals())


def inject_good_wrapper():
    if ORGANISM == "homo_sapiens" or ORGANISM == "mus_musculus":
        return injectable_ensembl_download_simplify
    else:
        return injectable_ensembl_download


rule DLrefAssembly:
    """
    Download genome sequence in fasta from Ensembl database
    """
    message:
        """--- Download reference genome ---"""
    input:
        REFFOL+"download_genomes_fastqscreen.done"
    params:
        organism = ORGANISM,
        version = RELEASE,
        composante = ASSEMBLY, # can be 'dna' or 'cdna' for genome or transcriptome
        extension = "fasta"
    output:
        outdir = directory(REFFOL+RELEASE+"/"+ORGANISM+"/assembly"),
        log = touch(REFFOL+RELEASE+"/"+ORGANISM+"/download_assembly.done")
    shell:
        inject(inject_good_wrapper(), locals())


if EXTANNOT == "gff":
    EXTANNOT = "gff3" # ensembl database keyword 'gff3' required
else:
    EXTANNOT


rule DLrefAnnotation:
    """
    Download reference annotation in gtf or gff from Ensembl database
        Note: Checksums are also downloaded.
        You can check if the files are correctly downloaded with md5sum or sha256 according to the key type
    """
    message:
        """--- Download reference annotation ---"""
    input:
        REFFOL+RELEASE+"/"+ORGANISM+"/download_assembly.done"
    params:
        organism = ORGANISM,
        version = RELEASE,
        extension = EXTANNOT
    output:
        outdir = directory(REFFOL+RELEASE+"/"+ORGANISM+"/annotation"),
        log = touch(REFFOL+RELEASE+"/"+ORGANISM+"/download_annotation.done")
    shell:
        inject(injectable_ensembl_download, locals())


rule unzipREF:
    """
    Unzip all references files from Ensembl before pipeline execution
    """
    message:
        """--- Unzip reference files ---"""
    conda:
        "../envs/pre-alignment/pigz.yml"
    input:
        REFFOL+RELEASE+"/"+ORGANISM+"/download_annotation.done"
    params:
        REFFOL+RELEASE+"/"+ORGANISM+""
    output:
        REFFOL+RELEASE+"/"+ORGANISM+"/unzip.done"
    shell:
        "pigz -d {params}/*/*.gz "
        "2>{output}"


rule Check_inputs:
    """
    Verify if files are ready to use
    Files checked are all fasta, fastq and gtf/gff annotation files
    """
    message:
        """--- Check input files ---"""
    conda:
        "../envs/stats/pandas.yml"
    input:
        check = REFFOL+RELEASE+"/"+ORGANISM+"/unzip.done",
        data = INPFOL,
        metadata = METADATA,
        metacontrast = METACONTRAST
    params:
        analysis = ANALYSIS,
        references = REFFOL+RELEASE+"/"+ORGANISM
    output:
        RESFOL+"check_input_files.log"
    script:
        "../src/check_input.py"


rule ExtractGeneID:
    """
    Extract informations of reference genome
    """
    message:
        """--- Biotype on reference genome ---"""
    conda:
        "../envs/stats/pandas.yml"
    log:
        REFFOL+RELEASE+"/"+ORGANISM+"/annotation/extractID.log"
    input:
        RESFOL+"check_input_files.log"
    params:
        pathto = REFFOL+RELEASE+"/"+ORGANISM+"/annotation/",
        extensionAnnotation = EXTANNOT,
    output:
        outfile = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene.tsv"
    script:
        "../src/ExtractIDfromGFF.py"


rule Biotype:
    """
    Extract informations on biotype of reference genome
    """
    message:
        """--- Biotype on reference genome ---"""
    conda:
        "../envs/stats/pandas.yml"
    input:
        REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene.tsv"
    params:
        inputbiodb = "'Ensembl Gene ID'",
        outputbiodb = "'Gene Symbol, Ensembl Gene Info'",
        organismcode = ORGANISMCODE,
        extensionAnnotation = EXTANNOT,
        outfile = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene_biotype.tsv"
    output:
        done = RESFOL+"biotype_of_genome_"+RELEASE+"/"+ORGANISM+".done"
    shell:
        "chmod +x translatome/src/biodbnet_db2db.py && "
        "python translatome/src/biodbnet_db2db.py "
        "--input_file {input} "
        "--identifier {params.inputbiodb} "
        "--other_id {params.outputbiodb} "
        "--organism {params.organismcode} "
        "--chunck 200 "
        "--biodb {params.outfile} "
        "1> {output.done}"


#######################################################################
# Author: RipollJ
# Contact: julie.ripoll87@gmail.com
# Created: 2021-02-22
# License: CeCILL
# Update: 2023-09-01
