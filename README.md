# Translatome

Workflow for polysome profiling analysis

- [Translatome](#translatome)
  - [Context](#context)
  - [LICENSE](#license)
  - [Contributions](#contributions)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Clone workflow into desired working directory](#clone-workflow-into-desired-working-directory)
    - [Installation of Snakemake through MiniConda](#installation-of-snakemake-through-miniconda)
    - [Installation of Snakemake through Docker](#installation-of-snakemake-through-docker)
    - [Installation of wrapper system package](#installation-of-wrapper-system-package)
  - [Documentation](#documentation)
  - [Configuration](#configuration)
    - [Edit the configuration file and complete with your project information](#edit-the-configuration-file-and-complete-with-your-project-information)
    - [Special requirements: Metadata and Metacontrast files](#special-requirements-metadata-and-metacontrast-files)
  - [Execution](#execution)
    - [Run the rule's documentation](#run-the-rules-documentation)
    - [Run the snakefile in dry-run mode](#run-the-snakefile-in-dry-run-mode)
    - [Steps](#steps)

---

## Context

Translation is the second phase of protein synthesis that involves converting the information on the mRNA sequence to an amino acid sequence.
This process requires mRNA, ribosomes, amino acids, tRNA, and other factors.

Translatomics is the study of all open reading frames (ORFs) that are being actively translated in a cell or organism.
This collection of ORFs is referred to as the translatome.
Characterizing a cell's translatome can give insight into the array of biological pathways that are active in the cell.

Polysome profiling is a technique that infer the translation status of a specific mRNA by analyzing the behavior of polysomes, i.e. to the group of ribosomes bound to an mRNA, that is present during the elongation phase of translation.
Polysome profiling consists in the quantification of mRNA abundance in polysome fractionation by sucrose density gradient centrifugation.
Based on the assumption that ribosome density is highly correlated with protein production, the translational efficiency of individual genes can be calculated.
Translational efficiency is calculated by comparing two RNA-seq profiles, one for transcription (all mRNA) and one for translation (mRNA with ribosomes bound).

We provide here a Snakemake workflow to analyze Polysome profiling (POL-SEQ). This workflow is composed of seven steps:

- Download of reference files and Input checking
- Quality checking
- Filtering of RNA-seq data
- Alignment on reference assembly file (genome or transcriptome)
- Read counts
- Statistics
- GSEA enrichment analysis

This workflow is provided as is, without guaranty if any modification is made.

Several packages and softwares are used here, please cite them or cite this GitHub Project.

[<small>[top↑]</small>](#)

## LICENSE

[LICENSE](https://github.com/Translatome/Translatome/blob/master/LICENSE.md) terms are in agreement with CeCILL License.
See: http://www.cecill.info/licences.fr.html

## Contributions

See the [contributing](https://github.com/Translatome/Translatome/blob/master/CONTRIBUTING.md) file and please follows the recommandation in the [code of conduct](https://github.com/Translatome/Translatome/blob/master/CODE_OF_CONDUCT.md)

[<small>[top↑]</small>](#)

## Requirements

This workflow requires:

- Python version >= 3.10
- Snakemake version >= 7.14.0
- Conda and/or Mamba
- Download from distant databases for steps 00 and 06

Most of the programs are installed using conda environments.

Tools used in this workflow:

- Control quality: FastQC
- Contamination control: FastqScreen
- Quality trimming: Fastp / Cutadapt
- Alignment on reference: STAR
- Converter: Samtools
- Counting: featureCounts / HTSeq-count
- Statistics: DESeq2 / limma-voom
- Enrichment: gprofileR and goatools
- Reports: multiQC and Rmarkdown

[<small>[top↑]</small>](#)

## Installation

### Clone workflow into desired working directory

Copy the Translatome workflow in the desired folder [path/to/workdir]:

```shell
git clone https://github.com/Translatome/Translatome.git [path/to/workdir]
```

### Installation of Snakemake through MiniConda

For installation of Conda, please refer to the [official documentation](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)
and for Mamba, please also refer to the [official documentation](https://mamba.readthedocs.io/en/latest/installation.html).

Here, this is an example:

1. Download the installer miniconda: [https://conda.io/miniconda.html](https://conda.io/miniconda.html)

2. In your Terminal window, run:

```shell
bash Miniconda3-latest-Linux-x86_64.sh
```

3. Follow the prompts on the installer screens.

4. Install Snakemake in a dedicated environment

```shell
conda install -c bioconda -c conda-forge -n snakemake snakemake python=3.10
conda activate snakemake
```

### Installation of Snakemake through Docker

The official docker image can be found at: [https://hub.docker.com/r/snakemake/snakemake](https://hub.docker.com/r/snakemake/snakemake)

This image contains conda, mamba and snakemake.

For more informations on docker image, read the [official documentation](https://docs.docker.com/get-started/overview/).

### Installation of wrapper system package

In order to install the wrappers on your computer, go to the translatome/ folder and use the:

```shell
# require pip for python 3
pip install .
```

[<small>[top↑]</small>](#)


## Documentation

Complete documentation is available at [Translatome](https://jripoll.lirmm.net/translatome)

The developer documentation was produced using sphinx.

[<small>[top↑]</small>](#)


## Configuration

### Edit the configuration file and complete with your project information

This workflow is managed by a global configuration defined in a yaml configuration file. Several example configuration files are available in the "configs" folder.

```shell
vim configs/config_PE.yml
```

The configuration options will change depending on whether the data is kinetic or not, and whether the data is paired end or single end.

These configuration files are divided into several sections, where each option available to the user is specified. Before each run of the workflow, the user must check that the options are correct according to the data.

The first option in this workflow is the choice of analysis type: pol-seq (classical POL-seq analysis combined with RNA-seq) or kinetic (kinetic POL-seq analysis combined with kinetic RNA-seq).

For more informations on the options available in the configuration file, see the [documentation](https://github.com/Translatome/Translatome/blob/master/docs/source/config_management.md).

### Special requirements: Metadata and Metacontrast files

This workflow requires **Metadata** and **Metacontrast** files. These files contain all the information needed for the statistical and enrichment parts of the analysis.

Examples are provided in the "ressources" folder.

**Metadata**

This file is required for the first part of the analysis, rule _Diff_gene_expression_ in _05_stats.smk_, where the simple comparisons are made based on sample identification.

This file should contains 4 columns:

| sample_name  | fraction | treatment | replicate |
| ------------ | -------- | --------- | --------- |
| F1_Treat1_R1 | F1       | Treat1    | R1        |
| F1_Treat1_R2 | F1       | Treat1    | R2        |
| F1_Treat1_R3 | F1       | Treat1    | R3        |
| F1_Treat2_R1 | F1       | Treat2    | R1        |
| F1_Treat2_R2 | F1       | Treat2    | R2        |
| ...          | ...      | ...       | ...       |

Nomenclature:

- F1 refers to the Cytoplasmic or Total fraction
- F2 to the Polysomal fraction.

In the case of kinetic analysis, a "time" column is added to the metadata file. See [example](https://github.com/Translatome/Translatome/blob/master/ressources/metadata_kinetic.txt).

**Metacontrast**

This file is required for the first and second part of the anlaysis, rules _Diff_gene_expression_ and _Categorization_ in _05_stats.smk_, where the simple comparisons are compared in a double comparison.

This file should contains 4 columns:

| fraction | treatment1 | fraction_bis | treatment2 |
| -------- | ---------- | ------------ | ---------- |
| F1       | Treat1     | F1           | Treat2     |
| F2       | Treat1     | F2           | Treat2     |

The direction of comparison of the fractions is fixed, first cytoplasmic/total, then polysomal.

In the case of a treatment/control comparison, "treatment1" corresponds to the treatment and Treatment 2 corresponds to the control.

If the comparison is between two treatments (without a control), the user must select the treatment to be considered as the control in order to assign it to the "treatment2" column in the table.

By default, the control should always be assigned to the "treatment2" column and the comparison will always be made by "treatment1" vs. "treatment2".

You can add as many contrasts as you like within the required format. The separator is a tab and the file format is ".txt".

In the case of a kinetic analysis, two columns are required for the time to be compared, "time" and "time_bis". Time_bis can be the same time or a different time according to your comparisons. See [example](https://github.com/Translatome/Translatome/blob/master/ressources/metaconstrast_kinetic.txt).

[<small>[top↑]</small>](#)

## Execution

### Run the rule's documentation

Here, this run shows the documentation associated to the rules using the "-l" option:

```shell
snakemake -s translatome/workflows/00_download_references.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
```

### Run the snakefile in dry-run mode

Here, this run tests the execution in dry mode:

```shell
snakemake -s translatome/workflows/00_download_references.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 1 -p -r -k -n
```

Remove the -n option to run the workflow.

To run the workflow on a cluster, please refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

### Steps

Carry out the steps according to the numbering of the snakefiles:

- 00_download_references.smk
- 01_quality_check.smk
- 02_clean_data.smk
- 03_mapping.smk
- 04_counts.smk
- 05_Stats.smk
- 06_enrichment.smk

After the _00_download_references.smk_ step, complete the references files section in the configuration file with the name of your downloaded file before running the rest of the workflow.

Take a look at all the log files before moving on to the next steps.

[<small>[top↑]</small>](#)


Note: This pipeline can be tested with these datasets:
- BioProject [PRJNA741225](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA741225) 75bp paired_end data, [Therizols et al. (2022)](https://doi.org/10.1038/s41467-021-27847-8), *Homo sapiens*, Illumina NextSeq 500
- BioProject [PRJNA397005](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA397005) 100bp single_end data, [Bernabo et al. (2017)](https://doi.org/10.1016/j.celrep.2017.10.010), *Mus muculus*, Illumina HiSeq 2000

