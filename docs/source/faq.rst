Frequently Asked Questions
==========================

**Q. What should I do if I get an error message saying that the wrapper paths are not available ?**

If you get an error message when running a snakefile saying that the access paths to the wrappers are not available, please run the following command in the root of the folder:

.. code-block:: shell

    pip install .

**Q. How do I activate my snakemake environment ?**

If you get the error:

.. code-block:: shell

     bash: snakemake: command not found

then run the command:

.. code-block:: shell

    conda activate name_env_snakemake


**Q. Why does the 02_clean_data.smk snakefile take a long time to run ?**

The *02_clean_data.smk* snakefile can take a long time to run when Cutadapt is selected, as it checks the list of supplied adapters one by one. The longer the list, the longer it will take to run. To make sure it is running, go to the *results/trimming* folder and follow the `*.done` files which list the trimming information for each sample.

**Q. The pipeline didn't consider my second read file. What's the problem ?**

This may be due to the fact that the paired option in the *config_PE.yaml* file being set to "NO". If this is the case, please change it to "YES". Don't forget to also check the other options to match the paired-end data. It is also possible that the configuration file has been reversed, so make that sure you select the configuration file *config_PE.yaml* for a non-kinetic pair-end analysis and *config_PE_kinetic.yml* for a kinetic pair-end analysis.

**Q. I can't run 00_download_references.smk or 06_enrichment.smk on a Compute Cluster, what's wrong ?**

Snakefiles that require an internet connection (API queries, reference file downloads, etc.) cannot run on computer clusters that do not have an internet connection. Run these snakefiles on your local machine.

**Q. Why aren't the images displayed in the HTML report ?**

One of the main reasons why report images are not displayed is that the pipeline folder tree structure is not respected. The tree must be structured as it was originally defined (see the git repository structure and pipeline documentation) and the folder/file names must remain unchanged. If you only want to share the report, please save it in pdf format (this will fix the images in the report).

**Q. Why do I get warnings when I run the general report rule ?**

The presence of warnings during the execution of the general report is normal. They should be ignored. When you open the HTML report, the images should still be there. 
Exmple of normal warning :

.. code-block:: shell

    [WARNING] Could not fetch resource img/LogoLIRMM_sansfond.png
    [WARNING] Could not fetch resource img/logo_um_2022_rouge_RVB.png


**Q. Snakemake can't find my input files, what should I do ?**

If snakemake cannot find the input files, this may be due to an error in the input paths. Check the configuration file. 
It is also important not to change the pipeline folder tree, nor the folder/file names.

**Q. I get empty enrichment files or I don't get any files for certain samples. Is this normal ?**

In some cases, you may receive empty files or files with a "No enrichment" line. These files indicate that there are no enrichment results for the entered genes. The same applies to missing files such as enrichment plots. As these are based on the output tables, the graphs will not be generated if no results are found.

**Q. Can I run this pipeline on less than 3 replicates ?**

No, it is designed to run with a minimum of 3 replicates. In addition, the statistical methods used can only be used with a minimum of 3 replicates.

**Q. I get a low percentage count with feature count or htseq counts.**

This may be due to a wrong choice of strandness, if you don't know this information (which must be provided by the sequencing platform), please test the three strandness options with Feature Counts. If the percentage of strandedness is high with the "unstranded" option and the "stranded" option, then the recommended choice is "stranded". If the percentage of strandness is high with the "unstranded" option and the "reversely stranded" option, then the recommended choice is "reversely stranded". If the percentage of strandedness is high with the "unstranded" option but is lower with the "stranded" and "reversely stranded" options, then the recommended choice is "unstranded".

**Q. When should kinetic analysis be used ?**

The kinetic analysis is adapted to POL-seq and RNA-seq data collected over several time points. The pipeline can be run on an unlimited number of time points.

**Q. How can I change the calculation resources that require a task to be executed ?**

It is possible to adjust the calculation resources by changing the number of threads using the "threads" option in the configuration file. The "-j" option is used to set the maximum number of CPUs used to run jobs in parallel.

**Q. I get a parser error when executing snakefile 06_enrichment.smk. What should I do ?**

If you get a parser error, run the same snakemake command again. After the second run, the error should no longer appear.

**Q. An error occurs when checking the metacontrast or metadata, even though I have respected the required format. What should I do ?**

Check the separator, it should be a tabulation and that the column names are as requested. The encoding of the file may cause problems. In this case, go back to the example files available on the repository and enter your information.
