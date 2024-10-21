import os
import sys


# Input data
INPFOL = config["inpfol"] # specify path of inputs
RESFOL = config["resfol"] # specify pathway for outputs
REFFOL = config["ref_fol"] # specify path to references files
METADATA = config["metadata"] # specify path of metadata file
METACONTRAST = config["metacontrast"] # specify path of metacontrast file

# Params
ANALYSIS = config["analysis"]
COUNT_TOOL = config["tools"]["count"]
STATS_TOOL = config["tools"]["stats"]
ENRICHMENT = config["tools"]["enrichment"]
ORGANISM = config["organism"] # organism name
RELEASE = config['release'] # release of the organism specified in config yaml file
PVAL = config["pval"] # desired threshold
LOGFC = config["logFC"] # desired threshold
FILTER_METHOD = config["filtering_method"] # method choosen for filtering
REPLICAT = config["replicat"] # number of replicates


#######################################################################

if int(REPLICAT) < 3:
    print("""
    WARNING: replicates should be at least >= 3, you cannot continue with this pipeline.
    """)
else:
    def get_all_inputs():
        """
        function to get all input files according to conditions defined in configuration file
        """
        inputs = {}

        if ENRICHMENT == "None":
            inputs["rep"] = "report/"+config["report"]["filename"]+".html"

        inputs["logplt"] = RESFOL+"logplot/SessionLogLogPlots.RData"

        return inputs


    rule all:
        """
        Statistics
        ==========

        Requirements
        ------------
        - snakemake with python 3.10 min. & conda / miniconda or mamba
        - R: deseq2 /limma /ggplot2
        - pandas
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
        >snakemake -s translatome/workflows/05_stats.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 1 -p -r -k -n
        
        # print doc
        >snakemake -s translatome/workflows/05_stats.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
        ``
        
        """
        input: 
            **get_all_inputs()


    #######################################################################

    print("------------------------------------------------------")
    print("Pipeline will be launch for:", ANALYSIS, " analysis")
    print("Tool choosen for statistical analysis:", STATS_TOOL)
    print("Available tools: deseq2 or limma")
    print("------------------------------------------------------")

    #######################################################################

    match STATS_TOOL:
        case "deseq2":
            def analysis_dependent_deseq2_script():
                if ANALYSIS == "kinetic":
                    return "../src/deseq2_kinetic.R"
                else:
                    return "../src/deseq2.R"


            rule Diff_gene_expression_deseq2:
                """
                DESeq2 analysis
                Aim of the rule: perform differential expression analysis
                """
                message:
                    """--- Perform differential analysis with DESeq2 ---"""
                conda: 
                    "../envs/stats/deseq2.yml"
                log:
                    RESFOL+"differential_expression_"+STATS_TOOL+"/deseq2.log"
                input:
                    data = RESFOL+"Counts_"+COUNT_TOOL+"/All_counts.tsv",
                    metadata = METADATA,
                    metacontrast = METACONTRAST
                params:
                    pval = PVAL,
                    logFC = LOGFC,
                    count_tool = COUNT_TOOL,
                    filtering_method = FILTER_METHOD # can be edger, deseq or default
                output:
                    directory(RESFOL+"differential_expression_"+STATS_TOOL+"/")
                script:
                    analysis_dependent_deseq2_script()


        case "limma":
            def analysis_dependent_limma_script():
                if ANALYSIS == "kinetic":
                    return "../src/limma_voom_kinetic.R"  
                else: 
                    return "../src/limma_voom.R"


            rule Diff_gene_expression_limma:
                """
                Limma analysis
                Aim of the rule: perform differential expression analysis
                """
                message:
                    """--- Perform differential analysis with limma ---"""
                conda: 
                    "../envs/stats/limma.yml"
                log:
                    RESFOL+"differential_expression_"+STATS_TOOL+"/limma.log"
                input:
                    data = RESFOL+"Counts_"+COUNT_TOOL+"/All_counts.tsv",
                    metadata = METADATA,
                    metacontrast = METACONTRAST
                params:
                    pval = PVAL,
                    logFC = LOGFC,
                    count_tool = COUNT_TOOL,
                    replicat = REPLICAT if ANALYSIS == "kinetic" else []
                output:
                    directory(RESFOL+"differential_expression_"+STATS_TOOL+"/")
                script:
                    analysis_dependent_limma_script()


    #######################################################################

    rule Categorization:
        """
        Categorization Transcription/Translation/Both and Biodbnet API call
        Aim of the rule: categorization of genes according to transcription or translation status and annotation of genes
        """
        message:
            """--- Categorize data and enrich IDs with biodbnet ---"""
        conda: 
            "../envs/stats/pandas.yml"
        log:
            RESFOL+"categorization_"+STATS_TOOL+"/explore_results.log"
        input:
            input_path = RESFOL+"differential_expression_"+STATS_TOOL+"/",
            gene_names = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene_biotype.tsv",
            metacontrast = METACONTRAST
        params:
            analysis = ANALYSIS,
            inputbiodb = 'Ensembl Gene ID',
            outputbiodb = 'Ensembl Gene Info, Gene Symbol', # can be 'Ensembl Gene Info, Gene Symbol' or 'Gene Info'
            config_tool = STATS_TOOL,
            organismcode = config["organismcode"],
            pval_threshold = PVAL,
            logFC_threshold = LOGFC
        output:
            directory(RESFOL+"categorization_"+STATS_TOOL+"/")
        script:
            "../src/categorize_results.py"


    rule logPLot:
        """
        Aim of the rule: perform log-log plots
        ## Warnings: with Windows WSL2 terminal, this can export all plots in one Rplots.pdf file
        """
        message:
            """--- Perform log-log plots ---"""
        conda: 
            "../envs/stats/ggplot2.yml"
        log:
            RESFOL+"logplot/logplots.log"
        input:
            input_path = RESFOL+"categorization_"+STATS_TOOL+"/"
        params:
            logFC = LOGFC,
            stat_tool = STATS_TOOL
        output:
            output_path = directory(RESFOL+"logplot/"),
            sessioninfo = RESFOL+"logplot/SessionLogLogPlots.RData"
        script:
            "../src/log-log_plots.R"


    #######################################################################

    if ENRICHMENT == "None":
        rule clean:
            """
            Aim of the rule: Clean Rplot in root folder.
            """
            input:
                RESFOL+"logplot/SessionLogLogPlots.RData"
            output:
                RESFOL+"cleaning.log"
            shell:
                "rm Rplots*.pdf && echo 'Cleaning success' > {output}"


        def report_script():
            if ANALYSIS == "kinetic":
                return "../../report/src/report_pol_seq_kinetic.Rmd"
            else:
                return "../../report/src/report_pol_seq.Rmd"


        rule general_report:
            """
            Aim of the rule: General report on the statistical part, 
            including the main graphs and tables of the analysis.
            """
            conda:
                "../envs/stats/rmarkdown.yml"
            input:
                metacontrast = METACONTRAST,
                mds_raw_data = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_raw_data.png",
                biotype = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene_biotype.tsv",
                cat_diff_genes = RESFOL+"categorization_"+STATS_TOOL+"/",
                logplot = RESFOL+"logplot/",
                diff_genes = RESFOL+"differential_expression_"+STATS_TOOL+"/",
                check = RESFOL+"cleaning.log"
            params:
                toolstat = STATS_TOOL,
                context = config["report"]["context"],
                comments = config["report"]["comments"],
                toolenrich = ENRICHMENT,
                aknow = config["report"]["acknowledgements"],
                author = config["report"]["author"],
                filtering_method = FILTER_METHOD 
            output:
                "report/"+config["report"]["filename"]+".html"
            script:
                report_script()


#######################################################################
# Author: "Julie Ripoll & Céline Mandier"
# Contact: julie.ripoll87@gmail.com & celine.mandier@gmail.com
# Created: 2021-12-09
# License: CeCILL
# Update: 2023-04-03
