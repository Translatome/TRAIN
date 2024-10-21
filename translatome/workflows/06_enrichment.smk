import os
import sys


# Input folder
INPFOL = config["inpfol"] # specify path of inputs
RESFOL = config["resfol"] # specify pathway for outputs
REFFOL = config["ref_fol"] # specify path to references files
METADATA = config["metadata"] # specify path of metadata file
METACONTRAST = config["metacontrast"] # specify path of metacontrast file

# params
ANALYSIS = config["analysis"]
COUNT_TOOL = config["tools"]["count"]
STATS_TOOL = config["tools"]["stats"]
FILTER_METHOD = config["filtering_method"] # method choosen for filtering of the data
ENRICHMENT = config["tools"]["enrichment"]
PVAL = config["gpro"]["pval_threshold"] # desired threshold
CORRECTION = config["gpro"]["correction"]
ORGANISM = config["organism"] # organism specified in config yaml file
ORGANISMCODE = config["organismcode"]
RELEASE = config['release'] # release of the organism specified in config yaml file
COMPLEX = config["gpro"]["complex"]
GO_SOURCES = config["gpro"]["go_sources"]
PATHWAY_SOURCES = config["gpro"]["pathway"]
TERM_NUMBER = config["gpro"]["term_nb"]

# wildcards
TREATMENT, = glob_wildcards(RESFOL+"categorization_"+STATS_TOOL+"/Categorization_{treat}.tsv")
CATEGORIES = ['Both_mRNA_UP', 'Both_mRNA_DOWN',
              'Transcription_UP', 'Transcription_DOWN',
              'Translation_UP', 'Translation_DOWN',
              'Divergent_UPDN', 'Divergent_DNUP']


#######################################################################

if ENRICHMENT == "None":
    print("------------------------------------------------------------------------------------------------------")
    print("WARNING: Pipeline will be launch with:", ENRICHMENT, " enrichment option.")
    print("Check your config file, if you want enrichment use 'gprofiler' keyword in the tools section.")
    print("------------------------------------------------------------------------------------------------------")
else:
    print("------------------------------------------------------------------------------------------------------")
    print("Pipeline will be launch for:", ANALYSIS, " analysis.")
    print("Wildcards used in this snakefile depend on the execution of 05_stats.smk snakefile.")
    print("------------------------------------------------------------------------------------------------------")

    rule all:
        """
        Enrichment step
        ===============

        Requirements
        ------------
        - snakemake with python 3.10 min. & conda / miniconda or mamba
        - R: gprofiler
        - pandas
        - goatools
        - jupyter notebook
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
        >snakemake -s translatome/workflows/06_enrichment.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -j 1 -p -r -k -n
        
        # print doc
        >snakemake -s translatome/workflows/06_enrichment.smk --use-conda --conda-frontend mamba --configfile configs/config_PE.yml -l
        ``
        
        """
        input:
            clean = RESFOL+"cleaning.log"

    #######################################################################

    rule gprofiler:
        """
        Aim of the rule: Enrichment of genes from the different categories
            Use gprofiler2 API, requires a connection to the web
            Sources used are: 
            - pathways: KEGG, WIKIPATHWAY and REACTOME
            - GO: BP, MF & CC
            - Protein complex: CORUM
            Change the value of term_number according to nb of terms wanted
            in plots.
        """
        conda:
            "../envs/stats/gprofiler.yml"
        wildcard_constraints:
            treat =str()
        log:
            RESFOL+"enrichment_GSEA/{treat}/gpro.log"
        input:
            data = RESFOL+"categorization_"+STATS_TOOL+"/Categorization_{treat}.tsv",
            genome = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene_biotype.tsv",
            gpro_names = "ressources/tables_resources/gprofiler_organisms.tsv"
        params:
            correction = CORRECTION,
            organism = ORGANISM,
            pval_threshold = PVAL,
            go_sources = GO_SOURCES,
            pathways = PATHWAY_SOURCES,
            complex_sources = COMPLEX,
            term_number = TERM_NUMBER,
            analysis = ANALYSIS
        output:
            main_folder = directory(RESFOL+"enrichment_GSEA/{treat}/"),
            go_output = directory(RESFOL+"enrichment_GSEA/{treat}/Gene_ontology/"),
            pathway = directory(RESFOL+"enrichment_GSEA/{treat}/Pathway/"),
            corum_output = directory(RESFOL+"enrichment_GSEA/{treat}/Corum/"),
            session_info = RESFOL+"enrichment_GSEA/{treat}/Session_gprofiler.RData"
        script:
            '../src/gprofiler.R'


    rule enrichment_summary_gprofiler:
        """
        Aim of the rule: Summary of the number of enrichment terms for each 
        category and source (KEGG, WIKIPATHWAY, REACTOME, GO: BP, MF & CC, and CORUM)
        """
        conda:
            "../envs/stats/ggplot2.yml"
        wildcard_constraints:
            treat =str()
        log:
            RESFOL+"enrichment_GSEA/{treat}/summary/gprofiler_summary.log"
        input:
            input_gpro = RESFOL+"enrichment_GSEA/{treat}/",
            gpro_names = "ressources/tables_resources/gprofiler_organisms.tsv",
            check = RESFOL+"enrichment_GSEA/{treat}/Session_gprofiler.RData"
        params:
            gpro_version = RESFOL+"enrichment_GSEA/gprofiler_api_version_info.txt",
            organism = ORGANISM,
            analysis = ANALYSIS
        output:
            output_gpro = RESFOL+"enrichment_GSEA/{treat}/summary/gprofiler_summary.csv",
            session_info = RESFOL+"enrichment_GSEA/{treat}/summary/Session_gprofiler_summary.RData"
        script:
            '../src/summary_gprofiler.R'


    rule ref_goatools:
        """
        Aim of the rule: Download reference files for hierarchy
            WARNINGS: GO hierarchy is updated daily than gprofiler
            databases is updated once per year. Differences can appear
            due to obsolence of some GO ids.
        """
        conda:
            "../envs/stats/goatools.yml"
        log:
            RESFOL+"enrichment_GSEA/goatools_reference_files/go_ref_hierarchy.log"
        input:
            check = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene_biotype.tsv"
        params:
            path_ref = RESFOL+"enrichment_GSEA/goatools_reference_files/",
            organismcode = ORGANISMCODE
        output:
            ref_hierarchy = RESFOL+"enrichment_GSEA/goatools_reference_files/BP_MF_CC_concise.txt",
            gene2go = RESFOL+"enrichment_GSEA/goatools_reference_files/gene2go",
            goobo = RESFOL+"enrichment_GSEA/goatools_reference_files/go-basic.obo"
        script:
            '../src/goatools_references.py'


    rule goatools:
        """
        Aim of the rule: Add hierarchy information to Gene Ontologies
            WARNINGS: GO hierarchy is updated daily than gprofiler
            databases is updated once per year. Differences can appear
            due to obsolence of some GO ids.
        """
        conda:
            "../envs/stats/goatools.yml"
        wildcard_constraints:
            treat =str()
        log:
            RESFOL+"enrichment_GSEA/{treat}/Gene_ontology_hierarchy/go_hierarchy_{cat}.log"
        input:
            ref_hierarchy = RESFOL+"enrichment_GSEA/goatools_reference_files/BP_MF_CC_concise.txt",
            gene2go = RESFOL+"enrichment_GSEA/goatools_reference_files/gene2go",
            goobo = RESFOL+"enrichment_GSEA/goatools_reference_files/go-basic.obo",
            check = RESFOL+"enrichment_GSEA/{treat}/summary/Session_gprofiler_summary.RData"
        params:
            gene_names = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene_biotype.tsv",
            table = RESFOL+"enrichment_GSEA/{treat}/Gene_ontology/GO_table_{cat}.tsv",
            path_ref = RESFOL+"enrichment_GSEA/{treat}/Gene_ontology_hierarchy",
            organismcode = ORGANISMCODE
        output:
            hierarchy_table = RESFOL+"enrichment_GSEA/{treat}/Gene_ontology_hierarchy/GO_table_{cat}_complete_hierarchy.tsv",
            depth_summary = RESFOL+"enrichment_GSEA/{treat}/Gene_ontology_hierarchy/summary_depth_{cat}_gene_ontology.tsv",
            level_summary = RESFOL+"enrichment_GSEA/{treat}/Gene_ontology_hierarchy/summary_level_{cat}_gene_ontology.tsv"
        script:
            '../src/goatools.py'


    rule enrichment_summary_goatools:
        """
        Aim of the rule: Summary of the number of terms according to the depth 
        and level of the hierarchy information for the Gene Ontologies
        """
        conda:
            "../envs/stats/ggplot2.yml"
        wildcard_constraints:
            treat =str()
        log:
            RESFOL+"enrichment_GSEA/{treat}/summary/goatools_summary.log"
        input:
            check = expand(RESFOL+"enrichment_GSEA/{{treat}}/Gene_ontology_hierarchy/GO_table_{cat}_complete_hierarchy.tsv", cat = CATEGORIES) # double braces required
        params:
            hierarchy_table = RESFOL+"enrichment_GSEA/{treat}/Gene_ontology_hierarchy/",
            analysis = ANALYSIS
        output:
            output_goa_depth = RESFOL+"enrichment_GSEA/{treat}/summary/goatools_summary_depth.csv",
            output_goa_level = RESFOL+"enrichment_GSEA/{treat}/summary/goatools_summary_level.csv",
            output_goa_depth_graph = RESFOL+"enrichment_GSEA/{treat}/summary/goatools_summary_depth_graph.png",
            output_goa_level_graph = RESFOL+"enrichment_GSEA/{treat}/summary/goatools_summary_level_graph.png",
            output_goa_depth_graph2 = RESFOL+"enrichment_GSEA/{treat}/summary/goatools_summary_depth_graph.pdf",
            output_goa_level_graph2 = RESFOL+"enrichment_GSEA/{treat}/summary/goatools_summary_level_graph.pdf",
            session_info = RESFOL+"enrichment_GSEA/{treat}/summary/Session_goatools_summary.RData"
        script:
            '../src/summary_goatools.R'


    #######################################################################

    def general_report_script():
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
        wildcard_constraints:
            treat =str()
        input:
            check = expand(RESFOL+"enrichment_GSEA/{treat}/summary/Session_goatools_summary.RData", treat = TREATMENT),
            metacontrast = METACONTRAST,
            mds_raw_data = RESFOL+"explore_counts_"+COUNT_TOOL+"/plotMDS_raw_data.png",
            biotype = REFFOL+RELEASE+"/"+ORGANISM+"/list_id_gene_biotype.tsv",
            cat_diff_genes = RESFOL+"categorization_"+STATS_TOOL+"/",
            logplot = RESFOL+"logplot/",
            diff_genes = RESFOL+"differential_expression_"+STATS_TOOL+"/"
        params:
            toolstat = STATS_TOOL,
            context = config["report"]["context"],
            comments = config["report"]["comments"],
            toolenrich = ENRICHMENT,
            aknow = config["report"]["acknowledgements"],
            author = config["report"]["author"],
            enrichment = RESFOL+"enrichment_GSEA/",
            filtering_method = FILTER_METHOD 
        output:
            "report/"+config["report"]["filename"]+".html"
        script:
            general_report_script()


    rule clean:
        """
        Aim of the rule: Clean Rplot in root folder.
        """
        input:
            check = "report/"+config["report"]["filename"]+".html"
        output:
            RESFOL+"cleaning.log"
        shell:
            "rm Rplots*.pdf && echo 'Cleaning success' > {output}"


#######################################################################
# Author: "Julie Ripoll, Céline Mandier"
# Contact: julie.ripoll87@gmail.com & celine.mandier@gmail.fr
# Created: 2021-12-09
# License: CeCILL
# Update: 2023-07-24
