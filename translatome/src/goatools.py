#!/usr/bin/env python
# coding: utf-8

"""
Gene Ontology hierarchy
=======================

Aim
---
Filter complete DAG hierarchy by GO terms found with gprofiler

Requirements
------------
- snakemake with python 3 & conda / miniconda or mamba
- pandas
- wget
- goatools
- certifi # for secure ssl

Except snakemake & conda, programs can be installed using environment yaml files

Execution in snakemake rule
---------------------------

 ::

    script:
        "../src/goatools.py"

"""


import os
import wget
import pandas as pd
from contextlib import redirect_stdout
import goatools as goatools


__authors__ = "Julie Ripoll & Céline Mandier"
__created__ = "2023-02-22"
__license__ = "CeCILL"
__version__ = "0.1.0"


if __name__ == '__main__':
    # Use scripts from goatools
    # but it's not a module
    from goatools.obo_parser import GODag
    from goatools.rpt.rpt_lev_depth import RptLevDepth
    from goatools.cli import wr_hierarchy

    # All snakemake values

    LOG_FILE = snakemake.log[0]

    ## input files
    GENE2GO = snakemake.input["gene2go"]
    GOOBO = snakemake.input["goobo"]
    HIERA_REF = snakemake.input["ref_hierarchy"]
    INPUT_FILE = snakemake.params["table"]

    ## params
    GENENAMES = snakemake.params["gene_names"]
    ORGANISMCODE = snakemake.params["organismcode"]
    PATH_REF = snakemake.params["path_ref"]

    ## outputs
    HIERA_TAB = snakemake.output["hierarchy_table"]
    DEPTH_TAB = snakemake.output["depth_summary"]
    LEVEL_TAB = snakemake.output["level_summary"]


    #####################################################################

    with open(LOG_FILE, 'w') as f:
        with redirect_stdout(f):

            # Import gprofiler files
            data = pd.read_csv(INPUT_FILE, sep="\t")

            if 'x' in data.columns[0]:
                print("This category file is empty, pass")
                data_res_wo_hier = pd.DataFrame()
                depth_table = pd.DataFrame()
                level_table = pd.DataFrame()
                data_res_wo_hier.to_csv(HIERA_TAB, sep="\t", index=None)
                depth_table.to_csv(DEPTH_TAB, sep="\t", index=True)
                level_table.to_csv(LEVEL_TAB, sep="\t", index=True)
            else:
                print(data.head())

                ## keep only GO ids
                term_ids = pd.DataFrame(data["term_id"])
                term_ids.to_csv(PATH_REF+"/terms.txt")
                print(term_ids)

                ## parse complete DAG correctly
                os.system('bash translatome/src/go_format.sh '
                          +HIERA_REF+' '+PATH_REF+
                        '/Complete_GO_table_hierarchy.txt')

                ## read GO hierar
                data_hier = pd.read_csv(
                    PATH_REF+"/Complete_GO_table_hierarchy.txt", sep="\t", header=None)
                data_hier = data_hier.drop(columns={2, 3, 4, 5, 8})
                data_hier = data_hier.drop_duplicates()
                data_hier.columns = ["term_id", "GO", "Level", "Depth"]
                
                ## filter on term ids
                data_filt = data_hier[data_hier["term_id"].isin(term_ids["term_id"])]               
                print("----------------------------------------")
                print("Number of terms output from gprofiler:", len(data))
                print("Number of terms output from goatools:", len(data_filt))
                print("Terms with an obsolete term id:", len(data)-len(data_filt))

                #####################################################################

                # Obsolete terms, differences should be due to obsolete ids
                ### https://wiki.geneontology.org/index.php/Obsoleting_an_Existing_Ontology_Term
                data_res_wo_hier = data_filt.merge(
                    data, on = ["term_id"], how = "outer")
                data_res_wo_hier.fillna("obsolete_terms", inplace=True)
                print(data_res_wo_hier)

                #####################################################################

                # summary of depth and level

                ## depth
                depth_table = pd.DataFrame(data_res_wo_hier.groupby(
                    "GO")["Depth"].value_counts().sort_index(level=['GO', 'Depth']))
                depth_table.to_csv(DEPTH_TAB, sep="\t", index=True)

                ## level
                level_table = pd.DataFrame(data_res_wo_hier.groupby(
                    "GO")["Level"].value_counts().sort_index(level=['GO', 'Level']))
                level_table.to_csv(LEVEL_TAB, sep="\t", index=True)

                #####################################################################

                # Add gene symbol info for intersections
                
                ## Import biotype file with gene symbol (require step 00_download_references)
                names = pd.read_csv(GENENAMES, sep="\t")
                
                # Initialize an empty dataframe for results storage
                table_id_symbol = pd.DataFrame(columns=["term_id", "Gene_Symbol"])

                # Iterate by length
                for i in range(len(data_res_wo_hier)):
                    # Intersection column in list of identifier
                    inters = data_res_wo_hier.loc[i, 'intersection'].split(',')
                    dt_inters = pd.DataFrame({"identifier": inters})
                    # Merge with gene symbol
                    result_inter = pd.merge(dt_inters, names.loc[:, ["identifier","Gene Symbol"]], on="identifier")
                    # Create a line for term i
                    line_id_ensembl = pd.DataFrame({
                        "term_id": [data_res_wo_hier.loc[i, 'term_id']],
                        "Gene_Symbol": [",".join(result_inter['Gene Symbol'])]
                    })
                    # Concat dataframes
                    table_id_symbol = pd.concat([table_id_symbol, line_id_ensembl], ignore_index=True)

                # Merge goatools results with the df containing the gene symbol column
                table_final = pd.merge(data_res_wo_hier, table_id_symbol, on="term_id")

                ## Export results
                table_final.to_csv(HIERA_TAB, sep="\t", index=None)


#####################################################################
# END
