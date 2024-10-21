#!/usr/bin/env python
# coding: utf-8

"""
Gene Ontology hierarchy
=======================

Aim
---
Download complete DAG hierarchy for latest version of gene ontology

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
        "../src/goatools_references.py"

"""


import os
import wget
import pandas as pd
from contextlib import redirect_stdout
import goatools as goatools


__authors__ = "Julie Ripoll & Fati Chen"
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

    ## params
    PATH_REF = snakemake.params["path_ref"]
    ORGANISMCODE = snakemake.params["organismcode"]

    ## outputs
    HIERA_REF = snakemake.output["ref_hierarchy"]
    GENE2GO = snakemake.output["gene2go"]
    GOOBO = snakemake.output["goobo"]


    #####################################################################

    with open(LOG_FILE, 'w') as f:
        with redirect_stdout(f):

            # Prepare hierarchy files

            ## import latest version of Gene Ontology files
            print("Download of go-basic.obo and gene2go files for analysis")

            if not os.path.exists(PATH_REF):
                os.makedirs(PATH_REF)

            ## go basic
            if not os.path.isfile(PATH_REF+"/go-basic.obo"):
                wget.download('http://purl.obolibrary.org/obo/go/go-basic.obo',
                                  out = (PATH_REF+"/go-basic.obo"))

            ## gene2go
            ### readme
            if not os.path.isfile(PATH_REF+"/README_gene2go"):
                wget.download("https://ftp.ncbi.nih.gov/gene/DATA/README",
                          out = (PATH_REF+"/README_gene2go"))
            ### file
            if not os.path.isfile(PATH_REF+"/gene2go.gz"):
                wget.download("https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz",
                                    out = (PATH_REF+"/gene2go.gz"))
            ### decompress
            os.system('gzip -f -d '+PATH_REF+'/gene2go.gz')

            ## initialize GO DAG object
            obodag = GODag(PATH_REF+"/go-basic.obo")

            ### initialize reporter class
            rptobj = RptLevDepth(obodag)

            ### depth level
            print("Depth and Level of Gene Ontology hierarchy downloaded")
            print(rptobj.write_summary_cnts_all())

            # Create hierarchy file for gene ontology
            
            ## use class from goatools wr_hierarchy.py
            objcli = wr_hierarchy.WrHierCli([
                    "BP", "MF", "CC",
                    '--gene2go='+GENE2GO,
                    '--dag='+GOOBO,
                    '--taxid='+ORGANISMCODE,
                    '--dash_len=17',
                    '--concise',
                    '-o '+HIERA_REF,
                ])
            
            ## required code from goatools wr_hierarchy.py
            fouts_txt = objcli.get_fouts()
            print(fouts_txt)
            if fouts_txt:
                for fout_txt in fouts_txt:
                    objcli.wrtxt_hier(fout_txt)
                print(fouts_txt)


#####################################################################
# END
