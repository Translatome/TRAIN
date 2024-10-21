#!/usr/bin/env python
# coding: utf-8

"""
Extract Identifier from GFF/GTF annotation file
===============================================

Requirements
------------
- snakemake with python 3 & conda / miniconda or mamba
- pandas
- BCBio

Except snakemake & conda, programs can be installed using environment yaml files

Execution in snakemake rule
---------------------------

 ::

    script:
        "../src/ExtractIDfromGFF.py"

"""


import os
import pandas as pd
from glob2 import iglob
from contextlib import redirect_stdout


__authors__ = "Julie Ripoll"
__created__ = "2023-05-12"
__license__ = "CeCILL"
__version__ = "0.1.0"


def attrParser(attribute, separator=" "):
    """
    Function to parse attributes in GFF/GTF file
    """
    attribute = attribute.strip()
    attrDict = {}
    buffer = ''
    key = ''
    quoted = False
    for char in attribute:
        if char == '"':
            quoted = not quoted
        elif not quoted and char == ';':
            attrDict[key] = buffer
            buffer = ''
            key = ''
        elif not quoted and char == separator:
            if buffer != '':
                key = buffer
                buffer = ''
        else:
            buffer = buffer + char
    return attrDict


if __name__ == '__main__':
    # All snakemake values
    
    LOG_FILE = snakemake.log[0]

    # input files
    INPUT_PATH = snakemake.params["pathto"]
    EXTENSION = snakemake.params["extensionAnnotation"]

    # output files
    OUTPUT = snakemake.output["outfile"]

    #####################################################################

    with open(LOG_FILE, 'w') as f:
        with redirect_stdout(f):

            annotation_file = pd.read_csv(next(iglob(INPUT_PATH+"*."+EXTENSION)),
                            sep = '\t',
                            comment = '#',
                            header = None,
                            names = ["seqid","source", "type", "start", "end", "score","strand","phase","attributes"],
                            low_memory = False)
            print("Annotation file downloaded:")
            print(annotation_file)

            # "=" pour gff et " " pour gtf
            if EXTENSION == "gtf":
                acc = annotation_file['attributes'].apply(attrParser, separator = " ")
                lastCol = pd.DataFrame(acc.tolist())
            elif EXTENSION == "gff3":
                acc = annotation_file['attributes'].apply(attrParser, separator = "=")
                lastCol = pd.DataFrame(acc.tolist())
            else:
                print("Other formats than GFF3 and GTF are not supported here, please create the output file with the genome identifier yourself for the next step")

            # Concat gff/gtf and parsed attributes
            annotation_file = pd.concat([annotation_file, lastCol], sort = False, axis = 1)
            print("Annotation file with parsed attributes:")
            print(annotation_file)

            # extract identifier
            genes = pd.DataFrame(annotation_file["gene_id"])
            genes.rename(columns = {"gene_id": "identifier"}, inplace = True)
            genes.drop_duplicates(inplace = True)
            print("Genes identifier")
            print(genes)

            # export file
            genes.to_csv(OUTPUT, sep = "\t", header = True, index = None)


#####################################################################
# END
