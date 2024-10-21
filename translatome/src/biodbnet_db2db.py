#!/usr/bin/env python
# coding: utf-8

"""
bioDBnet db2db
==============

biological DataBase network is an application integrating a vast number of biological databases including Gene, UniProt, Ensembl, GO, Affy, RefSeq etc. The databases are created by downloading data from various public resources. They are formatted and maintained in a relational structure at the Advanced Biomedical Computing Center

https://biodbnet-abcc.ncifcrf.gov/

Here, we use the db2db offer. db2db handles all the conversions from one database identifier to another

Requirements
------------
- snakemake (only if used in the workflow)
- conda / miniconda or mamba for creation of the environment
- python 3
- argparse
- pandas
- urllib

Except snakemake & conda, programs can be installed using environment yaml files

Execution in snakemake rules
----------------------------

 ::

    shell:
        "../src/biodbnet_db2db.py "
        "--input_file {input} "
        "--identifier {params.identifier} "
        "--other_id {params.other_id} "
        "--organism {params.organismcode} "
        "--chunck {params.chunck} "
        "--biodb {output.biodb} "
        "1> {log}"

"""


import argparse
import pandas as pd
import urllib
from io import BytesIO


__author__ = "Julie Ripoll"
__created__ = "2021-11-24"
__license__ = "CeCILL"
__version__ = "0.1.0"


def get_args():
    """
    CLI with argparse
    """
    parser = argparse.ArgumentParser(description="""
    
    bioDBnet API request for db2db
    
    Requirements
    ------------
    - python3
    - argparse
    - urllib
    - pandas

    Add rigths
    ----------
    chmod +x biodbnet_db2db.py

    Examples
    --------
    biodbnet_db2db.py \
    --input_file path_to_file/file.tsv \
    --identifier "Ensembl Gene ID" \
    --other_id "KEEG Gene ID, Gene Info" \
    --organism 9606 \
    --chunck 399 \
    --biodb path_to_file/outfile.tsv \
    1> path_to_file/log.txt

    """)
    parser.add_argument('--input_file', '-inp', type=str,
                        required=True, help="FILE Input file with identifier to convert, required tab separator")
    parser.add_argument('--identifier', '-ide', type=str,
                        required=True, help='STRING type of identifier used during mapping can be "Ensembl Gene ID" or "Gene ID" etc')
    parser.add_argument('--other_id', '-oid', type=str,
                        required=True, help='STRING type of results from the db2db database wanted e.g. "KEGG Gene ID,Gene ID,UniProt Accession,KEGG Pathway Info, Ensembl Gene Info"')
    parser.add_argument('--organism', '-org', type=int,
                        required=True, help='INT value for organism, example Human = 9606')
    parser.add_argument('--chunck', '-chu', type=int,
                        required=True, help='INT value for chunck list of IDs in the url, max for ensembl ID = 399')
    parser.add_argument('--biodb', '-bdb', type=str,
                        required=True, help='FILE output file with bioDBnet conversion of IDs')

    args = parser.parse_args()
    return args


def importIDlist(inputfile):
    """
    Import of dataframe using pandas
    Data's identifier should be specified as a column titled 'identifier'
    """
    data = pd.read_csv(inputfile, sep ="\t", header = 0)
    try:
        if "identifier" in data.columns:
            print(data)
    except NameError:
        raise NameError("""Oops!  Please renamed your input label column as: identifier """)
    return data


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def createURL(genes):
    """
    Function to create url
    """
    return urllib.parse.unquote(
        hostname + urllib.parse.urlencode({'method': 'db2db',
                                           'format': 'row',
                                           'input': INPUTDB,
                                           'inputValues': genes,
                                           'outputs': OUTPUTSDB,
                                           'taxonId': ORGANISM}))


def chunkURLCall(lst, n):
    """
    Function to request url and export as dataframe
        chunks multiple api calls and returns a single dataframe
    """
    acc = []
    for chunk in chunks(lst, n):
        acc.append(pd.read_json(createURL(",".join(chunk))))
    return pd.concat(acc, ignore_index=True)


if __name__ == '__main__':
    
    args = get_args()
    print(args)

    # Import data
    data = importIDlist(args.input_file)

    #####################################################################

    # API bioDBnet using urllib python package

    # Fix parameters for query
    ORGANISM = args.organism
    INPUTDB = args.identifier
    OUTPUTSDB = args.other_id

    # biodbnet url
    hostname = 'https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?'

    #####################################################################

    # Request

    # Filter null ids
    data_nonull = data[data["identifier"].notnull()]
    # No more than 500 ID in the url !!!
    df_final = chunkURLCall(data_nonull["identifier"], args.chunck)
    df_final = df_final.rename(columns={"InputValue": "identifier"})
    print(df_final)

    # Merge biodbnet ids with test df
    df_merge = pd.merge(data, df_final, on = ["identifier"])
    print(df_merge.head())

    # Parse Gene Info column if exists
    if 'Ensembl Gene Info' in df_merge.columns:
        lastCol = df_merge['Ensembl Gene Info'].str.findall("\[([^:]*): *([^\]]+)\]").apply(dict).apply(pd.Series)
        res = pd.concat([df_merge, lastCol], sort = False, axis = 1)
        res.drop('Ensembl Gene Info', axis = 1, inplace = True)
        res.to_csv(args.biodb, sep = "\t", header = True, index = False)
    elif 'Gene Info' in df_merge.columns:
        lastCol = df_merge['Gene Info'].str.findall("\[([^:]*): *([^\]]+)\]").apply(dict).apply(pd.Series)
        res = pd.concat([df_merge, lastCol], sort = False, axis = 1)
        res.drop('Gene Info', axis = 1, inplace = True)
        res.to_csv(args.biodb, sep = "\t", header = True, index = False)
    else:
        df_merge.to_csv(args.biodb, sep = "\t", header = True, index = False)


#####################################################################
# END
