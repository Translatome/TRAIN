#!/usr/bin/env python
# coding: utf-8

"""
Merge counts obtained from featureCounts
========================================

Requirements
------------
- snakemake with python 3 & conda / miniconda or mamba
- pandas
- searborn
- glob

Except snakemake & conda, programs can be installed using environment yaml files

Execution in snakemake rule
---------------------------

 ::

    script:
        "../src/merge_counts_feature.py"

"""

import os
import pandas as pd
import glob
import seaborn as sns
from wrappers.format_checker import is_gff3_ext


__author__ = "Julie Ripoll"
__created__ = "2021-11-24"
__license__ = "CeCILL"
__version__ = "0.1.0"


def get_id(path):
    """
    Separate path and basename
    """
    return os.path.basename(path).rsplit(".", 2)[0]


def get_file_df_gff(path):
    """
    Get featurecounts file and prepare it for use
        Counts were performed using a GFF annotation file
    """
    name = get_id(path)
    test = pd.read_csv(path,
                        sep='\t',
                        comment='#',
                        names=['Geneid','Chr','Start','End','Strand','Length','Counts'],
                        index_col=['Geneid']) # error if RefSeq
    test.drop(test.head(1).index, inplace=True)
    test[name] = test['Counts']
    return test[name]


def get_file_df_gtf(path):
    """
    Get featurecounts file and prepare it for use
        Counts were performed using a GTF annotation file
    """
    name = get_id(path)
    test = pd.read_csv(path, sep='\t', 
                    comment='#', 
                    names=['Geneid','Chr','Start','End','Strand','Length','Gene_name', 'Counts'],
                    index_col=['Geneid'])
    test.drop(test.head(1).index, inplace=True)
    test[name] = test['Counts']
    return test[name]


if __name__ == '__main__':
    # All snakemake declarations
    
    INPUT_PATH = snakemake.params["inputdir"]
    INPUT_ANNOT = snakemake.params["annotation"]
    OUTPUT_COUNTS = snakemake.output["counts"]
    OUTPUT_HEATMAP = snakemake.output["heatmap"]

    #####################################################################
    # Import counts and annotation files

    try:
        path = str(INPUT_PATH)
        all_filenames = glob.glob(path+'*.tsv')
        print(all_filenames)
    except FileNotFoundError:
        raise FileNotFoundError(f"""Oops!  Please provide a path to the folder containing count files in '.tsv' format.  Try again...
                                Provided: {path}*.tsv""")

    annot = str(INPUT_ANNOT)

    ## Combine all files in the list
    if is_gff3_ext(annot):
        series = [get_file_df_gff(f) for f in all_filenames]
        print(series)
    else:
        series = [get_file_df_gtf(f) for f in all_filenames]
        print(series)

    for f in series:
        print("Check if index is unique")
        print(f.index.is_unique)

    ## Combine all last columns
    combined_tsv = pd.concat([f for f in series if f.size > 0], axis=1)
    print("Combined counts")
    print(combined_tsv.head())
    print(type(combined_tsv))

    ## Export to csv
    combined_tsv.to_csv(OUTPUT_COUNTS, encoding='utf-8-sig', sep ="\t")

    #####################################################################

    # Heatmap of counts
    
    df = combined_tsv.astype(int)

    ## Plot distribution
    df.plot(figsize=(20,20), legend=False)

    ## delete empty rows
    dfWTNA = df.loc[(df!=0).any(axis=1)] 
    
    ## check counts 
    if len(dfWTNA) > 2:
        print(dfWTNA)
    else:
        raise ValueError(f"""Warning: all counts are zero, check params of Feature Counts""")

    ## na omit
    dfna = dfWTNA.fillna("0")
    print(dfna)

    ## heatmap
    fig = sns.clustermap(dfna.corr(), figsize=(12,12))

    ## save heatmap
    fig.savefig(OUTPUT_HEATMAP)

#####################################################################
# END
