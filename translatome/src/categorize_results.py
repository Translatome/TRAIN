#!/usr/bin/env python
# coding: utf-8

"""
Categorization translation/transcription/both
=============================================

Requirements
------------
- snakemake with python 3 & conda / miniconda or mamba
- pandas
- urllib
- certifi # for secure ssl

Except snakemake & conda, programs can be installed using environment yaml files

Execution in snakemake rule
---------------------------

 ::

    script:
        "../src/explore_results.py"

"""

import pandas as pd
import math
import urllib
import urllib.request
import urllib.parse
from contextlib import redirect_stdout


__authors__ = "Julie Ripoll and Céline Mandier"
__created__ = "2022-06-20"
__license__ = "CeCILL"
__version__ = "0.1.0"


def importFiles(
        path,
        file1,
        file2,
        extension,
        separator,
        decimalvalue
    ):
    """
    Import & merge & rename results from statistical analysis
        translation -> y / transcription -> x

    Args:
        path (_str_): path to files
        file1 (_object_): file _ results from transcription analysis
        file2 (_object_): file _ results from translation analysis
        extension (_str_): ".tsv"
        separator (_str_): tabulation
        decimalvalue (_str_): point

    Returns:
        _object_: pandas dataframe

    """
    transcription = pd.read_csv(path+"/"+str(file1)+extension,
                                sep=separator,
                                comment='#',
                                decimal=decimalvalue)
    translation = pd.read_csv(path+"/"+str(file2)+extension,
                              sep=separator,
                              comment='#',
                              decimal=decimalvalue)
    mergedata = pd.merge(transcription,
                         translation,
                         how="outer",
                         on="identifier",
                         suffixes=("_transcription", "_translation"))
    print("Lenght of transcription : ", len(transcription))
    print("Lenght of translation : ", len(translation))
    print("Lenght of merge dataframe : ", len(mergedata))
    mergedata.index = mergedata["identifier"]
    mergedata.drop("identifier", inplace=True, axis="columns")
    return mergedata


def categorization(data, pval_th, logFC_th):
    """
    Categorize df according to translation, transcription and both, with up or down logFC.
    """

    # define test
    def is_significative(attr):
        return pd.Series(attr) <= float(pval_th)

    def is_up(attr):
        return pd.Series(attr) >= float(logFC_th)

    def is_down(attr):
        return pd.Series(attr) <= -float(logFC_th)

    # copy df to avoid errors
    res = data.copy()
    if COUNT_TOOL == "deseq2":
        # categorize if transcription or translation is significative
        translation = {
            'logFC': res["log2FoldChange_translation"], 'pval': res["padj_translation"]}
        transcription = {
            'logFC': res["log2FoldChange_transcription"], 'pval': res["padj_transcription"]}
    elif COUNT_TOOL == "limma":
        # categorize if transcription or translation is significative
        translation = {
            'logFC': res["logFC_translation"], 'pval': res["adj.P.Val_translation"]}
        transcription = {
            'logFC': res["logFC_transcription"], 'pval': res["adj.P.Val_transcription"]}
    else:
        print(
            "Please add information on the tool used for statistics in the configuration file")

    # apply test
    res["IsTranslation"] = is_significative(translation['pval']) & (
        is_up(translation['logFC']) | is_down(translation['logFC']))
    res["IsTranscription"] = is_significative(transcription['pval']) & (
        is_up(transcription['logFC']) | is_down(transcription['logFC']))
    # categorization by type
    res["Both_mRNA"] = (res["IsTranslation"] == True) & (
        res["IsTranscription"] == True)
    res["Translation_only"] = (res["IsTranslation"] == True) & (
        res["IsTranscription"] == False)
    res["Transcription_only"] = (res["IsTranslation"] == False) & (
        res["IsTranscription"] == True)
    res["Not_significant"] = (res["IsTranslation"] == False) & (
        res["IsTranscription"] == False)
    # Summary
    res.loc[(res["Both_mRNA"] & is_up(translation['logFC']) & is_up(
        transcription['logFC'])), "Categories"] = "Both_mRNA_UP"
    res.loc[(res["Both_mRNA"] & is_down(translation['logFC']) & is_down(
        transcription['logFC'])), "Categories"] = "Both_mRNA_DOWN"
    res.loc[(res["Translation_only"] & is_up(
        translation['logFC'])), "Categories"] = "Translation_UP"
    res.loc[(res["Translation_only"] & is_down(
        translation['logFC'])), "Categories"] = "Translation_DOWN"
    res.loc[(res["Transcription_only"] & is_up(
        transcription['logFC'])), "Categories"] = "Transcription_UP"
    res.loc[(res["Transcription_only"] & is_down(
        transcription['logFC'])), "Categories"] = "Transcription_DOWN"
    res.loc[(res["Both_mRNA"] & (is_up(translation['logFC']) & is_down(
        transcription['logFC']))), "Categories"] = "Divergent_UPDN"
    res.loc[(res["Both_mRNA"] & (is_down(translation['logFC']) & is_up(
        transcription['logFC']))), "Categories"] = "Divergent_DNUP"
    res.loc[res["Not_significant"],
            "Categories"] = "Not_significant"
    return res


if __name__ == '__main__':
    # All snakemake values

    LOG_FILE = snakemake.log[0]

    # input files
    INPUT_PATH = snakemake.input["input_path"]
    GENENAMES = snakemake.input["gene_names"]
    METACONTRAST = snakemake.input["metacontrast"]

    # fix parameters
    ANALYSIS = snakemake.params["analysis"]
    ORGANISM = snakemake.params["organismcode"]
    INPUTDB = snakemake.params["inputbiodb"]
    OUTPUTSDB = snakemake.params["outputbiodb"]
    COUNT_TOOL = snakemake.params["config_tool"]
    PVAL_THRES = snakemake.params["pval_threshold"]
    LOGF_THRES = snakemake.params["logFC_threshold"]

    # output files
    OUTPUT_PATH = snakemake.output[0]

    # biodbnet url
    hostname = 'https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?'

    #####################################################################

    with open(LOG_FILE, 'w') as f:
        with redirect_stdout(f):

            # Import files

            meta_contrast = pd.read_csv(METACONTRAST, sep="\t")

            if ANALYSIS == "kinetic":
                files_names = meta_contrast.apply(lambda d: "_".join(d.array[:3]) + "_vs_" + "_".join(
                    d.array[3:]), axis=1).groupby(lambda x: math.floor(x/2)).apply(lambda d: list(d.array))
            else:
                files_names = meta_contrast.apply(lambda d: "_".join(d.array[:2]) + "_vs_" + "_".join(
                    d.array[2:]), axis=1).groupby(lambda x: math.floor(x/2)).apply(lambda d: list(d.array))

            files_names = files_names.to_list()

            files_df = pd.DataFrame(files_names, columns=["file1", "file2"])
            print(files_df)

            try:
                for row in files_df.index:

                    name1 = files_df['file1'][row]
                    print(name1)
                    name2 = files_df['file2'][row]
                    print(name2)

                    data = importFiles(INPUT_PATH, name1,
                                       name2, ".tsv", "\t", ".")
                    print(data.head())

                    # check index is unique
                    print("Index is unique ?")
                    print(data.index.is_unique)

                    # Use function to categorize
                    test = categorization(data, PVAL_THRES, LOGF_THRES)
                    print(test.head())

                    # check categories
                    print(test['Categories'].describe())

                    # check number of genes involved in translation processes
                    print(test["IsTranslation"].describe())

                    # check number of genes involved in transcription processes
                    print(test["IsTranscription"].describe())

                    # check number of genes involved in both processes
                    print(len(test[test["Both_mRNA"] == True]))

                    # check number of genes involved in translation only
                    print(len(test[test["Translation_only"] == True]))

                    # check number of genes involved in transcription only
                    print(len(test[test["Transcription_only"] == True]))

                    # check number of non significative genes
                    print(len(test[test["Not_significant"] == True]))

                    # Annotation of identifier with biodb databases
                    # import of gene_biotype file
                    names = pd.read_csv(GENENAMES, sep="\t")

                    # merge biodbnet ids with test df
                    df_merge = pd.merge(test, names, on=["identifier"])
                    df_merge.index = test.index

                    ###################################################################
                    # Export (+filter unused column)

                    # drop columns
                    if COUNT_TOOL == "deseq2":
                        searchfor = ["lfcSE_", "baseMean_", "stat_", "pvalue_"]
                        indexNames = df_merge.loc[:, df_merge.columns.str.contains(
                            '|'.join(searchfor))]
                        df_merge.drop(indexNames, inplace=True, axis=1)
                        print(df_merge.head())
                    elif COUNT_TOOL == "limma":
                        searchfor = ["AveExpr_", "t_", "B_", "P.Value_"]
                        indexNames = df_merge.loc[:, df_merge.columns.str.contains(
                            '|'.join(searchfor))]
                        df_merge.drop(indexNames, inplace=True, axis=1)
                        print(df_merge.head())
                    else:
                        print(
                            "Please add information on the tool used for statistics in the configuration file")

                    # build name
                    step = files_df['file1'][row].split("_")

                    # save results
                    if ANALYSIS == "kinetic":
                        df_merge.to_csv(OUTPUT_PATH+"/Categorization_"+step[1]+"vs"+step[5]+"_"+step[2]+"vs"+step[6]+".tsv",
                                        sep="\t",  # or ";"
                                        header=True,
                                        index=False)
                    else:
                        df_merge.to_csv(OUTPUT_PATH+"/Categorization_"+step[1]+"vs"+step[4]+".tsv",
                                        sep="\t",  # or ";"
                                        header=True,
                                        index=False)

                    # summary of categories
                    cat_vec = ['Not_significant',
                               'Total_sig_genes',
                               'Both_mRNA_UP',
                               'Both_mRNA_DOWN',
                               'Transcription_UP',
                               'Transcription_DOWN',
                               'Translation_UP',
                               'Translation_DOWN',
                               'Divergent_UPDN',
                               'Divergent_DNUP']
                    df_merge["Categories"] = pd.Categorical(
                        df_merge["Categories"], categories=cat_vec)
                    resTab = df_merge["Categories"].value_counts()
                    resTab.loc['Total_sig_genes'] = resTab.drop(
                        'Not_significant', axis=0).sum(axis=0)
                    resTab = resTab.reindex(cat_vec)
                    print(resTab)
                    # save summary
                    if ANALYSIS == "kinetic":
                        resTab.to_csv(OUTPUT_PATH+"/Summary_categories_"+step[1]+"vs"+step[5]+"_"+step[2]+"vs"+step[6]+".tsv",
                                      sep="\t",
                                      header=True,
                                      index=True)
                    else:
                        resTab.to_csv(OUTPUT_PATH+"/Summary_categories_"+step[1]+"vs"+step[4]+".tsv",
                                      sep="\t",
                                      header=True,
                                      index=True)

                    print("----------------------------------------------------")

            except FileNotFoundError:
                raise FileNotFoundError(f"""
                Oops!  That was no valid name of file.
                ----------------
                Names searched:
                {print(files_df)}
                ----------------
                Info: Input files should be tabulated with "." as decimal
                """)


#####################################################################
# END
