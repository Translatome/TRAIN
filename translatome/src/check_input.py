"""
Check input files format
========================

Requirements
------------
- snakemake with python 3 & conda / miniconda or mamba
- pandas
- glob
- wrappers.format_checker

Except snakemake & conda, programs can be installed using environment yaml files

Execution in snakemake rule
---------------------------

 ::

    script:
        "../src/check_input.py"

"""

import os
import glob
from contextlib import redirect_stdout

import pandas as pd

import wrappers.format_checker as fcheck


__author__ = "Julie Ripoll"
__created__ = "2023-04-25"
__license__ = "CeCILL"
__version__ = "0.1.0"


def check_input(
    INPUT_DATA, INPUT_METADATA, INPUT_METACONTRAST, INPUT_REF, ANALYSIS, OUTPUT_FILE
):
    """
    Select functions to execute according to the type of analysis, i.e. "pol-seq" or "kinetic"
    """
    with open(OUTPUT_FILE, "w") as f:
        with redirect_stdout(f):
            check_fastq_files(INPUT_DATA)

            if ANALYSIS == "pol-seq":
                check_input_metadata(INPUT_METADATA)

                check_input_metacontrast(INPUT_METACONTRAST)
            elif ANALYSIS == "kinetic":
                check_kinetic_input_metadata(INPUT_METADATA)

                check_kinetic_input_metacontrast(INPUT_METACONTRAST)

            check_reference_files(INPUT_REF)


def check_reference_files(INPUT_REF):
    """
    Check if fasta and GFF/GTF downloaded files are not corrupted.
    """
    if os.path.exists(INPUT_REF):
        # import all files with glob
        all_filenames = glob.glob(str(INPUT_REF) + "/*/*")
        # imported files
        print("Import of all files in: " + INPUT_REF + "/*/*")
        print(all_filenames)
        print("--------------------------------------------")

        success_message = """{path} seems to be well formatted and ready to use
        Please include the name of the desired file in your config.yml file, in the appropriate section
        ------------------------------------------------------------------------------------------------"""

        # check format, here should be fasta files
        for path in all_filenames:
            if fcheck.is_fasta_ext(path):
                # check the contents of the FASTA file
                if fcheck.is_fasta(path):
                    print(success_message.format(path=path))
                else:
                    raise ValueError(
                        f"""Error please fix this:
                        This fasta file is potentially corrupted: {path}"""
                    )
            elif fcheck.is_gff3_ext(path):
                # check the contents of the GFF file
                if fcheck.is_gff3_gtf(path):
                    print(success_message.format(path=path))
                else:
                    raise ValueError(
                        f"""Error please fix this:
                        This GFF file is potentially corrupted: {path}"""
                    )
            elif fcheck.is_gtf_ext(path):
                # check the contents of the GTF file
                if fcheck.is_gff3_gtf(path):
                    print(success_message.format(path=path))
                else:
                    raise ValueError(
                        f"""Error please fix this:
                        This GTF file is potentially corrupted: {path}"""
                    )
            else:
                print(f"{path}: is not useful for the workflow")
                print(
                    "------------------------------------------------------------------------------------------------"
                )
    else:
        raise FileNotFoundError(
            f"""Error please fix this:
            Please provide a path to the folder containing references files. Provided: {INPUT_REF}/* """
        )


def check_kinetic_input_metacontrast(INPUT_METACONTRAST):
    """
    Check if the metacontrast file is well formatted for the selected analysis, here "kinetic".
    """
    try:
        fcheck.raise_on_missing_file(INPUT_METACONTRAST)

        # read file with pandas
        metacontrast = pd.read_csv(INPUT_METACONTRAST, sep="\t", header=0)
        # check file format conditions
        try:
            fcheck.raise_bad_metacontrast_kinetic(metacontrast)
        except ValueError as error:
            print(repr(error))
        else:
            print(
                """The Metacontrast file seems to be well formatted and ready to use
                Note:
                - cytoplasmic or total RNA-seq should be the first element of your metacontrast file after the title
                - 'control' should be in the 'treatment2' column
                Example:
                --------
                fraction     treatment1  time  fraction_bis  treatment2  time_bis
                Cytoplasmic  treated     D1    Cytoplasmic   control     D1
                Polysomal    treated     D1    Polysomal     control     D1
                Cytoplasmic  treated     D2    Cytoplasmic   control     D2
                Polysomal    treated     D2    Polysomal     control     D2"""
            )
            print("--------------------------------------------------------")
    except FileNotFoundError as error:
        print(repr(error))


def check_kinetic_input_metadata(INPUT_METADATA):
    """
    Check if the metadata file is well formatted for the selected analysis, here "kinetic".
    """
    try:
        # Check metadata file
        fcheck.raise_on_missing_file(INPUT_METADATA)
        # read file with pandas
        metadata = pd.read_csv(INPUT_METADATA, sep="\t", header=0)
        # check file format conditions

        try:
            fcheck.raise_on_metadata_kinetic(metadata)
        except ValueError as error:
            print(repr(error))
        else:
            print("The Metadata file is well formatted and ready to use")
            print("----------------------------------------------------")
    except FileNotFoundError as error:
        print(repr(error))


def check_input_metacontrast(INPUT_METACONTRAST):
    """
    Check if the metacontrast file is well formatted for the selected analysis, here "pol-seq".
    """
    try:
        # Check metacontrast file
        fcheck.raise_on_missing_file(INPUT_METACONTRAST)
        # read file with pandas
        metacontrast = pd.read_csv(INPUT_METACONTRAST, sep="\t", header=0)
        # check file format conditions
        try:
            fcheck.raise_bad_metacontrast_format(metacontrast)
        except ValueError as error:
            print(repr(error))
        else:
            print("The Metacontrast file is well formatted and ready to use")
            print("--------------------------------------------------------")
    except FileNotFoundError as error:
        print(repr(error))


def check_input_metadata(INPUT_METADATA):
    """
    Check if the metadata file is well formatted for the selected analysis, here "pol-seq".
    """
    try:
        fcheck.raise_on_missing_file(INPUT_METADATA)

        metadata = pd.read_csv(INPUT_METADATA, sep="\t", header=0)

        try:
            # check file format conditions
            fcheck.raise_bad_metadata_format(metadata)
        except ValueError as error:
            print(repr(error))
        else:
            print("The Metadata file is well formatted and ready to use")
            print("----------------------------------------------------")
    except FileNotFoundError as error:
        print(repr(error))


def check_fastq_files(INPUT_DATA):
    """
    Check if fastq files are not corrupted.
    """
    if os.path.exists(INPUT_DATA):
        # import all files with glob
        path = str(INPUT_DATA)
        all_filenames = glob.glob(path + "*")
        # imported files
        print("Import of all files in: " + INPUT_DATA)
        print(all_filenames)
        # check format, here should be fastq files
        for path in all_filenames:
            if fcheck.is_fastq_ext(path):
                # check the contents of the file
                if fcheck.is_fastq(path):
                    print(f"{path} seems to be well formatted and ready to use")
                    print("---------------------------------------------------")
                else:
                    raise ValueError(
                        f"""Error please fix this:
                        This fastq file is potentially corrupted: {path}"""
                    )
    else:
        raise FileNotFoundError(
            f"""Error please fix this:
            Please provide a path to the folder containing fastq files. Provided: {INPUT_DATA}* """
        )


if __name__ == "__main__":
    # All snakemake values

    ## input files
    INPUT_DATA = snakemake.input["data"]
    INPUT_METADATA = snakemake.input["metadata"]
    INPUT_METACONTRAST = snakemake.input["metacontrast"]
    INPUT_REF = snakemake.params["references"]

    ## params
    ANALYSIS = snakemake.params["analysis"]

    ## output files
    OUTPUT_FILE = snakemake.output[0]

    #####################################################################

    check_input(
        INPUT_DATA, INPUT_METADATA, INPUT_METACONTRAST, INPUT_REF, ANALYSIS, OUTPUT_FILE
    )


#####################################################################
# END
