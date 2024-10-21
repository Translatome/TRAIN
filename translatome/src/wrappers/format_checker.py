"""
Utils for snakemake wrapper light-weight system
===============================================

"""

import os
import gzip
import math
import pandas as pd
from pandas.api.types import is_numeric_dtype


__author__ = "Julie Ripoll"
__created__ = "2020-05"
__license__ = "CeCILL"
__version__ = "0.1.0"
__vscript__ = "0.1.0"


def extensions(path, extension):
    """
    type: (path: str, extension: str | list | tuple) -> bool

    Checks if `path` matches the extension

    Parameters
    ----------
        path : str
            the path to the file

        extension : str, list, type
            possible accepted extension


    Returns
    -------
    `True` if `path` has `extension`, `False` otherwise

    Examples
    --------
        >>> extensions("file.pdf", ".pdf")
        True

        >>> extensions("file.jpg", (".jpg", ".jpeg"))
        True

        >>> extensions("file.pdf", [".png", ".jpeg",".jpg", ".gif"])
        False

    Notes
    -----
        Python > 3.5 :
            extensions(path: str, extension: str | list | tuple) -> bool
    
    """
    if isinstance(extension, (str, tuple)):
        if path.endswith(extension):
            return True
        return False
    elif isinstance(extension, list):
        if path.endswith(tuple(extension)):
            return True
        return False
    return False


def is_fasta_ext(input_fasta):
    """
    Boolean for extension of input file

    Require extensions()

    Returns
    -------
    `True` if `input` has `fasta extension`, `False` otherwise

    Examples
    --------
        >>> is_fasta_ext("test.fa")
        True

        >>> is_fasta_ext("test.fa.gz")
        True

        >>> is_fasta_ext("test.fq")
        False
    
    """
    return extensions(
        input_fasta,
        [".fa", ".fna", ".fas", ".fasta", ".fa.gz", ".fna.gz", ".fas.gz", ".fasta.gz"],
    )


def is_fastq_ext(input_fastq):
    """
    Boolean for extension of input file

    Require extensions()

    Returns
    -------
    `True` if `input` has `fastq extension`, `False` otherwise

    Examples
    --------
        >>> is_fastq_ext("test.fq")
        True

        >>> is_fastq_ext("test.fastq.gz")
        True

        >>> is_fastq_ext("test.fasta")
        False
    
    """
    return extensions(input_fastq, (".fq", ".fastq", ".fq.gz", ".fastq.gz"))


def is_gff3_ext(annotation_path):
    """
    Boolean for extension of input file

    Require extensions()

    Returns
    -------
    `True` if `input` has `gff extension`, `False` otherwise

    Examples
    --------
        >>> is_gff3_ext("test.gff")
        True

        >>> is_gff3_ext("test.gff3.gz")
        True

        >>> is_gff3_ext("test.gtf")
        False
    
    """
    return extensions(annotation_path, (".gff", ".gff.gz", ".gff3.gz", ".gff3"))


def is_gtf_ext(annotation_path):
    """
    Boolean for extension of input file

    Require extensions()

    Returns
    -------
    `True` if `input` has `gff extension`, `False` otherwise

    Examples
    --------
        >>> is_gtf_ext("test.gtf")
        True

        >>> is_gtf_ext("test.gtf.gz")
        True

        >>> is_gtf_ext("test.gff")
        False
    
    """
    return extensions(annotation_path, (".gtf", ".gtf.gz"))


def is_fasta(input_fasta: str):
    """
    Boolean for fasta file format
    
    Returns
    -------
    `True` if fasta input file is not empty, `False` otherwise.
    
    """
    if input_fasta.endswith("gz"):
        with gzip.open(input_fasta, "rb") as f:
            first_line = str(f.readline())
            return first_line.startswith("b'>")
    else:
        with open(input_fasta) as f:
            first_line = str(f.readline())
            return first_line.startswith(">")


def is_fastq(input_fastq) -> bool:
    """
    Boolean for fastq file format
    
    Returns
    -------
    `True` if fastq input file is not empty, `False` otherwise.

    Notes
    -----
    Requires gzip
    
    """
    if input_fastq.endswith("gz"):
        with gzip.open(input_fastq, "rb") as f:
            first_line = str(f.readline())
            return first_line.startswith("b'@")
    else:
        with open(input_fastq) as f:
            first_line = str(f.readline())
            return first_line.startswith("@")


def is_gff3_gtf(annotation_path):
    """
    Boolean for fasta file format
    
    Returns
    -------
    `True` if input GFF3/GTF annotation file is not empty and contains 9 elements, `False` otherwise.
    
    """
    if annotation_path.endswith("gz"):
        with gzip.open(annotation_path, "rb") as f:
            lines = f.readlines()
            for line in lines:
                if not str(line).startswith("b'#"):
                    check = str(line).split("\\t")
                    if len(check) != 9:
                        return False
            return True
    else:
        with open(annotation_path) as f:
            lines = f.readlines()
            for line in lines:
                if not str(line).startswith("#"):
                    check = str(line).split("\t")
                    if len(check) != 9:
                        return False
            return True


def raise_on_missing_file(path):
    """
    File checking
    
    Raises
    ------
    FileNotFoundError if file does not exists
    
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"""Error please fix this:
        Meta file not found: {path}"""
        )


def raise_bad_metacontrast_format(df: pd.DataFrame):
    """
    Metacontrast format checking
    
    Raises
    ------
    ValueError if df does not correspond to metacontrast format
    
    """
    expected_cols = ["fraction", "treatment1", "fraction_bis", "treatment2"]

    if not set(expected_cols).issubset(df.columns):
        # check if all items in set x are present in set y
        raise ValueError(
            f"""Error please fix this:
            Metacontrast file columns are incorrect, should be:
            {expected_cols}
            your file:
            {df.columns}"""
        )
    if len(df) < 2:
        raise ValueError(
            f"""Error please fix this:
            Metacontrast file should contain at least 2 contrasts:
            {df}"""
        )
    if not (df["fraction"] == df["fraction_bis"]).all():
        raise ValueError(
            """Error please fix this:
            Metacontrast fractions should be the same in 'fraction' and 'fraction_bis' columns."""
        )


def raise_bad_metadata_format(df: pd.DataFrame):
    """
    Metadata format checking
    
    Raises
    ------
    ValueError if df does not respect metadata format
    
    """
    expected_cols = ["sample_name", "fraction", "treatment", "replicate"]
    if not set(expected_cols).issubset(df.columns):
        raise ValueError(
            f"""Error please fix this:
            Metadata file columns are incorrect, should be:
            {expected_cols}
            found:
            {df.columns}"""
        )
    if is_numeric_dtype(df["replicate"]):
        raise ValueError(
            f"""Error please fix this:
            Metadata replicates should not be of the numeric type. Consider changing to R1, R2, etc:
            {df}"""
        )


def raise_bad_metacontrast_kinetic(df: pd.DataFrame):
    """
    Kinetic metacontrast file checking
    
    Raises
    ------
    ValueError if df does not respect metacontrast format
    
    """
    expected_cols = [
        "fraction",
        "treatment1",
        "time",
        "fraction_bis",
        "treatment2",
        "time_bis",
    ]

    if not set(expected_cols).issubset(df.columns):
        # check if all items in set x are present in set y
        raise ValueError(
            f"""Metacontrast file columns are incorrect, should be:
            {expected_cols}
            found:
            {df.columns}"""
        )

    if not (df["fraction"] == df["fraction_bis"]).all():
        raise ValueError(
            """Metacontrast fractions should be the same in 'fraction' and 'fraction_bis' columns."""
        )
    if is_numeric_dtype(df["time"]):
        raise ValueError(
            f"""Metacontrast time should not be of the numeric type. Consider changing to D1, D2, etc:
            {df}"""
        )
    if is_numeric_dtype(df["time_bis"]):
        raise ValueError(
            f"""Metacontrast time_bis should not be of the numeric type. Consider changing to D1, D2, etc:
            {df}"""
        )

    def test_fraction(file) -> pd.DataFrame:
        check_file = (
            file.apply(lambda d: "_".join(d.array[:1]), axis=1)
            .groupby(lambda x: math.floor(x / 2))
            .apply(lambda d: list(d.array))
        )
        check_file = check_file.to_list()
        return pd.DataFrame(check_file).nunique().eq(1)

    if not test_fraction(df).values.any():
        raise ValueError(
            """
              Error, your contrast should be paired by couple of lines.
              
              Note that cytoplasmic or total RNA-seq should be the first element of your metacontrast after the title !
              
              Also, 'control' should be in the 'treatment2' column
              
              Example:
              --------
              fraction     treatment1  time  fraction_bis  treatment2  time_bis
              Cytoplasmic  treated     D1    Cytoplasmic   control     D1
              Polysomal    treated     D1    Polysomal     control     D1
              Cytoplasmic  treated     D2    Cytoplasmic   control     D2
              Polysomal    treated     D2    Polysomal     control     D2
              """
        )


def raise_on_metadata_kinetic(df: pd.DataFrame):
    """
    Kinetic metadata file checking
    
    Raises
    ------
    ValueError on incorrect metadata file format
    
    """

    expected_cols = ["sample_name", "fraction", "treatment", "time", "replicate"]

    if not set(expected_cols).issubset(df.columns):
        raise ValueError(
            f"""Error please fix this:
            Metadata file columns are incorrect, should be:
            {expected_cols}
            found:
            {df.columns}"""
        )
    if is_numeric_dtype(df["time"]):
        raise ValueError(
            f"""Error please fix this:
            Metadata time should not be of the numeric type. Consider changing to D1, D2, etc:
            {df}"""
        )
    if is_numeric_dtype(df["replicate"]):
        raise ValueError(
            f"""Error please fix this:
            Metadata replicates should not be of the numeric type. Consider changing to R1, R2, etc:
            {df}"""
        )


if __name__ == "__main__":
    import doctest

    doctest.testmod(verbose=True)
