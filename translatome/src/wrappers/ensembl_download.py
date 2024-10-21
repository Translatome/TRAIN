"""
Ensembl db download
===================

Creates url for download of references files (genome, transcriptome, annotation)

"""

__author__ = "Julie Ripoll"
__created__ = "2020-06"
__license__ = "CeCILL"
__version__ = "0.1.0"
__vscript__ = "0.1.0"


# import
from wrappers.wrapper_system import join_str, val_mappy, mappy


def pathify(*arr):
    """
    Join characters/words with a "/"
    
    Examples
    --------
        >>> pathify('a','b','c')
        'a/b/c'
        >>> pathify('a')
        'a'
        
    """
    return "/".join(arr)


def compact_organism(organism: str):
    r"""
    Parameters
    ----------
        organism: str
            firstname_lastname

    Returns
    -------
        flastname

    Examples
    --------
        >>> compact_organism("gambusia_affinis")
        'gaffinis'
        >>> compact_organism("homo_sapiens")
        'hsapiens'
        
    """
    parts = organism.split("_", 1)
    species, *last_name = parts
    last_name_str = "".join(last_name)
    return f"{species[0]}{last_name_str}"


def injectable_ensembl_download(input, output, params, **kwargs):
    r"""
        function to create download link in snakemake rule
        
    Parameters
    ----------
        version: int
            release version
        extension: str
            extension of wanted file can be "fasta", "gtf", "gff3", "bed", "embl", "mysql"
        composante: str
            type of downloaded data can be "cdna", "dna", "cds", "dna_index", "ncrna", "pep"
        organism: str
            wanted organism
        output: str
            path for output files

    Examples
    --------
        output = outdir
        params = Params(fromdict={"organism": "mus_musculus", "version": "100", "extension": "gtf", "composante": ""})
        
        injectable_ensembl_download_simplify(output, params)
        'wget -P {output.outdir} ftp://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz'

    
    """

    path_version = f"release-{params.version}"

    base = "ftp://ftp.ensembl.org/pub"

    mappy(params, "extension", ("fasta", "gtf", "gff3", "bed", "embl", "mysql"))

    def wgets(formats, default):
        extended_base = pathify(base, path_version, formats)

        if formats == "fasta":
            composante = mappy(
                params,
                "composante",
                ("cdna", "dna", "cds", "dna_index", "ncrna", "pep"),
            )
            return f"wget -P {output.outdir} {pathify(extended_base, params.organism, composante, '*' )}"
        elif formats == "gff3" or formats == "gtf":
            return join_str(
                f"wget -P {output.outdir} {pathify(extended_base, params.organism, f'*{params.version}.{formats}.gz' )} &&",
                f"wget -P {output.outdir} {pathify(extended_base, params.organism, 'CHECKSUMS' )} &&",
                f"wget -P {output.outdir} {pathify(extended_base, params.organism, 'README' )}",
            )
        else:
            return f"wget -P {output.outdir} {pathify(extended_base, params.organism, '*')}"

    return mappy(params, "extension", wgets)


def injectable_ensembl_download_simplify(input, output, params, **kwargs):
    r"""
        function to create download link in snakemake rule
    
    Simplifed version for human and mus musculus to reduce download in case of dna
    
    See Also
    --------
    injectable_ensembl_download: generic function for download of Ensembl files
    
    """

    path_version = f"release-{params.version}"

    base = "ftp://ftp.ensembl.org/pub"

    mappy(params, "extension", ("fasta", "gtf", "gff3", "bed", "embl", "mysql"))

    def wgets(formats, default):
        extended_base = pathify(base, path_version, formats)

        if formats == "fasta":
            composante = mappy(
                params,
                "composante",
                ("cdna", "dna", "cds", "dna_index", "ncrna", "pep"),
            )
            if composante == "dna":
                return f"wget -P {output.outdir} {pathify(extended_base, params.organism, composante, '*dna.primary_assembly*' )}"
            else:
                return f"wget -P {output.outdir} {pathify(extended_base, params.organism, composante, '*' )}"
        elif formats == "gff3" or formats == "gtf":
            return join_str(
                f"wget -P {output.outdir} {pathify(extended_base, params.organism, f'*{params.version}.{formats}.gz' )} &&",
                f"wget -P {output.outdir} {pathify(extended_base, params.organism, 'CHECKSUMS' )} &&",
                f"wget -P {output.outdir} {pathify(extended_base, params.organism, 'README' )}",
            )
        else:
            return f"wget -P {output.outdir} {pathify(extended_base, params.organism, '*')}"

    return mappy(params, "extension", wgets)


if __name__ == "__main__":
    print(__doc__)
    import doctest

    doctest.testmod(verbose=True)
