"""
BBTools wrapper
===============

BB stands for Bestus Bioinformaticus. 
BBTools is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA sequence data. 


Citation
--------
Brian Bushnell and Jonathan Rood
BBMerge manuscript: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0185056
BBtools may be cited using the primary website: BBMap – Bushnell B. – sourceforge.net/projects/bbmap/


Documentation
-------------
https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/


"""

from wrappers.wrapper_system import join_str, mappy # for use in snakemake

__author__="Julie Ripoll"
__created__="2021-01-05"
__license__="CeCILL"
__version__="0.1.0"
__vwrapper__="0.1.0"
__vbbmap__="38.91"


def create_bbmerge(input, output, params, **kwargs):
    r"""
    BBMerge is designed to merge two overlapping paired reads into a single read. For example, a 2x150bp read pair with an insert size of 270bp would result in a single 270bp read. 

    Parameters
    ----------
    input: _dict
        | **inputR1**: object -- fastq file for R1
        | **inputR2**: object -- fastq file for R2

    params:
        **reads**: str -- '2m' for overlappping reads

    Returns
    -------
    output: _dict
    
        | **merged**: object -- File for merged reads.
        | **hist**: object -- Insert length histogram output file.
        | **done**: object -- gets the shell in a text file

    Notes
    -----
        others options were not provided here

    Examples
    --------
    >>> from wrappers.wrapper_system import _StubParams
    >>> from snakemake.io import Params

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_bbmerge(_StubParams("input"), _StubParams("output"), params)
    'bbmerge.sh in1={input.inputR1} in2={input.inputR2} out={output.merged} ihist={output.hist} reads={params.reads} 2>{output.done}'

    Specific cases

    >>> params = Params(fromdict={"paired": "YES", "reads": "read_params"})
    >>> create_bbmerge(_StubParams("input"), _StubParams("output"), params)
    'bbmerge.sh in1={input.inputR1} in2={input.inputR2} out={output.merged} ihist={output.hist} reads={params.reads} 2>{output.done}'

    """

    result = join_str(
        "bbmerge.sh",
        "in1={input.inputR1} in2={input.inputR2}",
        "out={output.merged}",
        "ihist={output.hist}",
        "reads={params.reads}",
        "2>{output.done}"
    )
    return result


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

