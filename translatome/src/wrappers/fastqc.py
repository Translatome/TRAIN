"""
FASTQC wrapper
==============

A quality control tool for high throughput sequence data


Citation
--------
FastQC v{version number} (Babraham Institute, Cambridge, UK)


Documentation
-------------
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


"""

from wrappers.wrapper_system import optional, join_str, mappy, val_mappy # for snakemake execution


__author__="Julie Ripoll"
__created__="2021-01-05"
__license__="CeCILL"
__version__="0.1.0"
__vwrapper__="0.1.0"
__vfastqc__="0.12.1"


def create_fastqc(input, output, params, **kwargs):
    r"""
    A quality control tool for high throughput sequence data

    Parameters
    ----------
    input: object
        **FReads**: directory where fastq files are present, can be passed as params is some case

    params: _dict
        Essential options:
            | **threads**: int -- default: THREADS
            | **realformat**: str -- extension of files specified in the config.yml, ex: fastq.gz
            | **create_folder**: bool -- create a directory for output (not mandatory)
        
        Other options:
            | **casava**: bool -- if casava output
            | **nano**: bool -- if nanopore output
            | **nofilter**: bool -- If running with --casava then don't remove read flagged by 
                casava as poor quality when performing the QC analysis.
            | **extract**: bool --  If set then the zipped output file will be uncompressed in
                the same directory after it has been created. If --delete is 
                also specified then the zip file will be removed after the 
                contents are unzipped.
            | **noextract**: bool -- don't uncompress output
            | **nogroup**: bool --  disable grouping of bases for reads >50bp
            | **min_lenght**: int -- Sets an artificial lower limit on the length of the sequence to be shown in the report
            | **dup_length**: int -- Sets a length to which the sequences will be truncated when defining them to be duplicates
            | **format**: str -- Bypasses the normal sequence file format detection and
                forces the program to use the specified format.
                Valid formats are bam,sam,bam_mapped,sam_mapped and fastq
            | **memory**: value -- Sets the base amount of memory, in Megabytes, used to process each file.  
                Defaults to 512MB
            | **svg**: bool -- Save the graphs in the report in SVG format
            | **contaminants**: object -- Specifies a non-default file which contains the list of
                contaminants to screen overrepresented sequences against.
                The file must contain sets of named contaminants in the
                form name[tab]sequence. Lines prefixed with a hash will be ignored.
            | **adapters**: object -- Specifies a non-default file which contains the list of
                adapters which will be explicity searched against the library.
                The file must contain sets of named adapters in the
                form name[tab]sequence. Lines prefixed with a hash will be ignored.
            | **limits**: object -- list file of set of criteria which will be used to determine the warn/error limits for the various modules
            | **kmers**: int -- length of kmers to look for in the kmer content module, must be between 2 and 10, default 7
            | **quiet**: bool -- Supress all progress messages
            | **dir**: object -- directory to be used for temporary files

    Returns
    -------
    output: _dict
    
        | **outQC**: object -- directory for reports
        | **done**: object -- gets the shell in a text file

    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_fastqc(input, output, params)
    'mkdir -p {output.outQC} && fastqc --outdir {output.outQC} [--casava] [--nano] [--nofilter] [--extract] [--noextract] [--nogroup] [--min_length <{{params.min_lenght}}, default=None>] [--dup_length <{{params.dup_length}}, default=None>] [--format <{{params.format}}, default=None>] [--memory <{{params.memory}}, default=None>] [--svg] [--threads <{{params.threads}}, default=None>] [--contaminants <{{params.contaminants}}, default=None>] [--adapters <{{params.adapters}}, default=None>] [--limits <{{params.limits}}, default=None>] [--kmers <{{params.kmers}}, default=None> ] [--quiet] [--dir <{{params.dir}}, default=None> ] [<{{params.FReads}}, default=None>/*.{params.realformat} ] [<{{input.FReads}}, default=None>/*.{params.realformat} ] 2>{output.done}'

    Specific cases

    >>> params = Params(fromdict={"FReads": "files", "realformat": "fastq.gz", "noextract": True, "create_folder": False})
    >>> create_fastqc(input, output, params)
    'fastqc --outdir {output.outQC} --noextract files/*.{params.realformat}  [<{{input.FReads}}, default=None>/*.{params.realformat} ] 2>{output.done}'

    >>> input = Params(fromdict={"FReads": "input"})
    >>> params = Params(fromdict={"realformat": "fastq.gz", "noextract": True, "create_folder": True, "threads": "3"})
    >>> create_fastqc(input, output, params)
    'mkdir -p {output.outQC} && fastqc --outdir {output.outQC} --noextract --threads 3 input/*.{params.realformat}  2>{output.done}'

    >>> input = Params(fromdict={"FReads": "input"})
    >>> params = Params(fromdict={"realformat": "fastq.gz", "noextract": True, "create_folder": False, "threads": "3"})
    >>> create_fastqc(input, output, params)
    'fastqc --outdir {output.outQC} --noextract --threads 3 input/*.{params.realformat}  2>{output.done}'
    
    """

    create_folder = mappy(params, 'create_folder',
                          (True, False), default=False)

    return join_str(
        (f"mkdir -p {output.outQC} &&" if create_folder else ""),
        "fastqc",
        f"--outdir {output.outQC}",
        optional(params, 'casava', "--casava"),
        optional(params, 'nano', "--nano"),
        optional(params, 'nofilter', "--nofilter"),
        optional(params, "extract", "--extract"),
        optional(params, 'noextract', "--noextract"),
        optional(params, 'nogroup', "--nogroup"),
        optional(params, "min_lenght", "--min_length {}"),
        optional(params, "dup_length", "--dup_length {}"),
        optional(params, 'format', "--format {}"),
        optional(params, "memory", "--memory {}"),
        optional(params, "svg", "--svg"),
        optional(params, "threads", "--threads {}"),
        optional(params, 'contaminants', "--contaminants {}"),
        optional(params, 'adapters', "--adapters {}"),
        optional(params, 'limits', "--limits {}"),
        optional(params, 'kmers', "--kmers {} "),
        optional(params, 'quiet', "--quiet"),
        optional(params, 'dir', "--dir {} "),
        optional(params, 'FReads', "{}/*.{{params.realformat}} "),
        optional(input, 'FReads', "{}/*.{{params.realformat}} "),
        f"2>{output.done}"
    )


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

