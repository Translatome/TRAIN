"""
FastqScreen wrapper
===================

FastQ Screen is a simple application which allows you to search a large sequence dataset against a panel of different genomes to determine from where the sequences in your data originate.


Citation
--------
Wingett SW and Andrews S. FastQ Screen: A tool for multi-genome mapping and quality control [version 2; referees: 4 approved]. F1000Research 2018, 7:1338 (https://doi.org/10.12688/f1000research.15931.2)


Documentation
-------------
https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html


"""


from wrappers.wrapper_system import optional, join_str, mappy # for snakemake execution


__author__ = "Julie Ripoll"
__created__ = "2021-01-05"
__license__ = "CeCILL"
__version__ = "0.1.0"
__vwrapper__ = "0.1.0"
__vfqscreen__ = "0.14.0"


def create_FastqScreenDownload(output, params, **kwargs):
    r"""
    Download references genome from FastqScreen database (fasta format)

    Parameters
    ----------
    input: object
        input file to be screened

    params: _dict
        | **threads**: int -- the number of CPU cores available for parallelization
        | **get_genomes**: bool -- require to download genomes from fastqscreen database

    Returns
    -------
    output: _dict
        | **outQC**: object -- repository for downloaded genomes
        | **done**: object -- gets the shell error output

    Example
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> output = _StubParams("output")

    >>> params = Params(fromdict={"get_genomes": True, "threads": 3})
    >>> create_FastqScreenDownload(output, params)
    'fastq_screen --get_genomes  --threads 3 --outdir {output.outQC} 2>{output.done}'
    """
    
    return join_str(
        "fastq_screen",
        optional(params, "get_genomes", "--get_genomes "),
        optional(params, "threads", "--threads {}"),
        f"--outdir {output.outQC}",
        f"2>{output.done}"
    )


# Function FastqScreen
def create_FastqScreen(input, output, params, **kwargs):
    r"""
    Check contamination

    Parameters
    ----------
    input: _dict
        | **FReads**: object -- path to your input file to be screened
        | **genomes**: object -- path to reference genomes previously downloaded with get_genomes

    params: _dict
        Essential options:
            | **threads**: int -- the number of CPU cores available for parallelization
            | **real_format**: str -- extension of your files, can be fastq, fq or gzip version
            | **conf**: object -- configuration file downloaded with genomes
            | **bowtie/bowtie2/bwa**: str -- Specify extra parameters to be passed to Bowtie/Bowtie2/BWA. 
                These parameters should be quoted to clearly delimit Bowtie/Bowtie2/BWA parameters from fastq_screen parameters.
                You should not try to use this option to override the normal search or 
                reporting options for bowtie which are set automatically 
                but it might be useful to allow reads to be trimmed before alignment etc.

        Other options:
            | **aligner**: func -- Specify the aligner to use for the mapping.
                Valid arguments are 'bowtie', bowtie2' (default) or 'bwa'.
                Bowtie maps with parameters -k 2, Bowtie 2 with parameters -k 2 --very-fast-local and BWA with mem -a.
                Local aligners such as BWA or Bowtie2 will be better at detecting the origin of chimeric reads.
            | **bisulfite**: bool -- Process bisulfite libraries.
                The path to the bisulfite aligner (Bismark) may be specified in the configuration file.
                Bismark runs in non-directional mode.
                Either conventional or bisulfite libraries may be processed, but not both simultaneously.
                The --bisulfite option cannot be used in conjunction with --bwa.
            | **filter**: str -- Produce a FASTQ file containing reads mapping to specified genomes.
                Pass the argument a string of characters (0, 1, 2, 3, -), in which each character corresponds to a reference genome (in the order the reference genome occurs in the configuration file).
                Below gives an explanation of each character.
                    0: Read does not map
                    1: Read maps uniquely
                    2: Read multi-maps
                    3: Read maps (one or more times)
                    4: Passes filter 0 or filter 1
                    5: Passes filter 0 or filter 2
                    -: Do not apply filter to this genome
                Consider mapping to three genomes (A, B and C), the string '003' produces a file in which reads do not map to genomes A or B, but map (once or more) to genome C.
                The string '--1' would generate a file in which reads uniquely map to genome C.
                Whether reads map to genome A or B would be ignored.
                A read needs to pass all the genome filters to be considered valid (unless --pass specified).
                When --filter is used in conjuction with --tag, FASTQ files shall be mapped, tagged and then filtered.
                If the --tag option is not selected however, the input FASTQ file should have been previously tagged.
            | **force**: bool -- Do not terminate if output files already exist, instead overwrite the files.
            | **illumina13**: bool -- Assume that the quality values are in encoded in Illumina v1.3 format.
                Defaults to Sanger format if this flag is not specified.
            | **nohits**: bool -- Writes to a file the sequences that did not map to any of the specified genomes.
                This option is equivalent to specifying --tag --filter 0000 (number of zeros corresponds to the number of genomes screened).
                By default the whole input file will be mapped, unless overridden by --subset.
            | **passIffilter**: int -- Used in conjunction with --filter.
                By default all genome filters must be passed for a read to pass the --filter option.
                However, a minimum number of genome filters may be specified that a read needs pass to be considered to pass the --filter option.(--pass 1 effecitively acts as an OR boolean operator for the genome filters.)
            | **subset**: int -- Don't use the whole sequence file, but create a temporary dataset of this specified number of reads.
                The dataset created will be of approximately (within a factor of 2) of this size.
                If the real dataset is smaller than twice the specified size then the whole dataset will be used.
                Subsets will be taken evenly from throughout the whole original dataset.
                By Default FastQ Screen runs with this parameter set to 100000.
                To process an entire dataset however, adjust --subset to 0.
            | **tag**: bool -- Label each FASTQ read header with a tag listing to which genomes the read did, or did not align.
                The first read in the output FASTQ file will list the full genome names along with a score denoting whether the read did not align (0), aligned uniquely to the specified genome (1), or aligned more than once (2).
                In subsequent reads the genome names are omitted and only the score is printed, in the same order as the first line.
                This option results in the he whole file being processed unless overridden explicitly by the user with the --subset parameter
            | **top**: [int / int, int] -- Don't use the whole sequence file, create a temporary dataset of the specified number of reads taken from the top of the original file.
                It is also possible to specify the number of lines to skip before beginning the selection e.g. --top 100000,5000000 skips the first five million reads and selects the subsequent one hundred thousand reads.
                While this option is usually faster than comparable --subset operations, it does not prevent biases arising from non-uniform distribution of reads in the original FastQ file.
                This option should only be used when minimising processing time is of highest priority.
            | **quiet**: bool -- Suppress all progress reports on stderr and only report errors.

    Returns
    -------
    output: _dict
    
        | **outQC**: object -- repository for downloaded genomes
        | **done**: object -- gets the shell error output

    Examples
    --------

    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_FastqScreen(input, output, params)
    'mkdir -p {output.outQC} && fastq_screen [--aligner <{{params.aligner}}, default=None>] [--bisulfite] [--bowtie <{{params.bowtie}}, default=None>] [--bowtie2 <{{params.bowtie2}}, default=None>] [--bwa <{{params.bwa}}, default=None>] [--conf <{{params.conf}}, default=None>] [--filter <{{params.filter}}, default=None>] [--force] [--illumina1_3] [--nohits] [--pass <{{params.passIffilter}}, default=None>] [--subset <{{params.subset}}, default=None>] [--tag] [--top <{{params.top}}, default=None>] [--quiet] [--threads <{{params.threads}}, default=None>] --outdir {output.outQC} {input.FReads}/*.{params.realformat} 2>{output.done}'

    Specific case

    >>> params = Params(fromdict={"bwa": "bowtie_text", "realformat":"fastq", "threads": 3})
    >>> create_FastqScreen(input, output, params)
    'fastq_screen --bwa bowtie_text --threads 3 --outdir {output.outQC} {input.FReads}/*.fastq 2>{output.done}'
    
    """

    create_folder = mappy(params, 'create_folder',
                          (True, False), default=False)

    return join_str(
        (f"mkdir -p {output.outQC} &&" if create_folder else ""),
        "fastq_screen",
        optional(params, "aligner", "--aligner {}"),
        optional(params, "bisulfite", "--bisulfite"),
        optional(params, "bowtie", "--bowtie {}"),
        optional(params, "bowtie2", "--bowtie2 {}"),
        optional(params, "bwa", "--bwa {}"),
        optional(params, "conf", "--conf {}"),
        optional(params, "filter", "--filter {}"),
        optional(params, "force", "--force"),
        optional(params, "illumina13", "--illumina1_3"),
        optional(params, "nohits", "--nohits"),
        optional(params, "passIffilter", "--pass {}"),
        optional(params, "subset", "--subset {}"),
        optional(params, "tag", "--tag"),
        optional(params, "top", "--top {}"),
        optional(params, "quiet", "--quiet"),
        optional(params, "threads", "--threads {}"),
        f"--outdir {output.outQC}",
        f"{input.FReads}/*.{params.realformat}",
        f"2>{output.done}"
    )


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
