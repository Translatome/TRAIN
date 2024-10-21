"""
Samtools wrapper
================

Samtools is a suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:
- Samtools : Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format
- BCFtools : Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants
- HTSlib : A C library for reading/writing high-throughput sequencing data

Samtools and BCFtools both use HTSlib internally, but these source packages contain their own copies of htslib so they can be built independently.


Citation
--------
Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9


Documentation
-------------
http://www.htslib.org/doc/samtools.html


"""

from wrappers.wrapper_system import optional, join_str # for use in snakemake


__author__="Julie Ripoll"
__created__="2023-03-16"
__license__="CeCILL"
__version__="0.1.0"
__vwrapper__="0.1.0"
__vsamtools__="1.13"


def create_samtools_sort(input, output, params, log, **kwargs):
    r"""
    Sort alignments by leftmost coordinates, or by read name

    Parameters
    ----------
    input: object
        Input SAM/BAM/CRAM files

    params: _dict
        | **threads**: int -- Number of additional threads to use, default: THREADS
        | **compression**: int -- Set compression level, from 0 (uncompressed) to 9 (best)
        | **uncompress**: bool -- Output uncompressed data (equivalent to -l 0)
        | **maxmem**: int -- Set maximum memory per thread; suffix K/M/G recognized [768M]
        | **minimiser**: bool -- Use minimiser for clustering unaligned/unplaced reads
        | **kmersize**: int -- Kmer size to use for minimiser [20]
        | **readname**: bool -- Sort by read name (not compatible with samtools index command)
        | **tagvalue**: tag -- Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
        | **prefixtmp**: str -- Write temporary files to PREFIX.nnnn.bam
        | **nopg**: bool -- do not add a PG line
        | **inputfmt**: value -- Specify a single input file format option in the form of OPTION or OPTION=VALUE
        | **outputfmt**: value -- Specify output format (SAM, BAM, CRAM)
        | **reference**: object -- Reference sequence FASTA FILE [null]
        | **writeindex**: bool -- Automatically index the output files [off]
        | **verbosity**: int -- Set level of verbosity

    Returns
    -------
    output: _dict
        | **out**: object -- Write final output to FILE rather than standard output. Specify output format (SAM, BAM, CRAM)
        | **console**: object -- gets the shell output

    Raises
    ------
    log: object
        gets the shell error output


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> log = _StubParams("log")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_samtools_sort(input, output, params, log)
    'samtools sort [-@ <{{params.threads}}, default=None>] [-l <{{params.compression}}, default=None>] [-u] [-m <{{params.maxmem}}, default=None>] [-M] [-K <{{params.kmersize}}, default=None>] [-n] [-t <{{params.tagvalue}}, default=None>] [-T <{{params.prefixtmp}}, default=None>] [--no-PG] [--input-fmt-option <{{params.inputfmt}}, default=None>] [--output-fmt-option <{{params.outputfmt}}, default=None>] [--reference <{{params.reference}}, default=None>] [--write-index <{{params.writeindex}}, default=None>] [--verbosity <{{params.verbosity}}, default=None>] -o {output.out} {input} 1>{output.console} 2>{log}'

    Specific cases

    >>> params = Params(fromdict={"maxmem": "50", "uncompress": "0", "threads": "3"})
    >>> create_samtools_sort(input, output, params, log)
    'samtools sort -@ 3 -u -m 50 -o {output.out} {input} 1>{output.console} 2>{log}'

    >>> params = Params(fromdict={"threads": "3"})
    >>> create_samtools_sort(input, output, params, log)
    'samtools sort -@ 3 -o {output.out} {input} 1>{output.console} 2>{log}'
    
    """

    result = join_str(
        "samtools sort",
        optional(params, "threads", "-@ {}"),
        optional(params, "compression", "-l {}"),
        optional(params, "uncompress", "-u"),
        optional(params, "maxmem", "-m {}"),
        optional(params, "minimiser", "-M"),
        optional(params, "kmersize", "-K {}"),
        optional(params, "readname", "-n"),
        optional(params, "tagvalue", "-t {}"),
        optional(params, "prefixtmp", "-T {}"),
        optional(params, "nopg", "--no-PG"),
        optional(params, "inputfmt", "--input-fmt-option {}"),
        optional(params, "outputfmt", "--output-fmt-option {}"),
        optional(params, "reference", "--reference {}"),
        optional(params, "writeindex", "--write-index {}"),
        optional(params, "verbosity", "--verbosity {}"),
        f"-o {output.out}",
        f"{input}",
        f"1>{output.console}",
        f"2>{log}"
    )
    return result



def create_samtools_index(input, output, params, log, **kwargs):
    r"""
    Index coordinate-sorted BGZIP-compressed SAM, BAM or CRAM files for fast random access.

    Parameters
    ----------
    input: object
        | Sorted Input SAM/BAM/CRAM files

    params: _dict
        | **threads**: int -- Number of additional threads to use, default: THREADS
        | **bai**: bool -- Generate BAI-format index for BAM files [default]
        | **csi**: bool -- Generate CSI-format index for BAM files
        | **min**: int -- Set minimum interval size for CSI indices to 2^INT [14]

    Returns
    -------
    output: object
        indexed file

    Raises
    ------
    log: object
        gets the shell error output


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> log = _StubParams("log")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_samtools_index(input, output, params, log)
    'samtools index [-b] [-c] [-m <{{params.min}}, default=None>] [-@ <{{params.threads}}, default=None>] {input} {output} 2>{log}'

    Specific case

    >>> params = Params(fromdict={"bai": True, "threads": "3"})
    >>> create_samtools_index(input, output, params, log)
    'samtools index -b -@ 3 {input} {output} 2>{log}'
    
    """

    result = join_str(
        "samtools index",
        optional(params, "bai", "-b"),
        optional(params, "csi", "-c"),
        optional(params, "min", "-m {}"),
        optional(params, "threads", "-@ {}"),
        f"{input}",
        f"{output}",
        f"2>{log}"
    )
    return result



def create_samtools_flagstat(input, output, params, **kwargs):
    r"""
    Does a full pass through the input file to calculate and print statistics to stdout. 

    Parameters
    ----------
    input: _dict
        | **bam**: object -- Sorted input SAM/BAM/CRAM files
        | **check**: object -- Indexed input SAM/BAM/CRAM files

    params: _dict
        | **threads**: int -- Number of additional threads to use, default: THREADS
        | **inputfmt**: value -- Specify a single input file format option in the form of OPTION or OPTION=VALUE
        | **outformat**: value -- Specify output format (json, tsv)
        | **verbose**: int -- Set level of verbosity

    Returns
    -------
    output: object
        gets the shell output


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_samtools_flagstat(input, output, params)
    'samtools flagstat [-@ <{{params.threads}}, default=None>] -O {params.outformat} [--input-fmt-option <{{params.inputfmt}}, default=None>] [--verbosity <{{params.verbosity}}, default=None>] {input.bam} 1>{output}'

    Specific case

    >>> params = Params(fromdict={"threads": "3", "outformat": "tsv"})
    >>> create_samtools_flagstat(input, output, params)
    'samtools flagstat -@ 3 -O tsv {input.bam} 1>{output}'
    
    """
    result = join_str(
        "samtools flagstat",
        optional(params, "threads", "-@ {}"),
        f"-O {params.outformat}",
        optional(params, "inputfmt", "--input-fmt-option {}"),
        optional(params, "verbosity", "--verbosity {}"),
        f"{input.bam}",
        f"1>{output}"
    )
    return result



def create_samtools_stats(input, output, params, **kwargs):
    r"""
    samtools stats collects statistics from BAM files and outputs in a text format. 

    Parameters
    ----------
    input: _dict
        | **bamfile**: object -- Sorted input SAM/BAM/CRAM files
        | **check**: object -- output from flagstat or other previous rule

    params: _dict
        | **threads**: int -- Number of additional threads to use, default: THREADS
        | **coverage**: [int,int,int] -- Coverage distribution min,max,step [1,1000,1]
        | **removedups**: bool -- Exclude from statistics reads marked as duplicates
        | **customindex**: object -- Use a customized index file
        | **requiredflag**: str/int -- Required flag, 0 for unset. See also `samtools flags` [0]
        | **filteringflag**: str/int -- Filtering flag, 0 for unset. See also `samtools flags` [0]
        | **GCdepth**: float -- the size of GC-depth bins (decreasing bin size increases memory requirement) [2e4]
        | **insertsize**: int --  Maximum insert size [8000]
        | **id**: str -- Include only listed read group or sample name
        | **readlength**: int -- Include in the statistics only reads with the given read length [-1]
        | **mostinserts**: float -- Report only the main part of inserts [0.99]
        | **splitprefix**: str -- Path or string prefix for filepaths output by -S (default is input filename)
        | **trimquality**: int -- The BWA trimming parameter [0]
        | **refseq**: object -- Reference sequence (required for GC-depth and mismatches-per-cycle calculation).
        | **split**: tag -- Also write statistics to separate files split by tagged field.
        | **targetregions**: object -- Do stats in these regions only. Tab-delimited file chr,from,to, 1-based, inclusive.
        | **sparse**: bool -- Suppress outputting IS rows where there are no insertions.
        | **removeoverlaps**: bool -- Remove overlaps of paired-end reads from coverage and base count computations.
        | **covthreshold**: int -- Only bases with coverage above this value will be included in the target percentage computation [0]
        | **inputfmt**: value -- Specify a single input file format option in the form of OPTION or OPTION=VALUE
        | **reference**: object -- Reference sequence FASTA FILE [null]
        | **verbosity**: int -- Set level of verbosity

    Returns
    -------
    output: object
        TSV file with statistics


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_samtools_stats(input, output, params)
    'samtools stats  [--coverage <{{params.coverage}}, default=None>] [--remove-dups] [--customized-index-file <{{params.customindex}}, default=None>] [--required-flag <{{params.requiredflag}}, default=None>] [--filtering-flag <{{params.filteringflag}}, default=None>] [--GC-depth <{{params.GCdepth}}, default=None>] [--insert-size <{{params.insertsize}}, default=None>] [--id <{{params.id}}, default=None>] [--read-length <{{params.readlength}}, default=None>] [--most-inserts <{{params.mostinserts}}, default=None>] [--split-prefix <{{params.splitprefix}}, default=None>] [--trim-quality <{{params.trimquality}}, default=None>] [--ref-seq <{{params.refseq}}, default=None>] [--split <{{params.split}}, default=None>] [--target-regions <{{params.targetregions}}, default=None>] [--sparse] [--remove-overlaps] [--cov-threshold <{{params.covthreshold}}, default=None>] [--input-fmt-option <{{params.inputfmt}}, default=None>] [--reference <{{params.reference}}, default=None>] [--verbosity <{{params.verbosity}}, default=None>] [-@ <{{params.threads}}, default=None>] {input} > {output}'

    Specific case

    >>> params = Params(fromdict={"threads": "3", "coverage": "1,100000,1"})
    >>> create_samtools_stats(input, output, params)
    'samtools stats  --coverage 1,100000,1 -@ 3 {input} > {output}'
    
    """
    result = join_str(
        "samtools stats ",
        optional(params, "coverage", "--coverage {}"),
        optional(params, "removedups", "--remove-dups"),
        optional(params, "customindex", "--customized-index-file {}"),
        optional(params, "requiredflag", "--required-flag {}"),
        optional(params, "filteringflag", "--filtering-flag {}"),
        optional(params, "GCdepth", "--GC-depth {}"),
        optional(params, "insertsize", "--insert-size {}"),
        optional(params, "id", "--id {}"),
        optional(params, "readlength", "--read-length {}"),
        optional(params, "mostinserts", "--most-inserts {}"),
        optional(params, "splitprefix", "--split-prefix {}"),
        optional(params, "trimquality", "--trim-quality {}"),
        optional(params, "refseq", "--ref-seq {}"),
        optional(params, "split", "--split {}"),
        optional(params, "targetregions", "--target-regions {}"),
        optional(params, "sparse", "--sparse"),
        optional(params, "removeoverlaps", "--remove-overlaps"),
        optional(params, "covthreshold", "--cov-threshold {}"),
        optional(params, "inputfmt", "--input-fmt-option {}"),
        optional(params, "reference", "--reference {}"),
        optional(params, "verbosity", "--verbosity {}"),
        optional(params, "threads", "-@ {}"),
        f"{input}",
        f"> {output}"
    )
    return result



def create_samtools_view(input, output, params, **kwargs):
    r"""
    With no options or regions specified, prints all alignments in the specified input alignment file (in SAM, BAM, or CRAM format) to standard output in SAM format (with no header).

    Parameters
    ----------
    input: _dict
        | **bamfile**: object -- Input SAM/BAM/CRAM files
        | **check**: object -- output from previous rule

    params: _dict
        | **threads**: int -- Number of additional threads to use, default: THREADS
        | **bam**: bool -- Output BAM
        | **cram**: bool -- Output CRAM (requires -T)
        | **fast**: bool -- Use fast BAM compression (implies --bam)
        | **uncompressed**: bool -- Uncompressed BAM output (implies --bam)
        | **withheader**: bool -- Include header in SAM output
        | **headeronly**: bool -- Print SAM header only (no alignments)
        | **noheader**: bool -- Print SAM alignment records only [default]
        | **count**: bool -- Print only the count of matching records
        | **unoutput**: object -- Output reads not selected by filters to FILE
        | **faireference**: object -- File listing reference names and lengths
        | **useindex**: bool -- Use index and multi-region iterator for regions
        | **regionfile**: object -- Use index to include only reads overlapping FILE
        | **customizedindex**: bool -- Expect extra index file argument after <in.bam>
        | **targetfile**: object -- overlap (BED) regions in FILE
        | **readgroup**: str -- are in read group STR
        | **readgroupfile**: object -- are in a read group listed in FILE
        | **qnamefile**: object -- whose read name is listed in FILE
        | **tag**: str1[:str2] -- have a tag STR1 (with associated value STR2)
        | **tagfile**: [str:file] -- have a tag STR whose value is listed in FILE
        | **minMQ**: int -- have mapping quality >= INT
        | **library**: str -- ...are in library STR
        | **minqlen**: int -- ...cover >= INT query bases (as measured via CIGAR)
        | **expr**: str -- ...match the filter expression STR
        | **requireflags**: flag -- ...have all of the FLAGs present
        | **excludeflags**: flag -- ...have none of the FLAGs present
        | **FLAG**: flag -- EXCLUDE reads with all of the FLAGs present
        | **subsample**: float -- Keep only FLOAT fraction of templates/read pairs
        | **subsampleseed**: int Influence WHICH reads are kept in subsampling [0]
        | **addflags**: flag -- Add FLAGs to reads
        | **removeflags**: flag -- Remove FLAGs from reads
        | **removetag**: str -- Strip tag STR from reads (option may be repeated)
        | **removeB**: bool -- Collapse the backward CIGAR operation
        | **noPG**: bool -- Do not add a PG line
        | **inputfmt**: value -- Specify a single input file format option in the form of OPTION or OPTION=VALUE
        | **outputfmt**: value -- Specify output format (SAM, BAM, CRAM)
        | **outputfmtoption**: value -- Specify a single output file format option in the form of OPTION or OPTION=VALUE
        | **reference**: object -- Reference sequence FASTA FILE [null]
        | **writeindex**: bool -- Automatically index the output files [off]
        | **verbosity**: int Set level of verbosity

    Returns
    -------
    output: object
        Write output to FILE [standard output]


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_samtools_view(input, output, params)
    'samtools view [--bam] [--cram] [--fast] [--uncompressed] [--with-header] [--header-only] [--no-header] [--count] [--unoutput <{{params.unoutput}}, default=None>] [--fai-reference <{{params.faireference}}, default=None>] [--use-index] [--region[s]-file <{{params.regionfile}}, default=None>] [--customized-index <{{params.customizedindex}}, default=None>] [--target[s]-file <{{params.targetfile}}, default=None>] [--read-group <{{params.readgroup}}, default=None>] [--read-group-file <{{params.readgroupfile}}, default=None>] [--qname-file <{{params.qnamefile}}, default=None>] [--tag <{{params.tag}}, default=None>] [--tag-file <{{params.tagfile}}, default=None>] [--min-MQ <{{params.minMQ}}, default=None>] [--library <{{params.library}}, default=None>] [--min-qlen <{{params.minqlen}}, default=None>] [--expr <{{params.expr}}, default=None>] [--require-flags <{{params.requireflags}}, default=None>] [--excl[ude]-flags <{{params.excludeflags}}, default=None>] [-G <{{params.FLAG}}, default=None>] [--subsample <{{params.subsample}}, default=None>] [--subsample-seed <{{params.subsampleseed}}, default=None>] [--add-flags <{{params.addflags}}, default=None>] [--remove-flags <{{params.removeflags}}, default=None>] [--remove-tag <{{params.removetag}}, default=None>] [--remove-B ] [--no-PG] [--input-fmt-option <{{params.inputfmt}}, default=None>] [--output-fmt <{{params.outputfmt}}, default=None>] [--output-fmt-option <{{params.outputfmtoption}}, default=None>] [--reference <{{params.reference}}, default=None>] [-@ <{{params.threads}}, default=None>] [--write-index <{{params.writeindex}}, default=None>] [--verbosity <{{params.verbosity}}, default=None>] {input} -o {output}'

    Specific case

    >>> params = Params(fromdict={"bam": True})
    >>> create_samtools_view(input, output, params)
    'samtools view --bam {input} -o {output}'
    
    """
    result = join_str(
        "samtools view",
        optional(params, "bam", "--bam"),
        optional(params, "cram", "--cram"),
        optional(params, "fast", "--fast"),
        optional(params, "uncompressed", "--uncompressed"),
        optional(params, "withheader", "--with-header"),
        optional(params, "headeronly", "--header-only"),
        optional(params, "noheader", "--no-header"),
        optional(params, "count", "--count"),
        optional(params, "unoutput", "--unoutput {}"),
        optional(params, "faireference", "--fai-reference {}"),
        optional(params, "useindex", "--use-index"),
        optional(params, "regionfile", "--region[s]-file {}"),
        optional(params, "customizedindex", "--customized-index {}"),
        optional(params, "targetfile", "--target[s]-file {}"),
        optional(params, "readgroup", "--read-group {}"),
        optional(params, "readgroupfile", "--read-group-file {}"),
        optional(params, "qnamefile", "--qname-file {}"),
        optional(params, "tag", "--tag {}"),
        optional(params, "tagfile", "--tag-file {}"),
        optional(params, "minMQ", "--min-MQ {}"),
        optional(params, "library", "--library {}"),
        optional(params, "minqlen", "--min-qlen {}"),
        optional(params, "expr", "--expr {}"),
        optional(params, "requireflags", "--require-flags {}"),
        optional(params, "excludeflags", "--excl[ude]-flags {}"),
        optional(params, "FLAG", "-G {}"),
        optional(params, "subsample", "--subsample {}"),
        optional(params, "subsampleseed", "--subsample-seed {}"),
        optional(params, "addflags", "--add-flags {}"),
        optional(params, "removeflags", "--remove-flags {}"),
        optional(params, "removetag", "--remove-tag {}"),
        optional(params, "removeB", "--remove-B "),
        optional(params, "noPG", "--no-PG"),
        optional(params, "inputfmt", "--input-fmt-option {}"),
        optional(params, "outputfmt", "--output-fmt {}"),
        optional(params, "outputfmtoption", "--output-fmt-option {}"), 
        optional(params, "reference", "--reference {}"),
        optional(params, "threads", "-@ {}"),
        optional(params, "writeindex", "--write-index {}"),
        optional(params, "verbosity", "--verbosity {}"),
        f"{input}",
        f"-o {output}"
    )
    return result


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

