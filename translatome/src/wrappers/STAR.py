"""
STAR wrapper
============

Spliced Transcripts Alignment to a Reference (STAR) is a fast RNA-seq read mapper, with support for splice-junction and fusion read detection.


Citation
--------
Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635


Documentation
-------------
https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf


"""


from wrappers.wrapper_system import optional, join_str, mappy # for use in snakemake
from wrappers.format_checker import is_gff3_ext


__author__="Julie Ripoll"
__created__="2020-05"
__license__="CeCILL"
__version__="0.1.0"
__vwrapper__="0.1.0"
__vstar__="2.7.9a"


def create_STARIndex(input, output, params, **kwargs):
    """
    Create index of reference genome/transcriptome for STAR

    Parameters
    ----------
    input: _dict
        | **faref**: object -- fasta reference file
        | **annot**: object -- GFF or GTF annotation file
        | **check_file**: object -- for rule consistency (default: None)

    params:
        | **create_folder**: create output dir if required (default: False)
        | **threads**: int -- number of threads to run STAR, default: THREADS from snakemake config file
        | **readlength**: int -- read length -1
        | **genomeSAindexNbases**: int -- length (bases) of the SA pre-indexing string. Typically between 10 and 15. 
            Longer strings will use much more memory, but allow faster searches. 
            For small genomes, the parameter --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1).
        |**sjdbFileChrStartEnd**: path -- path to the files with genomic coordinates (chr <tab> start <tab> end <tab> strand) 
            for the splice junction introns. Multiple files can be supplied wand will be concatenated.

    Returns
    -------
    output: path
        **refg**: directory for index


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> input = Params(fromdict={"annot": "test.gff", "faref": "test.fa"})
    >>> params = _StubParams("params")
    >>> create_STARIndex(input, output, params)
    'mkdir -p {output.refg} && STAR --runMode genomeGenerate --genomeDir {output.refg} --genomeFastaFiles test.fa --sjdbGTFfile test.gff [--runThreadN <{{params.threads}}, default=None>] --sjdbGTFtagExonParentTranscript Parent [--sjdbOverhang <{{params.readlength}}, default=None>] [--genomeSAindexNbases <{{params.genomeSAindexNbases}}, default=None>]'

    Specific cases

    >>> input = Params(fromdict={"annot": "test.gff", "faref": "test.fa"})
    >>> params = Params(fromdict={"create_folder": True})
    >>> create_STARIndex(input, output, params)
    'mkdir -p {output.refg} && STAR --runMode genomeGenerate --genomeDir {output.refg} --genomeFastaFiles test.fa --sjdbGTFfile test.gff --sjdbGTFtagExonParentTranscript Parent'

    >>> input = Params(fromdict={"annot": "test.gtf", "faref": "test.fa"})
    >>> params = Params(fromdict={"create_folder": True})
    >>> create_STARIndex(input, output, params)
    'mkdir -p {output.refg} && STAR --runMode genomeGenerate --genomeDir {output.refg} --genomeFastaFiles test.fa --sjdbGTFfile test.gtf'

    >>> input = Params(fromdict={"annot": "test.gtf", "faref": "test.fa"})
    >>> params = Params(fromdict={"create_folder": False})
    >>> create_STARIndex(input, output, params)
    'STAR --runMode genomeGenerate --genomeDir {output.refg} --genomeFastaFiles test.fa --sjdbGTFfile test.gtf'
    

    """

    create_folder = mappy(params, 'create_folder',
                          (True, False), default=False)
 
    return join_str(
        (f"mkdir -p {output.refg} &&" if create_folder else ""),
        "STAR",
        "--runMode genomeGenerate",
        f"--genomeDir {output.refg}",
        f"--genomeFastaFiles {input.faref}",
        f"--sjdbGTFfile {input.annot}",
        optional(params, "threads", "--runThreadN {}"),
        ("--sjdbGTFtagExonParentTranscript Parent" if is_gff3_ext(input.annot) else ""),
        optional(params, 'readlength', "--sjdbOverhang {}"),
        optional(params, 'genomeSAindexNbases', "--genomeSAindexNbases {}"),
        optional(input, 'sjdbFileChrStartEnd', "--sjdbFileChrStartEnd {}"))


def create_STAR(input, params, log, **kwargs):
    """
    Align on reference genome of selected isoforms

    Parameters
    ----------
    input: _dict
        | **fareads**: object -- fasta or fastq files to align
        | **refg**: path to index (e.g.: directory("../star-index/"))

    params: _dict
        Essential options:
            | **paired**: ["YES"/ "NO"] -- can be "YES" for paired-end data and "NO" for single-end
            | **annot**: object -- path to GFF or GTF annotation file
            | **outfilename**: str -- output filename
            | **threads**: int -- default: THREADS
            | **readFilesCommand**: [zcat | bzcat | -]  -- to uncompress, default in snakefile: zcat
                default STAR -

        Manage by wrapper:
            | **sjdbGTFtagExonParentTranscript**: ["Parent"] -- for GFF3 formatted annotations
                default: transcript id, works for GTF files

        Other options:
            | **outFilterMismatchNmax**: int -- alignment will be output only if it has no more mismatches than this value.
                default 10
            | **seedSearchStartLmax**: real -- seedSearchStartLmax normalized to read length (sum of mates' lengths for paired-end reads)
                default: 1.0
            | **seedPerReadNmax**: int>0 -- max number of seeds per read
                default: 1000
            | **seedPerWindowNmax**: int>0 -- max number of seeds per window
                default: 50
            | **alignTranscriptsPerReadNmax**: int>0 -- max number of different alignments per read to consider
                default: 10000
            | **alignTranscriptsPerWindowNmax**: int>0 -- max number of transcripts per window
                default: 100

        STAR has many other options, they have not been provided in this wrapper but may be added in future version.

    Returns
    -------
    output: object
        gets the shell in a text file

    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> log = _StubParams("log")
    >>> create_STAR(input, params, log)
    "STAR <mappy('{params.paired}', {{'YES': '--readFilesIn {{input.fareads1}} {{input.fareads2}}', 'NO': '--readFilesIn {{input.fareads}}'}})> --genomeDir {input.refg} --sjdbGTFfile {params.annot} --outFileNamePrefix {params.outfilename} [--readFilesCommand <{{params.readFilesCommand}}, default=None>] [--runThreadN <{{params.threads}}, default=None>] [-alignments <{{params.output_bam}}, default=None>] [--outFilterMismatchNmax <{{params.outFilterMismatchNmax}}, default=None>] [--seedSearchStartLmax <{{params.seedSearchStartLmax}}, default=None>] [--seedPerReadNmax <{{params.seedPerReadNmax}}, default=None>] [--seedPerWindowNmax <{{params.seedPerWindowNmax}}, default=None>] [--alignTranscriptsPerReadNmax <{{params.alignTranscriptsPerReadNmax}}, default=None>] [--alignTranscriptsPerWindowNmax <{{params.alignTranscriptsPerWindowNmax}}, default=None>] 1>{log}"

    Specific cases

    >>> input = Params(fromdict={"refg": "path/star_index/"})
    >>> params = Params(fromdict={"annot": "test.gff", "threads": "10", "paired": "YES", "outfilename": "path/{sample}"})
    >>> log = _StubParams("log")
    >>> create_STAR(input, params, log)
    'STAR --readFilesIn {input.fareads1} {input.fareads2} --genomeDir path/star_index/ --sjdbGTFfile test.gff --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix path/{sample} --runThreadN 10 1>{log}'

    >>> input = Params(fromdict={"refg": "path/star_index/"})
    >>> params = Params(fromdict={"annot": "test.gff", "threads": "10", "paired": "YES", "readFilesCommand": "zcat", "outfilename": "path/{sample}"})
    >>> log = _StubParams("log")
    >>> create_STAR(input, params, log)
    'STAR --readFilesIn {input.fareads1} {input.fareads2} --genomeDir path/star_index/ --sjdbGTFfile test.gff --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix path/{sample} --readFilesCommand zcat --runThreadN 10 1>{log}'

    >>> input = Params(fromdict={"refg": "path/star_index/"})
    >>> params = Params(fromdict={"annot": "test.gff", "threads": "10", "paired": "NO", "outfilename": "path/{sample}"})
    >>> log = _StubParams("log")
    >>> create_STAR(input, params, log)
    'STAR --readFilesIn {input.fareads} --genomeDir path/star_index/ --sjdbGTFfile test.gff --sjdbGTFtagExonParentTranscript Parent --outFileNamePrefix path/{sample} --runThreadN 10 1>{log}'

    >>> input = Params(fromdict={"refg": "path/star_index/"})
    >>> params = Params(fromdict={"annot": "test.gtf", "threads": "10", "paired": "NO", "outfilename": "path/{sample}"})
    >>> log = _StubParams("log")
    >>> create_STAR(input, params, log)
    'STAR --readFilesIn {input.fareads} --genomeDir path/star_index/ --sjdbGTFfile test.gtf --outFileNamePrefix path/{sample} --runThreadN 10 1>{log}'
    
    """

    fareadss = mappy(params, "paired", {
        "YES": "--readFilesIn {input.fareads1} {input.fareads2}",
        "NO": "--readFilesIn {input.fareads}"
    })

    result = join_str(
        "STAR",
        fareadss,
        f"--genomeDir {input.refg}",
        f"--sjdbGTFfile {params.annot}",
        ("--sjdbGTFtagExonParentTranscript Parent" if is_gff3_ext(params.annot) else ""),
        f"--outFileNamePrefix {params.outfilename}",
        optional(params, "readFilesCommand", "--readFilesCommand {}"),
        optional(params, "threads", "--runThreadN {}"),
        # optional(params, 'outSAMunmapped', "--outSAMunmapped {}"), # deprecated since STAR v2.7
        optional(params, 'output_bam', "-alignments {}"),
        optional(params, 'outFilterMismatchNmax', "--outFilterMismatchNmax {}"),
        optional(params, "seedSearchStartLmax", "--seedSearchStartLmax {}"),
        optional(params, "seedPerReadNmax", "--seedPerReadNmax {}"),
        optional(params, "seedPerWindowNmax", "--seedPerWindowNmax {}"),
        optional(params, "alignTranscriptsPerReadNmax", "--alignTranscriptsPerReadNmax {}"),
        optional(params, "alignTranscriptsPerWindowNmax", "--alignTranscriptsPerWindowNmax {}"),
        f"1>{log}"
    )
    return result


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose = True)
