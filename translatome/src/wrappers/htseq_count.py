"""
HTSeqCount wrapper
==================

This script takes one or more alignment files in SAM/BAM format and a feature file in GFF format and calculates for each feature the number of reads mapping to it. 


Citation
--------
Simon Anders, Paul Theodor Pyl, Wolfgang Huber HTSeq — A Python framework to work with high-throughput sequencing data Bioinformatics (2014)


Documentation
-------------
https://htseq.readthedocs.io/en/master/htseqcount.html


"""


from wrappers.wrapper_system import optional, join_str # for snakemake execution
from wrappers.format_checker import is_gff3_ext


__author__="Julie Ripoll & Celine Mandier"
__created__="2023-03-17"
__license__="CeCILL"
__version__="0.1.0"
__vwrapper__ ="0.1.0"
__vhtseqcounts__="0.13.5"


def create_htseqcount(input, params, output, log, **kwargs):
    r"""
    Given a file with aligned sequencing reads and a list of genomic features, htseq-count counts how many reads map to each feature.

    Parameters
    ----------
    input: object
        SAM/BAM alignment files

    params:
        Essential options:
            | **annotation**: object -- GTF/GFF file for annotation
            | **threads**: int -- Number of parallel CPU processes to use (default: 1), default in snakemake configuration file: THREADS
            | **stranded**: str -- Whether the data is from a strand-specific assay. 
                Specify 'yes', 'no', or 'reverse' (default: yes). 'reverse' means 'yes' with reversed strand interpretation
            | **featuretype**: str -- Feature type (3rd column in GTF file) to be used, all features of other type are ignored 
                (default, suitable for Ensembl GTF files: exon)
            | **mode**: str -- Mode to handle reads overlapping more than one feature 
                (choices: union, intersection-strict, intersection-nonempty; default: union)
            | **supplementaryalignments**: str -- Whether to score supplementary alignments (0x800 flag) can be 'score' or 'ignore' (default: 'score')

        Manage by wrapper:
            | **attributes**: str -- GTF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id). 
                All feature of the right type (see -t option) within the same GTF attribute will be added together. 
                The typical way of using this option is to count all exonic reads from each gene and add the exons 
                but other uses are possible as well.
            | **additional-attr**: str -- Additional feature attributes (default: none, suitable for Ensembl GTF files: gene_name). 
                Use multiple times for more than one additional attribute. 
                These attributes are only used as annotations in the output, while the determination of how the counts 
                are added together is done based on option -i.
        
        Other options:
            | **order**: str -- Sorting order of <alignment_file>, 'pos' or 'name' (default: name). 
                Paired-end sequencing data must be sorted either by position or by read name, and the sorting order must be specified. 
                Ignored for single-end data.
            | **maxreadsinbuffer**: int -- When <alignment_file> is paired end sorted by position, allow only so many reads to stay in memory 
                until the mates are found (raising this number will use more memory). 
                Has no effect for single end or paired end sorted by name
            | **minaqual**: int -- Skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10). 
                MAPQ is the 5th column of a SAM/BAM file and its usage depends on the software used to map the reads.
            | **addchrinfo**: bool -- Store information about the chromosome of each feature as an additional attribute 
                (e.g. column in the TSV output file).
            | **nonunique**: str -- Whether and how to score reads that are not uniquely aligned or ambiguously assigned to features 
                (choices: none, all, fraction, random; default: none)
            | **secondaryalignments**: str -- Whether to score secondary alignments (0x100 flag) can be 'score' or 'ignore' (default: 'score')
            | **samout**: object -- Write out all SAM alignment records into SAM/BAM files (one per input file needed), annotating each line with its feature assignment (as an optional field with tag 'XF'). See the -p option to use BAM instead of SAM.
            | **samoutformat**: [sam / bam] -- Format to use with the --samout / -o option. Choices : sam or bam (default: sam)
            | **delimiter**: [TAB] -- Column delimiter in output (default: TAB).
            | **appendoutput**: bool -- Append counts output. 
                This option is useful if you have already creates a TSV/CSV/similar file with a header for your samples 
                (with additional columns for the feature name and any additionl attributes) and want to fill in the rest of the file.
            | **featurequery**: str -- Restrict to features descibed in this expression. 
                Currently supports a single kind of expression: attribute == "one attr" to restrict the GFF to a single gene or transcript, 
                e.g. --feature-query 'gene_name == "ACTB"' - notice the single quotes around the argument of this option and the double quotes around the gene name. Broader queries might become available in the future.
            | **quiet**: bool -- Suppress progress report

    Returns
    -------
    output: object
        TSV file with counts

    Raises
    ------
    log: object
        gets the shell error output

    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")
    >>> log = _StubParams("log")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_htseqcount(input, params, output, log)
    "htseq-count [-r <{{params.order}}, default=None>] [--max-reads-in-buffer <{{params.maxreadsinbuffer}}, default=None>] [-s <{{params.stranded}}, default=None>] [-a <{{params.minaqual}}, default=None>] [-t <{{params.featuretype}}, default=None>] -i 'gene_id' --additional-attr 'gene_name' [--add-chromosome-info] [-m <{{params.mode}}, default=None>] [--nonunique <{{params.nonunique}}, default=None>] [--secondary-alignments <{{params.secondaryalignments}}, default=None>] [--supplementary-alignments <{{params.supplementaryalignments}}, default=None>] [-o <{{params.samout}}, default=None>] [-p <{{params.samoutformat}}, default=None>] [-d <{{params.delimiter}}, default=None>] [--append-output] [--feature-query <{{params.featurequery}}, default=None>] [-q] [-n <{{params.threads}}, default=None>] -c {output} {input} {params.annotation} 2>{log}"

    Specific cases

    >>> params = Params(fromdict={"stranded": "yes", "featuretype": "exon", "mode": "intersection-nonempty", "supplementaryalignments": "ignore", "threads": "1", "annotation": "test.gtf"})
    >>> create_htseqcount(input, params, output, log)
    "htseq-count -s yes -t exon -i 'gene_id' --additional-attr 'gene_name' -m intersection-nonempty --supplementary-alignments ignore -n 1 -c {output} {input} test.gtf 2>{log}"

    >>> params = Params(fromdict={"stranded": "no", "featuretype": "exon", "mode": "intersection-nonempty", "supplementaryalignments": "ignore", "annotation": "test.gff"})
    >>> create_htseqcount(input, params, output, log)
    "htseq-count -s no -t exon -i 'gene' -m intersection-nonempty --supplementary-alignments ignore -c {output} {input} test.gff 2>{log}"
    
    """


    result = join_str(
        "htseq-count",
        optional(params, "order", "-r {}"),
        optional(params, "maxreadsinbuffer", "--max-reads-in-buffer {}"),
        optional(params, "stranded", "-s {}"),
        optional(params, "minaqual", "-a {}"),
        optional(params, "featuretype", "-t {}"),
        ("-i 'gene'" if is_gff3_ext(params.annotation) else "-i 'gene_id'"),
        ("--additional-attr 'gene_name'" if is_gff3_ext(params.annotation) == False else ""),
        optional(params, "addchrinfo", "--add-chromosome-info"),
        optional(params, "mode", "-m {}"),
        optional(params, "nonunique", "--nonunique {}"),
        optional(params, "secondaryalignments", "--secondary-alignments {}"),
        optional(params, "supplementaryalignments", "--supplementary-alignments {}"),
        optional(params, "samout", "-o {}"),
        optional(params, "samoutformat", "-p {}"),
        optional(params, "delimiter", "-d {}"),
        optional(params, "appendoutput", "--append-output"),
        optional(params, "featurequery", "--feature-query {}"),
        optional(params, "quiet", "-q"),
        optional(params, "threads", "-n {}"),
        f"-c {output}", 
        f"{input} {params.annotation}",
        f"2>{log}"
    )
    return result


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
