"""
featureCounts wrapper
=====================

featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. 
It can be used to count both RNA-seq and genomic DNA-seq reads. It is available in the SourceForge Subread package or the Bioconductor Rsubread package. 


Citation
--------
Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014
Liao Y, Smyth GK and Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108, 2013


Documentation
-------------
https://subread.sourceforge.net/
https://janis.readthedocs.io/en/latest/tools/bioinformatics/subread/featurecounts.html


"""

__author__="Julie Ripoll & Celine Mandier" 
__created__="2023-03-20"
__license__="CeCILL"
__version__="0.1.0"
__vwrapper__ ="0.1.0"
__vsubread__="2.0.1"


from wrappers.wrapper_system import optional, join_str
from wrappers.format_checker import is_gff3_ext


def create_featurecounts(input, params, output, log, **kwargs):
    r"""
    Given a file with aligned sequencing reads and a list of genomic features, featureCounts counts how many reads map to each feature.

    Parameters
    ----------
    input: object
        SAM/BAM alignment files
            A list of SAM or BAM format files. They can be either name or location sorted.
            If no files provided, <stdin> input is expected.
            Location-sorted paired-end reads are automatically sorted by read names.

    params: _dict
        Essential options:
            | **annotation**: object -- GTF/GFF file for annotation. See -F option for more format information.
                Inbuilt annotations (SAF format) is available in 'annotation' directory of the package.
                Gzipped file is also accepted.
            | **threads**: int -- Number of the threads. 1 by default. Default in snakemake configuration file: THREADS
            | **formatREF**: [GTF / GFF / SAF] -- Specify format of the provided annotation file.
                Acceptable formats include 'GTF' (or compatible GFF format) and 'SAF'. 'GTF' by default.
                For SAF format, please refer to Users Guide.
            | **featuretype**: str -- Specify feature type(s) in a GTF annotation.
                If multiple types are provided, they should be separated by ',' with no space in between. 'exon' by default.
                Rows in the annotation with a matched feature will be extracted and used for read mapping.
            | **strand**: str -- Perform strand-specific read counting.
                A single integer value (applied to all input files) or a string of comma-separated values
                (applied to each corresponding input file) should be provided.
                Possible values include: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
                Default value is 0 (ie. unstranded read counting carried out for all input files).

        Manage by wrapper:
            | **attributetype**: str -- Specify attribute type in GTF annotation. 'gene_id' by  default.
                Meta-features used for read counting will be  extracted from annotation using the provided value.
            | **extraattributes**: str -- Extract extra attribute types from the provided GTF annotation
                and include them in the counting output.
                These attribute types will not be used to group features.
                If more than one attribute type is provided they should be separated by comma.

        Other options:
            | **chralias**: str -- Provide a chromosome name alias file to match chr names in annotation with those in the reads.
                This should be a two-column comma-delimited text file.
                Its first column should include chr names in the annotation
                and its second column should include chr names in the reads.
                Chr names are case sensitive. No column header should be included in the file.
            | **countelmt**: bool -- Perform read counting at feature level
                (eg. counting reads for exons rather than genes).
            | **overlapmetafeatures**: bool -- Assign reads to all their overlapping meta-features
                (or features if -f is specified).
            | **minoverlap**: int -- Minimum number of overlapping bases in a read that is required for read assignment.
                1 by default.
                Number of overlapping bases is counted from both reads if paired end.
                If a negative value is provided, then a gap of up to specified size will be allowed between read and the feature 
                that the read is assigned to.
            | **fracoverlap**: float -- Minimum fraction of overlapping bases in a read that is required for read assignment.
                Value should be within range [0,1]. 0 by default.
                Number of overlapping bases is counted from both reads if paired end.
                Both this option and '--minOverlap' option need to be satisfied for read assignment.
            | **fracoverlapfeature**: float -- Minimum fraction of overlapping bases in a feature that is required for read assignment.
                Value should be within range [0,1]. 0 by default.
            | **largestoverlap**: bool -- Assign reads to a meta-feature/feature
                that has the largest number of overlapping bases.
            | **nonoverlap**: int -- Maximum number of non-overlapping bases in a read (or a read pair)
                that is allowed when being assigned to a feature.
                No limit is set by default.
            | **nonoverlapfeature**: int -- Maximum number of non-overlapping bases in a feature
                that is allowed in read assignment.
                No limit is set by default.
            | **readextension5**: int -- Reads are extended upstream by <int> bases from their 5' end.
            | **readextension3**: int -- Reads are extended upstream by <int> bases from their 3' end.
            | **multimap**: bool -- Multi-mapping reads will also be counted.
                For a multi-mapping read, all its reported alignments will be counted.
                The 'NH' tag in BAM/SAM input is used to detect  multi-mapping reads.
            | **fraction**: bool -- Assign fractional counts to features.
                This option must be used together with '-M' or '-O' or both.
                When '-M' is specified, each reported alignment from a multi-mapping read (identified via 'NH' tag)
                will carry a fractional count of 1/x, instead of 1 (one),
                where x is the total number of alignments reported for the same read.
                When '-O' is specified, each overlapping feature will receive a fractional count of 1/y,
                where y is the total number of features overlapping with the read.
                When both '-M' and '-O' are specified, each alignment will carry a fractional count of 1/(x*y).
            | **minqual**: str -- The minimum mapping quality score a read must satisfy in order to be counted.
                For paired-end reads, at least one end should satisfy this criteria. 0 by default.
            | **splitonly**: bool -- Count split alignments only (ie. alignments with CIGAR string containing 'N').
                An example of split alignments is exon-spanning reads in RNA-seq data.
            | **nonsplitonly**: bool -- If specified, only non-split alignments (CIGAR strings do not contain letter 'N') will be counted.
                All the other alignments will be ignored.
            | **primary: bool -- Count primary alignments only.
                Primary alignments are identified using bit 0x100 in SAM/BAM FLAG field.
            | **ignoredup**: bool -- Ignore duplicate reads in read counting.
                Duplicate reads are identified using bit Ox400 in BAM/SAM FLAG field.
                The whole read pair is ignored if one of the reads is a duplicate read for paired end data.
            | **junction**: str -- Count number of reads supporting each exon-exon junction.
                Junctions were identified from those exon-spanning reads in the input (containing 'N' in CIGAR string).
                Counting results are saved to a file named '<output_file>.jcounts'
            | **formatrefseq**: object -- Provide the name of a FASTA-format file
                that contains the reference sequences used in read mapping that produced 
                the provided SAM/BAM files. This optional argument can be used with '-J' option to improve read counting for junctions.
            | **countfragment**: bool -- If specified, fragments (or templates) will be counted instead of reads.
                This option is only applicable for paired-end reads; single-end reads are always counted as reads.
            | **countpairs**: bool -- Only count read pairs that have both ends aligned.
            | **validitydist**: bool -- Check validity of paired-end distance when counting read pairs.
                Use -d and -D to set thresholds.
            | **minifraglength**: int -- Minimum fragment/template length, 50 by default.
            | **maxfraglength**: int -- Maximum fragment/template length, 600 by default.
            | **chimerism**: bool -- Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same 
                chromosome but on different strands.
            | **donotsort**: bool -- Do not sort reads in BAM/SAM input.
                Note that reads from the same pair are required to be located next to each other in the input.
            | **byreadgroup**: bool -- Assign reads by read group.
                "RG" tag is required to be present in the input BAM/SAM files.
            | **longread**: bool -- count long reads such as Nanopore and PacBio reads.
                Long read counting can only run in one thread and only reads (not read-pairs) can be counted.
                There is no limitation on the number of 'M' operations allowed in a CIGAR string in long read counting.
            | **formatassign**: str -- Output detailed assignment results for each read or read-pair.
                Results are saved to a file that is in one of the following formats: CORE, SAM and BAM.
                See Users Guide for more info about these formats.
            | **pathassignresult**: path -- Specify a directory to save the detailed assignment results.
                If unspecified, the directory where counting results are saved is used.
            | **tmpdir**: str -- Directory under which intermediate files are saved (later removed).
                By default, intermediate files will be saved to the directory specified in '-o' argument.
            | **maxmop**: int -- Maximum number of 'M' operations allowed in a CIGAR string. 10 by default.
                Both 'X' and '=' are treated as 'M' and adjacent 'M' operations are merged in the CIGAR string.
            | **verbose**: bool -- Output verbose information for debugging, such as unmatched chromosome/contig names.

    Returns
    -------
    output: object
        Name of output file including read counts.
            A separate file including summary statistics of counting results is also included in the output ('<string>.summary').
            Both files are in tab delimited format.

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
    >>> create_featurecounts(input, params, output, log)
    "featureCounts [-F <{{params.formatREF}}, default=None>] [-t <{{params.featuretype}}, default=None>] -g 'gene_id' --extraAttributes 'gene_name' [-A <{{params.chralias}}, default=None>] [-f] [-O] [--minOverlap <{{params.minoverlap}}, default=None>] [--fracOverlap <{{params.fracoverlap}}, default=None>] [--fracOverlapFeature <{{params.fracoverlapfeature}}, default=None>] [--largestOverlap] [--nonOverlap <{{params.nonoverlap}}, default=None>] [--nonOverlapFeature <{{params.nonoverlapfeature}}, default=None>] [--readExtension5 <{{params.readextension5}}, default=None>] [--readExtension3 <{{params.readextension3}}, default=None>] [-M] [--fraction] [-Q <{{params.minqual}}, default=None>] [--splitOnly] [--nonSplitOnly] [--primary] [--ignoreDup] [-s <{{params.strand}}, default=None>] [-J <{{params.junction}}, default=None>] [-G <{{params.formatrefseq}}, default=None>] [-B] [-P] [-d <{{params.minifraglength}}, default=None>] [-D <{{params.maxfraglength}}, default=None>] [-C] [--donotsort] [--byReadGroup] [-L] [-R <{{params.formatassign}}, default=None>] [--Rpath <{{params.pathassignresult}}, default=None>] [--tmpDir <{{params.tmpdir}}, default=None>] [--maxMOp <{{params.maxmop}}, default=None>] [--verbose] [-T <{{params.threads}}, default=None>] -a {params.annotation} -o {output} {input} 2>{log}"

    Specific cases

    >>> params = Params(fromdict={"stranded": "yes", "featuretype": "exon", "mode": "intersection-nonempty", "idattr": "gene_id", "additionalattr": "gene_name", "supplementaryalignments": "ignore", "paired": "YES", "threads": "1", "annotation": "file.gtf"})
    >>> create_featurecounts(input, params, output, log)
    "featureCounts -t exon -g 'gene_id' --extraAttributes 'gene_name' -p -T 1 -a file.gtf -o {output} {input} 2>{log}"

    >>> params = Params(fromdict={"stranded": "no", "featuretype": "exon", "mode": "intersection-nonempty", "idattr": "gene_id", "additionalattr": "gene_name", "supplementaryalignments": "ignore", "paired": "NO", "annotation": "file.gff"})
    >>> create_featurecounts(input, params, output, log)
    "featureCounts -t exon -g 'gene' -a file.gff -o {output} {input} 2>{log}"
    
    >>> params = Params(fromdict={"stranded": "yes", "featuretype": "exon", "mode": "intersection-nonempty", "idattr": "gene_id", "additionalattr": "gene_name", "paired": "YES", "annotation": "file.gff"})
    >>> create_featurecounts(input, params, output, log)
    "featureCounts -t exon -g 'gene' -p -a file.gff -o {output} {input} 2>{log}"
    
    """
    result = join_str(
        "featureCounts",
        optional(params, "formatREF", "-F {}"),
        optional(params, "featuretype", "-t {}"),
        ("-g 'gene'" if is_gff3_ext(params.annotation) else "-g 'gene_id'"),
        ("--extraAttributes 'gene_name'" if is_gff3_ext(params.annotation) == False else ""),
        optional(params, "chralias", "-A {}"),
        optional(params, "countelmt", "-f"),
        optional(params, "overlapmetafeatures", "-O"),
        optional(params, "minoverlap", "--minOverlap {}"),
        optional(params, "fracoverlap", "--fracOverlap {}"),
        optional(params, "fracoverlapfeature", "--fracOverlapFeature {}"),
        optional(params, "largestoverlap", "--largestOverlap"),
        optional(params, "nonoverlap", "--nonOverlap {}"),
        optional(params, "nonoverlapfeature", "--nonOverlapFeature {}"),
        optional(params, "readextension5", "--readExtension5 {}"),
        optional(params, "readextension3", "--readExtension3 {}"),
        optional(params, "multimap", "-M"),
        optional(params, "fraction", "--fraction"),
        optional(params, "minqual", "-Q {}"),
        optional(params, "splitonly", "--splitOnly"),
        optional(params, "nonsplitonly", "--nonSplitOnly"),
        optional(params, "primary", "--primary"),
        optional(params, "ignoredup", "--ignoreDup"),
        optional(params, "strand", "-s {}"),
        optional(params, "junction", "-J {}"),
        optional(params, "formatrefseq", "-G {}"),
        ("-p" if params.paired  == "YES" else ""),
        optional(params, "countpairs", "-B"),
        optional(params, "validitydist", "-P"),
        optional(params, "minifraglength", "-d {}"),
        optional(params, "maxfraglength", "-D {}"),
        optional(params, "chimerism", "-C"),
        optional(params, "donotsort ", "--donotsort"),
        optional(params, "byreadgroup", "--byReadGroup"),
        optional(params, "longread", "-L"),
        optional(params, "formatassign", "-R {}"),
        optional(params, "pathassignresult", "--Rpath {}"),
        optional(params, "tmpdir", "--tmpDir {}"),
        optional(params, "maxmop", "--maxMOp {}"),
        optional(params, "verbose", "--verbose"),
        optional(params, "threads", "-T {}"),
        f"-a {params.annotation}",
        f"-o {output}",
        f"{input}",
        f"2>{log}"
    )
    return result


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
