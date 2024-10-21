"""
Cutadapt wrapper
================

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.


Citation
--------
Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1):10-12, May 2011. DOI: http://dx.doi.org/10.14806/ej.17.1.200


Documentation
-------------
https://cutadapt.readthedocs.io/en/stable/index.html


"""

__author__="Julie Ripoll"
__created__="2021-01-05"
__license__="CeCILL"
__version__="0.1.0"
__vwrapper__="0.1.0"
__vcutadapt__="4.4"


from wrappers.wrapper_system import optional, join_str, mappy # for snakemake execution


def create_cutadapt(input, output, params, **kwargs):
    r"""
    Cleaning your data

    Parameters
    ----------
    input: object
        **fareads**: fasta or fastq file(s)

    params:
        Essential options:
            | **paired**: str -- can be paired-end "YES" or single-end "NO" (mandatory)
            
            | **strandtrim**: bool -- require if adpaters to trim (mandatory)
            | **strand**: int -- can be 3, 5 or 35 according to strand to trim, require strandtrim at True (mandatory)
            
            | **threads**: int -- cores value
            
            | **action**: str -- {trim,retain,mask,lowercase,none}
                What to do if a match was found. trim: trim adapter and up- or downstream
                sequence; retain: trim, but retain adapter; mask: replace with 'N' characters;
                lowercase: convert to lowercase; none: leave unchanged. Default: trim
            | **adapter1**: object -- Sequence of an adapter ligated to the read
            | **adapter2**: object -- adapter to be removed from second read in a pair, require paired as "YES" and strandtrim as True
            
            | **cut**: bool -- require if cut of bases
            | **cutlength1**: int -- Remove bases from each read (first read only if paired).
                If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming. Require cut as True
            | **cutlength2**: int -- Remove LENGTH bases from second read in a pair, require  paired as "YES" and cut as True

            | **debug**: bool -- print debugging information
        
            | **trimn**: bool -- Trim N's on ends of reads
            | **maxn**: count -- Discard reads with more than COUNT 'N' bases.
                If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
            | **polya**: bool -- Trim poly-A tails
            
                If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
            | **length**: int -- Shorten reads to LENGTH
            | **minimumlength**: int -- Discard trimmed reads that are shorter than LENGTH.
                Reads that are too short even before adapter removal are also discarded.
                In colorspace, an initial primer is not counted.
                Default: 0
            | **maximumlength**: int -- Discard trimmed reads that are longer than LENGTH.
                Reads that are too long even before adapter removal are also discarded.
                In colorspace, an initial primer is not counted. 
                Default: no limit
            
            | **qualitycutoff**: int -- Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal.
                Applied to both reads if data is paired.
                If one value is given, only the 3' end is trimmed.
        
        Others options:
            | **discardtrimmed**: bool -- Discard reads that contain an adapter. Use also -O to avoid discarding too many
                randomly matching reads.
            | **discarduntrimmed**: bool -- Discard reads that do not contain an adapter.
            | **discardcasava**: bool -- Discard reads that did not pass CASAVA filtering (header has :Y:).
            | **errorrate**: bool -- Maximum allowed error rate (no. of errors divided by the length of the matching region).
                Default: 0.1
            | **fasta**: bool -- Output FASTA to standard output even on FASTQ input.
            | **infofile**: object -- write information about each read and its adapter matches to file
            | **interleaved**: bool -- Read and write interleaved paired-end reads
            | **json**: object -- Dump report in JSON format to FILE
            | **lengthtag**: tag -- Search for TAG followed by a decimal number in the description field of the read.
                Replace the decimal number with the correct length of the trimmed read.
                For example, use --length-tag 'length=' to correct fields like 'length=123'.
            | **matchreadwildcards**: bool -- Interpret IUPAC wildcards in reads.
                Default: False
            | **maxexpectederrors**: int -- Discard reads whose expected number of errors (computed from quality values)
                exceeds ERRORS.
            | **nextseqtrim**: bool -- NextSeq-specific quality trimming (each read).
                Trims also dark cycles appearing as high-quality G bases (EXPERIMENTAL).
            | **noindels**: bool -- AAllow only mismatches in alignments.
                Default: allow both mismatches and indels
            | **nomatchadapterwildcards**: bool -- Do not interpret IUPAC wildcards in adapters
            | **overlap**: int -- If the overlap between the read and the adapter is shorter than MINLENGTH, the read is not modified.
                Reduces the no. of bases trimmed due to random adapter matches.
                Default: 3
            | **pairadapters**: bool -- Treat adapters given with -a/-A etc. as pairs. Either both or none are removed
                from each read pair.
            | **pairfilter**: str -- {any,both,first}
                Which of the reads in a paired-end read have to match the filtering criterion in
                order for the pair to be filtered. Default: any
            | **qualitybase**: int -- Assume that quality values in FASTQ are encoded as ascii(quality + QUALITY_BASE).
                This needs to be set to 64 for some old Illumina FASTQ files.
                Default: 33
            | **report**: str -- {full,minimal}
                Which type of report to print: 'full' or 'minimal'. Default: full
            | **rename**: object -- Rename reads using TEMPLATE containing variables such as {id}, {adapter_name}
                etc. (see documentation)
            | **restfile**: object -- When the adapter matches in the middle of a read, write the rest (after the
                adapter) to FILE.
            | **revcomp: bool -- Check both the read and its reverse complement for adapter matches. If match is
                on reverse-complemented version, output that one. Default: check only read
            | **stripsuffix**: str -- Remove this suffix from read names if present
            | **prefix: str -- add prefix to read names
                Use {name} to insert the name of the matching adapter.
            | **suffix**: str -- add suffix to read names
                Use {name} to insert the name of the matching adapter.
            | **times**: int -- Remove up to COUNT adapters from each read.
                Default: 1
            | **tooshortoutput**: object -- write readobject
            | **tooshortpairedoutput**: object -- Write second read in a pair to this file if pair is too short.
                Use together with --too-short-output.
            | **toolongoutput**: object -- write reads that are too long to file
                Default: discard reads
            | **toolongpairedoutput**: object -- Write second read in a pair to this file if pair is too long.
                Use together with --too-long-output.
            | **untrimmedoutput**: object -- Write reads that do not contain the adapter to FILE
            | **untrimmedpairedoutput**: object -- Use this option together with --untrimmed-output when trimming paired-end reads,
                Default: output to same file
            | **wildcardfile**: object -- When the adapter has N wildcard bases, write adapter bases matching wildcard
                positions to FILE. (Inaccurate with indels.)
            | **zerocap**: bool -- Change negative quality values to zero.

        .. deprecated::
            Deprecated options:
                | **fformat**: str -- Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. Ignored when reading csfasta/qual files.
                    Default: auto-detect from file name extension.
                | **maskadapter**: bool -- Mask adapters with 'N' characters instead of trimming them
                | **notrim**: bool -- Match and redirect reads to output/untrimmed-output as usual, but do not remove adapters
                | **colorspace**: bool -- Enable colorspace mode: 
                    Also trim the color that is adjacent to the found adapter.
                | **doubleencode**: bool -- Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).
                | **stripf3**: bool -- strip the _F3 suffix of read names
                | **maq**: bool -- maq color. This enables -c, -d, -t, --strip-f3 and -y '/1'.
                | **bwa**: bool -- bwa color. This enables -c, -d, -t, --strip-f3 and -y '/1'.
                | **nozerocap**: bool -- Do not change negative quality values to zero in colorspace data.
                    By default, they are since many tools have problems with negative qualities.

    Returns
    -------
    output: _dict
    
        | **cutreads**: object -- Write trimmed reads to FILE.
            FASTQ or FASTA format is chosen depending on input.
            The summary report is sent to standard output. 
            Use '{name}' in FILE to demultiplex reads into multiple files. Default: write to standard output
        | **done**: object -- gets the shell output (for summary report)


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> input = _StubParams("input")
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_cutadapt(input, output, params)
    'cutadapt <mappy(\'{params.paired}\', {{\'YES\': "<mappy(\'{{params.strand}}\', {{{{\'3\': \'-a {{{{params.adapter1}}}} -A {{{{params.adapter2}}}}\', \'5\': \'-g {{{{params.adapter1}}}} -G {{{{params.adapter2}}}}\', \'35\': \'-b {{{{params.adapter1}}}} -B {{{{params.adapter2}}}}\'}}}})>", \'NO\': "<mappy(\'{{params.strand}}\', {{{{\'3\': \'-a {{{{params.adapter1}}}}\', \'5\': \'-g {{{{params.adapter1}}}}\', \'35\': \'-b {{{{params.adapter1}}}}\'}}}})>"}})> [--cores] [--debug] [--action <{{params.action}}, default=None>] <mappy(\'{params.paired}\', {{\'YES\': \'-u {{params.cutlength1}} -U {{params.cutlength2}}\', \'NO\': \'-u {{params.cutlength1}}\'}})> <mappy(\'{params.paired}\', {{\'YES\': \'-u {{params.morelength1}} -U {{params.morelength2}}\', \'NO\': \'-u {{params.morelength1}}\'}})> [--quality-cutoff <{{params.qualitycutoff}}, default=None>] [--length <{{params.length}}, default=None>] [--minimum-length <{{params.minimumlength}}, default=None>] [--maximum-length <{{params.maximumlength}}, default=None>] [--trim-n] [--max-n <{{params.maxn}}, default=None>] [--poly-a] [--discard-trimmed <{{params.discardtrimmed}}, default=None>] [--discard-untrimmed] [--discard-casava] [---error-rate] [--fasta] [--info-file <{{params.infofile}}, default=None>] [--interleaved] [--json <{{params.json}}, default=None>] [--length-tag <{{params.lengthtag}}, default=None>] [--match-read-wildcards] [--max-expected-errors <{{params.maxexpectederrors}}, default=None>] [--nextseq-trim] [--no-indels] [--no-match-adapter-wildcards] [--overlap <{{params.overlap}}, default=None>] [--pair-adapters] [--pair-filter <{{params.pairfilter}}, default=None>] [--quality-base <{{params.qualitybase}}, default=None>] [--rename <{{params.rename}}, default=None>] [--report <{{params.report}}, default=None>] [--rest-file <{{params.restfile}}, default=None>] [--revcomp] [--strip-suffix <{{params.stripsuffix}}, default=None>] [--prefix <{{params.prefix}}, default=None>] [--suffix <{{params.suffix}}, default=None>] [--times <{{params.times}}, default=None>] [--too-short-output <{{params.tooshortoutput}}, default=None>] [--too-short-paired-output <{{params.tooshortpairedoutput}}, default=None>] [--too-long-output <{{params.toolongoutput}}, default=None>] [--too-long-paired-output <{{params.toolongpairedoutput}}, default=None>] [--untrimmed-output <{{params.untrimmedoutput}}, default=None>] [--untrimmed-paired-output <{{params.untrimmedpairedoutput}}, default=None>] [--wildcard-file <{{params.wildcardfile}}, default=None>] [--zero-cap] <mappy(\'{params.paired}\', {{\'YES\': \'-o {{output.cutreads1}} -p {{output.cutreads2}} \', \'NO\': \'-o {{output.cutreadse}} \'}})> <mappy(\'{params.paired}\', {{\'YES\': \'{{input.fareads1}} {{input.fareads2}}\', \'NO\': \'{{input.fareadse}}\'}})> 1>{output.infos}'

    Specific cases

    >>> params = Params(fromdict={"cut": True, "paired": "YES", "strandtrim": True, "strand": "35"})
    >>> create_cutadapt(input, output, params)
    'cutadapt -b {params.adapter1} -B {params.adapter2} -u {params.cutlength1} -U {params.cutlength2} -u {params.morelength1} -U {params.morelength2} -o {output.cutreads1} -p {output.cutreads2}  {input.fareads1} {input.fareads2} 1>{output.infos}'

    >>> params = Params(fromdict={"cut": True, "paired": "NO", "strandtrim": True, "strand": "5"})
    >>> create_cutadapt(input, output, params)
    'cutadapt -g {params.adapter1} -u {params.cutlength1} -u {params.morelength1} -o {output.cutreadse}  {input.fareadse} 1>{output.infos}'

    >>> params = Params(fromdict={"cut": False, "paired": "NO", "strandtrim": True, "strand": "3"})
    >>> create_cutadapt(input, output, params)
    'cutadapt -a {params.adapter1} -o {output.cutreadse}  {input.fareadse} 1>{output.infos}'
    
    """
    
    fareads = mappy(params, "paired", {"YES": f"{{input.fareads1}} {{input.fareads2}}", "NO": f"{{input.fareadse}}"})

    cutreads = mappy(params, "paired", {"YES": "-o {output.cutreads1} -p {output.cutreads2} ", "NO": "-o {output.cutreadse} "})

    cut = mappy(params, "cut", (True, False), default = False)

    cutlengths = mappy(params, "paired", {"YES": "-u {params.cutlength1} -U {params.cutlength2}" if cut else "", "NO": "-u {params.cutlength1}" if cut else ""})

    morelengths = mappy(params, "paired", {"YES": "-u {params.morelength1} -U {params.morelength2}" if cut else "", "NO": "-u {params.morelength1}" if cut else ""})

    strandtrim = mappy(params, "strandtrim", (True, False), default = False)

    strandsense = mappy(params, "paired", {"YES": mappy(params, "strand", 
                                                            {"3": "-a {params.adapter1} -A {params.adapter2}" if strandtrim else "",
                                                            "5": "-g {params.adapter1} -G {params.adapter2}" if strandtrim else "",
                                                            "35": "-b {params.adapter1} -B {params.adapter2}" if strandtrim else ""}),
                                            "NO": mappy(params, 'strand', 
                                                            {"3": "-a {params.adapter1}" if strandtrim else "",
                                                            "5": "-g {params.adapter1}" if strandtrim else "",
                                                            "35": "-b {params.adapter1}" if strandtrim else ""}) } )

    return join_str(
        "cutadapt",
        # essential options
        strandsense,
        optional(params, 'threads', "--cores"),
        optional(params, 'debug', "--debug"),
        optional(params, 'action', "--action {}"),
        cutlengths,
        morelengths,
        optional(params, 'qualitycutoff', "--quality-cutoff {}"),
        optional(params, 'length', "--length {}"),
        optional(params, 'minimumlength', "--minimum-length {}"),
        optional(params, 'maximumlength', "--maximum-length {}"),
        optional(params, 'trimn', "--trim-n"),
        optional(params, 'maxn', "--max-n {}"),
        optional(params, 'polya', "--poly-a"),
        # other options
        optional(params, 'discardtrimmed', "--discard-trimmed {}"),
        optional(params, 'discarduntrimmed', "--discard-untrimmed"),
        optional(params, 'discardcasava', "--discard-casava"),
        optional(params, 'errorrate', "---error-rate"),
        optional(params, 'fasta', "--fasta"),
        optional(params, 'infofile', "--info-file {}"),
        optional(params, 'interleaved', "--interleaved"),
        optional(params, 'json', "--json {}"),
        optional(params, 'lengthtag', "--length-tag {}"),
        optional(params, 'matchreadwildcards', "--match-read-wildcards"),
        optional(params, 'maxexpectederrors', "--max-expected-errors {}"),
        optional(params, 'nextseqtrim', "--nextseq-trim"),
        optional(params, 'noindels', "--no-indels"),
        optional(params, 'nomatchadapterwildcards', "--no-match-adapter-wildcards"),
        optional(params, 'overlap', "--overlap {}"),
        optional(params, 'pairadapters', "--pair-adapters"),
        optional(params, 'pairfilter', "--pair-filter {}"),
        optional(params, 'qualitybase', "--quality-base {}"),
        optional(params, 'rename', "--rename {}"),
        optional(params, 'report', "--report {}"),
        optional(params, 'restfile', "--rest-file {}"),
        optional(params, 'revcomp', "--revcomp"),
        optional(params, 'stripsuffix', "--strip-suffix {}"),
        optional(params, 'prefix', "--prefix {}"),
        optional(params, 'suffix', "--suffix {}"),
        optional(params, 'times', "--times {}"),
        optional(params, 'tooshortoutput', "--too-short-output {}"),
        optional(params, 'tooshortpairedoutput', "--too-short-paired-output {}"),
        optional(params, 'toolongoutput', "--too-long-output {}"),
        optional(params, 'toolongpairedoutput', "--too-long-paired-output {}"),
        optional(params, 'untrimmedoutput', "--untrimmed-output {}"),
        optional(params, 'untrimmedpairedoutput', "--untrimmed-paired-output {}"),
        optional(params, 'wildcardfile', "--wildcard-file {}"),
        optional(params, 'zerocap', "--zero-cap"),
        # input / output
        cutreads,
        fareads,
        ## log-error output
        f"1>{output.infos}"
        )



if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
