"""
fastp  wrapper
==============

A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance.


Citation
--------
Shifu Chen. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, Volume 34, Issue 17, 01 September 2018, Pages i884-i890 DOI: https://doi.org/10.1093/bioinformatics/bty560


Documentation
-------------
https://github.com/OpenGene/fastp
https://manpages.debian.org/testing/fastp/fastp.1.en.html


"""


from wrappers.wrapper_system import optional, join_str,  mappy # for snakemake execution


__author__ = "Julie Ripoll & Celine Mandier"
__created__ = "2023-02-21"
__license__ = "CeCILL"
__version__ = "0.1.0"
__vwrapper__ = "0.1.0"
__vfastp__ = "0.23.2"


def create_fastp(input, output, params, log, **kwargs):
    r"""
    Cleaning your data 

    Parameters
    ----------
    input: object
        **fareads**: fasta or fastq file(s)

    params: _dict
        Essential options:
            | **threads**: int -- worker thread number, default is 3 (int [=3]), default in snakemake configuration file: THREADS
            | **paired**: ["YES" / "NO"] -- can be "YES" for paired-end data and "NO" for single-end
            | **qualityphred**: int -- the quality value that a base is qualified. 
                Default 15 means phred quality >=Q15 is qualified. (int [=15])
            | **lengthrequired**: int -- reads shorter than length_required will be discarded, default is 15. (int [=15])
            | **trimpolyx**: bool -- enable polyX trimming in 3' ends.
            | **polyxminlen**: int -- the minimum length to detect polyX in the read tail. 10 by default. (int [=10])

        Manage by wrapper:
            | **detectadapterforpe**: bool -- by default, the auto-detection for adapter is for SE data input only, 
                turn on this option to enable it for PE data.

        Other options:
            | **unpaired1**: object -- for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. 
                Default is to discard it. (string [=])
            | **unpaired2**: object -- for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. 
                If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])
            | **failedout**: object -- specify the file to store reads that cannot pass the filters. (string [=])
            | **dontoverwrite**: bool -- don't overwrite existing files. 
                Overwritting is allowed by default.
            | **disablequalityfiltering**: bool -- quality filtering is enabled by default. 
                If this option is specified, quality filtering is disabled
            | **nbaselimit**: int -- if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
            | **unqualifiedpercentlimit**: int -- how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
            | **averagequal**: int -- if one read's average quality score <avg_qual, then this read/pair is discarded. 
                Default 0 means no requirement (int [=0])
            | **disablelengthfiltering**: bool -- length filtering is enabled by default. If this option is specified, length filtering is disabled
            | **lengthlimit**: int -- reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
            | **lowcomplexityfilter**: bool -- enable low complexity filter. 
                The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1])
            | **complexitythreshold**: int -- the threshold for low complexity filter (0~100). 
                Default is 30, which means 30% complexity is required. (int [=30])
            | **disableadaptertrimming**: bool -- adapter trimming is enabled by default. 
                If this option is specified, adapter trimming is disabled
            | **adaptersequence**: string -- the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. 
                For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
            | **adaptersequencer2**: string -- the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. 
                If not specified, it will be the same as <adapter_sequence> (string [=auto])
            | **adapterfasta**: object -- specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
            | **cutfront**: bool -- move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, 
                stop otherwise.
            | **cuttail**: bool -- move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, 
                stop otherwise.
            | **cutright**: bool -- move a sliding window from front to tail, if meet one window with mean quality < threshold, 
                drop the bases in the window and the right part, and then stop.
            | **cutwindowsize**: int -- the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
            | **cutmeanquality**: int -- the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. 
                Range: 1~36 default: 20 (Q20) (int [=20])
            | **cutfrontwindowsize**: int -- the window size option of cut_front, default to cut_window_size if not specified (int [=4])
            | **cutfrontmeanquality**: int -- the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
            | **cuttailwindowsize**: int -- the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
            | **cuttailmeanquality**: int -- the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
            | **cutrightwindowsize**: int -- the window size option of cut_right, default to cut_window_size if not specified (int [=4])
            | **cutrightmeanquality**: int -- the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
            | **correction**: bool -- enable base correction in overlapped regions (only for PE data), default is disabled
            | **overlaplenrequire**: int -- the minimum length to detect overlapped region of PE reads. 
                This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
            | **overlapdifflimit**: int -- the maximum number of mismatched bases to detect overlapped region of PE reads. 
                This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
            | **overlapdiffpercentlimit**: int -- the maximum percentage of mismatched bases to detect overlapped region of PE reads. 
                This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
            | **trimfront1**: int -- trimming how many bases in front for read1, default is 0 (int [=0])
            | **trimfront2**: int -- trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
            | **trimtail1**: int -- trimming how many bases in tail for read1, default is 0 (int [=0])
            | **trimtail2**: int -- trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
            | **maxlen1**: int -- if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. 
                Default 0 means no limitation (int [=0])
            | **maxlen2**: int -- if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. 
                Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
            | **trimpolyg**: bool -- force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
            | **polygminlen**: int -- the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
            | **disabletrimpolyg**: bool -- disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
            | **umi**: bool -- enable unique molecular identifier (UMI) preprocessing
            | **umiloc**: string -- specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
            | **umilen**: int -- if the UMI is in read1/read2, its length should be provided (int [=0])
            | **umiprefix**: string -- if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). 
                No prefix by default (string [=])
            | **umiskip**: int -- if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])
            | **split**: int -- split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to 
                output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
            | **splitbylines**: int -- split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added 
                to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
            | **splitprefixdigits**: int -- the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 
                0001.xxx, 0 to disable padding (int [=4])
            | **overrepresentationanalysis**: bool -- enable overrepresented sequence analysis.
            | **overrepresentationsampling**: int -- one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis 
                (1~10000), smaller is slower, default is 20. (int [=20])
            | **merge**: bool -- for paired-end input, merge each pair of reads into a single read if they are overlapped. 
                The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
            | **mergeout**: string -- in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged 
                output (string [=])
            | **includeunmerged**: bool -- in the merging mode, write the unmerged or unpaired reads to the file specified by --merge. 
                Disabled by default.
            | **dedup**: bool -- enable deduplication to drop the duplicated reads/pairs
            | **dupcalcaccuracy**: int -- accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G). 
                Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
            | **dontevalduplication: bool -- don't evaluate duplication rate to save time and use less memory.
            | **reporttitle**: str -- should be quoted with ' or ", default is "fastp report" (string [=fastp report])
            | **compression**: int -- compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
            | **phred64**: bool -- indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
            | **stdin**: bool -- input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
            | **stdout**: bool -- stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. 
                Disabled by default.
            | **interleavedin**: bool -- indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.
            | **reads2process**: int -- specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
            | **fixmgiid**: bool -- the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.
            | **filterindex1**: object -- specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
            | **filterindex2**: object -- specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
            | **filterindexthreshold**: int -- the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])
            | **verbose**: bool -- output verbose log information (i.e. when every 1M reads are processed).

    Returns
    -------
    output: _dict
    
        | **cutreads**: object -- cleaned fasta or fastq files
        | **json**: object -- reports in json
        | **html**: object -- reports in html
        | **done**: object -- gets the shell output

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
    >>> create_fastp(input, output, params, log)
    "fastp [--unpaired1 <{{params.unpaired1}}, default=None>] [--unpaired2 <{{params.unpaired2}}, default=None>] [--failed_out <{{params.failedout}}, default=None>] [--dont_overwrite] [--disable_quality_filtering] [--n_base_limit <{{params.nbaselimit}}, default=None>] [--qualified_quality_phred <{{params.qualityphred}}, default=None>] [--unqualified_percent_limit <{{params.unqualifiedpercentlimit}}, default=None>] [--average_qual <{{params.averagequal}}, default=None>] [--disable_length_filtering] [--length_required <{{params.lengthrequired}}, default=None>] [--length_limit <{{params.lengthlimit}}, default=None>] [--low_complexity_filter] [--complexity_threshold] [--disable_adapter_trimming] [--adapter_sequence <{{params.adaptersequence}}, default=None>] [--adapter_sequence_r2 <{{params.adaptersequencer2}}, default=None>] [--adapter_fasta <{{params.adapterfasta}}, default=None>] [--cut_front] [--cut_tail] [--cut_right] [--cut_window_size <{{params.cutwindowsize}}, default=None>] [--cut_mean_quality <{{params.cutmeanquality}}, default=None>] [--cut_front_window_size <{{params.cutfrontwindowsize }}, default=None>] [--cut_front_mean_quality <{{params.cutfrontmeanquality}}, default=None>] [--cut_tail_window_size <{{params.cuttailwindowsize }}, default=None>] [--cut_tail_mean_quality <{{params.cuttailmeanquality}}, default=None>] [--cut_righ_window_size <{{params.cutrightwindowsize }}, default=None>] [--cut_righ_mean_quality <{{params.cutrightmeanquality}}, default=None>] [--correction] [--overlap_len_require <{{params.overlaplenrequire}}, default=None>] [--overlap_diff_limit <{{params.overlapdifflimit}}, default=None>] [--overlap_diff_percent_limit <{{params.overlapdiffpercentlimit}}, default=None>] [--trim_front1 <{{params.trimfront1}}, default=None>] [--trim_front2 <{{params.trimfront2}}, default=None>] [--trim_tail1 <{{params.trimtail1}}, default=None>] [--trim_tail2 <{{params.trimtail2}}, default=None>] [--trim_poly_g] [--poly_g_min_len <{{params.polygminlen}}, default=None>] [--disable_trim_poly_g] [--trim_poly_x] [--poly_x_min_len <{{params.polyxminlen}}, default=None>] [--umi] [--umi_loc <{{params.umiloc}}, default=None>] [--umi_len <{{params.umilen}}, default=None>] [--umi_prefix <{{params.umiprefix}}, default=None>] [--umi_skip <{{params.umiskip}}, default=None>] [--split <{{params.split}}, default=None>] [--split_by_lines <{{params.splitbylines}}, default=None>] [--split_prefix_digits <{{params.splitprefixdigits}}, default=None>] [--overrepresentation_analysis] [--overrepresentation_sampling <{{params.overrepresentationsampling}}, default=None>] [merge] [--merge_out <{{params.mergeout}}, default=None>] [--include_unmerged] [--dedup] [--dup_calc_accuracy <{{params.dupcalcaccuracy}}, default=None>] [--dont_eval_duplication] [--report_title <{{params.reporttitle}}, default=None>] [--compression <{{params.compression}}, default=None>] [--phred64] [--stdin] [--stdout] [--interleaved_in] [--reads_to_process <{{params.reads2process}}, default=None>] [--fix_mgi_id] [--filter_by_index1 <{{params.filterindex1}}, default=None>] [--filter_by_index2 <{{params.filterindex2}}, default=None>] [--filter_by_index_threshold <{{params.filterindexthreshold}}, default=None>] [--verbose] [--thread <{{params.threads}}, default=None>] <mappy('{params.paired}', {{'YES': '--out1 {{output.cutreads1}} --out2 {{output.cutreads2}}', 'NO': '-o {{output.cutreadse}}'}})> <mappy('{params.paired}', {{'YES': '--in1 {{input.fareads1}} --in2 {{input.fareads2}}', 'NO': '-i {{input.fareadse}}'}})> --json {output.json} --html {output.html} 2>{log}"

    Specific cases

    >>> params = Params(fromdict={"paired": "YES", "adapterfasta": "adapters.file", "threads": "3"})
    >>> create_fastp(input, output, params, log)
    'fastp --detect_adapter_for_pe --adapter_fasta adapters.file --thread 3 --out1 {output.cutreads1} --out2 {output.cutreads2} --in1 {input.fareads1} --in2 {input.fareads2} --json {output.json} --html {output.html} 2>{log}'

    >>> params = Params(fromdict={"paired": "NO", "threads": "3"})
    >>> create_fastp(input, output, params, log)
    'fastp --thread 3 -o {output.cutreadse} -i {input.fareadse} --json {output.json} --html {output.html} 2>{log}'

    >>> params = Params(fromdict={"paired": "NO", "trimpolyg": True})
    >>> create_fastp(input, output, params, log)
    'fastp --trim_poly_g -o {output.cutreadse} -i {input.fareadse} --json {output.json} --html {output.html} 2>{log}'

    >>> params = Params(fromdict={"paired": "NO", "trimpolyx": True, "polyxminlen": "10", "lengthrequired": "91", "qualityphred": "33"})
    >>> create_fastp(input, output, params, log)
    'fastp --qualified_quality_phred 33 --length_required 91 --trim_poly_x --poly_x_min_len 10 -o {output.cutreadse} -i {input.fareadse} --json {output.json} --html {output.html} 2>{log}'


    """

    fareads = mappy(params, "paired", {"YES": f"--in1 {{input.fareads1}} --in2 {{input.fareads2}}", "NO": f"-i {{input.fareadse}}"})

    cutreads = mappy(params, "paired", {"YES": f"--out1 {{output.cutreads1}} --out2 {{output.cutreads2}}", "NO": f"-o {{output.cutreadse}}"})

    return join_str(
        "fastp",
        optional(params, 'unpaired1', "--unpaired1 {}"),
        optional(params, 'unpaired2', "--unpaired2 {}"),
        optional(params, 'failedout', "--failed_out {}"),
        optional(params, 'dontoverwrite', "--dont_overwrite"),
        optional(params, 'disablequalityfiltering', "--disable_quality_filtering"),
        optional(params, 'nbaselimit', "--n_base_limit {}"),
        optional(params, 'qualityphred', "--qualified_quality_phred {}"),
        optional(params, 'unqualifiedpercentlimit', "--unqualified_percent_limit {}"),
        optional(params, 'averagequal', "--average_qual {}"),
        optional(params, 'disablelengthfiltering', "--disable_length_filtering"),
        optional(params, 'lengthrequired', "--length_required {}"),
        optional(params, 'lengthlimit', "--length_limit {}"),
        optional(params, 'lowcomplexityfilter', "--low_complexity_filter"),
        optional(params, 'complexitythreshold', "--complexity_threshold"),
        optional(params, 'disableadaptertrimming', "--disable_adapter_trimming"),
        optional(params, 'adaptersequence', "--adapter_sequence {}"),
        optional(params, 'adaptersequencer2', "--adapter_sequence_r2 {}"),
        ("--detect_adapter_for_pe" if params.paired == "YES" else ""),
        optional(params, 'adapterfasta', "--adapter_fasta {}"),
        optional(params, 'cutfront', "--cut_front"),
        optional(params, 'cuttail', "--cut_tail"),
        optional(params, 'cutright', "--cut_right"),
        optional(params, 'cutwindowsize', "--cut_window_size {}"),
        optional(params, 'cutmeanquality', "--cut_mean_quality {}"),
        optional(params, 'cutfrontwindowsize ', "--cut_front_window_size {}"),
        optional(params, 'cutfrontmeanquality', "--cut_front_mean_quality {}"),
        optional(params, 'cuttailwindowsize ', "--cut_tail_window_size {}"),
        optional(params, 'cuttailmeanquality', "--cut_tail_mean_quality {}"),
        optional(params, 'cutrightwindowsize ', "--cut_righ_window_size {}"),
        optional(params, 'cutrightmeanquality', "--cut_righ_mean_quality {}"),
        optional(params, 'correction', "--correction"),
        optional(params, 'overlaplenrequire', "--overlap_len_require {}"),
        optional(params, 'overlapdifflimit', "--overlap_diff_limit {}"),
        optional(params, 'overlapdiffpercentlimit', "--overlap_diff_percent_limit {}"),
        optional(params, 'trimfront1', "--trim_front1 {}"),
        optional(params, 'trimfront2', "--trim_front2 {}"),
        optional(params, 'trimtail1', "--trim_tail1 {}"),
        optional(params, 'trimtail2', "--trim_tail2 {}"),
        optional(params, 'trimpolyg', "--trim_poly_g"),
        optional(params, 'polygminlen', "--poly_g_min_len {}"),
        optional(params, 'disabletrimpolyg', "--disable_trim_poly_g"),
        optional(params, 'trimpolyx', "--trim_poly_x"),
        optional(params, 'polyxminlen', "--poly_x_min_len {}"),
        optional(params, 'umi', "--umi"),
        optional(params, 'umiloc', "--umi_loc {}"),
        optional(params, 'umilen', "--umi_len {}"),
        optional(params, 'umiprefix', "--umi_prefix {}"),
        optional(params, 'umiskip', "--umi_skip {}"),
        optional(params, 'split', "--split {}"),
        optional(params, 'splitbylines', "--split_by_lines {}"),
        optional(params, 'splitprefixdigits', "--split_prefix_digits {}"),
        optional(params, 'overrepresentationanalysis', "--overrepresentation_analysis"),
        optional(params, 'overrepresentationsampling', "--overrepresentation_sampling {}"),
        optional(params, 'merge', "merge"),
        optional(params, 'mergeout', "--merge_out {}"),
        optional(params, 'includeunmerged ', "--include_unmerged"),
        optional(params, 'dedup ', "--dedup"),
        optional(params, 'dupcalcaccuracy', "--dup_calc_accuracy {}"),
        optional(params, 'dontevalduplication', "--dont_eval_duplication"),
        optional(params, 'reporttitle', "--report_title {}"),
        optional(params, 'compression', "--compression {}"),
        optional(params, 'phred64', "--phred64"),
        optional(params, 'stdin', "--stdin"),
        optional(params, 'stdout', "--stdout"),
        optional(params, 'interleavedin', "--interleaved_in"),
        optional(params, 'reads2process', "--reads_to_process {}"),
        optional(params, 'fixmgiid', "--fix_mgi_id"),
        optional(params, 'filterindex1', "--filter_by_index1 {}"),
        optional(params, 'filterindex2', "--filter_by_index2 {}"),
        optional(params, 'filterindexthreshold', "--filter_by_index_threshold {}"),
        optional(params, 'verbose', "--verbose"),
        optional(params, "threads", "--thread {}"),
        cutreads,
        fareads,
        f"--json {output.json}",
        f"--html {output.html}",
        f"2>{log}"
    )


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
