""" 
MultiQC wrapper
===============

MultiQC, a tool to create a single report visualising output from multiple tools across many samples, enabling global trends and biases to be quickly identified.


Citation
--------
Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics (2016). doi: 10.1093/bioinformatics/btw354


Documentation
-------------
https://multiqc.info/docs/


"""


from wrappers.wrapper_system import optional, join_str, mappy # for snakemake execution


__author__ = "Julie Ripoll"
__created__ = "2021-01-05"
__license__ = "CeCILL"
__version__ = "0.1.0"
__vwrapper__ = "0.1.0"
__vmultiqc__ = "1.14"


def create_multiqc(output, params, **kwargs):
    r"""
    MultiQC aggregates results from bioinformatics analyses across many samples into a single report.

    Parameters
    ----------
    input: object
        path to the raw input files or a check file for rule consistency

    params: _dict
        Essential options:
            | **module**: str -- Use only this module. Can specify multiple times
            | **title**: str -- Report title. Printed as page header, used for filename if not otherwise specified.
            | **filename**: str -- Report filename. Use 'stdout' to print to standard out.
            | **dirs**: object -- Prepend directory to sample names, specify *.extension of files required
            | **fullnames**: bool -- Do not clean the sample names (leave as full file name)
            | **zipdatadir**: bool -- Compress the data directory.
            | **export**: bool -- Export plots as static images in addition to the report
            | **interactive**: bool -- Interactive js plots
            | **verbose**: bool -- Increase output verbosity.

        Other options:
            | **viewtags**: bool -- View the available tags and which modules they load
            | **force**: bool -- Overwrite any existing reports
            | **config**: path -- Specific config file to load, after those in MultiQC dir / home dir / working dir.
            | **clconfig**: str -- Specify MultiQC config YAML on the command line
            | **ignore**: str -- Ignore analysis files (glob expression)
            | **ignoresamples**: str -- Ignore sample names (glob expression)
            | **ignoresymlinks**: bool -- Ignore symlinked directories and files
            | **filelist**: object -- Supply a file containing a list of file paths to be searched, one per row
            | **exclude**: str -- Do not use this module. Can specify multiple times.
            | **tag**: str -- Use only modules which tagged with this keyword, eg. RNA
            | **viewtags**: bool -- View the available tags and which modules they load
            | **dirsdepth**: int -- Prepend <int> directories to sample names. Negative number to take from start of path.
            | **fn_as_s_names**: bool -- Use the log filename as the sample name 
            | **replacenames**: object -- TSV file to rename sample names during report generation (PATH)
            | **comment**: str -- Custom comment, will be printed at the top of the report.
            | **template**: [default|default_dev|geo|simple|sections] Report template to use.
            | **samplesnames**: object -- File containing alternative sample names
            | **samplefilters**: object -- TSV file containing show/hide patterns for the report (PATH)
            | **custom_css**: object -- Custom CSS file to add to the final report (PATH)
            | **flat**: bool -- Use only flat plots (static images)
            | **datadir**: bool -- Force the parsed data directory to be created.
            | **nodatadir**: bool -- Prevent the parsed data directory from being created.
            | **dataformat**: [tsv|yaml|json] -- Output parsed data in a different format. Default: tsv
            | **noreport**: bool -- Do not generate a report, only export data and plots
            | **pdf**: bool -- Creates PDF report with 'simple' template. Requires Pandoc to be installed
            | **quiet**: bool -- Only show log warnings
            | **lint**: bool -- Use strict linting (validation) to help code development
            | **profileruntime**: bool -- Add analysis of how long MultiQC takes to run to the report
            | **nomegaqc**: bool -- Don't upload generated report to MegaQC, even if MegaQC options are found
            | **noansi**: bool -- Disable coloured log output

    Returns
    -------
    output: _dict
    
        | **outQC**: folder -- create report in the specified output directory
        | **done**: object -- gets the shell in a text file


    Examples
    --------
    >>> from snakemake.io import Params
    >>> from wrappers.wrapper_system import _StubParams
    >>> output = _StubParams("output")

    All optional values should be present

    >>> params = _StubParams("params")
    >>> create_multiqc(output, params)
    'mkdir -p {output.outQC} && multiqc [-f] [-c <{{params.config}}, default=None>] [--cl-config <{{params.clconfig}}, default=None>] [-n <{{params.filename}}, default=None>] --outdir {output.outQC} [-x <{{params.ignore}}, default=None>] [--ignore-samples <{{params.ignoresamples}}, default=None>] [--ignore-symlinks] [-l <{{params.filelist}}, default=None>] [-m <{{params.module}}, default=None>] [-e <{{params.exclude}}, default=None>] [--tag <{{params.tag}}, default=None>] [--view-tags] [-d <{{params.dirs}}, default=None>] [-dd <{{params.dirsdepth}}, default=None>] [-s] [--fn_as_s_name] [--replace-names <{{params.replacenames}}, default=None>] [-i <{{params.title}}, default=None>] [-b <{{params.comment}}, default=None>] [-t <{{params.template}}, default=None>] [--sample-names <{{params.samplesnames}}, default=None>] [--sample-filters <{{params.samplefilters}}, default=None>] [--custom-css-file <{{params.custom_css}}, default=None>] [-fp] [-ip] [-p] [--data-dir] [--no-data-dir] [-k <{{params.dataformat}}, default=None>] [-z] [--no-report] [--pdf] [-v] [-q] [--lint] [--profile-runtime] [--no-megaqc-upload] [--no-ansi] 2>{output.done}'

    Specific case

    >>> params = Params(fromdict={"config": "config_tools.yml", "module": "'fastqc'", "title": "'Report of FQC on raw data'", "filename": "'report_fastqc.html'", "dirs": "test/results/", "fullnames":True, "viewtags": True, "zipdatadir":True, "export": True, "verbose":True, "interactive": True})
    >>> create_multiqc(output, params)
    "multiqc -c config_tools.yml -n 'report_fastqc.html' --outdir {output.outQC} -m 'fastqc' --view-tags -d test/results/ -s -i 'Report of FQC on raw data' -ip -p -z -v 2>{output.done}"
    
    """

    create_folder = mappy(params, 'create_folder',
                          (True, False), default=False)
    
    return join_str(
        (f"mkdir -p {output.outQC} &&" if create_folder else ""),
        "multiqc",
        ## Main options
        optional(params, "force ", "-f"),
        optional(params, "config", "-c {}"),
        optional(params, "clconfig", "--cl-config {}"),
        optional(params, "filename", "-n {}"),
        f"--outdir {output.outQC}",
        optional(params, "ignore", "-x {}"),
        optional(params, "ignoresamples", "--ignore-samples {}"),
        optional(params, "ignoresymlinks", "--ignore-symlinks"),
        optional(params, "filelist", "-l {}"),
        ## Choose modules to run
        optional(params, "module", "-m {}"),
        optional(params, "exclude", "-e {}"),
        optional(params, "tag", "--tag {}"),
        optional(params, "viewtags", "--view-tags"),
        ## Sample hangling
        optional(params, "dirs", "-d {}"),
        optional(params, "dirsdepth", "-dd {}"),
        optional(params, "fullnames", "-s"),
        optional(params, "fn_as_s_names", "--fn_as_s_name"),
        optional(params, "replacenames", "--replace-names {}"),
        ## Report customisation
        optional(params, "title", "-i {}"),
        optional(params, "comment", "-b {}"),
        optional(params, "template", "-t {}"),
        optional(params, "samplesnames", "--sample-names {}"),
        optional(params, "samplefilters", "--sample-filters {}"),
        optional(params, "custom_css", "--custom-css-file {}"),
        ## Output files
        optional(params, "flat", "-fp"),
        optional(params, "interactive", "-ip"),
        optional(params, "export", "-p"),
        optional(params, "datadir", "--data-dir"),
        optional(params, "nodatadir", "--no-data-dir"),
        optional(params, "dataformat", "-k {}"),
        optional(params, "zipdatadir", "-z"),
        optional(params, "noreport", "--no-report"),
        optional(params, "pdf", "--pdf"),
        ## MultiQC behaviour
        optional(params, "verbose", "-v"),
        optional(params, "quiet", "-q"),
        optional(params, "lint", "--lint"),
        optional(params, "profileruntime", "--profile-runtime"),
        optional(params, "nomegaqc", "--no-megaqc-upload"),
        optional(params, "noansi", "--no-ansi"),
        ## log-error output
        f"2>{output.done}"
    )


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
