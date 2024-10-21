Metafiles
=========

There are two files that need to be filled in to ensure that the information about the samples are correct:

    * the **Metadata** file: sample name, fraction type, replicate number (and for kinetic time).
    * the **Metacontrast** file: indicates the desired single and double comparisons.

These files contain all the information needed for the statistical and enrichment parts of the analysis.

Examples are provided in the *ressources/* folder.


Metadata
--------

This file is required for the first part of the analysis, rule *Diff_gene_expression* in *05_stats.smk*, where the simple comparisons are made based on sample identification.

This file should contains 4 columns:

============   ========  =========  =========
sample_name    fraction  treatment  replicate
============   ========  =========  =========
F1_Treat1_R1   F1        Treat1     R1
F1_Treat1_R2   F1        Treat1     R2
F1_Treat1_R3   F1        Treat1     R3
F1_Treat2_R1   F1        Treat2     R1
F1_Treat2_R2   F1        Treat2     R2
...            ...       ...        ...
============   ========  =========  =========

Nomenclature:

- F1 refers to the Cytoplasmic or Total fraction
- F2 to the Polysomal fraction.

In the case of kinetic analysis, a "time" column is added to the metadata file. See `example <https://github.com/Translatome/TRAIN/blob/main/ressources/metadata_kinetic.txt>`_.


Metacontrast
------------

This file is required for the first and second part of the anlaysis, rules *Diff_gene_expression* and *Categorization* in *05_stats.smk*, where the simple comparisons are compared in a double comparison.

This file should contains 4 columns:

========  ==========  ============   ==========
fraction  treatment1  fraction_bis   treatment2
========  ==========  ============   ==========
F1        Treat1      F1             Treat2
F2        Treat1      F2             Treat2
========  ==========  ============   ==========

The direction of comparison of the fractions is fixed, first cytoplasmic/total, then polysomal.

In the case of a treatment/control comparison, "treatment1" corresponds to the treatment and "treatment2" corresponds to the control.

If the comparison is between two treatments (without a control), the user must select the treatment to be considered as the control in order to assign it to the "treatment2" column in the table.

By default, the control should always be assigned to the "treatment2" column and the comparison will always be made by "treatment1" vs. "treatment2".

You can add as many contrasts as you like within the required format. The separator is a tab and the file format is ".txt".

In the case of a kinetic analysis, two columns are required for the time to be compared, "time" and "time_bis". Time_bis can be the same time or a different time according to your comparisons. See `example <https://github.com/Translatome/TRAIN/blob/main/ressources/metacontrast_kinetic.txt>`_.
