Additional programs
^^^^^^^^^^^^^^^^^^^

When you ``conda install`` or ``pip install`` flair, the following helper programs will be in your $PATH.
The options of each are generated from the program itself, so they are always what the program accepts.

annotate_aaseq_with_uniprot
===========================

.. include:: cli/annotate_aaseq_with_uniprot.rst

This is for aaseq predictions from FLAIR predictProductivity.

To obtain a UniProt reference file:
 - First go to: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/
     - README file has mapping of your organism name to the proteome directory
 - Navigate to your organism’s directory
     - Human files are in: https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/Eukaryota/UP000005640/
 - Download the fasta and additional.fasta files
     - Human: download UP000005640_9606.fasta.gz (primary isoforms) and UP000005640_9606_additional.fasta.gz (alternative isoforms)
 - Cat together zipped files
     - Ex: cat UP000005640_9606.fasta.gz UP000005640_9606_additional.fasta.gz > UniProt_human_withadditional_100925.fasta.gz

**Output**

A fasta file similar to the input but with aaidXX changed to UniProt ID when possible

bed_to_gtf
==========

.. include:: cli/bed_to_gtf.rst

Convert isoforms in BED12 format to a GTF, written to standard output.

diff_iso_usage
==============

.. include:: cli/diff_iso_usage.rst

Calculates the usage of each isoform as a fraction of the total expression of the gene and compares this between samples.

Identifies and calculates significance of alternative isoform usage between two
samples using Fisher's exact tests. The more differentially used the isoforms are
between the first and second condition, the lower the p-value.

Output file format columns are as follows, where the two sample names are the ones
given on the command line, so the direction of ``delta_PSI`` can be read from the file:

 - gene name
 - isoform name
 - p-value
 - this isoform's count in each of the two samples
 - the gene's other isoforms' counts in each of the two samples
 - each sample's PSI, and the difference between them

diffsplice_fishers_exact
========================

.. include:: cli/diffsplice_fishers_exact.rst

Identifies and calculates the significance of alternative splicing
events between two samples without replicates using Fisher's exact
tests, from a flair diffsplice ``tsv`` of alternative splicing calls for
one splicing event type.

**Output**

The output file contains the original columns with an additional column
containing the p-values appended.

fasta_seq_lengths
=================

.. include:: cli/fasta_seq_lengths.rst

flair_partition
===============

.. include:: cli/flair_partition.rst

Defines non-overlapping regions across all of the input files, which is how FLAIR
splits work across threads.

gtf_to_bed
==========

.. include:: cli/gtf_to_bed.rst

Convert a GTF to a BED12 file.

identify_gene_isoform
=====================

.. include:: cli/identify_gene_isoform.rst

Identifies the most likely gene id associated with each isoform and renames the
isoform.

junctions_from_sam
==================

.. include:: cli/junctions_from_sam.rst

mark_intron_retention
=====================

.. include:: cli/mark_intron_retention.rst

Assumes the bed has the correct strand information.

**Outputs**

 - the input BED with a ``retains_intron`` column added, 1 when another isoform of the
   locus has an exon where this one has an intron, 0 otherwise.  A header line names the
   columns, so the file is a TSV rather than a BED.
 - a TSV of the retained introns, with columns ``chrom``, ``start``, ``end`` and
   ``strand``.  An intron appears once however many isoforms retain it.

plot_isoform_usage
==================

.. include:: cli/plot_isoform_usage.rst

Visualization script for FLAIR isoform structures and the percent usage
of each isoform in each sample for a given gene. If you supply the
isoforms.bed file from running ``predictProductivity``, then isoforms
will be filled according to the predicted productivity (solid for
``PRO``, hatched for ``PTC``, faded for ``NGO`` or ``NST``). The gene
name supplied should correspond to a gene name in your isoform file and
counts file.

The script will produce two images, one of the isoform models and another of the usage proportions.

The most highly expressed isoforms across all the samples will be plotted.

The minor isoforms are aggregated into a gray bar. You can toggle min_reads or
color_palette to plot more isoforms.

**Outputs**

 - gene_name_isoforms.png of isoform structures
 - gene_name_usage.png of isoform usage by sample

For example:

.. figure:: img/toy_diu_isoforms.svg

.. figure:: img/toy_diu_usage.svg

spliceevents_outliers
=====================

.. include:: cli/spliceevents_outliers.rst

Finds outlier samples of the splicing events from flair spliceevents: samples
whose PSI of an event is far from the median of all samples', at least 0.1 from
it and more than 3 deviations (half the interquartile range), of events where at
least 5 samples have ``--event_support`` reads.  It reads
``OUTPUT.spliceevents.counts.tsv`` and ``OUTPUT.spliceevents.totals.tsv``, and
calls outliers of either PSI: ``--psi junc``, of the event's reads, which is more
sensitive, or ``--psi tot``, of the reads covering the event's locus, which finds
strong overall events.

**Outputs**

 - ``OUTPUT.spliceevents.outliers.tsv``, every outlier
 - ``OUTPUT.spliceevents.outliers.filtered.tsv``, the outliers more than 6
   deviations from the median that a simpler event of the same junctions doesn't
   explain as well: a complex exon skip by an exon skip, and an alternative
   splice site by an exon skip or a complex exon skip

.. note::

   Starts and ends (tss and tts) are left out of the filtered outliers on
   purpose.  Their events, from read ends, give far too many false positive
   outliers for now.  They are still in ``OUTPUT.spliceevents.outliers.tsv``.
