.. _file-formats-label:

FLAIR file formats
==================

Every file FLAIR reads and writes is one of five things: a BED, a GTF, a FASTA, a
BAM or VCF of the standard form, or a tab separated table.  This page gives the
columns of each.  :ref:`modules-label` gives the files each subcommand writes and
what they are for, and :ref:`ids-label` gives the identifiers in them.

Conventions
-----------

- A table's first line is its column names.  A reader of a table whose columns are
  fixed checks that line, and reports what it found against what it expected.
- A column of several values is comma separated.  The exception is ``source_ids``
  of the ``flair combine`` maps, separated by ``;`` and a space, since its values
  contain commas.
- An unknown or inapplicable value is an empty field, never 0 or ``NA``.  A PSI is
  empty when there were too few reads to take a fraction from, which is not the
  same as a PSI of 0.
- Any text file FLAIR reads may be gzipped, and ``--compress`` gzips what it
  writes; see :ref:`compression-label`.
- BED coordinates are zero based and half open, as BED is defined.  GTF and VCF
  coordinates are one based and inclusive, as those formats are defined.
- The event names and coordinate columns of ``flair spliceevents`` are BED
  coordinates written 5' to 3', so the first number is the larger one on the minus
  strand: ``chr12:25206932-25205246(-)``.  ``flair diffsplice`` writes a junction
  the same way.
- A sample column is named by the sample id from the manifest, so the names are
  arbitrary and the columns are read by position after the named ones.

Input files
-----------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - File
     - What it must be
   * - genome
     - a FASTA, read by position through its ``.fai`` index.  A compressed genome
       must be bgzip compressed and indexed
   * - annotation
     - a GTF with ``gene_id`` and ``transcript_id`` attributes, given with ``-f``
   * - reads
     - a BAM of the reads aligned to the genome, coordinate sorted, with its index
       beside it; see :ref:`alignment-label`
   * - ``--junction_bed``
     - splice junctions as BED6 to BED9 with the supporting read count in the score
       column, as intron-prospector writes
   * - ``--junction_star``
     - splice junctions as STAR's ``SJ.out.tab``.  The same evidence as
       ``--junction_bed`` in another form, not an additional kind
   * - ``--region_bed``
     - a BED of the regions to work in; only the first three columns are read
   * - ``--isoform_bed``
     - a FLAIR isoform BED, the 25 columns below
   * - ``--manifest``
     - a sample per line, tab separated; see below
   * - phased VCF
     - a VCF with genotypes, of any caller, as ``flair isoalleles`` takes; its
       reads are a haplotagged BAM

The manifest
~~~~~~~~~~~~

``flair quantify``, ``flair combine`` and ``flair spliceevents`` take a manifest of
four fields, the sample id, its condition, its batch and its BAM::

   sample1	ctl	batch1	mydata/sample1.bam
   sample2	test	batch1	mydata/sample2.bam

``flair spliceevents`` also takes the two field form, a sample id and its BAM, and
then writes no conditions beside its counts.  Every line has the same number of
fields; a manifest of both forms is an error.  Lines starting with ``#`` and blank
lines are skipped, and a relative path is read against the directory FLAIR runs in.

The sample id, the condition and the batch are written into a counts matrix's
columns, so none of them may contain an underscore: a matrix of a FLAIR before 3.1
named each column ``sample_condition_batch``, and those names are split on the
underscore when such a matrix is read.

The FLAIR isoform BED
---------------------

``flair transcriptome``, ``flair combine``, ``flair fusion`` and ``flair
isoalleles`` write their isoforms as a BED12 with 13 columns of their own.  It is
read by ``flair quantify``, ``flair diffsplice``, ``flair spliceevents``,
``flair isoalleles`` and ``bed_to_gtf``.

.. list-table::
   :header-rows: 1
   :widths: 8 22 70

   * - Column
     - Name
     - Contents
   * - 1 to 12
     - the BED12 columns
     - the isoform's position and exons.  ``name`` is its FLT id, ``score`` is
       ``read_support`` clamped to the 0 to 1000 a genome browser shades, and
       ``itemRgb`` colors an isoform by its productivity, or red for a fusion and
       purple for a readthrough
   * - 13
     - ``gene_id``
     - the FLG id of its gene
   * - 14
     - ``ref_transcript_id``
     - the transcript id of the annotation, where the isoform matches one
   * - 15
     - ``ref_gene_mappings``
     - the gene ids of the annotation it overlaps, or the novel region FLAIR named
   * - 16
     - ``read_support``
     - the full length reads assigned to it
   * - 17
     - ``frac_support``
     - its share of its gene's reads, spliced and single exon together, to four
       decimal places
   * - 18
     - ``productivity``
     - ``PRO`` productive, ``PTC`` premature termination codon, ``NGO`` no start
       codon, ``NST`` a start codon with no stop
   * - 19
     - ``transcript_class``
     - ``basic``, ``readthrough`` for adjacent genes, or ``fusion``
   * - 20
     - ``fused_genes``
     - the FLG ids a fusion joins, in transcript order
   * - 21
     - ``pos_in_fusion``
     - which fragment of the fusion this row is, counting from 0
   * - 22
     - ``samples``
     - the samples the isoform was found in, which ``flair combine`` fills
   * - 23
     - ``source_isoform``
     - for an allele, the isoform it was split from
   * - 24
     - ``allele_group``
     - for an allele, its group, as ``H2-S1``
   * - 25
     - ``aaseq_id``
     - the FLP id of its protein, empty where there is no coding prediction

A fusion isoform is written as one row per fragment, each with the same name and a
``pos_in_fusion`` of its own.

The same columns are written as a bigBed, ``.bb``, whose field definitions are
``src/flair/flair_bed.as``.  ``bed_to_gtf`` converts the BED to a GTF, and
``gtf_to_bed`` converts an annotation GTF to a BED; see :ref:`scripts-label`.

Counts matrices
---------------

A counts matrix is the table the differential tests read, and ``flair quantify``,
``flair isoalleles`` and ``flair spliceevents`` all write one.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Contents
   * - ``gene_id``
     - what the rows are grouped by and tested within: the gene for ``flair
       quantify``, the isoform for ``flair isoalleles``, the event for ``flair
       spliceevents``
   * - ``isoform_id``
     - what the row counts within that group: an isoform, an allele
       ``isoform|group``, or an alternative of an event
   * - one per sample
     - the count in that sample, named by its sample id

``<output>.isoforms.cpm.tsv`` has the same columns with the counts as counts per
million.  A matrix of an earlier FLAIR joined the two ids into one ``isoform_gene``
column; that form is rejected rather than split, since no split is right for a gene
id containing an underscore.

Beside the matrix, ``<output>.isoforms.sample_info.tsv`` names the columns:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Contents
   * - ``sample_id``
     - names one counts column, in the order of those columns
   * - ``condition``
     - the condition the manifest assigned it
   * - ``batch``
     - the batch the manifest assigned it

Read maps and id maps
---------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Columns
   * - a read map, written with ``--generate_map``
     - ``isoform_id``, the isoform or ``isoform|group``, and ``reads``, its read
       names
   * - ``<output>.isoforms.counts.tsv`` of ``flair transcriptome``
     - ``isoform_id`` and ``full_length_reads``, that sample's own counts
   * - ``<output>.combined.isoforms.map.tsv``, ``.combined.genes.map.tsv``
     - ``new_id`` and ``source_ids``, each written ``samples:id``.  ``flair
       combine`` numbers afresh, so these are the only link to the per-sample ids
   * - ``<output>.aaseq.tsv``
     - ``aaseq_id`` and ``aaseq``, each protein sequence once, however many
       isoforms translate to it

The isoform FASTA holds one record per isoform, named by its FLT id, or
``FLT00000001|H1`` for an allele of ``flair isoalleles``.

flair spliceevents tables
-------------------------

Three files share one set of columns and differ only in what the sample columns
contain: ``.spliceevents.counts.tsv``, the reads of the alternative,
``.spliceevents.psi_junc.tsv``, its share of the event's reads, and
``.spliceevents.psi_tot.tsv``, its share of the reads covering the event's locus.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Contents
   * - ``event``
     - the event, ``<type>|<chrom>:<start>-<end>(<strand>)|<gene>``
   * - ``alternative``
     - the alternative the row counts, or ``other`` for the event's reads of none
       of them
   * - ``event_type``
     - the kind of event
   * - ``gene``, ``gene_name``
     - the gene the event was found in, and its name
   * - ``additional_genes``, ``additional_gene_names``
     - other genes sharing the event's junctions
   * - ``junctions_included``, ``junctions_excluded``
     - the junctions of the alternative, and the junctions it skips
   * - ``outer_junctions``
     - the junctions flanking the event
   * - ``exons``
     - the exons the alternative includes, or the position of a start or end event
   * - ``novel_junctions``
     - its junctions that are not annotated
   * - one per sample
     - the count or the PSI

``.spliceevents.totals.tsv`` has the gene columns and a ``total`` column, two rows
per event: ``event_reads``, which PSIjunc divides by, and ``locus_reads``, which
PSItot divides by.  ``.spliceevents.bed`` has a row per alternative, named
``event|alternative``, and ``.spliceevents.regions.bed`` is the three column BED of
the regions the events were called in.

``spliceevents_outliers`` writes the samples whose PSI is far from the median of
the others; see :ref:`spliceevents-outliers-label`.

flair isoalleles tables
-----------------------

``flair isoalleles`` writes ``<output>.isoalleles.<table>.tsv``:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Table
     - Columns
   * - ``sequences``
     - ``gene_id``, ``isoform_id``, ``allele_group``, ``variants``,
       ``terminal_exon_somatic_added``, ``productivity``, ``productivity_note``,
       ``aaseq_id``, ``protein_change``, ``somatic_protein_change``, ``cds``,
       ``phase_sets``
   * - ``variants``
     - ``gene_id``, ``isoform_ids``, ``variant``, ``kind``, ``haplotype_or_sets``,
       ``detail``
   * - ``somatic_sets``
     - ``gene_id``, ``somatic_set``, ``variants``, ``reads``
   * - ``excluded``
     - ``gene_id``, ``variant``, ``normal_reads``, ``isoform_ids``, ``reason``
   * - ``warnings``
     - ``gene_id``, ``isoform_id``, ``warning``, ``detail``
   * - ``terminal_sites``
     - ``gene_id``, ``isoform_id``, ``variant``, ``allele_group``,
       ``reads_reaching``, ``reads_carrying``, ``decision``
   * - ``genes``
     - ``gene_id``, ``sample``, ``spliced_reads``, ``assigned_spliced_reads``,
       ``assigned_share``, ``assigned_reads``

A variant is written ``chrom:pos:ref>alt``, the position one based as a VCF gives
it.  A CDS is written ``chrom:start-end``, zero based and half open as a BED gives
it.  A group is a haplotype, ``H1``, ``H2`` or ``untagged``, and in a tumor and
normal pair a somatic set with it, ``H2-S1``, or ``H2-S?`` where the set is
unknown.

flair diffexp and flair diffsplice tables
-----------------------------------------

The tables of these two are written by DESeq2, DRIMSeq and the FLAIR t-test, so
their columns are those tools', and :ref:`diffexp-label` and :ref:`diffsplice-label`
give them with the files they belong to.  Both read a counts matrix and the sample
info beside it, which is how a comparison gets its conditions.
