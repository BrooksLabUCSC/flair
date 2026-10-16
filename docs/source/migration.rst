.. _migration-label:

Migrating from FLAIR 3.0
========================

What a command line written for FLAIR 3.0 has to change. Options and outputs that are
not listed here are unchanged.

Everywhere
----------

- Option names are snake_case. Run-together names such as ``--maxintronlen`` and
  camelCase names such as ``--conditionA`` are gone; each one is listed under its
  command below.
- ``-z``, ``--compress`` gzips the text files a run writes, adding ``.gz`` to their
  names, and bgzips FASTA so it can be indexed. It is a global option, so it comes
  before the subcommand: ``flair -z transcriptome ...``. Compressed input is read
  whether or not it is given.
- The TSVs FLAIR writes are named ``.tsv``, not ``.txt``: ``flair.read.map.txt`` is
  ``flair.isoform.read.map.tsv``.
- Genes, isoforms and protein sequences carry FLAIR identifiers, and the gene is a
  column of its own rather than part of the isoform name. See :ref:`ids-label`.
- ``--quiet`` is gone. Logging is set by the global ``--log-level``, ``--log-stderr``,
  ``--log-conf`` and ``--log-debug``.

Commands that moved
-------------------

- ``flair correct`` and ``flair collapse`` are gone. ``flair transcriptome`` does both,
  starting from a genome aligned BAM, with the junction evidence given by ``--gtf``,
  ``--junction_bed`` or ``--junction_star``.
- ``flair variants`` is now ``flair alleles``, with a different interface.
- ``flair isoalleles`` and ``flair variantquant`` are new.

flair align
-----------

**Renamed**

- ``--nvrna`` is ``--native_rna``
- ``--filtertype`` is ``--filter_type``
- ``--maxintronlen`` is ``--max_intron_len``
- ``--minfragmentsize`` is ``--min_fragment_size``

**Gone**

- ``-f``, ``--gtf``
- ``--remove_internal_priming``, ``--intprimingthreshold``, ``--intprimingfracAs``
- ``--remove_singleexon``

flair transcriptome
-------------------

**Renamed**

- ``--genomealignedbam`` is ``-b``, ``--genome_aligned_bam``
- ``--parallelmode`` is ``--parallel_mode``
- ``--keep_sup`` is ``--keep_supplementary``
- ``--noaligntoannot`` is ``--no_align_to_annot``
- ``--junction_tab`` is ``--junction_star``; the old name still works and warns

**Defaults flipped**

- ``--stringent`` is now the default; ``--no_stringent`` turns it off
- ``--check_splice`` is now the default; ``--no_check_splice`` turns it off

**Split**

- ``-s``, ``--support`` became three options: ``--sjc_support`` for spliced isoforms
  (default 1), ``--single_exon_support`` for single exon isoforms (default 3), and
  ``--frac_support``, the fraction of the gene locus an isoform must make up
  (default 0.05)

**Gone**

- ``-n``, ``--no_redundant``
- ``--predictCDS``

**New**

- ``--sample_name`` is required, and is recorded in the output files
- ``--trust_strand``, ``--trust_ends``, ``--normalize_ends``
- ``--generate_map``, ``--fusion_breakpoints``
- ``--allow_paralogs``, which is not implemented and is an error rather than a no-op

**Outputs**

- ``flair.isoforms.bed``, ``flair.isoforms.gtf``, ``flair.isoforms.fa`` as before
- ``flair.read.map.txt`` is ``flair.isoform.read.map.tsv``, written with
  ``--generate_map``

flair quantify
--------------

**Renamed**

- ``-r``, ``--reads_manifest`` is ``--manifest``
- ``--isoformbed`` is ``--isoform_bed``

**Changed**

- The manifest names each sample's reads as a sorted, indexed BAM aligned to the
  genome, where it named a reads file before. Its first three columns, sample id,
  condition and batch, are unchanged.
- The isoforms are given as ``--isoform_bed`` with ``-g``, ``--genome``. The isoform
  FASTA option ``-i``, ``--isoforms`` is gone.

**Gone**

- ``--output_bam``, ``--sample_id_only``, ``--temp_dir``
- ``--quality``, ``--stringent``, ``--check_splice``

**New**

- ``--normalize_ends``

**Outputs**

- ``<output>.counts.tsv`` has a ``gene_id`` and an ``isoform_id`` column, in place of
  the single ``ids`` column that joined the two.
- Each sample's condition and batch are in ``<output>.sample_info.tsv`` rather than
  packed into the counts matrix column names, which now hold the sample id alone.
  ``flair diffexp`` and ``flair diffsplice`` read that file.
- Read maps are ``<output>.<sample>.read.map.tsv``.

flair combine
-------------

**Renamed**

- ``--output_prefix`` is ``-o``, ``--output``
- ``--endwindow`` is ``--end_window``
- ``--minpercentusage`` is ``--min_percent_usage``

**Default flipped**

- Single exon isoforms are kept. ``--include_se`` is gone; ``--remove_single_exon``
  drops them.

**Gone**

- ``-c``, ``--convert_gtf``: the GTF is always written
- ``-f``, ``--filter``

**New**

- ``--prefixes`` and ``--isoform_beds`` name the transcriptomes to combine without a
  manifest file; they and ``--manifest`` are mutually exclusive
- ``-g``, ``--genome``
- ``--max_ends``, ``--min_reads``

**Changed**

- The manifest is one isoform BED path per line. The sample name, type, FASTA and read
  map columns are gone.

flair fusion
------------

**Renamed**

- ``--genomechimbam`` is ``-b``, ``--genome_aligned_bam``
- ``--maxloci`` is ``--max_loci``
- ``--minfragmentsize`` is ``--min_fragment_size``
- ``-s``, ``--support`` is ``--min_support``

**Defaults flipped**

- ``--stringent`` and ``--check_splice`` are the default; ``--no_stringent`` and
  ``--no_check_splice`` turn them off

**Gone**

- ``-r``, ``--reads``
- ``--transcriptchimbam``

**New**

- ``--sample_name`` is required
- ``--max_dist_to_tss``, ``--min_dist_between_breakpoints``
- ``--keep_intermediate``, ``--allow_paralogs``

flair diffexp
-------------

**Renamed**

- ``--out_dir`` is ``-o``, ``--output``
- ``-of``, ``--out_dir_force`` is ``--overwrite_output``
- ``-e``, ``--exp_thresh`` is ``--min_expression``

**New**

- ``--condition_a`` and ``--condition_b`` name the two conditions to compare. Without
  them the two are taken in sorted order rather than from the column order of the
  counts matrix.

**Outputs**

- ``isoforms_sig_exp_change_norm_by_gene.tsv``, isoforms whose usage changes, from a
  t-test on counts normalized by gene.

flair diffsplice
----------------

**Renamed**

- ``-i``, ``--isoforms`` is ``--isoform_bed``
- ``--outDir`` is ``-o``, ``--output``
- ``-of``, ``--out_dir_force`` is ``--overwrite_output``
- ``--conditionA`` and ``--conditionB`` are ``--condition_a`` and ``--condition_b``
- ``--drim1`` is ``--min_samps_gene_expr``
- ``--drim2`` is ``--min_samps_feature_expr``
- ``--drim3`` is ``--min_gene_expr``
- ``--drim4`` is ``--min_feature_expr``

**Outputs**

- Event names end in the strand, so ``inclusion_chr1:500-650`` is
  ``inclusion_chr1:500-650:+``. The same goes for the ``coordinate`` column. Without
  it, an event called on both strands at one locus appeared twice under one name.
- All four ``drimseq_<event>_A_v_B.tsv`` are always written. An event type that called
  nothing, or whose events DRIMSeq found nothing significant in, leaves a file of only
  its header.

flair alleles
-------------

``flair variants`` under a new name and a new interface: it takes the tumor and normal
BAMs and VCFs rather than a manifest and an isoform FASTA.

**Renamed**

- ``-o``, ``--output_prefix`` is ``-o``, ``--output``

**Gone**

- ``-m``, ``--manifest``
- ``-i``, ``--isoforms``
- ``-b``, ``--bedisoforms``

**New**

- ``--tumor_bam``, ``--normal_bam``, ``--vcf``, ``--normal_vcf``
- ``--isoform_bed``, ``--frac_support``, ``--read_support``
- ``--annotate_bams``, ``--generate_map``

flair isoalleles and flair variantquant
---------------------------------------

Both are new.

- ``flair isoalleles`` builds allele labelled isoforms from the output of
  ``flair alleles``; see :ref:`isoalleles-label`.
- ``flair variantquant`` counts the reads carrying each variant, per sample; see
  :ref:`variantquant-label`.
