Modules
^^^^^^^

Each FLAIR analysis is a sub-command of ``flair``.  See :ref:`running-label`
for the options that come before the sub-command and for how compression
works.

If you want to compare multiple samples, there are two primary ways of doing this:
 - Combine the fastq or fasta reads of all samples and run FLAIR transcriptome on all samples together (will generate the most comprehensive transcriptome)
 - Run FLAIR transcriptome on each sample separately (better for large sets of samples) and
   use FLAIR combine to merge results

.. figure:: img/flair_workflow_compartmentalized.svg


.. figure:: img/flair_workflow_compartmentalized.svg


.. _alignment-label:

Aligning reads to the genome
============================

FLAIR's modules take reads aligned to the genome as a sorted, indexed bam.  

Please see the protocols page with more details on pros and cons of different alignment options, but here is an explanation of minimap options if you chose to use minimap and want to customize your run.

The simplest protocol is to align with `minimap2 <https://github.com/lh3/minimap2>`__ and sort with samtools, in one line:

.. code:: sh

   minimap2 -ax splice -s 40 -G 350k -t 25 --MD --secondary=no genome.fa sample.fastq | samtools sort -o sample.genomealigned.bam - && samtools index sample.genomealigned.bam

The reads can be FASTA or FASTQ, gzipped or not. We recommend input with stranded reads without polyA tails trimmed, but unstranded and trimmed reads will work as well.

Options to consider:

``-ax splice``
    Spliced alignment of long reads, with SAM output. You can also use splice:hq with high-quality reads, although FLAIR has been tested with splice (optimized for lower-quality long-reads).
``-s 40``
    The minimum score of an alignment, roughly its length in matching bases: shorter
    pieces of a read are left unaligned.  For fusion detection, ``-s 40`` keeps shorter
    partner pieces, as a short 5' exon, and matches flair fusion's
    ``--min_fragment_size`` default.
``-G 350k``
    The longest intron.  Longer allows novel isoforms with long introns.  A fusion from
    a deletion between genes closer than this aligns as one spliced read, which flair
    fusion does not report: to find those, align with a smaller ``-G`` and lower flair
    fusion's ``--min_dist_between_breakpoints`` to match.
``--MD``
    Adds the MD tags that flair alleles needs to find the reads' mismatches.
``--secondary=no``
    Leaves out secondary alignments.  FLAIR uses each read's primary alignment, and
    flair fusion its supplementary ones.
``-t 20``
    Threads.
``-uf``
    For stranded reads, oriented 5' to 3' as their RNA: splice junctions are taken on
    the read's strand.  flair fusion needs stranded reads.
``--junc-bed annotation.bed``
    We do not recommend this option, as minimap2 tends to overcorrect to the submitted junctions

If you're using human sequences, the best reference genome is
`GCA_000001405.15_GRCh38_no_alt_analysis_set <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz>`__ as described in this
`helpful blog post by Heng Li <https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use>`__.
If your input sequences are Oxford nanopore reads, please use `Pychopper <https://github.com/epi2me-labs/pychopper>`__ before aligning them.


.. _transcriptome-label:

flair transcriptome
===================

This module generates a transcriptome of high confidence isoforms (bed, gtf, and fasta files) directly from a bam file of
reads aligned to the genome, as :ref:`aligning reads to the genome <alignment-label>` describes.
For when to count splicing events from reads with flair spliceevents instead, see
:ref:`choosing-analysis-label`.


**Outputs**

 - ``flair.isoforms.bed``
 - ``flair.isoforms.gtf``
 - ``flair.isoforms.fa``
 - ``flair.isoform.read.map.tsv``

Options
-------

.. include:: cli/transcriptome.rst


.. _fusions-label:

flair fusion
============

This identifies gene fusions and generates a fusion transcriptome. 
To incorporate this fusion transcriptome in downstream analysis, 
use flair combine to merge it with normal isoforms. 

Its input is a bam of all the reads aligned to the genome, as :ref:`aligning reads to the
genome <alignment-label>` describes, with ``-s 40`` to keep short partner pieces.

flair fusion requires stranded reads, each oriented 5' to 3' as its RNA.
The 5' and 3' genes of a fusion come from the order of a read's pieces, so a
reversed read gives the fusion backwards, which is then rejected as having its
5' gene in antisense.  Data that is poorly stranded, as some older PacBio
data, will miss fusions.

Each piece of a final fusion isoform gets its gene as flair transcriptome gives an
isoform its gene: by the annotated junctions and splice sites it shares with genes,
counting its splice sites at the fusion junctions, so a piece of only a gene's first
exon has its donor, and never by where reads end.  A piece with neither keeps the gene
of the preliminary fusion.  An isoform two of whose pieces are of one gene is left
out: flair fusion does not call rearrangements within a gene.

The fusion isoforms are built by flair transcriptome from the fusion reads, which are
often 5' truncated.  By default, as transcriptome's ``--filter nosubset``, an isoform
that is a part of another is dropped, and a truncated read counts for no isoform, so a
fusion can have many fewer final reads than preliminary ones: of a KIF5B--RET fusion's
123 reads, mostly 5' truncated, 6 were in its isoforms.  For the most recall, as for
reporting fusions, give ``--filter comprehensive``, every spliced subset isoform kept,
and ``--frac_support 0``, every isoform however few of the fusion's reads it has: that
fusion's isoforms then had 111 of its reads.

flair fusion is configured for human genomes: its paralog table and its filters of
immunoglobulin genes and mitochondrial fusions are human.  It has been tested with
GENCODE, RefSeq and Ensembl references, and recognises chromosomes as UCSC and GENCODE
(chr12), RefSeq (NC_000012.12) and Ensembl (12) name them, leaving out mitochondria and
unplaced scaffolds.

**Output**

sample.fusions.tsv
    Table of the fusions, most final reads first: the genes 5' to 3', those of its
    final isoforms, the reads of the preliminary fusion and of its final isoforms, the number of isoforms, the best
    supported breakpoint and isoform, and their reads.  A breakpoint is each junction
    between consecutive loci, in fusion order and comma-separated, as the last base of
    the 5' locus and the first base of the 3' one, 1-based, chrom:pos(strand)--chrom:pos(strand)
sample.fusions.isoforms.bed
    Bed file of fusion transcriptome, with a line for each locus of a fusion isoform and
    its position in the fusion in the pos_in_fusion column
sample.fusions.isoforms.gtf, sample.fusions.isoforms.bb
    The fusion transcriptome as GTF and bigBed
sample.fusions.isoforms.fa
    Fasta file of fusion transcriptome
sample.fusions.isoform.counts.tsv
    Read counts of the fusion isoforms, as flair transcriptome's isoform.counts.txt
sample.fusions.isoform.read.map.tsv
    read map of reads to fusion isoforms
sample.fusions.aaseq.tsv
    Protein sequences of the productive fusion isoforms, by the aaseq_id of the bed
sample.prelimfusions.bed
    The preliminary fusions' loci, with their reads and distances to a TSS
sample.genomealigned.chim.bam, sample.transcriptomealigned.chim.bam
    The chimeric reads of the genome alignment and of the alignment to annotated
    transcripts, reused by a later run with the same output name

With ``--keep_intermediate``, the synthetic fusion genome (sample-syntheticFusionGenome.fa),
its annotation, the reads aligned to it (sample.syntheticAligned.bam) and the other
intermediate files are kept as well, for viewing fusion isoforms on the synthetic genome.

Options
-------

.. include:: cli/fusion.rst


.. _combine-label:

flair combine
=============

Combines FLAIR transcriptomes with other FLAIR transcriptomes or annotation
transcriptomes to generate an accurate combined transcriptome. The transcriptomes to
combine are named in one of three mutually exclusive ways: --manifest, --prefixes or
--isoform_beds.

Manifest example, one bed file path per line (we suggest using absolute file paths to
point to your files though):

.. code:: text

    sample1.FLAIR.isoforms.bed
    sample2.FLAIR.isoforms.bed
    sample1.fusion.isoforms.bed
    sample2.fusion.isoforms.bed

To combine FLAIR transcriptomes with annotated transcripts, convert the annotation gtf
to a bed file with gtf_to_bed (see Additional Programs) and name that bed here.

Isoforms of the same splice junctions and strand from different samples are combined
into one.  We recommend giving flair combine an annotation GTF (``-f``).  Each combined
isoform then gets its gene, reference transcript, CDS and productivity from that
annotation, as flair transcriptome would give them, whatever annotation or FLAIR version
the samples were run with.  This lets samples run with different annotations be
combined.  The proteins are numbered anew, in ``combined.aaseq.tsv``.  Fusion isoforms
keep the samples' genes.  Without ``-f``, the samples' genes, CDS and productivity are
used, an isoform that samples gave different genes or reference transcripts is an
error, and isoforms have no protein ids, as each sample numbers its own.

Isoforms of alleles, from flair isoalleles, can be combined, but without their allele
information: the alleles of an isoform are combined into one isoform, with a warning.
Their support is per allele, so each allele is a smaller fraction of its gene's reads
than its isoform is, and ``--min_percent_usage`` should be set lower for them.  Either way, isoforms without an
annotated gene are grouped into novel genes as flair transcriptome groups them, and
isoform and gene ids are numbered in genome order, so they don't depend on the order
of the samples.

Flair combine will generate a counts file, but for the most accurate quantification, we recommend
running FLAIR quantify using all samples against the combined transcriptome

Options
-------

.. include:: cli/combine.rst


.. _quantify-label:


flair quantify
==============

Counts the reads of each isoform in each sample, from the samples' reads aligned to
the genome, for a set of isoforms such as ``flair combine``'s.  Each read is counted
for at most one isoform, one it is full-length for, so the counts of different
samples, and of the isoforms of a gene, can be compared.

**Samples**

``--manifest`` is a tab-separated file of a line per sample, of four fields: its id, its
condition, its batch, and its reads aligned to the genome as a sorted, indexed BAM::

   sample1	condition1	batch1	mydata/sample1.bam
   sample2	condition1	batch1	mydata/sample2.bam
   sample3	condition2	batch1	mydata/sample3.bam
   sample4	condition2	batch1	mydata/sample4.bam

Lines starting with ``#`` and blank lines are skipped.  A line of other than four
fields, fields separated by spaces, a sample id given twice, a missing or unindexed
BAM, and a BAM with none of the isoforms' chromosomes (as ``chr1`` against ``1``) are
errors.  Absolute paths are safest.

**Isoforms**

``--isoform_bed`` is a FLAIR isoform bed, the 25 columns ``flair transcriptome``,
``flair combine`` and ``flair fusion`` write: an isoform's gene id is its 13th column.
A gene's isoforms must be on one chromosome and strand, and an isoform's name used
once on a chromosome.  A set from an older ``flair combine`` can have a gene on both
strands: combine its samples again.  The pieces of a fusion isoform, sharing its name,
are counted as the fusion.

**How reads are counted**

A read counts for an isoform it is full-length for: aligned across all of the
isoform's splice junctions, with few mismatches and indels within 5 bases of each, and
into its first and last exons.  A read whose junctions are a contiguous run of a longer
isoform's is the shorter isoform's, not the longer one's, which it doesn't cover.

- A read whose junctions are exactly those of one isoform, aligned cleanly at its
  splice sites, is counted for it without realignment.
- Other spliced reads are aligned with minimap2 to the sequences of the isoforms of
  their locus and counted for the isoform their best alignment is full-length for.  A
  large deletion in an isoform's first or last exon, as reverse transcriptase skipping
  across double-stranded RNA leaves in a 3' UTR, is tolerated, unless it is on a GT-AG,
  GC-AG or AT-AC splice motif: a read spliced where the isoform isn't is another
  isoform's.
- An unspliced read is counted for the single-exon isoform it covers more than half
  of whose ends are nearest its own, and never for a spliced isoform, as in ``flair
  transcriptome``.  An unspliced read's strand isn't known, so a single-exon isoform of
  either strand can take it.
- A fusion's reads are its chimeric reads, those with a supplementary alignment or 40
  or more bases clipped at an end, aligned to its fused transcript, its pieces' exons
  in order: a read must be full-length across the fusion junctions too.  They aren't
  counted for the fusion's partner genes as well.

Genes that overlap are counted together, and each read is considered once, so a read
fitting the isoforms of two overlapping genes is counted for the one it fits best,
not for both.  With ``--trust_ends``, a read also needs its ends within 50 bases of a
spliced isoform's, for protocols giving full-length molecules.  With
``--normalize_ends``, each spliced isoform is given, for counting, the furthest ends
of its gene's isoforms with terminal exons like its own, so reads aren't lost to where
an isoform's ends were placed.

**Output**

``<output>.counts.tsv`` opens with the two id columns, the gene and the isoform, and
has a column of read counts per sample, named by its id, a row per isoform of a gene:

.. code:: text

   gene_id	isoform_id	sample1	sample2	sample3	sample4
   FLG00000001	FLT00000001	21	12	10	10
   FLG00000001	FLT00000002	7	6	7	15

Each row is one isoform of one gene, named by the two leading id columns described
in :ref:`ids-label`.  A fusion isoform's gene is the fusion's.

The condition and batch of each column are in ``<output>.sample_info.tsv``, written
beside the counts file, which ``flair diffexp`` and ``flair diffsplice`` read:

.. code:: text

   sample_id	condition	batch
   sample1	condition1	batch1
   sample2	condition1	batch1
   sample3	condition2	batch1
   sample4	condition2	batch1

Because these are columns of their own, the id, condition and batch fields may
contain any characters. A counts matrix from an earlier FLAIR, whose columns are
named ``sample_condition_batch``, is still read by taking the fields back out of the
column name.

With ``--tpm``, ``<output>.tpm.tsv`` has the counts as transcripts per million; a
sample with no reads counted has 0 throughout, with a warning.  With
``--generate_map``, ``<output>.<sample id>.read.map.tsv`` lists each sample's reads by
isoform, under a header, a row of an isoform's name and its reads' names,
comma-separated.  With ``--compress`` these are all written gzipped.

**Temporary files**

The reads realigned and minimap2's alignments are written to a directory of the
run's own in ``--temp_dir``, by default ``$TMPDIR`` or the system temporary directory,
which is best on a local disk, and removed when the run ends, even on failure, unless
``--keep_intermediate``.  ``-t`` sets how many batches of genes are counted at once,
and minimap2 is given the threads left idle as the run finishes.

Options
-------

.. include:: cli/quantify.rst



.. _diffexp-label:

flair diffexp
=============


The standard `conda` environment no long installed `R` and the required packages.
These maybe added do the environment as describe in :ref:`installing-label` 

This module performs differential *expression* and differential *usage* analyses between **exactly two** conditions with 
3 or more replicates. Name the two conditions with ``--condition_a`` and ``--condition_b``;
``--condition_a`` is the reference, so fold changes are reported for ``--condition_b``
relative to it. With neither given, the two conditions are used in sorted order, which
makes the control the reference when its name sorts first (eg ctl and test). It does so
by running these R packages:

 - `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__ on genes and isoforms. This tests for differential expression.
 - `DRIMSeq <https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html>`__ is used on isoforms only and tests for differential usage. This is done by testing if the ratio of isoforms changes between conditions.

If you do not have replicates you can use the `diff_iso_usage <#diffisoscript>`__ standalone script.

If you have more than two sample condtions, either split your counts matrix ahead of time or run DESeq2 and DRIMSeq yourself. 

**Outputs**

After the run, the output directory (``--output``) contains the following, where ``A``
is ``--condition_a`` and ``B`` is ``--condition_b``. Fold changes are reported for
``B`` relative to ``A``, so naming the conditions the other way round negates them
and renames these files.

 - ``genes_deseq2_A_v_B.tsv`` Filtered differential gene expression table.
 - ``genes_deseq2_QCplots_A_v_B.pdf`` QC plots, see the `DESeq2 manual <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`__ for details.
 - ``isoforms_deseq2_A_v_B.tsv`` Filtered differential isoform expression table.
 - ``isoforms_deseq2_QCplots_A_v_B.pdf`` QC plots
 - ``isoforms_drimseq_A_v_B.tsv`` Filtered differential isoform usage table
 - ``isoforms_sig_exp_change_norm_by_gene.tsv`` Isoforms whose usage changes, from a t-test on counts normalized by gene.
 - ``workdir`` Temporary files including unfiltered output files.


Options
-------

.. include:: cli/diffexp.rst


Notes
-----

DESeq2 and DRIMSeq are optimized for short read experiments and expect many reads for each expressed gene. Lower coverage (as expected when using long reads) will tend to result in false positives.

For instance, look at this counts table with two groups (s and v) of three samples each:

.. code:: text

    gene   s1    s2      s3      v1      v2      v3
       A    1     0       2       0       4       2
       B  100    99     101     100     104     102

Gene A has an average expression of 1 in group s, and 2 in group v but the total variation in read count is 0-4. The same variation is true for gene B, but it will not be considered differentially expressed.

Flair does not remove low count genes as long as they are expressed in all samples of at least one group so please be careful when interpreting results.

Results tables are filtered and reordered by p-value so that only p<0.05 differential genes/isoforms remain. Unfiltered tables can be found in ``workdir``

.. _diffsplice-label:

flair diffsplice
================

The standard `conda` environment no long installed `R` and the required packages.
These maybe added do the environment as describe in :ref:`installing-label` 

This module calls alternative splicing (AS) events from isoforms. Currently supports
the following AS events: 

 - intron retention (ir)
 - alternative 3’ splicing (alt3)
 - alternative 5’ splicing (alt5)
 - cassette exons (es)

If there are 3 or more samples per condition, then you can run with
``--test`` and DRIMSeq will be used to calculate differential usage of
the alternative splicing events between two conditions. See below for
more DRIMSeq-specific arguments. 

If conditions were sequenced without replicates, then the diffSplice output files can
be input to the `diffsplice_fishers_exact <#diffsplice_fishers>`__
script for statistical testing instead.

**Outputs**

After the run, the output directory (``--output``) contains the following tab separated files:

 - ``diffsplice.alt3.events.quant.tsv``
 - ``diffsplice.alt5.events.quant.tsv``
 - ``diffsplice.es.events.quant.tsv``
 - ``diffsplice.ir.events.quant.tsv``

If DRIMSeq was run, where ``A`` is ``--condition_a`` and ``B`` is ``--condition_b``:

 - ``drimseq_alt3_A_v_B.tsv``
 - ``drimseq_alt5_A_v_B.tsv``
 - ``drimseq_es_A_v_B.tsv``
 - ``drimseq_ir_A_v_B.tsv``
 - ``workdir`` Temporary files including unfiltered output files.

All four are written whichever events were called: an event type that called nothing, or
whose events DRIMSeq found nothing significant in, leaves a file of only its header.

Options
-------

.. include:: cli/diffsplice.rst


Notes
-----

Results tables are filtered and reordered by p-value so that only p<0.05 differential genes/isoforms remain. Unfiltered tables can be found in ``workdir``

For a complex splicing example, please note the 2 alternative 3’ SS, 3
intron retention, and 4 exon skipping events in the following set of
isoforms that ``flair diffSplice`` would call and the isoforms that are
considered to include or exclude the each event:

.. figure:: img/toy_isoforms_coord.svg

.. code::

   a3ss_feature_id       coordinate                    sample1 sample2 ... isoform_ids
   inclusion_chr1:80(+)   chr1:80-400(+)_chr1:80-450(+)     75.0    35.0    ... a,e
   exclusion_chr1:80(+)   chr1:80-400(+)_chr1:80-450(+)     3.0     13.0    ... c
   inclusion_chr1:500(+)  chr1:500-650(+)_chr1:500-700(+)   4.0     18.0    ... d
   exclusion_chr1:500(+)  chr1:500-650(+)_chr1:500-700(+)   70.0    17.0    ... e

.. code::

   ir_feature_id             coordinate        sample1 sample2 ... isoform_ids
   inclusion_chr1:500-650(+)  chr1:500-650(+)    46.0    13.0    ... g
   exclusion_chr1:500-650(+)  chr1:500-650(+)    4.0     18.0    ... d
   inclusion_chr1:500-700(+)  chr1:500-700(+)    46.0    13.0    ... g
   exclusion_chr1:500-700(+)  chr1:500-700(+)    70.0    17.0    ... e
   inclusion_chr1:250-450(+)  chr1:250-450(+)    50.0    31.0    ... d,g
   exclusion_chr1:250-450(+)  chr1:250-450(+)    80.0    17.0    ... b

.. code::

   es_feature_id             coordinate        sample1 sample2 ... isoform_ids
   inclusion_chr1:450-500(+)  chr1:450-500(+)    83.0    30.0    ... b,c
   exclusion_chr1:450-500(+)  chr1:450-500(+)    56.0    15.0    ... f
   inclusion_chr1:200-250(+)  chr1:200-250(+)    80.0    17.0    ... b
   exclusion_chr1:200-250(+)  chr1:200-250(+)    3.0     13.0    ... c
   inclusion_chr1:200-500(+)  chr1:200-500(+)    4.0     18.0    ... d
   exclusion_chr1:200-500(+)  chr1:200-500(+)    22.0    15.0    ... h
   inclusion_chr1:400-500(+)  chr1:400-500(+)    75.0    35.0    ... e,a
   exclusion_chr1:400-500(+)  chr1:400-500(+)    56.0    15.0    ... f

.. _spliceevents-label:

flair spliceevents
==================

For when to count splicing events from reads rather than build isoforms, see
:ref:`choosing-analysis-label`.

This module counts alternative splicing events in each sample directly from its
reads aligned to the genome, without building isoforms first, so it suits
fragmented or degraded libraries and datasets of variable quality, which flair
transcriptome handles less well.  Its events resemble flair diffsplice's, but they are
counted from reads, while diffsplice derives them from isoforms and their counts, so
the two are separate analyses, not steps of one.  Reads are corrected against the annotation, and
optionally short-read junctions, as flair transcriptome corrects them, and assigned to
genes as flair transcriptome assigns them: to annotated genes, a read through two genes
to both, and reads of no annotated gene, over all samples, to novel genes named for
their locus.  Unspliced reads count only as fragments, for starts and ends: one
whose 5' end lies in an annotated last exon, or within 100 bases before it, is a 3'
fragment, whose end counts as one of the gene's transcript ends, and one whose 3' end
lies in an annotated first exon is a 5' fragment, whose start counts as a start.  A
read overlapping several genes' exons is the fragment of the one whose exons share the
most bases with it, if it's at that gene's first or last exon; a tie leaves it out.
Without ``--trust_strand``, the genes of both strands count, as an unspliced read's
strand isn't known.  The first and last exons are those of ``--annot_basic``, if given.  Other unspliced reads are left out.  It finds these
events:

 - exon skipping (es), and complex exon skipping (ces): a skipping junction against
   every way the reads include the exons between its ends, several exons, or one
   exon with alternative splice sites
 - alternative 5' and 3' splice sites (alt5, alt3), and alternative 5' and 3'
   splice sites shared by a group of junctions (alt5ss, alt3ss)
 - intron retention (ir): a junction's reads against the reads with an exon across it,
   and in its ``other`` row the reads of junctions both of whose splice sites are
   within 100 bases of its own, which don't retain it either.  A junction with under
   ``--junc_frac_of_event`` of all the reads over it in every sample, those of every
   junction overlapping it too, is none, as a shifted copy of an intron or a junction
   of few reads nested in a long one, whose retained reads are the long intron's
 - alternative first and last exons (afe, ale): only exons whose junction reads,
   over all samples together, rarely have an exon before (or after) them count, as
   a truncated read makes an internal exon look like a first or last one
 - alternative transcription starts and ends (tss, tts)

Each event has a row for each of its alternatives, such as an exon's inclusion and
exclusion, and an ``other`` row for its reads that fit none of them.  An event is
found from its alternatives of at least ``--feature_support`` reads in some sample
that are alternative: their part of the event's reads at least
``--junc_frac_of_event`` in some sample, and at most 1 minus it in some sample, the
same one or another, so one used alone in a sample and not at all in another is.
How that part differs between samples is left to what reads the tables, as
spliceevents_outliers.  An event's main
alternative, too common to be alternative, is a row too when it has a tenth of the
event's reads in some sample, and the reads of the rest are its ``other`` row.

The events of splicing (es, ces, alt5, alt3, alt5ss, alt3ss, ir) are found among
reads linked by shared splice sites, whatever their genes, so overlapping genes and
readthroughs don't split an event's reads.  Each is then given the gene annotated
with most of its introns, those of all its alternatives, the ``other`` row's too:
each intron scores 2 for a gene annotated with it, or failing that 1 for each of its
splice sites the gene has, weighted by the intron's share of the reads of all the
event's introns, so a rare alternative counts little.  A tie goes to the gene of more annotated transcripts, so
a gene rather than a readthrough of it.  Genes still tied are all given, comma-separated.
An event sharing no intron or splice site with the annotation, as in a novel gene,
gets the gene of most of its reads.  First and last exons and starts and ends are
found per gene: a read's first exon and start with the gene it starts in, its last
exon and end with the gene it ends in.

An event is named ``TYPE|CHROM:START-END(STRAND)|GENE``, its locus, 1-based, being
its skipped exon for es, and otherwise the outer bounds of its introns, so it can be
pasted into a genome browser, and its gene, as the ``gene`` column, the genes
comma-separated if tied, so the events of two genes of one locus, as the ends of a gene
and a readthrough of it, are told apart.  Each row is one alternative of its event, named by
what sets it apart from the event's other rows: ``inc`` and ``exc`` for es and
``exc`` for ces, the exons
included of a ces's other rows, ``retained`` and ``spliced`` for ir, the introns of
each alternative of alt5, alt3, alt5ss, alt3ss, afe and ale, the position of each
end of tss and tts, and ``other``.  Introns are written as STAR writes them, from
their first base to their last, and exons from their first base to their last.
For example, an exon skip and an alternative 3' splice site::

   es|chr12:25215437-25215560(-)|ENSG00000133703.12    inc
   alt3|chr1:201-350(+)|G1                             201-300

The tables' coordinate columns use the same 1-based coordinates.  An event's PSI
is computed two ways: of the event's reads (PSIjunc), which is more sensitive, and of
all reads covering the event's locus (PSItot), which favours events that are strong
overall.  A sample needs ``--event_support`` reads of an event for a PSI.

With ``--region_bed``, only events inside its regions are reported: those whose
junctions and exons all lie in a region, whose reads all overlap it, and the start,
end, first and last exon events of genes whose annotated junctions (a novel gene's
reads' junctions) all lie in one, so a region around a gene's annotated exons, or only
its spliced part, reports them.
Reads overlapping a region are corrected over their whole length, wherever the region
ends.  The samples of a region are corrected in parallel.

Without ``--region_bed``, each chromosome is cut in intergenic gaps of the annotation
into pieces of about 5000 annotated exons, which are corrected and searched for events
in parallel.  A piece owns the reads starting in it, corrected over their whole length,
and reports the events starting in it, so an event close to a cut can miss the reads
of the piece before.  Reads are aligned to their piece's annotated transcripts, not
the whole chromosome's.  ``--partition no`` keeps chromosomes whole, slower but
without those losses, and ``--partition yes`` cuts ``--region_bed`` regions too.

The manifest is a tab-separated file of a sample name and the path to its sorted and
indexed bam file, a line per sample::

   sample1	sample1.genomealigned.bam
   sample2	sample2.genomealigned.bam

**Outputs**

 - ``OUTPUT.spliceevents.counts.tsv``, the reads of each event row in each sample,
   with the row's ``event``, ``alternative``, ``event_type``, ``gene`` and
   ``gene_name``, the annotation's name of the gene, empty if it has none,
   ``additional_genes`` and ``additional_gene_names``, the other annotated genes
   sharing an intron or splice site with a splicing event, best scoring first, to
   find all the events that could be a gene's, and its
   ``junctions_included``, ``junctions_excluded``, ``outer_junctions`` and
   ``exons``, and ``novel_junctions``, the row's introns that aren't annotated:
   those its reads splice, the skipping junctions of an exon skip's exclusion.  An
   event of an annotated gene can have novel alternatives, so these, rather than
   the gene, find novel splicing
 - ``OUTPUT.spliceevents.PSIjunc.tsv`` and ``OUTPUT.spliceevents.PSItot.tsv``, the
   rows' PSIs, empty where a sample has too few reads
 - ``OUTPUT.spliceevents.totals.tsv``, each event's reads (``event_reads``) and the
   reads covering its locus (``locus_reads``) in each sample, the PSIs' denominators
 - ``OUTPUT.spliceevents.bed``, the event rows for a genome browser, named
   ``EVENT|ALTERNATIVE`` and scored by their median PSIjunc

Outlier samples, whose PSI of an event is far from the other samples', can be found
from these tables with the `spliceevents_outliers <scripts.html#spliceevents-outliers>`__
script.

Options
-------

.. include:: cli/spliceevents.rst

.. _variantquant-label:

flair variantquant
==================

Counts, in each sample, the reads showing each SNV allele at the sites of a vcf,
from reads aligned to the genome.  It was built to quantify RNA editing at the sites
of a database such as REDIportal, and counts the variants called in the samples as
well, by haplotype too, for the allele bias of a somatic variant.

**Samples**

A sample is given by ``--bam``, named by ``--sample_name`` or else for its bam file,
or several by ``--manifest``, a tab-separated file of a line per sample, either of a
sample name and its bam file, with the sites from ``--vcf``::

   tumor	tumor.genomealigned.bam
   normal	normal.genomealigned.bam

or of a sample name, its bam file and the vcf of its own variants, on every line::

   tumor	tumor.genomealigned.bam	tumor.vcf.gz
   normal	normal.genomealigned.bam	normal.vcf.gz

whose sites are merged, every sample counted at all of them, so a site called in one
sample has the reads of the others too.  Lines starting with ``#`` and blank lines are
skipped; a vcf on some lines only, duplicate sample names and missing files are
errors.  The bam files are read whole, so they needn't be sorted or indexed, and
``-t`` of them are read at once, each in its own process.

**Sites**

The vcf may be plain or gzipped, and only its first five columns are needed, so a
database's list of sites can be given as a vcf of them.  Only SNVs are counted: an
indel or a substitution of more than one base, which a read's base can't show, is
left out with a warning.  Each alt of a multi-allelic site is an allele of its own.
Only records whose FILTER is ``PASS`` or ``.``, or that have no FILTER, are sites, so a
DeepVariant ``RefCall`` isn't; ``--keep_flags`` counts the records of other flags too,
as ``--keep_flags RefCall``, or ``all`` of them.  A site in more than one vcf has the
alts of all of them, and two vcfs giving it different REFs are an error.

The vcf must be of the genome's strand: a record's REF is the genome's base, so an
A-to-I edit of a minus strand gene is T>C, not A>G.  A record of the transcript's
strand is never counted, and one whose REF is neither the genome's base nor its
complement counts every read of the genome's base as showing its alt.  Given the
genome, ``-g``, variantquant checks each site's REF and leaves out those that aren't
the genome's base, with a warning of how many, how many are the complement, and a few
of them.  Of REDIportal's 15.6 million A-to-I sites, 41,442 are of the transcript's
strand and 247 of neither, all left out by ``-g``.

**Reads**

A read's primary and supplementary alignments are counted, not its secondary ones,
which would count its bases again, and only alignments of at least ``--quality``
mapping quality (0, all, by default).  A read's base at a site is found from its
cigar, so MD tags aren't needed; a site in a read's deletion or intron isn't covered by
it.  A read showing one of a site's alts is counted as showing it, and any other base
as not.

With ``--trust_strand``, as for oriented reads or direct RNA, each read is of its
alignment's strand, and each strand of reads at a site is counted in a row of its
own, so the reads of an antisense transcript aren't counted with the sense
transcript's.  With ``--haplotypes``, the reads of a haplotagged bam, as WhatsHap
haplotag writes it, are counted by their phase set and haplotype, their ``PS`` and
``HP`` tags, a row each, and the reads without an ``HP`` tag in a row of their own: at
a somatic variant, the haplotype carrying it shows its alt reads against the other's,
and the normal's haplotypes both show the reference.  A haplotype is comparable only
within a phase set and one phasing, so a run's samples should be one patient's,
haplotagged with its phased variants; samples sharing no phase set are warned of.

**Genes**

A site is named by the genes of the exons over it, of ``--isoform_bed``, ``--gtf`` or
both, and is of no gene without them.  An isoform's gene is its bed's 13th column, as
a flair isoform bed gives it.  An isoform without one, as of a BED12, given a GTF, is
of the GTF gene its exons overlap most on its strand, so a site in a novel exon of a
known gene is named by it; with no GTF it names no site.  Given a GTF alone, its
exons name the sites by their ``gene_id``.  A site of more than one gene is named by
the first, sorted.

**Outputs**

 - ``OUTPUT.var.counts.tsv``, a row per allele of a site, of each group of its reads,
   of the columns ``varpos``, the site, ``CHROM:POS`` 1-based, ``ref``, ``alt``,
   ``strand``, ``phase_set`` and ``haplotype``, ``.`` when the reads aren't split by
   them, and ``gene``, then a column per sample of ``UNMODIFIED;MODIFIED``: the reads
   covering the site not showing the alt, and those showing it, so the two add up to
   the reads covering the site whatever its other alts.  An allele's rows are written
   if a read of some sample shows the alt, or with ``--output_all``, each if some
   sample has ``--min_coverage`` reads of its group covering the site

**Speed**

Each bam is read once, whatever the number of sites, about 3 minutes for 15 million
reads with a few hundred sites and 36 with REDIportal's 15.6 million, and ``-t``
samples at once.  For 100,000 REDIportal sites across the genome, in a bam of 15.5
million reads, it takes 5 minutes, against 3.2 hours for ``samtools mpileup -l`` and
about 9 hours for pysam's ``pileup``, their counts the same.  ``samtools mpileup`` is
faster for sites packed in a small region, which it reads alone: 1 minute, against 3,
for 100,000 sites in 10 Mb.  The sites are loaded once: REDIportal's in 37 seconds
and about 2.2 GB from its vcf, 43 seconds with a flair isoform bed's genes, 57 with
GENCODE's, and 7 seconds more to check them against the genome with ``-g``.
Temporary files, a file per sample of its reads covering the sites, are written in a
directory of each run's own in ``--temp_dir``, by default ``$TMPDIR``, and removed
when it ends unless ``--keep_intermediate``.

Options
-------

.. include:: cli/variantquant.rst

.. _alleles-label:

flair alleles
=============

Calls alleles from variants and groups the reads that carry them, writing a vcf of
allele groups.  A normal bam and vcf may be given alongside the tumor ones to call
variants as somatic or not.

Options
-------

.. include:: cli/alleles.rst

.. _isoalleles-label:

flair isoalleles
================

Groups alleles by isoform, using the allele calls from flair alleles together with an
isoform bed, and predicts the protein each isoform allele produces.

Options
-------

.. include:: cli/isoalleles.rst
