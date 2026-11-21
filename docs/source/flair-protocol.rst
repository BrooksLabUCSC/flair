.. _flair-protocol-label:

FLAIR protocol
==============

- Align reads to the genome
- Obtain novel splice junction reference
- Choose an analysis: isoforms, or splicing events counted from reads
- For isoforms:

  - Generate transcriptomes, and detect gene fusions if wanted
  - Combine transcriptomes between samples
  - Quantify splicing relative to the combined transcriptome
  - Optionally, split each isoform's reads by haplotype and somatic variants
  - Test for differences with ``flair diffexp`` and ``flair diffsplice``

- For splicing events counted from reads: run ``flair spliceevents``, and test
  for differences with ``flair diffexp``

Aligning reads to the genome
----------------------------

FLAIR takes reads aligned to the genome as a sorted, indexed BAM, from any spliced
aligner.  Align each sample on its own.  Samples of fewer than 10M reads can instead
be aligned as merged replicates, for the coverage; FLAIR has been run on samples of
about 1M reads.

- minimap2, recommended: the aligner the other two wrap, run directly, with all of
  its options available

  - Raw reads in a BAM, as PacBio and ONT write them, keep their tags, such as UMIs
    and base modifications: ``samtools fastq -T '*'`` writes every tag into each
    read's FASTQ comment, and minimap2's ``-y`` copies the comments back into the
    alignments as tags
  - From a BAM of raw reads:

    - ``samtools fastq -T '*' SAMPLE.reads.bam | minimap2 -ax splice -s 40 -G 350k -t 25 --MD --secondary=no -y GENOME.fa - | samtools sort -o SAMPLE.genome_aligned.bam - && samtools index SAMPLE.genome_aligned.bam``
    - From a FASTA or FASTQ of reads, give it to minimap2 in place of ``-``,
      without ``samtools fastq`` and ``-y``
    - For a stranded library, add ``-uf``, which takes a read's splice junctions on
      the read's own strand

- pbmm2 and dorado align: the minimap2 wrappers of PacBio and ONT

  - They align a raw BAM directly, keeping its tags, UMIs and base modifications
  - They expose few of minimap2's options, and they change more often than
    minimap2 does

:ref:`alignment-label` explains the minimap2 options FLAIR depends on.

Obtaining a novel splice junction reference
-------------------------------------------

FLAIR reports only splice junctions that are in the annotation or in a file of
supported junctions: ``flair transcriptome`` and ``flair spliceevents`` leave out a
read with any other junction, which in a locus with novel splicing can be most of
its reads.  Give a junction file even with a well annotated transcriptome, or the
novel junctions of that locus are lost.  The junctions come either from short reads
of the same material or from the long reads themselves.

Use one junction file for every sample of a dataset.  Where that is not possible, as
in a dataset that grows over time, a file per sample works, and several files are
combined by concatenating them.

- Junctions from short reads

  - They are evidence independent of the long reads, so they do not share the
    artifacts of long-read alignment, and they place splice sites accurately
  - They require sequencing of their own, and short reads of other material than
    the long reads are less sensitive to the sample's own junctions
  - Align the short reads with STAR, which writes SJ.out.tab by default, and give
    that file to ``flair transcriptome`` with ``--junction_star``
  - Short reads aligned with another aligner: extract the junctions with
    `intron-prospector <https://github.com/diekhans/intron-prospector>`__, as for
    long reads below, and give them with ``--junction_bed``

- Junctions from the long reads

  - They are the junctions of the reads the transcripts are built from, and
    intron-prospector leaves out the alignment artifacts
  - They are not independent evidence, since they have the same biases as the reads
    they come from, and a junction of few reads is less likely to be kept
  - From the indexed BAM of the genome alignment:
    ``intron-prospector --genome-fasta=GENOME.fa --intron-bed6=SAMPLE.ip_junctions.bed -C 0.0 SAMPLE.genome_aligned.bam``

.. _choosing-analysis-label:

Choosing an analysis: isoforms or splicing events
-------------------------------------------------

FLAIR compares splicing between samples in two ways, and the choice follows from how
complete the reads are.

- Isoforms: ``flair transcriptome``, then ``flair combine``, ``flair quantify``, and
  ``flair diffsplice`` or ``flair diffexp``

  - Builds full-length isoforms from the reads of each sample, combines them into
    one transcriptome for all samples, and counts each sample's reads that match an
    isoform.  ``flair diffsplice`` then derives splicing events from the isoforms
    and their counts and tests them between conditions; ``flair diffexp`` tests the
    expression of genes and isoforms.
  - It reports which exons, splice sites and ends are used together in one
    transcript, and from that the protein the transcript encodes.  The events
    ``flair diffsplice`` reports name the isoforms they come from.
  - It depends on reads that cover most of a transcript.  A read is counted only if
    it matches an isoform, and ``flair quantify`` counts only unambiguous matches,
    so 55 to 90 percent of the reads are counted.  With degraded RNA, or library
    artifacts that truncate reads, fewer reads match, and the samples most affected
    lose the most.

- Splicing events counted from reads: ``flair spliceevents``

  - Finds splicing events (exon skipping, alternative splice sites, intron
    retention, alternative first and last exons, and transcript starts and ends)
    in each sample's reads, without building isoforms, and gives the reads of each
    alternative of an event and its PSI, percent spliced in, the share of the
    event's reads that have that alternative.
  - It recovers more read support for an event: a read counts toward every event it
    covers, whether or not it spans a whole transcript, so fragmented and degraded
    samples keep most of their reads, and samples of different quality can be
    compared.  It also runs on a few loci on their own, with ``--region`` or
    ``--region_bed``.
  - An event is harder to interpret than an isoform: it states nothing about the
    other exons of the transcript it comes from, events overlap, one read counts
    toward several, and an event is tied to no isoform or protein.  An event uses
    only what a short read spanning the same junctions would show.  Unspliced reads
    count only toward starts and ends, as fragments of an annotated first or last
    exon.
  - Its counts are tested between conditions by ``flair diffexp`` and
    ``diff_iso_usage``, as an isoform's usage is, and outlier samples are found
    with ``spliceevents_outliers``.

The two analyses are comparable only in their results: the events of ``flair
spliceevents`` resemble those of ``flair diffsplice``, which follows transcriptome,
combine and quantify, but they are counted from reads rather than derived from
isoforms, so the two are separate analyses, not steps of one.

- Use isoforms for reads that are mostly full length, and for transcript structures,
  proteins, or isoform and gene expression.
- Use splicing events counted from reads for reads that are fragmented or of uneven
  quality between samples, and for specific events or loci.
- The two can be used together, for example isoforms from the samples with the most
  complete reads and splicing events from all of them, or isoforms for short genes
  and splicing events for long genes whose reads rarely cover a whole transcript.

Generating transcriptomes
-------------------------

``flair transcriptome`` builds a sample's isoforms from its genome-aligned BAM.  It
replaces the correct and collapse modules of FLAIR 3.0: it builds transcripts one
cluster of overlapping reads at a time, in parallel and directly from the BAM, in
less memory and less time, and it reports novel genes and transcripts that those
modules did not.  It runs with or without an annotation GTF, which gives its isoforms
their genes, reference transcripts and coding sequences.  The commands below are a
starting point; :ref:`transcriptome-label` documents every option.

- With junctions from short reads:

  - ``flair transcriptome -g GENOME.fa -f ANNOTATION.gtf --junction_star SHORTREAD.SJ.out.tab -b SAMPLE.genome_aligned.bam --sample_name SAMPLE -o SAMPLE``

- With junctions from long reads:

  - ``flair transcriptome -g GENOME.fa -f ANNOTATION.gtf --junction_bed LONGREAD.ip_junctions.bed --junction_support 2 -b SAMPLE.genome_aligned.bam --sample_name SAMPLE -o SAMPLE``

- The sample name labels the sample's isoforms and reads in the outputs, and
  ``flair combine``'s columns, so give each sample its own
- For a few loci, add ``--region chr12:25,205,001-25,210,000``, repeated per
  locus, or ``--region_bed REGIONS.bed``.  The regions are expanded over the reads
  and annotated transcripts overlapping them, so the transcripts of the genes they
  cover are those of a whole-genome run.  The expansion can reach into a gene next
  to a region, whose transcripts are reported but may be built from only part of
  its reads

Detecting gene fusions
----------------------

To study gene fusions alongside normal isoforms, run ``flair fusion`` on the same
BAM file as ``flair transcriptome``, with the same sample name and output prefix.
Its fusion isoforms are then combined with the sample's normal isoforms and
quantified with them, below.  ``flair fusion`` needs stranded reads, and reads
aligned with ``-s 40``, as in the recommended minimap2 command, to keep short
partner pieces.

- ``flair fusion -g GENOME.fa -f ANNOTATION.gtf -b SAMPLE.genome_aligned.bam --sample_name SAMPLE -o SAMPLE``

This writes ``SAMPLE.fusions.isoforms.bed`` beside ``flair transcriptome``'s
``SAMPLE.isoforms.bed``, and ``SAMPLE.fusions.tsv``, the fusions found.

Use the same sample name for both runs.  The name only reaches ``flair combine``'s
outputs, its ``samples`` column and counts file, which then have one column per
sample rather than one for its isoforms and another for its fusions.  FLAIR
quantify names its columns from its own manifest.

For the most fusion isoforms, as for reporting fusions, ``flair fusion``'s
``--filter comprehensive --frac_support 0`` keeps the isoforms of a fusion's
5' truncated reads (see the ``flair fusion`` documentation).  ``flair combine`` filters
fusion isoforms as it does normal ones, at ``--frac_support`` of the
fusion's reads and ``--min_reads``, so it drops the rarest of them again unless
those are lowered, for all isoforms.  Report fusions from ``SAMPLE.fusions.tsv``.

Combine transcriptomes between samples
--------------------------------------

Samples are compared against one transcriptome that covers all of them: build a
transcriptome for each sample, or for each batch of samples, then combine the
transcriptomes with ``flair combine``.

- Make a manifest of the transcriptomes, one path per line; see
  :ref:`combine-label`.  With ``flair fusion``, list each sample's
  ``SAMPLE.fusions.isoforms.bed`` as well as its ``SAMPLE.isoforms.bed``
- Run ``flair combine`` (``-f ANNOTATION.gtf`` is optional, and gives each
  combined isoform its annotated gene and protein):

  - More stringent combination (keep only isoforms used at 5% or more of
    their locus in at least 1 sample, with at least 3 reads over all samples,
    the defaults)

    - ``flair combine --manifest MANIFEST.txt -g GENOME.fa -f ANNOTATION.gtf -o OUTPUT``

  - Less stringent combination (keeps anything supported in any file,
    just combines based on splice junctions and similar ends)

    - ``flair combine --manifest MANIFEST.txt -g GENOME.fa -f ANNOTATION.gtf -o OUTPUT --frac_support 0 --min_reads 1``

Quantify splicing relative to the combined transcriptome
--------------------------------------------------------

Count every sample again on the combined transcriptome with ``flair quantify``, so
that the counts of all samples are of the same isoforms.  It counts a read for an
isoform only where the read matches that isoform and no other, with no expectation
maximization and no partial assignment, so each count is of reads that fit the
isoform alone.  That is why fewer reads are counted than a gene-level method counts,
55 to 90 percent of them.

- Make a new manifest of all samples, a line per sample with its id, condition,
  batch and its reads aligned to the genome as a sorted, indexed BAM,
  tab-separated; see :ref:`quantify-label`.  These are the BAM files the
  transcriptomes were built from.

- Run ``flair quantify``:

  - ``flair quantify --manifest MANIFEST.tsv -g GENOME.fa --isoform_bed COMBINED.isoforms.bed -o OUTPUT``
  - Add ``--normalize_ends`` unless your reads' ends are reliable
    (``--trust_ends``).  The counts matrix names each isoform's gene in its own
    column, as ``flair diffexp`` expects
  - Fusion isoforms in the combined transcriptome are counted too, from the same
    BAM files.  A read counted for a fusion is not counted again for an isoform
    of either partner gene
  - If performing variant calling and incorporate into isoforms with FLAIR isoalleles, add ``--generate_map`` to FLAIR 
    quantify for each sample's read map. The read maps must be made with the same isoform BED given to FLAIR isoalleles.

Isoforms by haplotype and somatic variants
------------------------------------------

``flair isoalleles`` splits each isoform's reads by haplotype, and in a tumor and
normal pair by the somatic variants they have, and gives the sequence and protein
of each group, such as an isoform's tumor reads on the haplotype with a KRAS
G12V and on the other haplotype.  It does not call or phase variants itself: it
takes the variants and haplotagged reads of upstream tools, and works with any
variant caller whose VCF has genotypes.  See the ``flair isoalleles`` documentation
for every option and output.

1. Build the isoforms and quantify every sample, as above, adding
   ``--generate_map`` to ``flair quantify`` for each sample's read map.  The read maps
   must be made with the same isoform BED given to ``flair isoalleles``.

2. Call variants from each sample's genome-aligned BAM.  Any caller writing a VCF
   with genotypes works.  These are the callers FLAIR is run with:

   - DeepVariant, which calls SNVs and indels, with a model for long-read RNA
     where one fits your data, such as its MAS-Seq model for PacBio Kinnex reads
   - longshot or longcallR, if you do not need indels: both call SNVs only from long
     reads.  For example
     ``longshot --bam SAMPLE.bam --ref GENOME.fa --out SAMPLE.vcf``

3. Phase the germline variants and haplotag the reads with WhatsHap, which is not
   in the conda environment; see :ref:`whatshap-install-label`.  In a pair,
   phase the normal's variants, using both samples' reads, and haplotag both BAMs
   with that phasing, so H1 and H2 mean the same in the tumor and the normal.
   ``--ignore-read-groups`` is needed when a BAM's read groups do not name the VCF's
   sample.

   - ``whatshap phase --ignore-read-groups --reference GENOME.fa -o NORMAL.phased.vcf.gz NORMAL.vcf.gz NORMAL.bam TUMOR.bam``
   - ``tabix -p vcf NORMAL.phased.vcf.gz``
   - ``whatshap haplotag --ignore-read-groups --reference GENOME.fa -o TUMOR.haplotagged.bam NORMAL.phased.vcf.gz TUMOR.bam``,
     then the same for the normal, and ``samtools index`` each

   For one sample, phase and haplotag its own variants and reads the same way.
   If you phase only some of the variants, as only the heterozygous ones, give FLAIR
   isoalleles all of them with ``--germline_vcf`` too, or the homozygous variants
   are left out of the sequences.

4. Run ``flair isoalleles``.

   - One sample: its reads by haplotype, H1, H2 or untagged

     - ``flair isoalleles --isoform_bed COMBINED.isoforms.bed -g GENOME.fa --bam SAMPLE.haplotagged.bam --read_map OUTPUT.isoforms.SAMPLE.read.map.tsv --phased_vcf SAMPLE.phased.vcf.gz --sample_name SAMPLE -o SAMPLE``

   - A tumor and normal pair: the tumor's reads by haplotype and somatic variants,
     such as H2-S1, with the normal's as a control.  Give the normal's BAM and read
     map, its phased VCF, and the somatic variants in one of three ways:

     - ``--tumor_vcf TUMOR.vcf.gz``, the tumor's calls from the same caller as the
       normal's: those the normal does not call are somatic, if at least 10 normal
       reads cover the site (``--min_normal_depth``)
     - ``--somatic_vcf SOMATIC.vcf.gz``, a somatic caller's calls, such as ClairS or
       DeepSomatic, taken as they are
     - ``--variant CHROM:POS:REF>ALT``, only the variants given, as a known driver
       that the normal also calls because tumor cells are in it
     - ``flair isoalleles --isoform_bed COMBINED.isoforms.bed -g GENOME.fa --bam TUMOR.haplotagged.bam --read_map OUTPUT.isoforms.TUMOR.read.map.tsv --sample_name TUMOR --normal_bam NORMAL.haplotagged.bam --normal_read_map OUTPUT.isoforms.NORMAL.read.map.tsv --normal_name NORMAL --phased_vcf NORMAL.phased.vcf.gz --tumor_vcf TUMOR.vcf.gz -o PAIR``

   ``OUTPUT.isoalleles.counts.tsv`` is a counts matrix of each isoform's groups, so
   ``diff_iso_usage`` compares their usage between the tumor and the normal.

Testing for differences between conditions
------------------------------------------

``flair quantify`` writes the conditions and batches of its manifest beside the counts
matrix, in ``OUTPUT.isoforms.sample_info.tsv``, which ``flair diffexp`` and
``flair diffsplice`` read.
Name the two conditions to compare; the first is the reference that changes are
measured against.  ``flair diffexp`` needs at least 3 replicates of each condition,
and ``flair diffsplice`` by default tests only events covered in 6 samples
(``--min_samps_gene_expr``), as with 3 of each.

- Differential gene and isoform expression, and isoform usage (DESeq2 and
  DRIMSeq):

  - ``flair diffexp --counts_matrix OUTPUT.isoforms.counts.tsv --condition_a CONTROL --condition_b TREATED -o DIFFEXP``

- Differential splicing events, derived from the isoforms and their counts, and
  tested with DRIMSeq when the conditions are named (or with ``--test``, which
  takes the only two conditions in sorted order); otherwise only the events'
  counts are written:

  - ``flair diffsplice --isoform_bed COMBINED.isoforms.bed --counts_matrix OUTPUT.isoforms.counts.tsv --condition_a CONTROL --condition_b TREATED -o DIFFSPLICE``

- Without replicates, compare isoform usage between two samples with Fisher's
  exact test, with p-values adjusted for the isoforms tested:

  - ``diff_iso_usage OUTPUT.isoforms.counts.tsv SAMPLE1 SAMPLE2 OUTPUT.usage.tsv``

Counting splicing events from reads
-----------------------------------

If you chose splicing events counted from reads, ``flair spliceevents`` takes the
genome-aligned BAM files directly, with a manifest of each sample's name and
BAM file, tab-separated.  As with ``flair transcriptome``, include a junction
reference.

- ``flair spliceevents -m MANIFEST.tsv -g GENOME.fa -f ANNOTATION.gtf --junction_bed LONGREAD.ip_junctions.bed -o OUTPUT``
- For a few loci, add ``--region chr12:25,205,001-25,210,000`` or
  ``--region_bed REGIONS.bed``; only the regions' reads and annotation are used,
  and an event reaching past a region is left out rather than counted from part
  of its reads.
- To find samples whose PSI of an event is far from the others':
  ``spliceevents_outliers -i OUTPUT -o OUTPUT``
- To test events between conditions, ``flair diffexp`` reads
  ``OUTPUT.spliceevents.counts.tsv`` and tests each event's alternatives as it
  does an isoform's usage (DRIMSeq only).  Write each sample's condition and
  batch beside it in ``OUTPUT.spliceevents.sample_info.tsv``, a header of
  ``sample_id``, ``condition`` and ``batch`` and a line per sample, in the order
  of its columns:

  - ``flair diffexp --counts_matrix OUTPUT.spliceevents.counts.tsv --condition_a CONTROL --condition_b TREATED -o DIFFEXP``
  - Without replicates:
    ``diff_iso_usage OUTPUT.spliceevents.counts.tsv SAMPLE1 SAMPLE2 OUTPUT.usage.tsv``
