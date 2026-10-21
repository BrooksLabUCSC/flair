.. _flair-protocol-label:

FLAIR protocol
==============

- Align reads to the genome
- Obtain novel splice junction reference
- Choose an analysis: isoforms, or splicing events counted from reads
- For isoforms:

  - Generate transcriptomes
  - Combine transcriptomes between samples
  - Quantify splicing relative to the combined transcriptome
  - Test for differences with FLAIR diffexp and diffsplice

- For splicing events counted from reads: run FLAIR spliceevents

Aligning reads to the genome
----------------------------

FLAIR is now very flexible about read alignment strategies; please use
whatever makes the most sense for your data. We now recommend aligning
each sample individually. The exception is if your samples have less
than 10M reads. In that case, you may want to combine replicates for
better coverage. However, we have tested and obtained good results with
~1M reads.

- Pbmm2 / Dorado align

  - These are the minimap2 wrappers developed by PacBio and ONT
  - Pros: will align raw BAM files directly and preserve tags, UMIs, and
    modification information
  - Cons: minimal flexibility in terms of alignment options, less stable
    than minimap (more frequent updates)

- Minimap2

  - This is the basal aligner used by all of these methods, which can
    also be run directly
  - Pros: maximum flexibility and control in terms of alignment options
  - Cons: this requires FASTA/FASTQ read files, so if your raw reads are
    provided as BAM files, you would need to convert them to
    FASTA/FASTQ, losing information stored in the BAM tags. Also, if
    your reads are already aligned, realigning is cumbersome.
  - Recommended alignment command:

    - ``minimap2 -ax splice -s 40 -G 350k -t 25 --MD --secondary=no GENOME.fa SAMPLE.fa | samtools view -hb - | samtools sort - > SAMPLE.genomealigned.bam; samtools index SAMPLE.genomealigned.bam``
    - If your library is stranded, you may want to add ``-uf`` to the
      minimap2 command (forces splice junctions to be on the same strand
      as the read)

Obtaining novel splice junction reference
-----------------------------------------

FLAIR requires some kind of splice junction reference. Even if you have
a well-annotated reference transcriptome, we highly recommend including
some reference for novel splice junctions. FLAIR will not report any
splice junctions not found in either the reference transcriptome or the
supported junctions file: both FLAIR transcriptome and FLAIR spliceevents
leave out reads with any other junction, which in a locus with novel
splicing can be most of its reads. You can either use short-read junctions (from
paired short-read sequencing or from some other sequencing of the same
tissue type) or calculate supported junctions directly from the
long-read data

- NOTE: We usually recommend that you use the same reference splice
  junctions file for all samples within a dataset, but when that is not
  possible (when the dataset is growing over time), using
  sample-specific files is okay.

  - For combining files between samples, concatenating the files into
    one is fine.

- Using short-read junctions:

  - Pros: provides orthogonal support, minimizes long-read specific
    biases/artifacts
  - Cons: requires additional sequencing, or if non-matched sequencing
    is used, could be less sensitive to sample-specific junctions
  - Usage:

    - Preferred: align short-reads with STAR. STAR will output a
      SJ.out.tab file by default; you can use this file directly with
      the –junction_tab option in FLAIR transcriptome
    - If your short-reads have already been aligned with another
      aligner, you can extract the junctions with the junctions_from_sam
      script included with FLAIR or using
      `intron-prospector <https://github.com/diekhans/intron-prospector>`__.

- Using long-read junctions:

  - Pros: detects and removes artifacts, direct junction detection from
    the same reads you’re assembling transcripts from (best match)
  - Cons: Potentially less sensitive to rare events with low read
    support, not orthogonal support (has the same lrRNA biases)
  - Usage:

    - After genomic alignment of your reads (have indexed BAM file), use
      intron-prospector to generate a junctions.bed file:
    - ``intron-prospector --genome-fasta=GENOME.fa --intron-bed6=SAMPLE.IPjunctions.bed -C 0.0 SAMPLE.genomealigned.bam``

.. _choosing-analysis-label:

Choosing an analysis: isoforms or splicing events
-------------------------------------------------

FLAIR can compare splicing between samples in two ways, and the choice
depends mostly on how complete your reads are.

- Isoforms: FLAIR transcriptome, then FLAIR combine, FLAIR quantify, and
  FLAIR diffsplice (or FLAIR diffexp)

  - Builds full-length isoforms from the reads of each sample, combines
    them into one transcriptome for all samples, and counts each sample's
    reads that match an isoform.  FLAIR diffsplice then derives splicing
    events from the isoforms and their counts, and tests them between
    conditions; FLAIR diffexp tests the expression of genes and isoforms.
  - Pros: makes full use of long reads: which exons, splice sites and ends
    are used together in one transcript, and so the proteins it can encode
    and novel transcript structures.  Isoforms are easy to interpret and to
    view, and the events FLAIR diffsplice reports are explained by them.
  - Cons: depends on reads that cover most of a transcript.  A read is
    counted only if it matches an isoform, and FLAIR quantify assigns only
    unambiguous matches, so you will likely recover 55-90% of your reads.
    With degraded RNA, or library artifacts that truncate reads, many fewer
    reads match, and the samples most affected lose the most.

- Splicing events counted from reads: FLAIR spliceevents

  - Finds splicing events (exon skipping, alternative splice sites, intron
    retention, alternative first and last exons, and transcript starts and
    ends) directly from each sample's reads, without building isoforms, and
    gives the reads of each alternative of an event and its PSI in every
    sample.
  - Pros: recovers more read support for events.  A read counts toward
    every event it covers, whether or not it spans a whole transcript, so
    fragmented and degraded samples keep most of their reads, and samples
    of different quality can be compared.  It also runs quickly on a few
    loci of interest, with ``--region_bed``.
  - Cons: harder to interpret.  An event says nothing about which other
    exons the transcript it comes from has, events overlap, one read counts
    toward several, and events aren't tied to an isoform or protein.  It
    doesn't take advantage of long-read technology: it only uses what a
    short read spanning the same junctions could show.  It has no
    differential test between conditions; its tables of counts and PSIs
    can be tested with your own statistics, and outlier samples can be found
    with the spliceevents_outliers script.  Unspliced reads only count toward
    starts and ends, as fragments of an annotated first or last exon.

Their outputs are only comparable at the end: FLAIR spliceevents' events
resemble those of FLAIR diffsplice, after transcriptome, combine and quantify,
but they are counted from reads rather than derived from isoforms, so the two
are separate analyses, not steps of one.

- Use isoforms when your reads are mostly full length and you want
  transcript structures, proteins, or isoform and gene expression.
- Use splice events counted from reads when your reads are fragmented or
  of uneven quality between samples, or when you want specific events or loci.
- The two can be used together, for example isoforms from the samples with
  the most complete reads and splice events from all of them, or isoforms
  for shorter genes and splice events just on tricky long genes

Generating transcriptomes
-------------------------

This is the core function of FLAIR. We have recently switched from using
two modules to assemble the transcriptome (correct and collapse) to
using FLAIR transcriptome. This single module allows parallelized
transcriptome building directly from the aligned BAM file, reducing
memory usage and increasing speed. FLAIR transcriptome is also more
sensitive to novel genes and transcripts than previous versions. FLAIR
transcriptome can be run with or without an annotated GTF file, but
having a good reference does improve performance. Here we provide a
standard recommended command, but FLAIR transcriptome has many options
that can be tuned for your application, see its dedicated page to learn
more.

- With short-read orthogonal support:

  - ``flair transcriptome -g GENOME.fa -f ANNOTATION.gtf --junction_tab SHORTREAD.SJ.out.tab -b SAMPLE.genomealigned.bam -o SAMPLE``

- With junctions from long reads:

  - ``flair transcriptome -g GENOME.fa -f ANNOTATION.gtf --junction_bed LONGREAD.IPjunctions.bed --junction_support 2 -b SAMPLE.genomealigned.bam -o SAMPLE``

Combine transcriptomes between samples
--------------------------------------

If you have a dataset composed of many samples that you want to compare,
you want to have a reference transcriptome that works equally well for
all samples. To achieve this, first generate a transcriptome for each
sample or batch of samples. Next, use FLAIR combine to combine the
transcriptomes.

- Make a manifest file pointing to your transcriptomes for all of your
  samples (see FLAIR combine documentation for more details)
- Run FLAIR combine (``-f ANNOTATION.gtf`` is optional, and gives each
  combined isoform its annotated gene and protein):

  - More stringent combination (keep only isoforms used at 5% or more of
    their locus in at least 1 sample, with at least 3 reads over all samples,
    the defaults)

    - ``flair combine --manifest MANIFEST.txt -g GENOME.fa -f ANNOTATION.gtf -o OUTPUT``

  - Less stringent combination (keeps anything supported in any file,
    just combines based on splice junctions and similar ends)

    - ``flair combine --manifest MANIFEST.txt -g GENOME.fa -f ANNOTATION.gtf -o OUTPUT --min_percent_usage 0 --min_reads 1``

Quantify splicing relative to the combined transcriptome
--------------------------------------------------------

For the most consistent transcript quantification, we recommend
re-quantifying all of your samples together on the combined
transcriptome using FLAIR quantify. FLAIR quantify is a stringent
quantification method. It assigns reads to unambiguous transcript
matches and does not support EM or partial support algorithms. This is
excellent for downstream approaches that require unambiguous
read-to-transcript matching and is a stringent quantification method,
but you will likely recover only 55-90% of your reads compared to
gene-level quantification.

- Make a new manifest file of all samples, a line per sample with its id,
  condition, batch, and its reads aligned to the genome, as a sorted and
  indexed BAM file, tab-separated (see FLAIR quantify documentation for
  more details). These are the same BAM files you generated transcriptomes
  from.

- Run FLAIR quantify:

  - ``flair quantify --manifest MANIFEST.tsv -g GENOME.fa --isoform_bed COMBINED.isoforms.bed -o OUTPUT``
  - Add ``--normalize_ends`` unless your reads' ends are reliable
    (``--trust_ends``).  The counts matrix names each isoform's gene in its own
    column, as FLAIR diffexp expects

Counting splicing events from reads
-----------------------------------

If you chose splicing events counted from reads, FLAIR spliceevents takes the
genome-aligned BAM files directly, with a manifest of each sample's name and
BAM file, tab-separated.  As with FLAIR transcriptome, include a junction
reference.

- ``flair spliceevents -m MANIFEST.tsv -g GENOME.fa -f ANNOTATION.gtf --junction_bed LONGREAD.IPjunctions.bed -o OUTPUT``
- For a few loci, add ``--region_bed REGIONS.bed``; only the regions' reads
  and annotation are used.
- To find samples whose PSI of an event is far from the others':
  ``spliceevents_outliers -i OUTPUT -o OUTPUT``
