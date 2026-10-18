# Major user-visible changes

## [v3.x.xx]
* General
  * fixed problems with some diffsplice_fishers_exact, and other
    auxiliary console script installs
  * flair now uses real subcommands, so `flair <subcommand> --help` and the
    subcommand list in `flair --help` describe what each one accepts.
  * Command line documentation is generated from the programs themselves, so it
    can no longer drift from what they accept.  flair variantquant, alleles and
    isoalleles are documented for the first time.
  * All subcommands now honor the logging options; previously only flair align did.
  * Bad input is reported as one line rather than a Python traceback: every
    installed command now reports errors this way, and the errors for mismatched
    annotation fasta and isoform names, alignments without AS or MD tags, and the
    option and input problems in junctions_from_sam and identify_vars say what to
    do about them.
  * `flair diffexp` gained `--condition_a` and `--condition_b`, matching
    `flair diffsplice`.  `condition_a` is the reference that fold changes are
    measured against.
  * The `diffexp` and `diffsplice` R dependencies are now part of
    `misc/flair_conda_env.yaml`, and `misc/flair_diffexp_conda_env.yaml` is gone.
    The BioConda package still does not carry them; `installing.rst` gives the
    `conda install` command to add them.
  * `flair fusion` rejects fusions between paralogs with a GENCODE annotation.
    The bundled paralog table has unversioned Ensembl gene ids, so the
    versioned GENCODE ids never matched it and the check did nothing.
  * `flair fusion` no longer rejects fusions of genes whose names share their
    first 4 letters, which threw out unrelated genes with a common stem, as the
    known MCF7 fusion BCAS4--BCAS3.  Genes of one name, as a gene's X and Y
    copies, are still rejected, and paralogs by the paralog table.
  * `flair fusion` gives each chimeric read piece the shortest gene sharing exons
    with the one it matched, so a readthrough annotation gives way to the gene it
    contains.  Reads of the MCF7 fusion MTAP--CDKN2B-AS1 were split between MTAP
    and a readthrough lncRNA starting with MTAP's exons, and the fusion was
    reported twice with part of its support each time.
  * `flair fusion` leaves out fold-back reads, chimeric reads with pieces on
    opposite strands within 1 kb of each other: cDNA folded back on itself, not
    a fusion.  They made fusions of genes with their antisense or neighbouring
    opposite-strand genes, as SLC7A5--SLC7A5-AS1 and EXOC7--ZACN in MCF7.
  * `flair fusion` accepts a non-canonical junction across a fusion breakpoint
    only if it carries at least half the reads crossing the breakpoint.  Any such
    junction with 2 reads used to do, so a few reads that happened to agree among
    many scattered ones made a fusion isoform, as TAOK3--RPS2 in MCF7 (3 of 21
    reads).  Canonical breakpoint junctions are unchanged.
  * `flair fusion --frac_support` sets the minimum fraction of a fusion's reads for
    a fusion isoform to be reported (default 0.05, as before).  A fusion with one
    dominant isoform hid its minor ones: in MCF7, BCAS4--BCAS3 has isoforms with
    a BCAS3 exon upstream of the usual junction exon, with 2% of its reads.
  * `flair fusion --quality` (default 10) leaves out chimeric reads with a piece
    below that MAPQ in the genome alignment, whose partner could as well be
    elsewhere, as the multi-copy subtelomeric DDX11L genes in MCF7.  It applies to
    the genome alignments only, not to reads aligned to annotated transcripts.
  * `flair fusion --min_dist_between_breakpoints` measured the gap between minus
    strand pieces from the far end of one to the far end of the other, so
    readthrough on the minus strand got through when the two pieces' lengths took
    it past the minimum, as TMEM64--LINC01030 in MCF7 (90 kb apart).
  * `flair fusion` assigns a chimeric read piece to a gene only when its aligned
    bases overlap the gene's exons.  A piece inside a gene's introns, or spliced
    across its exons without aligning to them, is now a non-genic locus: in MCF7,
    pieces inside RPS29's introns made it a hub with 26 partners.
  * `flair fusion` leaves out reads chimeric in the alignment to annotated
    transcripts but aligned to the genome in one piece, where the transcript
    alignment puts consecutive sequence on two overlapping or duplicated genes.
    Given all reads rather than only the chimeric ones, this removed the 14
    false fusions that the transcript alignment added in SeraCare and UPK3B,
    ANAPC1P2 and NPEPPSP1--TBC1D3L in MCF7.  A partner piece the genome
    alignment leaves clipped still counts (MYO6--SENP6 in MCF7).  The `-b` help
    now asks for all reads aligned with minimap2, and the
    `--min_dist_between_breakpoints` help explains finding fusions between
    close genes with a smaller minimap2 `-G`.
  * `flair fusion --max_dist_to_tss` now measures to the 5' gene's transcript
    starts upstream of a fusion's 5' end (or up to 500 bp downstream), not to any
    start: one downstream of the reads can't be their promoter.  The MCF7 calls are
    unchanged.  In SeraCare, whose reads are 1 kb fragments of longer spike-in RNAs,
    EGFR--SEPTIN14 had passed only through an EGFR start 1 kb downstream of its
    reads and is now rejected.  The help suggests a very high value or -1 for
    fragmented or 5' truncated reads: with -1, SeraCare finds all 16 spike-in
    fusions, ETV6--NTRK3 and EGFR--SEPTIN14 included.
  * `flair fusion --max_dist_to_tss` takes the nearer of two 5' ends of a fusion's
    reads: the most upstream one held by 2 or more reads, which with 5' truncated
    reads is the one near the promoter, and now also the median one.  In MCF7,
    ESR1--CCDC170's reads mostly start at an annotated ESR1 TSS, but 3 start 14 kb
    upstream, 19 kb from the nearest annotated start upstream of them.
  * `flair fusion --no_transcriptome_alignment` finds fusions from the chimeric
    reads of the genome alignment only, without aligning reads to annotated
    transcripts: about 5 minutes on MCF7 and SeraCare against 10 and 9, at a
    third of the memory, but it misses partner pieces the genome alignment leaves
    clipped, often short 5' exons (MYO6--SENP6 in MCF7 falls from 16 reads to 2).
    With the changes below it still finds all 19 MCF7 reference fusions.
  * `flair fusion --max_partners` (default 1000) sets the number of partners a
    gene can have in one role, 5' gene, 3' gene or between them, before its
    fusions are rejected as promiscuous.  It was fixed at 4, counted over all
    roles, which in MCF7 rejected SULF2, in the 20q13 amplicon, with all its
    partners: ARFGEF2--SULF2 (336 reads) and SULF2--PRICKLE2 were found only
    through the alignment to annotated transcripts, where SULF2 had fewer.  The
    help suggests lowering it for genomes without rearrangements.
  * `flair fusion` failed with an invalid GTF record on an annotation sorted by
    position: minus strand exons were put in genomic order by reversing their
    order in the file, which only works for GENCODE's transcript order.
  * The `-b` and `-o` help of `flair fusion` say that the chimeric read bams of an
    earlier run with the same `-o` are reused, whatever `-b` is.
  * `flair fusion -t` is passed on to flair transcriptome on the synthetic fusion
    genome, which used its own default of 12 threads whatever `-t` was.
  * `flair fusion` aligns to annotated transcripts only the reads that can show a
    fusion the genome alignment misses: unmapped, chimeric, or with
    `--min_fragment_size` or more clipped at an end: 11% of the reads in MCF7 and
    1.5% in SeraCare, where all were aligned.  The calls are identical, and the runs
    take 10 and 9 minutes against 24 and 37.  The genome bam is read once, for the
    chimeric reads and these together, where it was read three times.
  * `flair fusion` no longer writes OUTPUT.fusionreads.fa, the sequences of the
    final isoforms' reads: OUTPUT.fusion.isoform.read.map.txt names them, and
    OUTPUT.fusionreads.prelim.fa holds the sequences of every preliminary fusion read.
  * Without `--keep_intermediate`, `flair fusion` removes its intermediate files,
    where it removed only those of flair transcriptome: the annotated transcripts
    (1 GB with GENCODE), the preliminary fusion reads, and the synthetic fusion
    genome with its annotation, breakpoints and alignments.  It keeps the fusion
    isoforms, their read map, OUTPUT.prelimfusions.bed and the chimeric read bams,
    which a later run with the same `-o` reuses.
  * `flair fusion` writes OUTPUT.fusions.tsv, a table of the fusions, most final
    reads first: the genes, the reads of the preliminary fusion and of its final
    isoforms, the number of isoforms, and the best supported breakpoint and isoform
    with their reads.  A breakpoint is each junction between consecutive loci, as
    chrom:pos(strand)--chrom:pos(strand), comma-separated for 3 or more loci.
  * `flair fusion` found no fusions with a RefSeq or Ensembl reference: it took only
    chromosomes named chr... as nuclear, to leave out unplaced scaffolds, so rejected
    every fusion on NC_000012.12 or 12.  Chromosomes are now recognised as UCSC,
    GENCODE, RefSeq and Ensembl name them, still leaving out mitochondria and unplaced
    scaffolds, and non-genic loci by their chrom:start-end form rather than a chr
    prefix.  The docs say flair fusion is configured for human genomes.
  * `flair fusion` writes its isoforms as flair transcriptome does, all named
    OUTPUT.fusions.*: with the bed, fasta, GTF and bigBed, it now writes
    OUTPUT.fusions.isoform.counts.txt and OUTPUT.fusions.aaseq.tsv, the proteins
    of the productive isoforms, whose ids the bed's aaseq_id column now carries.
    The read map is renamed OUTPUT.fusions.isoform.read.map.txt (was
    OUTPUT.fusion.isoform.read.map.txt).
  * `flair fusion` drew a fusion piece the coding region doesn't reach as all
    coding, in the bed, GTF and bigBed: the 3' piece of an isoform whose stop
    codon is in its 5' piece, as 5 PTC isoforms in MCF7.
    Such a piece now has no coding region.
  * `flair fusion` gives each piece of a final fusion isoform its gene as flair
    transcriptome gives an isoform its gene, by the junctions and splice sites it
    shares with annotated genes, counting its splice sites at the fusion junctions,
    rather than taking the preliminary fusion's, which depended on where reads ended.
    A piece of only a first exon has its donor, and of genes a piece reads through,
    it has the one at its fusion junction.  A piece with no annotated junction or
    splice site gets the gene whose exons cover more than half of its own, and a
    piece in no gene's exons is a non-genic locus spanning the final pieces there,
    not the preliminary fusion's locus, which spanned all its reads and could be
    over 100 kb wider.  An isoform two of whose pieces are then of one gene is left
    out, as a rearrangement within a gene, which flair fusion doesn't call: a
    BPTF--BPTF call in SeraCare.  OUTPUT.fusions.tsv names fusions by their final
    isoforms' genes and non-genic loci.
  * `flair fusion` reads its GTF once, rather than four times, which takes 1 to 1.5
    minutes off a GENCODE run, but holds it in memory: 6 GB at its peak for SeraCare,
    from 3 GB.
  * `flair fusion` groups chimeric reads into preliminary fusions by their loci as
    well as their genes.  Reads of one fusion whose pieces were given different
    genes, by where the reads ended, were split into several preliminary fusions,
    whose synthetic contigs held the same sequence, and realigning divided the reads
    between them at random: 87 reads at MAPQ 0 in MCF7, now 22.  Preliminary fusions
    whose loci overlap, padded as their synthetic contigs are, are now one, named by
    the genes most of its reads were given.  The chimeric reads of the genome
    alignment and of the alignment to annotated transcripts are taken together,
    `--min_support` counting both, each read's pieces from one alignment: those of
    the alignment to transcripts if they find a gene the genome alignment missed, as
    a 3-gene fusion's short middle exon, and explain as much of the read.
  * `flair fusion` accepts a preliminary fusion whose reads mostly agree on pieces
    that abut at the fusion junction, though the rest overlap more: SMARCC1--EMCN
    in MCF7, where 14 of 36 reads' pieces abut and the others' SMARCC1 piece runs
    125 or 400 bases into the sequence their EMCN piece aligns.  Overlaps of more
    than 10 bases otherwise still reject a fusion: in MCF7 and SeraCare those are
    paralogs and pseudogenes the paralog table misses, template switches and RN7SL
    artifacts, and no reference fusions.
* `flair align` is removed.  It ran minimap2 and filtered its output, which no
  other module uses: flair transcriptome and flair fusion take a bam of all reads
  aligned to the genome.  Align with minimap2 as the docs show,
  `minimap2 -ax splice -s 40 -G 350k --MD --secondary=no`, with the options
  explained in the new docs section on aligning reads.  The tests align their
  reads with minimap2 directly, and
  `make test-alignments` checks their alignments, so a minimap2 that aligns them
  differently shows there first.
* flair transcriptome
  * Isoforms of the same coordinates on different sequences got the same isoform
    id: an isoform's name, by which ids are given, hashed its exons but not its
    chromosome.  In LRGASP WTC11 with spike-ins, single-exon ERCCs starting at 1
    shared ids (5 duplicated in PacBio, 4 in ONT).  The name now includes the
    chromosome, hashed so that names stay the same between runs.
  * About 3x faster: full WTC11 PacBio (3 million reads, GENCODE, 24 threads)
    went from 75 to 23 minutes and from 1,490 to 500 CPU minutes.  Setup before
    the partitions run is under 2 minutes (was 11); minimap2 threads stay within
    `--threads`, and the largest partitions run first.  Reads that match a
    first-pass isoform exactly are assigned without the final realignment, and
    reads no isoform can take are not realigned, which cuts the final
    realignment's minimap2 time 6-fold.  Reads are aligned to the annotation only
    when it can change their correction.
  * Isoform ends come from the reads assigned to each isoform in the final
    realignment, against transcripts padded so that reads running past an
    isoform's end are no longer rejected: more reads are assigned to isoforms.
    The densest cluster of read ends gives its mode.
  * Single-exon isoforms are reported at each clustered pair of read ends,
    whatever `--max_ends` is; single-exon isoforms that are 3' UTR fragments of a
    spliced transcript, and internally primed single-exon reads (an A-rich aligned
    end, or SQANTI3's 60% A after an untailed end), are removed.
  * A spliced isoform removed as a subset is reported after all when every
    isoform it is a subset of fails support and it has at least
    `--subset_backup_support` reads (default 10) covering its own sequence past
    what those isoforms contain; truncated reads of those isoforms are not its
    reads.  It keeps its own gene and reference transcript.
  * An isoform is a subset only of isoforms of its own gene, or of any gene's
    when its junctions reach past its gene's annotated junctions, as a fragment
    across a readthrough.  In WTC11, SRSF1's main transcript (1,042 reads) and
    BTG1-201 (88 reads) were removed as subsets of 1-read isoforms of
    overlapping readthrough genes, leaving their reads unassigned.
  * An isoform's unique sequence boundaries, past which its reads must reach,
    are those its own densest read ends reach past.  Isoforms retaining a
    terminal intron, whose reads had been discarded, are reported.
  * Isoforms are assigned to genes by junctions and splice sites, and novel genes
    are grouped by shared splicing.  Genes that match an isoform equally go to
    the one with more annotated transcripts, so a gene wins over a readthrough
    of it, which is usually annotated with one or two.
  * Semi-canonical junction motifs (GC-AG, AT-AC) are weak strand evidence:
    without an annotation, the strand of a read whose junctions all have them
    comes from gene identification.  `--trust_junctions` trusts every input
    junction's strand instead.
  * `--trust_strand` applies to spliced reads too; `--directRNA` declares direct
    RNA input (it implies `--trust_strand` and allows no large UTR deletions);
    `--total_rna` keeps single-exon reads that end at genomic A runs and needs no
    poly(A) tails (it requires `--trust_strand`).
  * Temporary files go in `--temp_dir`, by default $TMPDIR.
* flair combine
  * `-f/--gtf` gives each combined isoform its gene, reference transcript, CDS,
    productivity and protein from that annotation, as flair transcriptome would,
    whatever annotation or FLAIR version the samples were run with, and writes the
    proteins to `combined.aaseq.tsv`.  Samples whose transcriptomes gave the same
    junction chain different genes, which stopped combine, can be combined with it.
  * Without `--gtf`, isoforms that samples gave different genes or reference
    transcripts are an error that names the samples and their genes, and the
    samples' protein ids, which each sample numbers from 1, are dropped.
  * All the single-exon isoforms on a chromosome strand were combined into one
    isoform: they share the empty junction chain, whose end groups `--max_ends 1`
    merged.  Each group of single-exon isoforms with ends within `--end_window` is
    now its own isoform.  Isoforms on opposite strands are no longer combined.
  * Isoforms without an annotated gene are grouped into novel genes as flair
    transcriptome groups them, by shared splicing, rather than through whichever
    sample genes were seen first.  Isoform and gene ids are numbered in genome
    order, and the map files are sorted, so the output no longer depends on the
    order of the samples.
  * Isoforms of alleles, from flair isoalleles, are combined without their allele
    information, with a warning.
  * Fusion isoforms are combined only if they share their fusion junctions too,
    not only each piece's own junctions: two fusions of one 5' piece with 3'
    partners spliced at different acceptors, as ATP1A1 with ZFP64 and with a
    non-genic chr20 locus in MCF7, were taken for one isoform with two sets of
    genes, an error.
  * With `--gtf`, the pieces of fusion isoforms get their genes as flair fusion gives
    them, so samples that gave a fusion isoform different genes no longer stop
    combine, and their proteins, translated from flair fusion's coding region, are
    written to `combined.aaseq.tsv`, which left them out.
  * An isoform found in many samples no longer stops combine making
    `combined.isoforms.bb`: its sample list, like its genes, can be longer than 255
    characters, bedToBigBed's limit for a string.
* flair spliceevents
  * Is a flair subcommand, `flair spliceevents`, documented with the other modules;
    it was only runnable as `python flair_spliceevents.py`.  Its options are named
    as flair transcriptome's: `-f/--gtf` (was `--annot`), `--keep_supplementary`
    (`--keep_sup`) and `--no_align_to_annot` (`--noaligntoannot`).
    `--junc_support`, easily confused with `--junction_support`, is
    `--feature_support`: the reads a junction, splice site or end of an event needs.
    `--quality` sets the minimum mapping quality, which was fixed at 0.  Its defaults
    are flair transcriptome's: `--ss_window` 10 (was 15), `--junction_support` 2
    (was 1), `--quality` 1 and `--threads` 12 (was 4).
  * Corrects reads as flair transcriptome does, with the same code: the alignment to
    annotated transcripts (with transcriptome's reference, whose terminal exons reach
    their gene's furthest ends, and without aligning reads it couldn't change), then
    junction correction using the genome's splice motifs for strand.  It adds
    transcriptome's `--trust_strand` and `--trust_junctions`.  Single-exon reads are
    still left out.
  * With `--region_bed`, only events inside its regions are reported: those whose
    junctions and exons all lie in a region, and the start, end, first and last exon
    events of genes whose annotated junctions all lie in one.  A gene the region only
    clipped got events from only its reads overlapping the region.  Reads overlapping
    a region are corrected over their whole length: a region's annotation and
    junctions now cover the reads it fetches, so reads running past its end are not
    dropped.
  * A chromosome or region without annotated genes no longer stops the run when the
    regions' tables are combined.
  * A complex exon skip's (ces) inclusion counts every read including its exons.
    Reads counted before the first read skipping them were left out, so samples
    listed early in the manifest, or without skipping, had too low an inclusion,
    and some events were lost.
  * An event in only one sample with enough reads for a PSI is reported, if its PSI
    is at least `--junc_frac_of_event` from 0 and 1, as in a run of one sample.  It
    was left out, as its PSI couldn't differ from another sample's.
  * A complex exon skip's BED line of more than one included exon has the right
    number of blocks, which was always 3.  Its rows are written in the same order in
    every run.
  * Events and rows are named consistently, and every coordinate is 1-based.  An event
    is `TYPE|CHROM:START-END(STRAND)`, its locus being the skipped exon for es
    and otherwise the outer bounds of its introns (was, for example,
    `es-of-chr1:300-400`, without strand, or `alt3-relTo-chr1:200(+)-G1`).
    The tables have separate `event` and `alternative` columns (were `eventname`,
    the row joined to the event, and `event`), and `event_type` (was `eventtype`); a
    row's alternative is `inc`, `exc`, `retained`, `spliced` or `other`, the exons
    included of a ces row, or the introns, as STAR writes them, of an alternative
    splice site or first or last exon.  The coordinate columns are 1-based too, and
    AFE and ALE are afe and ale.  The BED lines are named `EVENT|ALTERNATIVE`.  spliceevents_outliers
    reads and writes the new columns.  Tables of earlier versions aren't read.
  * A tss or tts taken from the reads' starts (tss on +, tts on -) was placed one
    base before them.
  * `--output_read_ends` and its OUTPUT.spliceevents.readends.bed are gone.  Each
    sample's corrected reads are handed from correction to event calling as a line
    per junction chain with its reads' ends, rather than a line per read with its
    name, a sixth of the size.
  * Alternative first and last exons (afe, ale) are told from truncated reads over
    all samples together: an exon is a first exon only if at most 5% of the reads
    of its junctions have an exon before it (25% for a first exon of
    `--annot_basic`), and its reads are the event's only then.  An internal exon
    whose reads were all truncated in any one sample was taken for a first exon,
    and every truncated read of a gene was in its afe and ale totals.
  * Reads are given genes by flair transcriptome's gene identification, rather than
    spliceevents' own: by their junction chain, the junctions and splice sites they
    share with genes, then exon overlap, a gene's strand resolving a read of unknown
    strand.  Overlapping genes are told apart as transcriptome tells them, so on PC9
    reads of EEF1G and BSCL2 are no longer given an overlapping unnamed gene or a
    readthrough.  A read through two genes is given both, its events those of the
    gene `G1,G2`, and reads of no annotated gene, which were left out, are grouped
    over all samples into novel genes named for their locus, as transcriptome names
    them.
  * Events of splicing (es, ces, alt5, alt3, ir) are found among the junction
    chains linked by shared splice sites, not per gene, so overlapping genes and
    readthroughs no longer split an event's reads between genes.  A read's first
    exon and start (afe, tss) count with the gene it starts in, and its last exon
    and end (ale, tts) with the gene it ends in.  An event's PSItot is over all the
    strand's reads covering it.
  * Each splicing event is given its own gene: the one annotated with most of the
    introns of all its alternatives, each intron scoring 2 for a gene annotated with
    it, or failing that 1 for each of its splice sites the gene has, weighted by its
    share of the event's intron reads, a tie going to
    the gene of more annotated transcripts, or, sharing none, the gene of most of
    its reads.  Genes still tied are listed.  An event's gene had been
    the genes of all its reads, a list for 495 events of PC9's chr19.  The tables have a
    `novel_junctions` column, a row's introns that aren't annotated, a `gene_name`
    column, the annotation's name of the gene, and `additional_genes` and
    `additional_gene_names`, the other genes sharing an intron or splice site with a
    splicing event, as do the outliers.
  * Unspliced reads that are fragments of an annotated first or last exon count
    toward starts and ends (tss, tts): a read whose 5' end is in a last exon, or
    within 100 bases before it, adds its end to its gene's transcript ends, and one
    whose 3' end is in a first exon adds its start, with first and last exons of
    `--annot_basic` if given.  A read overlapping several genes' exons belongs to
    the one whose exons share the most bases with it, of its own strand with
    `--trust_strand`, else of either.  They were all left out.
  * Whole chromosomes are cut into pieces in intergenic gaps of the annotation, about
    5000 annotated exons each, and each piece's reads are corrected and its events
    found in parallel.  A piece owns the reads starting in it, however far they run,
    and reports the events starting in it, so an event near a cut can miss reads
    starting in the piece before.  Reads are aligned to the annotated transcripts of
    their piece rather than of their chromosome.  On PC9's chr19 this takes 8 minutes
    rather than 63, and 0.8 rather than 2.6 GB, and 95% of events are unchanged, the
    other 5% from the alignment to fewer transcripts: without it, 99.9% are.
    `--partition no` corrects whole chromosomes, and `--partition yes` cuts
    `--region_bed` regions, which are not cut by default.
  * The genome-wide annotation of first and last exons was pickled into every
    region's task, 1.7 seconds each, which kept the threads waiting on the main
    process.  A 1% sample of PC9 over the whole genome now takes under 5 minutes
    rather than 18.
  * `--annot_basic` first and last exons were matched by position on any
    chromosome, so a first exon's splice site counted as annotated if another
    chromosome had one at the same position.
  * spliceevents_outliers compares an outlier with those of the events sharing its
    junctions, not of its gene, and reads only tables with `novel_junctions`.
  * Event calling is faster: on PC9's chr19, 63 rather than 71 seconds.  Each read's
    chain was compared with every other chain of its gene, for a filter of
    truncated read ends that was no longer applied.
  * With `--region_bed`, only the annotation of the regions is loaded, from `-f` and
    `--annot_basic`: for one locus of GENCODE, in 10 seconds and 0.2 GB rather than
    75 seconds and 3.5 GB.
  * `--annot_basic` is read as a GTF, not converted to BED and cut per region with
    bedtools.  A gene id containing `_`, as in older GENCODE releases' `_PAR_Y`
    genes, lost its annotated first and last exons, as the gene was taken from the
    transcript name.
  * The samples of a region are corrected in parallel, within the `--threads` budget,
    not one after another: a run on one gene used a single thread.
  * Its outputs are named OUTPUT.spliceevents.*, not OUTPUT.diffsplice.*, which
    were easily taken for flair diffsplice's.  The counts and PSI tables have an
    `event` column, the event each row belongs to, and OUTPUT.spliceevents.totals.tsv
    has each event's reads and the reads covering its locus per sample, the
    denominators of PSIjunc and PSItot.
  * Outlier calling moved out of flair spliceevents into the console script
    `spliceevents_outliers`, which reads OUTPUT.spliceevents.counts.tsv and
    totals.tsv, so its outliers can be called of either PSI (`--psi junc` or `tot`)
    without rerunning spliceevents.  `--check_outliers` is removed.  A sample's PSI
    is used when it has at least `--event_support` reads, as in the PSI tables,
    rather than more than that.
* Incompatibles
  * Removed flair correct and collapse modules, the functionality is replaced
    by flair transcriptome.
  * Remove flair align options that are no longer need without flair correct.
  * Logging options (--log-level, --log-stderr, --log-conf, --log-debug) and
    --version must now precede the subcommand: `flair --log-debug align ...`.
  * Options renamed for one name per concept and one meaning per short option.
    Old names are not accepted.
    * everywhere: --out_dir, --output_prefix to --output; --bed_paths to
      --isoform_beds; --bedisoforms and diffsplice --isoforms to --isoform_bed;
      --out_dir_force to --overwrite_output; --norm_ends to --normalize_ends
    * align: --filtertype to --filter_type; --minfragmentsize to
      --min_fragment_size; --maxintronlen to --max_intron_len; --nvrna to
      --native_rna
    * transcriptome: --se_support to --single_exon_support; --keep_sup to
      --keep_supplementary
    * combine: --endwindow to --end_window; --minpercentusage to
      --min_percent_usage; --remove_se to --remove_single_exon
    * fusion: --support to --min_support; --maxloci to --max_loci;
      --max_dist_to_TSS to --max_dist_to_tss; --min_dist_between_bp to
      --min_dist_between_breakpoints
    * diffexp: --exp_thresh to --min_expression
    * diffsplice: --drim1 to --drim4 to --min_samps_gene_expr,
      --min_samps_feature_expr, --min_gene_expr, --min_feature_expr, after the
      DRIMSeq dmFilter arguments they are passed to; --conditionA and
      --conditionB to --condition_a and --condition_b
    * variantquant: --input_bam to --transcriptome_bam; --threshold to
      --min_coverage
    * alleles and isoalleles: --bam to --tumor_bam; --norm_bam to --normal_bam;
      --norm_vcf to --normal_vcf; --allele_read_map_norm and --iso_read_map_norm
      to --allele_read_map_normal and --iso_read_map_normal
  * Short options dropped where a letter meant two different things: -w, -p, -q,
    -e, -s, -i, -k, -m, -v and -of.  -b is now always the aligned BAM, -r the
    reads, -t threads, -g the genome, -f the annotation and -o the output.
  * Removed options that were accepted and never used: flair align --quiet and
    flair quantify --quality.
  * Removed `--allow_paralogs` from flair transcriptome and fusion, which only
    reported that it was not implemented.  A read aligning equally well to
    several transcripts, as to identical paralogs, is assigned to one of them, now
    the same one for every such read (by transcript name) rather than whichever
    minimap2 reported first, which could split a paralog pair's reads between
    them.
  * flair align now requires exactly one of --genome and --mm_index, rather than
    failing later when neither was given.
  * junctions_from_sam, fasta_seq_lengths, mark_intron_retention,
    identify_annotated_gene and diffsplice_fishers_exact now use standard option
    parsing and accept --help; bed_to_gtf --noCDS is now --no_cds.
  * annotate_aaseq_with_uniprot is now installed as a command.
  * `flair diffexp` no longer takes the two conditions from the first and last
    column of the counts matrix.  With `--condition_a` and `--condition_b` left
    out, the two conditions are used in sorted order.  For a counts matrix whose
    first and last columns are not the conditions in sorted order, the reference
    changes, which flips the sign of the reported fold changes and renames the
    output files `prefix_X_v_Y.tsv`.  Results are not comparable across this
    change; name the conditions to get a specific direction.
  * `flair diffexp` previously took the two conditions from the first and last
    column even when they were the same condition, as in a column order A,B,B,A,
    which filtered one condition twice and ignored the other.  Such runs gave
    wrong results and now give correct ones.
  * `diff_iso_usage` output columns are named after the two samples given rather
    than sample1 and sample2.
  * `flair quantify` writes `<output>.sample_info.tsv` beside the counts file,
    giving the condition and batch of each counts column, and the counts columns
    are now named after the sample alone rather than `sample_condition_batch`.
    `flair diffexp` and `flair diffsplice` read that file.  A counts matrix from
    an earlier FLAIR is still read by taking the fields out of the column names.
  * Because condition and batch are columns of their own, the id, condition and
    batch fields of the quantify manifest may now contain underscores.
  * `flair quantify --sample_id_only` is gone; the counts columns always name the
    sample alone, which is what it asked for.
  * `flair transcriptome --trust_ends` is removed: ends now come from the reads
    assigned to each isoform.
  * `flair transcriptome --filter bysupport` is replaced by `--filter <N>X`, which
    keeps a subset isoform with more than N times its supersets' reads;
    `bysupport` was `1.2X`.  The default is still `nosubset`.
  * `flair transcriptome --normalize_ends` now reports, for each spliced isoform,
    the ends it shares with the isoforms of its gene whose terminal splice sites
    are near its own.
  * `flair transcriptome --ss_window` defaults to 10 (was 15).
  * `count_sam_transcripts` no longer takes `--quality`, `--fusion_dist` or
    `--end_norm_dist`.
    
## [v3.0.0] 2025-11-31
* General
  * Bug fixes since v3.0.0b1

## [v3.0.0b1] 2025-11-27 (beta 1)
* General
  * Removed dependencies on pandas and rpy2.
  * Added proper logging throughout
* FLAIR align
  * Cleaned up output files - now outputs unfiltered BAM and filtered BED
* FLAIR correct
  * Changed orthogonal junction inputs - now specify whether your file is 
  a bed (--junction_bed) or whether it comes from STAR short-read RNA alignment (--junction_tab)
  * Now support detecting junctions directly from long read data: 
  run intronProspector to generate junction bed, input that via --junction_bed, 
  specify desired junction read support with --junction_support
  * Removed -g genome option, increases speed
  * Now produces sorted output files
* FLAIR collapse
  * cleaned up output files, now only gives: isoforms.bed, isoforms.gtf, isoforms.fa, isoform.read.map.txt
  * (experimental) added filter for removing reads with internal priming. 
  Options: --remove_internal_priming, --intprimingthreshold, --intprimingfracAs
  * allow CDS prediction in collapse with --predictCDS
  * removed genomic range option (too fragile). For parallelization, run FLAIR transcriptome
  * made gtf (annotation) parsing more robust, especially for unsorted files 
  * improved fractional support filtering (with --support < 1)
  * improved isoform haplotyping through longshot - this is being deprecated though, please use FLAIR variants
* FLAIR quantify
  * improved fraction of reads assigned to isoforms and recovered (better handling of ambiguous alignments)
  * improved processing of multiple samples so running requires less available memory and is faster
* FLAIR diffexp and diffsplice
  * recoded directly in R to improve performance
* New modules
  * FLAIR transcriptome
    * Combines the functions of correct and collapse
    * Runs directly from an aligned BAM file
    * Performs more effective parallelization (specify in --parallelmode)
  * FLAIR fusion
    * Detects gene fusions and fusion isoforms
    * Fusion detection accuracy is comparable to JAFFAL and CTAT-lr-fusion
    * Fusion isoform detection is much more accurate
  * FLAIR combine
    * Allows combining transcriptomes generated from different samples
    * Fusion isoforms can be combined with collapsed isoforms to form a full transcriptome
    * Manipulate filtering of single exon isoforms with --include_se
  * FLAIR variants
    * Allows detection of variant-aware transcripts through 
    identification of read clusters with shared variants
    * Uses variants detected either from WGS or from lr-RNA-seq. 
    We recommend Longshot.
    * Good for identifying splicing of specific variants
    * Can be used in conjunction with FLAIR diffexp or diff_iso_usage 
    to identify changes in variant-aware transcripts between samples/groups
* Added FLAIR protocol to improve guidance on running flair


## [v2.2.1] 2025-04-06
* Fixed case where --quality=0 default didn't work, as it tested for greater
  than the cutoff.
* Converted all Python code to use 4-space indentation.  Previously some code
  used tabs, other used 4-spaces.
* Added missing dependencies for plot_isoform_usage.


## [v2.2.0] 2025-05-06
* Returned diffexp and diffsplice as standard modules.  The BioConda
  environment does not include the dependencies for these modules
  and required software does not run on Apple Silicon (ARM64) systems.
* The flair combine functionality is now to a module.  The
  `flair_combine` program is run with `flair combine`.
* Changed default MAPQ minimum quality score to 0. This allows more reads to
  be used in identifying isoforms, which tends to improve the overall models
  with out adversely affecting the accuracy.
* GitHub releases include the Conda YAML files for building FLAIR
  environments.  Useful if the BioConda release has not been manually
  reviewed.
* The FLAIR Docker now includes all dependencies to run diffexp and diffsplice.
* Reorganized the installation documentation.
* Fixed flair_quantify --output_bam crash
* Other bug fixes

## [v2.1.2] 2025-04-17
* Address issue getting BioConda to work
* Bug fixes for collapse command line parsing


## [v2.1.1] 2025-04-10
* converted all programs to use console scripts to allow BioConda to work


## [v2.1.0] 2025-03-27
* Numerous bug fixes.
* Removed support for PSL format.
* Remove `flair 123` to run multiple modules at once.
* Compatibility with Python 3.12 
* Compatibility with Apple ARM64 systems.
* Deprecated `diffExp` and `diffSplice`, they will be removed in a future release.
  Lets us know if you use this functionality.  Their dependencies are no longer
  part of the conda package, they can be added the conda environment with
  `misc/flair_diffexp_conda_env.yaml`.
