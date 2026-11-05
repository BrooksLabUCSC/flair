import os
from flair import resolve_deprecated_option
import shutil
import pysam
import logging
from flair.partition_runner import PartitionRunner, run_flair_partition
from flair import SeqRange
from flair.junction_correct import junction_corrector_factory
from collections import namedtuple
from flair.io_utils import make_temp_dir, read_manifest_rows, check_genome_fasta
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import annot_data_from_gtf
from flair.spliceevents_events import thresholds_of
from flair.spliceevents_genes import AnnotatedSplicing, gene_junctions, novel_genes_of_chains, process_region_events
from flair.spliceevents_reads import (annotated_terminal_exons, annotated_terminal_splice_sites,
                                      read_chains_path, read_read_chains, sample_read_chains, write_read_chains)
from flair.spliceevents_regions import read_windows, region_key
from flair.region_spec import REGION_FORMAT_HELP, parse_regions, write_regions_bed
from flair.spliceevents_tables import EventTables, combine_region_tables


# the input items, reads of all the samples and annotated transcripts, a partition is
# grouped up to: enough that a task is worth its fork and its annotation subsetting,
# few enough to keep the threads fed.  On A498's chr19, 328k reads and 18k transcripts
# in 1140 clusters, a task per cluster costs 9% more CPU than 16 partitions
PARTITION_ITEMS = 20000

# shared by the --region_bed and --region help: what restricting the run to regions does
REGION_RESTRICT_HELP = (
    'If no region is given, all chromosomes are used.  A region is expanded over the reads and '
    'annotated transcripts that overlap it, until nothing more overlaps, and its reads are corrected '
    'over their whole length.  '
    'The events reported are those inside a region: all their junctions and exons, or, for the start, '
    'end, first and last exon events, which count all of a gene\'s reads, all its annotated junctions.  '
    'An event reaching past a region is left out, not reported from part of its reads.  The regions '
    'reads were corrected in are written to <output>.spliceevents.regions.bed')


def add_subparser(subparsers):
    desc = ('counts alternative splicing events, and their PSIs, of each sample directly from bam '
            'files of long-read RNA-seq alignments, for data too fragmented or variable in '
            'quality for flair transcriptome')
    parser = subparsers.add_parser('spliceevents', help='Count alternative splicing events from aligned reads',
                                   description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-m', '--manifest', required=True,
                          help='tab-separated manifest of sample name and path to its sorted and indexed bam '
                               'file aligned to the genome, a line per sample')
    required.add_argument('-g', '--genome', required=True, type=str,
                          help='FASTA of reference genome, can be minimap2 indexed')
    required.add_argument('-f', '--gtf', dest='annot', metavar='GTF', required=True,
                          help='GTF annotation file, used for correcting reads and assigning them to genes')
    parser.add_argument('-o', '--output', default='flair',
                        help='output file name base (default: %(default)s)')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='number of threads to run with (default: %(default)s)')
    parser.add_argument('--annot_basic', default='',
                        help='GTF annotation file of full-length transcripts, such as GENCODE basic, whose '
                             'first and last exons may also be used as internal ones in some reads and '
                             'still count as first and last exons, and which unspliced reads are fragments of, '
                             'rather than those of -f.  A comprehensive annotation includes truncated '
                             'transcripts, whose ends are not real first and last exons')
    parser.add_argument('--junction_star', help='splice junctions in STAR SJ.out.tab format, as STAR writes '
                                                'when it aligns reads')
    parser.add_argument('--junction_tab', help='deprecated name for --junction_star')
    parser.add_argument('--junction_bed', help='splice junctions as BED6 to BED9 with the number of supporting '
                                               'reads in the score column, as intron-prospector writes from a BAM '
                                               'of either short or long reads.  An alternative format to '
                                               '--junction_star, not a different kind of evidence')
    parser.add_argument('--junction_support', type=int, default=2,
                        help='minimum number of supporting reads to keep a junction, from either --junction_star '
                             'or --junction_bed; for bed the score column holds the count (default: %(default)s)')
    parser.add_argument('--ss_window', type=int, default=10,
                        help='window size for correcting splice sites (default: %(default)s)')
    parser.add_argument('--partition_items', type=int, default=PARTITION_ITEMS,
                        help='run the reads in partitions of at least this many input items, reads of all '
                             'the samples and annotated transcripts (default: %(default)s), or of a whole '
                             'chromosome with 0.  A partition is one or more clusters of items that overlap, '
                             'and a cluster is corrected and searched for events on its own, so this sets how '
                             'much work a task is and changes no result')
    parser.add_argument('--region_bed',
                        help='bed file of the regions to report events in, such as the genes of '
                             'interest.  ' + REGION_RESTRICT_HELP)
    parser.add_argument('--region', dest='regions', action='append', default=[],
                        help='region to report events in, with everything --region_bed says.  '
                             + REGION_FORMAT_HELP)
    parser.add_argument('--feature_support', type=int, default=2,
                        help='minimum number of reads of a junction, splice site or end of an event in some '
                             'sample (default: %(default)s)')
    parser.add_argument('--event_support', type=int, default=20,
                        help='minimum number of reads of an event in a sample for its PSI to be calculated '
                             '(default: %(default)s)')
    parser.add_argument('--event_frac_of_tot', type=float, default=0.1,
                        help='minimum (total reads in splicing event)/(total coverage over locus) to call event '
                             '(default: %(default)s)')
    parser.add_argument('--junc_frac_of_event', type=float, default=0.05,
                        help="how alternative a junction, splice site or end must be to make an event: its "
                             "reads' part of the event's reads at least this in some sample, and at most 1 minus "
                             "it in some sample.  Set to 0 for max recall (default: %(default)s)")
    parser.add_argument('--quality', type=int, default=1,
                        help='minimum mapping quality of genomic alignments (default: %(default)s)')
    parser.add_argument('--trust_strand', action='store_true',
                        help='trust the stranding of the input reads and do not attempt strand correction: '
                             'reads keep the strand of their alignment, and spliced reads are corrected only '
                             'with splice junctions on that strand, and are only assigned to transcripts '
                             'on that strand')
    parser.add_argument('--trust_junctions', action='store_true',
                        help='trust the strand of every input splice junction for read correction, without '
                             'checking splice motifs.  By default, an unannotated junction without a GT-AG '
                             'motif, such as GC-AG or AT-AC, is only weak evidence of strand: a read stranded '
                             'only by such junctions gets its strand from gene identification, or is dropped')
    parser.add_argument('--keep_supplementary', action='store_true',
                        help='keep supplementary alignments')
    parser.add_argument('--no_align_to_annot', action='store_true',
                        help="skip the initial alignment to the annotated sequences and correct reads only "
                             "from the genomic alignment. Slightly faster but less accurate when the annotation "
                             "is good")
    parser.add_argument('--keep_intermediate', action='store_true',
                        help='keep intermediate and temporary files for debugging')
    parser.set_defaults(entry=spliceevents_cmd)


def spliceevents_cmd(args):
    resolve_deprecated_option(args, deprecated="junction_tab", current="junction_star")
    check_genome_fasta(args.genome)
    spliceevents(args)


def _correct_cluster_sample(args, cluster, work_dir, bamfile, genome, junction_corrector, gtf_data,
                            terminal_exons, reports):
    """One sample's read chains in a cluster of overlapping reads and annotated
    transcripts.  Reads are aligned to the cluster's annotated transcripts, not the
    partition's, so a chain is the same however the clusters are grouped into
    partitions"""
    # annot_data_from_gtf keeps every transcript of the GtfData it is given, so the
    # annotation is cut to the cluster here: it is what reads are aligned to, and the
    # partition's would make a read's alignment depend on how many clusters it holds
    annots = annot_data_from_gtf(gtf_data.subset_for_region(cluster.name, cluster.start, cluster.end), cluster)
    if len(annots.transcripts) == 0:
        return []
    cluster_terminal = [exon for exon in terminal_exons.get(cluster.name, ())
                        if exon.start < cluster.end and cluster.start < exon.end]
    windows = read_windows(_overlapping_reports(cluster, reports))
    if len(windows.ranges) == 0:
        return []
    # the junctions of the cluster, as the partition's would make a read's correction
    # depend on how many clusters the partition holds
    cluster_juncs = junction_corrector.subset_for_region(cluster.name, cluster.start, cluster.end)
    return sample_read_chains(args, cluster, work_dir, bamfile, genome, cluster_juncs, annots,
                              cluster_terminal, windows)


def _correct_region_sample(*, partition, gtf_data, junction_corrector, args, item, terminal_exons,
                           region_reports, region_clusters):
    """Correct one sample's reads in a partition, a cluster at a time, for _run_region.
    A task of its own, so the samples of a partition are corrected in parallel"""
    sample, bamfile = item
    path = read_chains_path(partition, sample)
    if os.path.exists(path):
        return
    region = partition.region
    reports = region_reports[region_key(region)]
    chains = []
    with pysam.FastaFile(args.genome) as genome:
        for num, cluster in enumerate(region_clusters[region_key(region)]):
            chains += _correct_cluster_sample(args, cluster, f'{path}.{num}.correction', bamfile, genome,
                                              junction_corrector, gtf_data, terminal_exons, reports)
    write_read_chains(path, region.name, chains)


def _run_region(*, partition, gtf_data, junction_corrector, args, allsamples, region_reports,
                annot_afe_ss, annot_ale_ss):
    "Find the events of a region in its samples' corrected reads, from _correct_region_sample"
    region = partition.region
    annots = annot_data_from_gtf(gtf_data, region)
    samples = [sample for sample, _ in allsamples]
    if len(annots.transcripts) == 0:
        # no annotated genes, so no events: empty tables, as every partition's are combined
        with EventTables(partition.file_prefix, samples, args.event_support):
            return
    sample_chains = [read_read_chains(read_chains_path(partition, sample)) for sample in samples]
    novel_gene = novel_genes_of_chains(region.name, sample_chains)
    process_region_events(partition.file_prefix, region.name, samples, sample_chains, novel_gene, thresholds_of(args),
                          annot_afe_ss.get(region.name, {}), annot_ale_ss.get(region.name, {}),
                          region_reports[region_key(region)], gene_junctions(sample_chains, novel_gene, annots),
                          AnnotatedSplicing(annots))
    if not args.keep_intermediate:
        for sample in samples:
            os.remove(read_chains_path(partition, sample))


def read_manifest(path):
    "a manifest's samples, (name, bam file) each (io_utils.read_manifest_rows)"
    return [tuple(row) for row in read_manifest_rows(path, (2,), 'a sample name and its sorted, indexed bam file',
                                                     {1: 'bam file'})]


def load_gtf(path, regions=None):
    "an annotation's transcripts and exons, of regions, or all"
    logging.info(f'loading {path}')
    return gtf_data_parser(path, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES, regions=regions)


def chromosomes(genome, names):
    return [SeqRange(chrom, 0, genome.get_reference_length(chrom)) for chrom in names]


# the partitions reads are corrected in, the report regions of each by its key, the
# clusters it is corrected a cluster at a time (by key), the number of input items in
# each, and the annotation, loaded for them
CorrectionPlan = namedtuple('CorrectionPlan', ('regions', 'reports', 'clusters', 'weights', 'gtf_data'))


def _overlapping_reports(region, report_regions):
    """The parts of the regions events are reported in that lie in a region, which are
    also the reads to fetch for it.  Clipping leaves every event test as it was, as no
    read crosses a cluster, so an event of its reads is inside it; without it a whole
    chromosome's reads are fetched for each of its partitions"""
    return [SeqRange(region.name, max(report.start, region.start), min(report.end, region.end))
            for report in report_regions
            if (report.name == region.name) and (report.start < region.end) and (region.start < report.end)]


def _partition_of(clusters, counts):
    "a partition of the clusters grouped into it, its region and its input items"
    return (SeqRange(clusters[0].name, clusters[0].start, clusters[-1].end), list(clusters), sum(counts))


def _work_units(clusters, counts, items):
    """The clusters grouped into partitions of at least items input items each, or a
    partition per chromosome when items is 0.  A cluster is never split and is
    corrected on its own, so the grouping changes no result; it only sets how much work
    a task is"""
    limit = float('inf') if items <= 0 else items
    units, group, group_counts = [], [], []
    for cluster, count in zip(clusters, counts):
        if group and ((cluster.name != group[0].name) or (sum(group_counts) >= limit)):
            units.append(_partition_of(group, group_counts))
            group, group_counts = [], []
        group.append(cluster)
        group_counts.append(count)
    if group:
        units.append(_partition_of(group, group_counts))
    return units


def _plan_of(units, report_regions, gtf_data):
    """The CorrectionPlan of work units, dropping one that reaches no report region: an
    expanded region can reach a gene whose events are all outside the regions"""
    kept = [(region, clusters, count) for region, clusters, count in units
            if len(_overlapping_reports(region, report_regions)) > 0]
    return CorrectionPlan([region for region, _, _ in kept],
                          {region_key(region): _overlapping_reports(region, report_regions)
                           for region, _, _ in kept},
                          {region_key(region): clusters for region, clusters, _ in kept},
                          [count for _, _, count in kept], gtf_data)


def plan_correction(args, report_regions, bamfiles, restricted, temp_dir):
    """The CorrectionPlan of the regions events are reported in: flair_partition's
    clusters of overlapping reads and annotated transcripts, grouped into partitions of
    about --partition_items items to run in parallel.  With --region/--region_bed the
    regions are expanded over the reads and transcripts overlapping them first, which is
    what partitioning over them does"""
    region_bed = (write_regions_bed(report_regions, os.path.join(temp_dir, 'report-regions.bed'))
                  if restricted else None)
    clusters, counts = run_flair_partition(bamfiles, args.annot, args.threads, region_bed=region_bed)
    units = _work_units(clusters, counts, args.partition_items)
    return _plan_of(units, report_regions, load_gtf(args.annot, clusters if restricted else None))


def load_terminal_exons(args, gtf_data, regions):
    """The first and last exons unspliced fragments lie in, by chromosome, of
    --annot_basic, or else of the annotation, and of --annot_basic their splice sites,
    by chromosome and strand, which mark first and last exons that may also be
    internal ones"""
    if not args.annot_basic:
        return annotated_terminal_exons(gtf_data.transcripts), {}, {}
    terminal_exons = annotated_terminal_exons(load_gtf(args.annot_basic, regions).transcripts)
    return (terminal_exons, *annotated_terminal_splice_sites(terminal_exons))


def spliceevents(args):
    temp_dir = make_temp_dir(args.output)
    logging.info('temp directory: %s', temp_dir)
    allsamples = read_manifest(args.manifest)
    bamfiles = [bamfile for _, bamfile in allsamples]

    regions = parse_regions(args.regions, args.region_bed)
    with pysam.FastaFile(args.genome) as genome:
        report_regions = regions if regions is not None else chromosomes(genome, genome.references)
    plan = plan_correction(args, report_regions, bamfiles, regions is not None, temp_dir)
    write_regions_bed(plan.regions, args.output + '.spliceevents.regions.bed')

    logging.info('building intron support database')
    junction_corrector = junction_corrector_factory(args.ss_window, args.junction_support, annot_gtf_data=plan.gtf_data,
                                                    intron_beds=args.junction_bed, star_sj_tabs=args.junction_star)
    terminal_exons, annot_afe_ss, annot_ale_ss = load_terminal_exons(args, plan.gtf_data,
                                                                     plan.regions if regions is not None else None)

    logging.info(f'correcting reads of {len(allsamples)} samples in {len(plan.regions)} partitions of '
                 f'{sum(len(c) for c in plan.clusters.values())} clusters, then finding events, '
                 f'with {args.threads} threads')
    runner = PartitionRunner(plan.regions, temp_dir, gtf_data=plan.gtf_data, junction_corrector=junction_corrector,
                             threads=args.threads, weights=plan.weights)
    runner.run(_correct_region_sample, items=allsamples, args=args, terminal_exons=terminal_exons,
               region_reports=plan.reports, region_clusters=plan.clusters)
    runner.run(_run_region, args=args, allsamples=allsamples, region_reports=plan.reports,
               annot_afe_ss=annot_afe_ss, annot_ale_ss=annot_ale_ss)

    combine_region_tables(args.output, [p.file_prefix for p in runner])
    if not args.keep_intermediate:
        shutil.rmtree(temp_dir)
