import os
from flair import FlairInputDataError
import shutil
import pysam
import logging
from flair.partition_runner import PartitionRunner
from flair import SeqRange
from flair.junction_correct import junction_corrector_factory
from collections import namedtuple
from flair.io_utils import make_temp_dir, read_manifest_rows
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import annot_data_from_gtf
from flair.pycbio.hgdata.bed import BedReader
from flair.spliceevents_events import thresholds_of
from flair.spliceevents_genes import AnnotatedSplicing, gene_junctions, novel_genes_of_chains, process_region_events
from flair.spliceevents_reads import (annotated_terminal_exons, annotated_terminal_splice_sites,
                                      get_juncs_single_sample, read_chains_path, read_read_chains)
from flair.spliceevents_regions import (combine_regions, correction_regions, intergenic_cuts, partition_regions,
                                        read_windows, region_key)
from flair.spliceevents_tables import EventTables, combine_region_tables


# with --partition, the annotated exons to put in each piece, so a chromosome is cut
# by its number of genes, and the same way whatever the number of threads
PARTITION_EXONS = 5000


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
                          help='FastA of reference genome, can be minimap2 indexed')
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
    parser.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    parser.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from long-read alignment with intron-prospector)')
    parser.add_argument('--junction_support', type=int, default=2,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support '
                             '(default: %(default)s)')
    parser.add_argument('--ss_window', type=int, default=10,
                        help='window size for correcting splice sites (default: %(default)s)')
    parser.add_argument('--partition', choices=('yes', 'no'), default=None,
                        help='whether to cut the regions reads are corrected and events found in, in intergenic '
                             'gaps of the annotation, into pieces of about %d annotated exons, each owning the reads '
                             'that start in it.  By default, yes for whole chromosomes and no for --region_bed '
                             'regions' % PARTITION_EXONS)
    parser.add_argument('--region_bed',
                        help='bed file of the regions to report events in, such as the genes of interest; if not '
                             'specified, all chromosomes are used.  Reads overlapping a region are corrected over '
                             'their whole length, and the events reported are those inside a region: all their '
                             'junctions and exons, or, for the start, end, first and last exon events, which count '
                             'all of a gene\'s reads, all its annotated junctions')
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
    if not os.path.exists(args.genome):
        raise FlairInputDataError(f'Genome file path does not exist: {args.genome}')
    spliceevents(args)


def _correct_region_sample(*, partition, gtf_data, junction_corrector, args, item, terminal_exons, region_reports,
                           read_starts):
    """Correct one sample's reads in a region, for _run_region.  A task of its own,
    so the samples of a region are corrected in parallel"""
    sample, bamfile = item
    path = read_chains_path(partition, sample)
    if os.path.exists(path):
        return
    annots = annot_data_from_gtf(gtf_data, partition.region)
    if len(annots.transcripts) > 0:
        with pysam.FastaFile(args.genome) as genome:
            region = partition.region
            region_terminal = [exon for exon in terminal_exons.get(region.name, ())
                               if exon.start < region.end and region.start < exon.end]
            windows = read_windows(read_starts.get(region_key(region)), region_reports[region_key(region)])
            get_juncs_single_sample(args, region, path, bamfile, genome, junction_corrector, annots, region_terminal,
                                    windows)


def _run_region(*, partition, gtf_data, junction_corrector, args, allsamples, region_reports, annot_afe_ss, annot_ale_ss,
                read_starts):
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
                          AnnotatedSplicing(annots), read_starts.get(region_key(region)))
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


# the regions reads are corrected in, the report regions of each by its key, the range
# each one's reads start in when they are cut (by key), their numbers of reads when
# known, and the annotation, loaded for them
CorrectionPlan = namedtuple('CorrectionPlan', ('regions', 'reports', 'read_starts', 'weights', 'gtf_data'))


def partitioned(args):
    "whether to cut the regions: --partition, by default yes for whole chromosomes and no for --region_bed"
    if args.partition is None:
        return not args.region_bed
    return args.partition == 'yes'


def plan_correction(args, genome, report_regions, bamfiles):
    """The CorrectionPlan of the regions events are reported in: cut in intergenic gaps
    of the annotation, by default for whole chromosomes, else extended over their
    reads, for --region_bed, or whole chromosomes as they are"""
    if partitioned(args):
        if args.region_bed:
            # the annotation of the regions' whole chromosomes, to cut them in their gaps
            gtf_data = load_gtf(args.annot, chromosomes(genome, sorted({r.name for r in report_regions})))
        else:
            gtf_data = load_gtf(args.annot)
        cuts = {region_key(region): intergenic_cuts(region, gtf_data.transcripts, PARTITION_EXONS) for region in report_regions}
        plan = CorrectionPlan([], {}, {}, [], gtf_data)
        for region, owned, report, nreads in partition_regions(report_regions, bamfiles, cuts, args.threads):
            plan.regions.append(region)
            plan.reports[region_key(region)] = [report]
            plan.read_starts[region_key(region)] = owned
            plan.weights.append(nreads)
        return plan
    if args.region_bed:
        regions, reports = correction_regions(report_regions, bamfiles, args.threads)
        # only the annotation of the regions reads are corrected in
        return CorrectionPlan(regions, reports, {}, None, load_gtf(args.annot, regions))
    return CorrectionPlan(report_regions, {region_key(region): [region] for region in report_regions}, {}, None,
                          load_gtf(args.annot))


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

    with pysam.FastaFile(args.genome) as genome:
        if args.region_bed:
            report_regions = combine_regions([SeqRange(bed.chrom, bed.chromStart, bed.chromEnd) for bed in BedReader(args.region_bed)])
        else:
            report_regions = chromosomes(genome, genome.references)
        plan = plan_correction(args, genome, report_regions, bamfiles)
    logging.info(f'input regions: {len(plan.regions)}')

    logging.info('building intron support database')
    junction_corrector = junction_corrector_factory(args.ss_window, args.junction_support, annot_gtf_data=plan.gtf_data,
                                                    intron_beds=args.junction_bed, star_sj_tabs=args.junction_tab)
    terminal_exons, annot_afe_ss, annot_ale_ss = load_terminal_exons(args, plan.gtf_data,
                                                                     plan.regions if args.region_bed else None)

    logging.info(f'correcting reads of {len(allsamples)} samples in {len(plan.regions)} regions, then finding events, '
                 f'with {args.threads} threads')
    runner = PartitionRunner(plan.regions, temp_dir, gtf_data=plan.gtf_data, junction_corrector=junction_corrector,
                             threads=args.threads, weights=plan.weights)
    runner.run(_correct_region_sample, items=allsamples, args=args, terminal_exons=terminal_exons,
               region_reports=plan.reports, read_starts=plan.read_starts)
    runner.run(_run_region, args=args, allsamples=allsamples, region_reports=plan.reports,
               annot_afe_ss=annot_afe_ss, annot_ale_ss=annot_ale_ss, read_starts=plan.read_starts)

    combine_region_tables(args.output, [p.file_prefix for p in runner], [sample for sample, _ in allsamples])
    if not args.keep_intermediate:
        shutil.rmtree(temp_dir)
