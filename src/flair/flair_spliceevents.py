import os
from bisect import bisect_left, bisect_right
import os.path as osp
from flair import FlairInputDataError
import shutil
import pysam
import logging
from flair.partition_runner import PartitionRunner, combine_temp_files_by_suffix
from flair import SeqRange
from statistics import median
from flair.junction_correct import junction_corrector_factory
from collections import namedtuple, Counter
from flair.io_utils import make_temp_dir
from flair.read_processing import generate_genomic_alignment_read_to_clipping_file
from flair.read_correction import filter_correct_group_reads
from flair.annotation_alignment import identify_good_match_to_annot
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import annot_data_from_gtf
from flair.pycbio.hgdata.bed import Bed, BedBlock, BedReader
from flair.isoform_data import Isoform, Junc
from flair.junction_correct import UNKNOWN_STRAND
from flair.gene_identification import get_gene_name_firstpass, get_gene_name_unknown_strand, novel_genes, group_by_splicing
from flair.interval_index import IntervalIndex


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
    parser.add_argument('--junc_frac_of_event', type=float, default=0.02,
                        help='minimum (reads with junction)/(total reads in splicing event). Set to 0 for max recall '
                             '(default: %(default)s)')
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


# the thresholds of event calling, from the command line
Thresholds = namedtuple('Thresholds', ('feature_support', 'event_frac_of_tot', 'junc_frac_of_event', 'event_support'))


def thresholds_of(args):
    return Thresholds(args.feature_support, args.event_frac_of_tot, args.junc_frac_of_event, args.event_support)


# Reads are counted per sample in Counters, sample -> reads, without the samples of none


def add_reads(counts, key, sample, numreads):
    "add a sample's reads to counts[key]"
    if key not in counts:
        counts[key] = Counter()
    key_counts = counts[key]
    key_counts[sample] = key_counts.get(sample, 0) + numreads


def add_counts(counts, key, more):
    "add more, reads per sample, to counts[key]"
    if key not in counts:
        counts[key] = Counter()
    counts[key].update(more)


def total(all_counts):
    "the sum of reads per sample"
    result = Counter()
    for counts in all_counts:
        result.update(counts)
    return result


def psi_varies(psis, junc_frac_of_event):
    """Is a junction, splice site or end of an event alternative, from its PSIs in the
    samples with enough reads of the event: differing across them by
    junc_frac_of_event, or, in the only sample with a PSI, at least junc_frac_of_event
    from both 0 and 1, so an event in one sample, as in a small dataset, is reported"""
    if len(psis) == 0:
        return False
    if len(psis) == 1:
        return junc_frac_of_event <= psis[0] <= 1 - junc_frac_of_event
    return max(psis) - min(psis) >= junc_frac_of_event


def sample_psis(counts, totals, samples, min_total):
    "the PSIs, counts over totals, of the samples with at least min_total reads in totals"
    return [counts[s] / totals[s] for s in samples if totals[s] >= min_total]


def junc_between(ss1, ss2):
    "the junction between two splice sites"
    return (min(ss1, ss2), max(ss1, ss2))


def lies_left(upstream, strand):
    """whether a feature lies left of another, at lower coordinates, on strand: if
    upstream (5') of it on the + strand, or downstream (3') of it on the - strand"""
    return strand == '+' if upstream else strand == '-'


def internal_exons(juncs):
    "the exons between a chain's junctions"
    return [(juncs[i][1], juncs[i + 1][0]) for i in range(len(juncs) - 1)]


def event_name(etype, chrom, start, end, strand):
    """an event's name: its type and outer bounds, start 0-based, written 1-based as a
    genome browser's locus.  Not its gene, which may be several, or decided only once
    the event is found"""
    return f'{etype}|{chrom}:{start + 1}-{end}({strand})'


def row_name(ename, alternative):
    "an event row's name, in BED records: the event's and the row's alternative"
    return f'{ename}|{alternative}'


def coord_ranges(ranges):
    "introns or exons, (start, end) with 0-based starts, written 1-based"
    return ','.join(f'{start + 1}-{end}' for start, end in sorted(ranges))


def span(ranges):
    "the outer bounds of introns or exons"
    return min(r[0] for r in ranges), max(r[1] for r in ranges)


def junction_blocks(junc, exons=()):
    "the BED blocks of a junction: 10 bases each side of it, and the exons included between its ends"
    return [(junc[0] - 10, junc[0]), *exons, (junc[1], junc[1] + 10)]


# the BED colors of each event type: first exons as starts, last exons as ends, and
# the sites sharing alternative ones as the alternative ones
EVENT_COLORS = {'es': '66,105,245', 'ces': '66,105,245', 'ir': '144,66,245',
                'alt3': '245,215,66', 'alt3ss': '245,215,66', 'alt5': '43,184,39', 'alt5ss': '43,184,39',
                'tss': '255,0,0', 'afe': '255,0,0', 'tts': '0,0,255', 'ale': '0,0,255'}


# an event's row: its reads per sample, its included and excluded junctions (or ends),
# the skipping junctions of an exon skip, its exons, and its BED blocks, a list per
# BED record
EventRow = namedtuple('EventRow', ('counts', 'inc_juncs', 'exc_juncs', 'outer_juncs', 'exons', 'blocks'))


class SplicingEvent:
    """An event: its rows, EventRows by alternative, its reads per sample (totjunc),
    the reads covering its locus (totoverlap), and the reads of its other row, those
    of the alternatives too rare or invariant to be rows"""
    def __init__(self, name, eventtype, chrom, strand, totjunc, gene=''):
        self.name = name
        self.eventtype = eventtype
        self.chrom = chrom
        self.strand = strand
        self.gene = gene
        self.totjunc = totjunc
        self.totoverlap = None
        self.rows = {}
        self.other = Counter()
        # all the introns of its alternatives, those of other too, for its gene
        self.introns = set()
        # the other annotated genes sharing an intron or splice site with it
        self.additional_genes = []

    def add(self, alternative, counts, inc_juncs, exc_juncs, outer_juncs=(), exons=(), blocks=()):
        self.rows[alternative] = EventRow(counts, set(inc_juncs), set(exc_juncs), set(outer_juncs), set(exons), list(blocks))

    def keep_other(self, other, min_reads):
        "count the reads of the event's other junctions or sites, if some sample has min_reads of them"
        if any(n >= min_reads for n in other.values()):
            self.other = other


class EventCaller:
    """What finding the events of a strand of a region needs: its chromosome, strand,
    samples, Thresholds and ReadCoverage, with the tests of an event's alternatives"""
    def __init__(self, chrom, strand, samples, thresholds, coverage):
        self.chrom = chrom
        self.strand = strand
        self.samples = samples
        self.thresholds = thresholds
        self.coverage = coverage

    def supported(self, counts, min_reads=None):
        "do some sample's reads reach min_reads, by default feature_support"
        if min_reads is None:
            min_reads = self.thresholds.feature_support
        return any(counts[s] >= min_reads for s in self.samples)

    def varies(self, counts, totals):
        "does the PSI of counts, of totals, vary over the samples with event_support reads (psi_varies)"
        return psi_varies(sample_psis(counts, totals, self.samples, self.thresholds.event_support),
                          self.thresholds.junc_frac_of_event)

    def some_psi_between(self, counts, totals):
        "is counts, of totals, between junc_frac_of_event and 1 minus it in some sample"
        frac = self.thresholds.junc_frac_of_event
        return any(frac <= counts[s] / totals[s] <= 1 - frac for s in self.samples if totals[s] > 0)

    def split(self, alternatives, totals, need_support=True):
        """the alternatives, key -> reads per sample, of an event, which vary and, with
        need_support, are supported, and the reads of the rest, its other row"""
        good, other = [], Counter()
        for key, counts in alternatives.items():
            if (not need_support or self.supported(counts)) and self.varies(counts, totals):
                good.append(key)
            else:
                other += counts
        return good, other

    def event(self, eventtype, start, end, totals):
        return SplicingEvent(event_name(eventtype, self.chrom, start, end, self.strand), eventtype, self.chrom,
                             self.strand, totals)

    def keep_other(self, event, other):
        event.keep_other(other, self.thresholds.feature_support)

    def covering(self, ref):
        "the reads per sample covering ref"
        return self.coverage.covering(ref)


TERMINAL_JUNCTION_DIST_FROM_ANNOT = 20
# the most of a first or last exon's junction reads, pooled over the samples, that
# have an exon before (or after) it: truncated reads can't add an exon, so a real
# first or last exon has none, while an internal exon whose reads in one sample are
# all truncated has them in the others
MAX_NONTERMINAL_FRACTION = 0.05
# for a first or last exon of --annot_basic, which may also be used as an internal one
MAX_NONTERMINAL_FRACTION_ANNOTATED = 0.25
MIN_TERMINAL_JUNCTION_SEPARATION = 100
MIN_TERMINAL_SS_FRACTION = 0.1
# read ends within this of the next are one start or end
END_GROUP_WINDOW = 100


# a cluster's splicing, reads per sample: of each junction, of each 3' splice site by
# its 5' site and each 5' site by its 3' site, of each pair of consecutive junctions,
# of each block, an internal exon or a read's first or last, and of each skipping
# junction by the junctions of the reads not skipping between its ends
ClusterSplicing = namedtuple('ClusterSplicing', ('juncs', 'ss5to3', 'ss3to5', 'junc_pairs', 'blocks', 'skipped'))


def extract_splicing_info(samples, cluster_chains, strand):
    "the ClusterSplicing of a cluster's chains, junctions -> [ReadEnds] per sample"
    info = ClusterSplicing({}, {}, {}, {}, {}, {})
    for sample, chains in zip(samples, cluster_chains):
        for juncs, ends in chains.items():
            numreads = len(ends)
            for junc in juncs:
                ss5, ss3 = junc if strand == '+' else junc[::-1]
                add_reads(info.ss5to3.setdefault(ss5, {}), ss3, sample, numreads)
                add_reads(info.ss3to5.setdefault(ss3, {}), ss5, sample, numreads)
                add_reads(info.juncs, junc, sample, numreads)
            for prev, nxt in zip(juncs, juncs[1:]):
                add_reads(info.blocks, (prev[1], nxt[0]), sample, numreads)
                add_reads(info.junc_pairs, (prev, nxt), sample, numreads)
            for read in ends:
                add_reads(info.blocks, (read.start, juncs[0][0]), sample, 1)
                add_reads(info.blocks, (juncs[-1][1], read.end), sample, 1)
    # once all the junctions are known: a chain counted before the skipping junction
    # was seen, as in a sample listed before the samples with skipping, was left out
    for sample, chains in zip(samples, cluster_chains):
        for juncs, ends in chains.items():
            add_skipped(info, juncs, sample, len(ends))
    return info


def add_skipped(info, juncs, sample, numreads):
    """count a chain's reads for each skipping junction from one of its junctions' start
    to another's end, by its junctions between them"""
    starts, ends = {j[0] for j in juncs}, {j[1] for j in juncs}
    for outer in info.juncs:
        if outer not in juncs and outer[0] in starts and outer[1] in ends:
            inner = tuple(j for j in juncs if outer[0] <= j[0] and j[1] <= outer[1])
            add_reads(info.skipped.setdefault(outer, {}), inner, sample, numreads)


def _flanking(outer, exon):
    "the junctions each side of exon inside the skipping junction outer"
    return (outer[0], exon[0]), (exon[1], outer[1])


def exon_skipping_events(caller, info):
    """Skipped exons (es): each internal exon, included by its flanking junctions and
    excluded by the junctions joining their outer ends.  Also the reads of each
    junction flanking an included exon by its skipping junction, which alt5 and alt3
    events count apart, and each event's junctions, which ces events don't repeat"""
    outers, inc, exc = {}, {}, {}
    for (prev, nxt), counts in info.junc_pairs.items():
        outer, exon = (prev[0], nxt[1]), (prev[1], nxt[0])
        outers.setdefault(exon, set()).add(outer)
        add_counts(inc, exon, counts)
        add_counts(exc, exon, info.juncs.get(outer, {}))
    events, es_juncs, seen = [], {}, set()
    for exon, exon_outers in outers.items():
        tot = inc[exon] + exc[exon]
        if not (caller.supported(inc[exon]) and caller.supported(exc[exon]) and caller.varies(inc[exon], tot)):
            continue
        # the skipping junctions whose inclusion is a large enough part of the event
        good = {outer for outer in exon_outers if caller.supported(info.junc_pairs[_flanking(outer, exon)])
                and caller.some_psi_between(info.junc_pairs[_flanking(outer, exon)], tot)}
        flanking = {j for outer in good for j in _flanking(outer, exon)}
        event = caller.event('es', *exon, tot)
        events.append(event)
        event.introns = {j for outer in exon_outers for j in _flanking(outer, exon)} | (exon_outers & info.juncs.keys())
        event.add('inc', inc[exon], flanking, (), good, [exon], [junction_blocks(outer, [exon]) for outer in good])
        event.add('exc', exc[exon], (), flanking, good, [exon])
        for outer in good:
            for j in _flanking(outer, exon):
                add_counts(es_juncs.setdefault(j, {}), outer, info.junc_pairs[_flanking(outer, exon)])
        ref = good or exon_outers
        event.totoverlap = caller.covering((max(j[0] for j in ref), min(j[1] for j in ref)))
        seen.add(frozenset(good | flanking))
    return events, es_juncs, seen


def complex_skipping_events(caller, info, es_juncs, seen):
    """Exons skipped together (ces): each skipping junction, and the junctions between
    its ends of the reads not skipping, unless an es event has the same junctions.
    Adds the reads of the junctions flanking the included exons to es_juncs"""
    events = []
    for outer, inner_counts in info.skipped.items():
        exc = info.juncs[outer]
        tot = total([exc, *inner_counts.values()])
        if not (caller.supported(exc) and caller.varies(exc, tot)):
            continue
        good, other = caller.split(inner_counts, tot)
        inner_juncs = set().union(*good)
        if not good or frozenset(inner_juncs | {outer}) in seen:
            continue
        event = caller.event('ces', *outer, tot)
        events.append(event)
        event.introns = {outer}.union(*inner_counts)
        exons = {exon for inner in good for exon in internal_exons(inner)}
        event.add('exc', exc, (), inner_juncs, [outer], exons, [junction_blocks(outer)])
        for inner in good:
            inner_exons = internal_exons(inner)
            event.add(coord_ranges(inner_exons), inner_counts[inner], inner, inner_juncs - set(inner), [outer], inner_exons,
                      [junction_blocks(outer, inner_exons)])
            for j in inner:
                add_counts(es_juncs.setdefault(j, {}), outer, inner_counts[inner])
        caller.keep_other(event, other)
        event.totoverlap = caller.covering(outer)
    return events


def _apart_from_skipping(caller, site, alt, alternatives, es_juncs):
    """is the junction of site and alt its own alternative: not flanking a skipped exon,
    or with reads apart from the skipping, over 10% of them and feature_support, in some
    sample"""
    junc = junc_between(site, alt)
    if junc not in es_juncs:
        return True
    skipping = total([es_juncs[junc][outer] for outer in (junc_between(site, other) for other in alternatives if other != alt)
                      if outer in es_juncs[junc]])
    counts = alternatives[alt]
    return any(counts[s] > 0 and (counts[s] - skipping[s]) / counts[s] > 0.1
               and counts[s] - skipping[s] > caller.thresholds.feature_support for s in caller.samples)


def alternative_site_events(caller, eventtype, site_alternatives, es_juncs, terminal_juncs):
    """Splice sites alternative to each other, spliced to one site (alt5 or alt3, by
    eventtype), from site_alternatives, each site's reads by the site it is spliced
    to.  A junction flanking a skipped exon counts only with reads apart from the
    skipping (es_juncs), and the junctions of a first or last exon event,
    terminal_juncs, aren't repeated.  Also each group of alternative sites, by the
    sites spliced to all of them, with their reads, for shared_site_events"""
    events, groups = [], {}
    for site, alternatives in site_alternatives.items():
        if len(alternatives) <= 1:
            continue
        tot = total(alternatives.values())
        good, other = caller.split(alternatives, tot, need_support=False)
        juncs = {junc_between(site, alt) for alt in good}
        if frozenset(juncs) in terminal_juncs:
            continue
        if len(good) > 1:
            groups.setdefault(frozenset(good), {})[site] = tot
        # only if alternative splicing, not counting the reads of junctions skipping exons
        own = [alt for alt in good if _apart_from_skipping(caller, site, alt, alternatives, es_juncs)]
        if len(own) <= 1 or not caller.supported(tot, caller.thresholds.event_support):
            continue
        # from the shared site to the furthest alternative one, so unique to the shared site
        start, end = span(juncs)
        event = caller.event(eventtype, start, end, tot)
        events.append(event)
        event.introns = {junc_between(site, alt) for alt in alternatives}
        for alt in good:
            junc = junc_between(site, alt)
            event.add(coord_ranges([junc]), alternatives[alt], [junc], juncs - {junc}, blocks=[[junc]])
        caller.keep_other(event, other)
        # all the alternative sites lie on the same side of the shared one
        event.totoverlap = caller.covering((site, min(good)) if site < min(good) else (max(good), site))
    return events, groups


def shared_site_events(caller, eventtype, groups):
    """Splice sites spliced alike to the same alternative sites (alt5ss or alt3ss, by
    eventtype), from the groups of alternative_site_events"""
    events = []
    for alts, site_totals in groups.items():
        tot = total(site_totals.values())
        if len(site_totals) <= 1 or not caller.supported(tot, caller.thresholds.event_support):
            continue
        good, other = caller.split(site_totals, tot, need_support=False)
        if len(good) <= 1:
            continue
        juncs = {junc_between(site, alt) for site in good for alt in alts}
        # unique: each site is in one group, and the event ends at its outermost site
        start, end = span(juncs)
        event = caller.event(eventtype, start, end, tot)
        events.append(event)
        event.introns = {junc_between(site, alt) for site in site_totals for alt in alts}
        for site in good:
            inc = {junc_between(site, alt) for alt in alts}
            # 10 bases of the exon, on its side of the site: upstream of a 5' site, downstream of a 3' one
            exon_left = lies_left(eventtype == 'alt5ss', caller.strand)
            block = (site - 10, site) if exon_left else (site, site + 10)
            event.add(coord_ranges(inc), site_totals[site], inc, juncs - inc, blocks=[[block]])
        caller.keep_other(event, other)
        event.totoverlap = caller.covering((max(good), min(alts)) if median(good) < median(alts) else (max(alts), min(good)))
    return events


def intron_retention_events(caller, info):
    "Retained introns (ir): each junction's reads, and those of the blocks containing it"
    events = []
    for junc, spliced in info.juncs.items():
        retained = total([counts for block, counts in info.blocks.items() if block[0] < junc[0] and junc[1] < block[1]])
        tot = spliced + retained
        if caller.supported(spliced) and caller.supported(retained) and caller.varies(retained, tot):
            event = caller.event('ir', *junc, tot)
            events.append(event)
            event.introns = {junc}
            event.add('spliced', spliced, [junc], (), blocks=[junction_blocks(junc)])
            event.add('retained', retained, (), [junc], blocks=[[junc]])
            event.totoverlap = caller.covering(junc)
    return events


def extract_terminal_info(samples, chains_by_sample, key, strand, terminal):
    """the junctions, first (afe) or last (ale) exons and read starts (or ends) of the
    chains of key, a 5' (or 3') gene: reads per sample of each junction, of each first
    (or last) exon's junction by its splice site, and (position, sample) per read.  Its
    fragments, under no junctions, have only their starts (or ends)"""
    juncs, terminal_exons, ends = {}, {}, []
    # a first exon's splice site, and the read's start, are upstream, a last exon's downstream
    take_start = lies_left(terminal == 'afe', strand)
    for sample, chains in zip(samples, chains_by_sample):
        for chain_juncs, reads in chains.get(key, {}).items():
            for junc in chain_juncs:
                add_reads(juncs, junc, sample, len(reads))
            if chain_juncs:
                junc, site = (chain_juncs[0], chain_juncs[0][0]) if take_start else (chain_juncs[-1], chain_juncs[-1][1])
                add_reads(terminal_exons.setdefault(site, {}), junc, sample, len(reads))
            ends.extend((read.start if take_start else read.end, sample) for read in reads)
    return juncs, terminal_exons, ends


def group_ends(ends, window):
    """read ends, (position, sample), clustered where they lie within window of the
    next, each cluster's reads per sample counted at its median"""
    grouped, cluster = {}, []
    for end in sorted(ends):
        if cluster and end[0] - cluster[-1][0] > window:
            _add_end_cluster(grouped, cluster)
            cluster = []
        cluster.append(end)
    if cluster:
        _add_end_cluster(grouped, cluster)
    return grouped


def _add_end_cluster(grouped, cluster):
    grouped.setdefault(int(median([pos for pos, _ in cluster])), Counter()).update([sample for _, sample in cluster])


def near_annotated_site(site, annot_sites):
    "is site within TERMINAL_JUNCTION_DIST_FROM_ANNOT of one of annot_sites, a sorted list"
    i = bisect_left(annot_sites, site - TERMINAL_JUNCTION_DIST_FROM_ANNOT)
    return i < len(annot_sites) and annot_sites[i] <= site + TERMINAL_JUNCTION_DIST_FROM_ANNOT


def is_terminal_site(junc_counts, juncs, min_read_support, annotated):
    """Is a splice site of first (or last) junctions really a first (or last) exon's:
    pooled over the samples, at least min_read_support reads have its junctions as
    their first (or last), and few of its junctions' reads have an exon before (or
    after) it, more being allowed of an annotated one"""
    terminal = sum(sum(counts.values()) for counts in junc_counts.values())
    all_reads = sum(sum(juncs[junc].values()) for junc in junc_counts)
    max_fraction = MAX_NONTERMINAL_FRACTION_ANNOTATED if annotated else MAX_NONTERMINAL_FRACTION
    return terminal >= min_read_support and all_reads - terminal <= max_fraction * all_reads


def terminal_exon_events(caller, eventtype, terminal_exons, juncs, annot_sites):
    """Alternative first (afe) or last (ale) exons of a gene, from terminal_exons, its
    reads' first or last junctions by their splice site.  Only the sites really of a
    first or last exon (is_terminal_site) are its alternatives, and only their reads
    are the event's: a truncated read's first junction, after an internal exon, is
    neither.  Also the event's junctions, which alt5 and alt3 events don't repeat"""
    site_counts = {site: total(junc_counts.values()) for site, junc_counts in terminal_exons.items()
                   if is_terminal_site(junc_counts, juncs, caller.thresholds.feature_support,
                                       near_annotated_site(site, annot_sites))}
    tot = total(site_counts.values())
    if len(site_counts) <= 1 or not caller.supported(tot, caller.thresholds.event_support):
        return [], set()
    good, other = caller.split(site_counts, tot)
    # the sites of a large enough part of the reads, far enough apart
    major = [site for site in good
             if any(site_counts[site][s] / tot[s] >= MIN_TERMINAL_SS_FRACTION for s in caller.samples if tot[s] > 0)]
    if len(major) <= 1 or max(major) - min(major) < MIN_TERMINAL_JUNCTION_SEPARATION:
        return [], set()
    event_juncs = {junc for site in good for junc, counts in terminal_exons[site].items()
                   if caller.supported(counts) and caller.some_psi_between(counts, tot)}
    # each first or last exon by its introns: those of the event's, or failing that, all its own
    site_juncs = {site: ({j for j in event_juncs if site in j} or set(terminal_exons[site])) for site in good}
    event = caller.event(eventtype, *span([j for js in site_juncs.values() for j in js]), tot)
    for site in good:
        # 100 bases of the exon, on its side of the splice site: upstream of a first exon's, downstream of a last's
        exon_left = lies_left(eventtype == 'afe', caller.strand)
        block = (site - 100, site) if exon_left else (site, site + 100)
        inc = {j for j in event_juncs if site in j}
        event.add(coord_ranges(site_juncs[site]), site_counts[site], inc, event_juncs - inc, blocks=[[block]])
    caller.keep_other(event, other)
    event.totoverlap = caller.covering((max(j[0] for j in event_juncs), min(j[1] for j in event_juncs)))
    return [event], {frozenset(event_juncs)}


def end_events(caller, eventtype, ends):
    "Alternative starts (tss) or ends (tts) of a gene, from its reads' ends grouped by group_ends"
    tot = total(ends.values())
    if not caller.supported(tot, caller.thresholds.event_support):
        return []
    good, other = caller.split(ends, tot)
    if len(good) <= 1:
        return []
    # the 1-based position of each end's base: the reads' starts are 0-based, on +
    # their tss and on - their tts
    from_starts = lies_left(eventtype == 'tss', caller.strand)
    positions = {end: end + 1 if from_starts else end for end in good}
    event = caller.event(eventtype, min(positions.values()) - 1, max(positions.values()), tot)
    for end, pos in positions.items():
        event.add(str(pos), ends[end], {pos}, set(positions.values()) - {pos}, blocks=[[(pos - 1, pos)]])
    caller.keep_other(event, other)
    return [event]


class ReadCoverage:
    """The reads of a strand of a region, per sample, by their spans rounded to 10
    bases, for the reads covering an event's locus, whatever their gene"""
    def __init__(self):
        self.spans = {}
        self.index = None

    def add(self, start, end, sample):
        add_reads(self.spans, (round(start, -1), round(end, -1)), sample, 1)

    def covering(self, ref):
        "reads per sample whose spans contain ref, (start, end)"
        if self.index is None:
            self.index = IntervalIndex()
            for read_span, counts in self.spans.items():
                self.index.add(read_span[0], read_span[1], (read_span, counts))
            self.index.build()
        # a span containing ref overlaps the range between its ends, which can be in either
        # order, as a terminal exon event's (its introns' last start, first end)
        low, high = min(ref), max(ref)
        return total([counts for read_span, counts in self.index.overlap(low, max(high, low + 1))
                      if read_span[0] < ref[0] and ref[1] < read_span[1]])


def junc_string(chrom, features):
    "junctions or exons, 1-based, or end positions, each with its chromosome"
    return ','.join(f'{chrom}:{f[0] + 1}-{f[1]}' if isinstance(f, tuple) else f'{chrom}:{f}' for f in sorted(features))


def is_expressed(event, samples, thresholds):
    """Is an event expressed enough to write: in some sample, its reads at least
    event_frac_of_tot of its locus reads, and some row, but an exon skip's exclusion,
    varying as a part of the locus reads"""
    locus = event.totoverlap
    if not any(locus[s] > 0 and thresholds.event_frac_of_tot <= event.totjunc[s] / locus[s] for s in samples):
        return False
    return any(psi_varies(sample_psis(row.counts, locus, samples, 1), thresholds.junc_frac_of_event)
               for alternative, row in event.rows.items() if event.eventtype not in ('es', 'ces') or alternative != 'exc')


class EventTables:
    """A region's tables of events, open for writing: each row's counts and two PSIs,
    over the event's reads (PSIjunc) and over the reads covering its locus (PSItot),
    with the row's introns that aren't annotated, its BED lines scored by its median
    PSIjunc, and the event's two totals, which spliceevents_outliers divides by"""
    def __init__(self, prefix, samples, event_support):
        self.samples = samples
        self.event_support = event_support
        self.files = {suffix: open(prefix + suffix, 'w') for suffix in REGION_OUTPUT_SUFFIXES}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        for fh in self.files.values():
            fh.close()

    def _line(self, suffix, columns):
        self.files[suffix].write('\t'.join(columns) + '\n')

    def _write_row(self, columns, counts, event):
        "a row's counts and PSIs; returns its PSIjuncs, of the samples with event_support reads"
        self._line('.spliceevents.counts.tsv', columns + [str(counts[s]) for s in self.samples])
        for suffix, totals in (('.spliceevents.PSIjunc.tsv', event.totjunc), ('.spliceevents.PSItot.tsv', event.totoverlap)):
            self._line(suffix, columns + [str(round(counts[s] / totals[s], 4)) if totals[s] >= self.event_support else ''
                                          for s in self.samples])
        return [counts[s] / event.totjunc[s] for s in self.samples if event.totjunc[s] >= self.event_support]

    def _write_beds(self, event, alternative, blocks_list, psis):
        "a row's BED records, scored by its median PSIjunc, none if no sample has a PSI"
        if not psis:
            return
        for blocks in blocks_list:
            start, end = blocks[0][0], blocks[-1][1]
            Bed(event.chrom, start, end, row_name(event.name, alternative), score=round(median(psis) * 100),
                strand=event.strand, thickStart=start, thickEnd=end, itemRgb=EVENT_COLORS[event.eventtype],
                blocks=[BedBlock(s, e) for s, e in blocks]).write(self.files['.spliceevents.bed'])

    def write(self, event, annotated):
        genes = [event.gene, annotated.gene_names(event.gene.split(',')),
                 ','.join(event.additional_genes), annotated.gene_names(event.additional_genes)]
        for alternative, row in event.rows.items():
            novel = {j for j in row_introns(event, alternative) if not annotated.is_annotated(j, event.strand)}
            columns = [event.name, alternative, event.eventtype, *genes,
                       *(junc_string(event.chrom, features)
                         for features in (row.inc_juncs, row.exc_juncs, row.outer_juncs, row.exons, novel))]
            self._write_beds(event, alternative, row.blocks, self._write_row(columns, row.counts, event))
        if any(event.other[s] > 0 for s in self.samples):
            self._write_row([event.name, 'other', event.eventtype, *genes, '', '', '', '', ''], event.other, event)
        for label, counts in (('event_reads', event.totjunc), ('locus_reads', event.totoverlap)):
            self._line('.spliceevents.totals.tsv', [event.name, event.eventtype, *genes, label] + [str(counts[s]) for s in self.samples])


def chain_isoforms(chrom, chains):
    """ReadChains as Isoforms, by (strand, junctions), with all their reads, spanning
    the median of their starts and of their ends"""
    isoforms = {}
    for chain in chains:
        key = (chain.strand, chain.juncs)
        if key not in isoforms:
            isoforms[key] = Isoform(chrom, chain.strand, tuple(Junc(*j) for j in chain.juncs))
        isoforms[key].reads.extend(chain.ends)
    for isoform in isoforms.values():
        isoform.start, isoform.end = int(median(isoform.starts)), int(median(isoform.ends))
    return isoforms


def _splicing_clusters(chain_isoforms):
    """The clusters of junction chains, (strand, junctions) keys of chain_isoforms,
    linked by shared splice sites (group_by_splicing), each named for its strand and
    order: each chain's cluster, and each cluster's strand.  Clusters whose junctions
    only overlap, as a gene nested in another's intron, are not merged: they share
    no event"""
    cluster_of, key_strand = {}, {}
    for strand in ('+', '-'):
        keys = sorted(k for k in chain_isoforms if k[0] == strand)
        for group in sorted(group_by_splicing(keys, chain_isoforms, merge_overlapping=False), key=lambda g: min(k[1][0][0] for k in g)):
            cluster = f'{strand}{len(key_strand)}'
            key_strand[cluster] = strand
            cluster_of.update((k, cluster) for k in group)
    return cluster_of, key_strand


def group_chains(chrom, samples, sample_chains, novel_gene):
    """Each sample's reads grouped three ways, each a list per sample of
    key -> junctions -> [ReadEnds]: by cluster of junction chains linked by shared
    splice sites, on one strand, for the events of splicing, which don't depend on
    genes; by each chain's 5'-most gene, for the first exons and starts; and by its
    3'-most gene, for the last exons and ends.  A gene's 5' fragments are with its
    starts and its 3' fragments with its ends, under no junctions, and nowhere else.
    Also each key's strand, each chain's genes, by (strand, junctions), a novel gene's
    if of none annotated, and each strand's ReadCoverage, of its spliced reads"""
    chain_genes = {}
    coverage = {'+': ReadCoverage(), '-': ReadCoverage()}
    for sample, chains in zip(samples, sample_chains):
        for chain in chains:
            if not chain.fragment:
                key = (chain.strand, chain.juncs)
                chain_genes[key] = chain.genes or (novel_gene[key],)
                for read in chain.ends:
                    coverage[chain.strand].add(read.start, read.end, sample)
    isoforms = chain_isoforms(chrom, (chain for chains in sample_chains for chain in chains if not chain.fragment))
    cluster_of, key_strand = _splicing_clusters(isoforms)
    return _chain_views(sample_chains, chain_genes, cluster_of, key_strand), key_strand, chain_genes, coverage


def _chain_views(sample_chains, chain_genes, cluster_of, key_strand):
    "group_chains' three groupings of each sample's chains, adding the genes' strands to key_strand"
    views = {'cluster': [], 'five': [], 'three': []}
    for chains in sample_chains:
        by = {'cluster': {}, 'five': {}, 'three': {}}
        for chain in chains:
            if chain.fragment:
                gene = chain.genes[0]
                by['five' if chain.fragment == 'first' else 'three'].setdefault(gene, {})[()] = chain.ends
                key_strand[gene] = chain.strand
                continue
            key = (chain.strand, chain.juncs)
            genes = chain_genes[key]
            for view, group in (('cluster', cluster_of[key]), ('five', genes[0]), ('three', genes[-1])):
                by[view].setdefault(group, {})[chain.juncs] = chain.ends
                key_strand[group] = chain.strand
        for view in views:
            views[view].append(by[view])
    return views


class AnnotatedSplicing:
    """The annotated introns and splice sites of each strand, and the genes they are
    annotated in, to give events their genes and tell their novel junctions"""
    def __init__(self, annots=None):
        self.intron_genes = {'+': {}, '-': {}}
        self.site_genes = {'+': {}, '-': {}}
        self.num_transcripts = {}
        self.names = {}
        if annots is not None:
            self.num_transcripts = annots.gene_to_num_transcripts
            self.names = annots.gene_to_name
            for gene, juncs in annots.gene_to_annot_juncs.items():
                strand = annots.gene_to_strand.get(gene)
                if strand not in self.intron_genes:
                    continue
                for j in juncs:
                    self.intron_genes[strand].setdefault((j.start, j.end), set()).add(gene)
                    for site in (j.start, j.end):
                        self.site_genes[strand].setdefault(site, set()).add(gene)

    def is_annotated(self, junc, strand):
        return junc in self.intron_genes[strand]

    def gene_names(self, genes):
        "the annotation's names of genes, a list, comma-separated, empty if it names none of them"
        names = [self.names.get(gene, '') for gene in genes]
        return ','.join(names) if any(names) else ''


def _top(scores):
    "the keys of the highest score, several on a tie"
    if not scores:
        return []
    best = max(scores.values())
    return sorted(key for key, score in scores.items() if score == best)


def row_introns(event, alternative):
    """the introns spliced by the reads of an event row: those it includes, or the
    skipping junctions of an exon skip's exclusion"""
    row = event.rows[alternative]
    if event.eventtype in ('es', 'ces') and alternative == 'exc':
        return set(row.outer_juncs)
    return {j for j in row.inc_juncs if isinstance(j, tuple)}


def event_gene_scores(introns, strand, annotated, junc_reads):
    """the genes sharing an intron or splice site with an event, on its strand, and their
    scores: 2 for an intron the gene has annotated, or failing that 1 for each of its
    splice sites the gene has, each intron weighted by its share of the reads of all
    the introns, junc_reads, so a rare alternative counts little"""
    total_reads = sum(junc_reads.get(junc, 0) for junc in introns)
    scores = Counter()
    for junc in introns:
        weight = junc_reads.get(junc, 0) / total_reads if total_reads > 0 else 1
        intron_genes = annotated.intron_genes[strand].get(junc, set())
        for gene in intron_genes:
            scores[gene] += 2 * weight
        for site in junc:
            for gene in annotated.site_genes[strand].get(site, set()) - intron_genes:
                scores[gene] += weight
    # genes of the same introns tie, whatever order their weights were added in
    return {gene: round(score, 9) for gene, score in scores.items()}


def annotated_event_genes(scores, annotated):
    """the genes of the highest event_gene_scores.  A tie goes to the gene of most
    annotated transcripts, so a gene and a readthrough of it give the gene; genes still
    tied are all given"""
    best = _top(scores)
    return _top({gene: annotated.num_transcripts.get(gene, 0) for gene in best})


def read_event_genes(introns, cluster_chains, chain_genes, strand):
    """the genes of most reads of an event's introns, by the genes the reads were given
    in correction, each gene of a readthrough counted: for an event of no annotated
    intron or splice site"""
    reads = Counter()
    for chains in cluster_chains:
        for juncs, ends in chains.items():
            if introns.intersection(juncs):
                for gene in chain_genes[(strand, juncs)]:
                    reads[gene] += len(ends)
    return _top(reads)


def _keys(view):
    return sorted(set().union(*[sample_groups.keys() for sample_groups in view]))


def gene_terminal_events(caller, view, gene, terminal, end_type, annot_sites):
    """A gene's first (afe) or last (ale) exon event, and start (tss) or end (tts)
    event, from its chains in view, and the first or last exon event's junctions"""
    juncs, terminal_exons, ends = extract_terminal_info(caller.samples, view, gene, caller.strand, terminal)
    events, event_juncs = terminal_exon_events(caller, terminal, terminal_exons, juncs, annot_sites)
    events += end_events(caller, end_type, group_ends(ends, END_GROUP_WINDOW))
    for event in events:
        event.gene = gene
    return events, event_juncs


def cluster_events(caller, cluster_chains, terminal_juncs, annotated, chain_genes):
    """The splicing events of a cluster of junction chains, each sample's, each event
    given the gene annotated with most of its introns, those of all its alternatives,
    and their splice sites, weighted by their reads (annotated_event_genes), or with
    none, that of most of its reads (read_event_genes); several only on a tie"""
    info = extract_splicing_info(caller.samples, cluster_chains, caller.strand)
    es, es_juncs, seen = exon_skipping_events(caller, info)
    ces = complex_skipping_events(caller, info, es_juncs, seen)
    alt3, alt3_groups = alternative_site_events(caller, 'alt3', info.ss5to3, es_juncs, terminal_juncs)
    alt5, alt5_groups = alternative_site_events(caller, 'alt5', info.ss3to5, es_juncs, terminal_juncs)
    events = (es + ces + alt3 + shared_site_events(caller, 'alt5ss', alt3_groups)
              + alt5 + shared_site_events(caller, 'alt3ss', alt5_groups) + intron_retention_events(caller, info))
    junc_reads = {junc: sum(counts.values()) for junc, counts in info.juncs.items()}
    for event in events:
        scores = event_gene_scores(event.introns, caller.strand, annotated, junc_reads)
        genes = annotated_event_genes(scores, annotated) or read_event_genes(event.introns, cluster_chains, chain_genes, caller.strand)
        event.gene = ','.join(genes)
        # the other genes sharing a splice site, best first
        event.additional_genes = sorted(set(scores) - set(genes), key=lambda gene: (-scores[gene], gene))
    return events


def process_region_events(temp_prefix, chrom, samples, sample_chains, novel_gene, thresholds, annot_afe_ss, annot_ale_ss,
                          report_regions, gene_juncs, annotated=None, owned=None):
    """A region's events, from each sample's ReadChains, written to its tables.  The
    events of splicing are found per cluster of junction chains linked by shared
    splice sites, without genes (cluster_events); the alternative first exons and
    starts per each chain's 5'-most gene, and the last exons and ends per its 3'-most
    gene, so a readthrough's reads count with the gene each end lies in.  An event's
    locus reads are all the strand's reads covering it.  Only the events inside
    report_regions, and given owned, starting in it, are written (event_in_regions)"""
    annotated = annotated or AnnotatedSplicing()
    views, key_strand, chain_genes, coverage = group_chains(chrom, samples, sample_chains, novel_gene)
    callers = {strand: EventCaller(chrom, strand, samples, thresholds, coverage[strand]) for strand in ('+', '-')}
    events = []  # (event, the span of its gene's junctions, for the regions)
    terminal_juncs = {'+': set(), '-': set()}
    for view, terminal, annot_ss, end_type in (('five', 'afe', annot_afe_ss, 'tss'), ('three', 'ale', annot_ale_ss, 'tts')):
        for gene in _keys(views[view]):
            caller = callers[key_strand[gene]]
            gene_events, event_juncs = gene_terminal_events(caller, views[view], gene, terminal, end_type,
                                                            annot_ss.get(caller.strand, []))
            terminal_juncs[caller.strand] |= event_juncs
            gene_span = gene_junction_span(gene_juncs.get(gene, ()))
            events.extend((event, gene_span) for event in gene_events)
    for cluster in sorted(_keys(views['cluster']), key=lambda c: int(c[1:])):
        caller = callers[key_strand[cluster]]
        cluster_chains = [chains.get(cluster, {}) for chains in views['cluster']]
        events.extend((event, (0, 0)) for event in cluster_events(caller, cluster_chains, terminal_juncs[caller.strand],
                                                                  annotated, chain_genes))
    reported = {}
    for event, gene_span in events:
        if event_in_regions(event, report_regions, gene_span, owned):
            if event.name in reported:
                raise AssertionError(f'two events named {event.name}, of {reported[event.name].gene} and {event.gene}')
            reported[event.name] = event
    with EventTables(temp_prefix, samples, thresholds.event_support) as tables:
        for event in reported.values():
            if event.totoverlap is None:
                event.totoverlap = event.totjunc
            if is_expressed(event, samples, thresholds):
                tables.write(event, annotated)


def correct_sample_reads(args, region, work_dir, bamfile_name, genome, annots, junction_corrector, read_starts=None):
    """The sample's reads in region, corrected as flair transcriptome corrects them,
    grouped by junction chain: aligned to the annotated transcripts, then their
    junctions corrected, the unspliced reads under no junctions.  Given read_starts,
    a SeqRange, only the reads aligned to start in it.  The alignment's
    intermediates are written in work_dir, removed after"""
    os.makedirs(work_dir, exist_ok=True)
    temp_prefix = osp.join(work_dir, 'reads')
    # the reads as samtools fasta writes them: primary alignments, in the read's own orientation
    with pysam.AlignmentFile(bamfile_name, 'rb') as bam, open(temp_prefix + '.reads.fasta', 'w') as out:
        for read in bam.fetch(region.name, region.start, region.end):
            if not (read.is_secondary or read.is_supplementary) \
                    and (read_starts is None or read_starts.start <= read.reference_start < read_starts.end):
                out.write(f'>{read.query_name}\n{read.get_forward_sequence()}\n')
    sj_to_ends = {}
    with pysam.AlignmentFile(bamfile_name, 'rb') as bam_file:
        generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, region.name, region.start, region.end)
        read_to_transcript = identify_good_match_to_annot(
            temp_prefix, region, annots, genome, bam_file, no_align_to_annot=args.no_align_to_annot,
            directRNA=False, check_splice=True, trust_strand=args.trust_strand, fusion_breakpoints=None)
        filter_correct_group_reads(bam_file=bam_file, region=region,
                                   read_to_annot_transcript=read_to_transcript,
                                   annots=annots, junction_corrector=junction_corrector,
                                   genome=genome,
                                   quality=args.quality, keep_sup=args.keep_supplementary,
                                   sj_to_ends=sj_to_ends,
                                   allow_outside_range=True,
                                   keep_single_exon=True, trust_strand=args.trust_strand,
                                   check_motifs=not args.trust_junctions, read_starts=read_starts)
    shutil.rmtree(work_dir)
    return sj_to_ends


# a read's corrected ends, all event calling needs of it beyond its junction chain
ReadEnds = namedtuple('ReadEnds', ('start', 'end'))


def read_chains_path(partition, sample):
    "the file of a sample's corrected reads in a partition, handed from its correction to the region's event calling"
    return f'{partition.file_prefix}_{sample}.read_chains.tsv'


# a sample's corrected junction chain: its annotated genes, 5' to 3', or () for none,
# its junctions, strand, and its reads' ends.  Or its unspliced reads that are
# fragments of one gene, with no junctions, and fragment 'last' for 3' fragments,
# 'first' for 5' (AnnotatedFragmentExons)
ReadChain = namedtuple('ReadChain', ('genes', 'juncs', 'strand', 'ends', 'fragment'), defaults=('',))


def write_read_chains(path, chrom, chains):
    """Write a sample's corrected reads, ReadChains, a line per junction chain: its
    annotated genes, chromosome, junctions, strand, the ends of each of its reads,
    0-based start-end pairs, and, for fragments, without junctions, which kind.  It's
    written whole or not at all, as a rerun with --keep_intermediate skips the samples
    whose file exists"""
    with open(path + '.tmp', 'w') as out:
        for chain in chains:
            out.write('\t'.join([','.join(chain.genes), chrom, ','.join(f'{s}-{e}' for s, e in chain.juncs), chain.strand,
                                 ','.join(f'{r.start}-{r.end}' for r in chain.ends), chain.fragment]) + '\n')
    os.replace(path + '.tmp', path)


def _ranges(field):
    return tuple(tuple(int(x) for x in pair.split('-')) for pair in field.split(',')) if field else ()


def read_read_chains(path):
    """a sample's corrected reads from write_read_chains, ReadChains; of an earlier
    version's files, without fragments, too"""
    chains = []
    with open(path) as fh:
        for line in fh:
            genes, _, juncs, strand, ends, *fragment = line.rstrip('\n').split('\t')
            chains.append(ReadChain(tuple(genes.split(',')) if genes else (), _ranges(juncs), strand,
                                    [ReadEnds(*pair) for pair in _ranges(ends)], *fragment))
    return chains


def annotated_genes(isoform, annots):
    """A corrected junction chain's annotated genes, 5' to 3', and its strand, as
    flair transcriptome finds them (gene_identification): () if it has none, and a
    strand of None for a chain of unknown strand that no gene gives one"""
    if isoform.strand == UNKNOWN_STRAND:
        genes, _, strand = get_gene_name_unknown_strand(isoform, annots)
        if genes is None:
            return (), None
    else:
        genes, _ = get_gene_name_firstpass(isoform, annots)
        if genes is None:
            return (), isoform.strand
        strand = annots.gene_to_strand[genes[0]]
    # a readthrough's genes in the order it reads them
    return tuple(sorted(genes, key=lambda g: annots.gene_to_exons[g][0][0], reverse=(strand == '-'))), strand


def fragment_chains(reads, fragment_exons, trust_strand):
    """the unspliced reads that are fragments of one gene (AnnotatedFragmentExons), as a ReadChain
    per gene and kind"""
    by_fragment = {}
    for read in reads:
        found = fragment_exons.fragment(read.start, read.end, read.strand if trust_strand else None)
        if found is not None:
            by_fragment.setdefault(found, []).append(ReadEnds(read.start, read.end))
    return [ReadChain((gene,), (), strand, ends, kind) for (kind, gene, strand), ends in sorted(by_fragment.items())]


def get_juncs_single_sample(args, region, path, bamfile_name, genome, junction_corrector, annots, terminal_exons,
                            read_starts=None):
    """correct a sample's reads in region, give their junction chains annotated genes,
    and write them to path, with the unspliced reads that are fragments of one gene's
    first or last exons, terminal_exons; the chains of no annotated gene get novel
    genes in event calling, over all samples"""
    sj_to_ends = correct_sample_reads(args, region, path + '.correction', bamfile_name, genome, annots, junction_corrector,
                                      read_starts)
    fragment_exons = AnnotatedFragmentExons(terminal_exons, annots)
    chains = []
    for (_, juncs), isoform in sj_to_ends.items():
        if len(juncs) > 0:
            genes, strand = annotated_genes(isoform, annots)
            if strand is None:
                logging.debug(f'reads dropped: strand unknown and not given by a gene: {juncs}')
                continue
            chains.append(ReadChain(genes, tuple((j.start, j.end) for j in juncs), strand, isoform.reads))
        else:
            chains.extend(fragment_chains(isoform.reads, fragment_exons, args.trust_strand))
    write_read_chains(path, region.name, chains)


def novel_genes_of_chains(chrom, sample_chains):
    """(strand, junctions) -> novel gene of the ReadChains of no annotated gene, of all
    samples together, so each sample's reads get the same novel genes, grouped and
    named as flair transcriptome does (novel_genes)"""
    novel = chain_isoforms(chrom, (chain for chains in sample_chains for chain in chains if not chain.genes))
    novel_gene = {}
    for strand in ('+', '-'):
        for gene, keys in novel_genes(sorted(k for k in novel if k[0] == strand), novel, chrom, strand):
            novel_gene.update((k, gene) for k in keys)
    return novel_gene


def gene_junctions(sample_chains, novel_gene, annots):
    """Each gene's junctions, for --region_bed: an annotated gene's annotated ones,
    a novel gene's those of its reads"""
    gene_juncs = {}
    for chains in sample_chains:
        for chain in chains:
            if chain.genes:
                for gene in chain.genes:
                    if gene not in gene_juncs:
                        gene_juncs[gene] = {(j.start, j.end) for j in annots.gene_to_annot_juncs.get(gene, ())}
            else:
                gene_juncs.setdefault(novel_gene[(chain.strand, chain.juncs)], set()).update(chain.juncs)
    return gene_juncs


# an exon's first or last exon by kind: 'first' or 'last', strand, start, end, gene
AnnotatedTerminalExon = namedtuple('AnnotatedTerminalExon', ('kind', 'strand', 'start', 'end', 'gene'))


def annotated_terminal_exons(transcripts):
    """chrom -> the first and last exons, AnnotatedTerminalExons, of the spliced transcripts, as
    of GENCODE basic, whose first and last exons are real ones"""
    exons_by_chrom = {}
    for trans in transcripts:
        exons = sorted((exon.start, exon.end) for exon in trans.exons)
        if len(exons) > 1:
            first, last = (exons[0], exons[-1]) if trans.strand == '+' else (exons[-1], exons[0])
            terminal = exons_by_chrom.setdefault(trans.chrom, set())
            terminal.add(AnnotatedTerminalExon('first', trans.strand, *first, trans.gene_id))
            terminal.add(AnnotatedTerminalExon('last', trans.strand, *last, trans.gene_id))
    return {chrom: sorted(exons) for chrom, exons in exons_by_chrom.items()}


def annotated_terminal_splice_sites(terminal_exons):
    """chrom -> strand -> the splice sites, sorted, of terminal_exons' first exons, and
    of their last exons, which mark alternative first and last exons.  By strand, not
    gene, so they don't depend on the reads' genes"""
    afe_ss, ale_ss = {}, {}
    for chrom, exons in terminal_exons.items():
        for sites in (afe_ss, ale_ss):
            sites[chrom] = {'+': set(), '-': set()}
        for exon in exons:
            # the splice site of a first exon is its 3' end, of a last exon its 5' end
            exon_left = lies_left(exon.kind == 'first', exon.strand)
            site = exon.end if exon_left else exon.start
            (afe_ss if exon.kind == 'first' else ale_ss)[chrom][exon.strand].add(site)
    return tuple({chrom: {strand: sorted(strand_sites) for strand, strand_sites in by_strand.items()}
                  for chrom, by_strand in sites.items()} for sites in (afe_ss, ale_ss))


# how far past its splice site, into the intron, the end of an unspliced read in a
# first or last exon may be, as flair transcriptome's --end_window
FRAGMENT_SPLICE_SITE_WINDOW = 100


class AnnotatedFragmentExons:
    """The annotated first and last exons of a region, and the exons of all its genes,
    to tell the unspliced reads that are fragments of a gene's transcripts: in a last
    exon, a 3' fragment, whose end is a transcript's end, or in a first exon, a 5'
    fragment, whose start is a transcript's start"""
    def __init__(self, terminal_exons, annots):
        self.terminal = IntervalIndex()
        for exon in terminal_exons:
            self.terminal.add(exon.start, exon.end, exon)
        # each gene's exons, merged, so the bases a read shares with a gene count once
        self.exons = IntervalIndex()
        for gene, exons in annots.gene_to_exons.items():
            strand = annots.gene_to_strand[gene]
            merged = []
            for start, end in sorted(exons):
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            for start, end in merged:
                self.exons.add(start, end, (gene, strand, start, end))

    def fragment(self, start, end, strand=None):
        """('last' or 'first', gene, strand) of an unspliced read, or None: of the genes
        whose exons overlap the read, the one sharing the most bases with it, if that's
        one gene, and the read's end at its splice site, its 5' end in a last exon or
        its 3' end in a first, lies in the exon or within FRAGMENT_SPLICE_SITE_WINDOW
        past the site, while its other end may run past the transcript's end.  Given
        strand, the read's trusted strand, only the genes on it; otherwise the genes of
        both strands, as an unspliced read's strand isn't known"""
        shared = {}
        for gene, gene_strand, exon_start, exon_end in self.exons.overlap(start, end):
            if strand in (None, gene_strand):
                key = (gene, gene_strand)
                shared[key] = shared.get(key, 0) + min(end, exon_end) - max(start, exon_start)
        if len(shared) == 0:
            return None
        most = max(shared.values())
        best = [key for key, bases in shared.items() if bases == most]
        if len(best) != 1:
            return None
        gene, gene_strand = best[0]
        found = set()
        for exon in self.terminal.overlap(start, end, slack=FRAGMENT_SPLICE_SITE_WINDOW):
            if exon.gene != gene:
                continue
            # a first exon is upstream of its splice site, a last exon downstream
            if lies_left(exon.kind == 'first', exon.strand):  # the splice site is the exon's end
                at_site = exon.start < end <= exon.end + FRAGMENT_SPLICE_SITE_WINDOW
            else:
                at_site = exon.start - FRAGMENT_SPLICE_SITE_WINDOW <= start < exon.end
            if at_site:
                found.add(exon.kind)
        if len(found) != 1:
            return None
        return found.pop(), gene, gene_strand


def _correct_region_sample(*, partition, gtf_data, junction_corrector, args, item, terminal_exons, read_starts):
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
            get_juncs_single_sample(args, region, path, bamfile, genome, junction_corrector, annots, region_terminal,
                                    read_starts.get(_region_key(region)))


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
                          region_reports[_region_key(region)], gene_junctions(sample_chains, novel_gene, annots),
                          AnnotatedSplicing(annots), read_starts.get(_region_key(region)))
    if not args.keep_intermediate:
        for sample in samples:
            os.remove(read_chains_path(partition, sample))


# the tables each region writes, combined into the run's
REGION_OUTPUT_SUFFIXES = ('.spliceevents.bed', '.spliceevents.counts.tsv', '.spliceevents.PSIjunc.tsv',
                          '.spliceevents.PSItot.tsv', '.spliceevents.totals.tsv')


def _region_key(region):
    return (region.name, region.start, region.end)


def gene_clusters(region, transcripts):
    """The spans of region's annotated genes, both strands, merged where they overlap,
    and the number of exons in each, sorted"""
    spans = sorted((max(t.start, region.start), min(t.end, region.end), len(t.exons)) for t in transcripts
                   if t.chrom == region.name and t.start < region.end and region.start < t.end)
    clusters = []
    for start, end, exons in spans:
        if clusters and start < clusters[-1][1]:
            clusters[-1][1] = max(clusters[-1][1], end)
            clusters[-1][2] += exons
        else:
            clusters.append([start, end, exons])
    return clusters


def intergenic_cuts(region, transcripts, min_exons=0):
    """Where to cut region into partitions: the middle of the gaps between its
    gene clusters, each partition with at least min_exons annotated exons, but
    the last, so 0 cuts at every gap"""
    clusters = gene_clusters(region, transcripts)
    cuts, exons = [], 0
    for (_, end, n), (start, _, _) in zip(clusters, clusters[1:]):
        exons += n
        if exons >= min_exons:
            cuts.append((end + start) // 2)
            exons = 0
    return cuts


def _read_spans(bamfiles, region):
    "the (start, end) of each mapped, not secondary, alignment overlapping region, in all bam files"
    for bamfile in bamfiles:
        with pysam.AlignmentFile(bamfile, 'rb') as bam:
            for read in bam.fetch(region.name, region.start, region.end):
                if read.is_mapped and not read.is_secondary:
                    yield read.reference_start, read.reference_end


def correction_regions(report_regions, bamfiles):
    """The regions to correct reads in, for the regions events are reported in: each
    report region extended over the reads overlapping it in any sample, so a read
    running past its end is corrected with the annotation and junctions it reaches,
    as anywhere else.  Extended regions that overlap are merged.  Returns the
    regions, and a map of each one's key to the report regions in it."""
    extended = []
    for region in report_regions:
        start, end = region.start, region.end
        for read_start, read_end in _read_spans(bamfiles, region):
            start, end = min(start, read_start), max(end, read_end)
        extended.append((SeqRange(region.name, start, end), region))
    merged = combine_regions([ext for ext, _ in extended])
    reports = {_region_key(region): [] for region in merged}
    for ext, region in extended:
        container = next(m for m in merged if m.name == ext.name and m.start <= ext.start and ext.end <= m.end)
        reports[_region_key(container)].append(region)
    return merged, reports


def partition_regions(report_regions, bamfiles, cuts):
    """The partitions of report regions cut at cuts, a dict of region key to sorted
    positions: each piece owns the reads aligned to start in it, the first also
    those starting before the region, and is corrected over their whole length.
    Reads crossing a cut are only the piece's they start in, so the events past a
    cut lose them.  Returns, for each piece with reads, its region to correct reads
    in, the range its reads start in, its report region, and its number of reads"""
    pieces = []
    for region in report_regions:
        bounds = [region.start] + cuts.get(_region_key(region), []) + [region.end]
        low, high = bounds[:-1], bounds[1:]
        reads = [0] * len(low)
        for read_start, read_end in _read_spans(bamfiles, region):
            i = max(bisect_right(bounds, read_start) - 1, 0)
            low[i], high[i] = min(low[i], read_start), max(high[i], read_end)
            reads[i] += 1
        for i in range(len(low)):
            if reads[i] > 0:
                owned = SeqRange(region.name, bounds[i] if i > 0 else 0, bounds[i + 1])
                pieces.append((SeqRange(region.name, low[i], high[i]), owned, region, reads[i]))
    return pieces


# events whose reads are all of a gene's: its starts and ends, first and last exons
GENE_WIDE_EVENTS = {'tss', 'tts', 'afe', 'ale'}


def _event_positions(event):
    "the splice sites and ends of an event's junctions, exons and ends"
    positions = []
    for row in event.rows.values():
        for feature in (*row.inc_juncs, *row.exc_juncs, *row.outer_juncs, *row.exons):
            positions.extend(feature if isinstance(feature, tuple) else (feature,))
    return positions


def gene_junction_span(juncs):
    """The span of a gene's junctions, its annotated ones, or a novel gene's of its
    reads, from the first one's start to the last one's end"""
    if not juncs:
        return (0, 0)
    return (min(j[0] for j in juncs), max(j[1] for j in juncs))


def event_in_regions(event, regions, gene_span, owned=None):
    """Is an event inside one of the regions events are reported in: all its
    junctions and exons, so every read reaching it overlaps the region and its
    counts are complete.  An event of all its gene's reads (GENE_WIDE_EVENTS) only
    when all the gene's annotated junctions are, gene_span (gene_junction_span),
    as a user choosing a region around a gene would think of it: its reads are
    spliced, so nearly all have a junction in the region and were fetched.  The
    gene's first and last exons may run past the region.  Given owned, the range
    a partition's reads start in, only an event starting in it, so each event is
    reported by one partition"""
    if event.eventtype in GENE_WIDE_EVENTS:
        low, high = gene_span
    else:
        positions = _event_positions(event)
        if not positions:
            return False
        low, high = min(positions), max(positions)
    if owned is not None and not owned.start <= low < owned.end:
        return False
    return any(region.name == event.chrom and region.start <= low and high <= region.end for region in regions)


def combine_regions(regions):
    "regions sorted, with those overlapping or touching merged"
    merged = []
    for region in sorted(regions):
        if merged and merged[-1].name == region.name and region.start <= merged[-1].end:
            merged[-1] = SeqRange(region.name, merged[-1].start, max(merged[-1].end, region.end))
        else:
            merged.append(SeqRange(region.name, region.start, region.end))
    return merged


# the columns of the counts and PSI tables before the samples': a row of an event,
# its event, and its junctions, and of the totals table
EVENT_GENE_COLUMNS = ['gene', 'gene_name', 'additional_genes', 'additional_gene_names']
EVENT_ROW_COLUMNS = ['event', 'alternative', 'event_type', *EVENT_GENE_COLUMNS, 'junctions_included', 'junctions_excluded',
                     'outer_junctions', 'exons', 'novel_junctions']
EVENT_TOTALS_COLUMNS = ['event', 'event_type', *EVENT_GENE_COLUMNS, 'total']


def write_with_header(output, prefixes, suffix, header):
    "the partitions' tables of one suffix, as OUTPUT + suffix under header"
    with open(output + suffix, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for prefix in prefixes:
            with open(prefix + suffix) as fh:
                shutil.copyfileobj(fh, out)


def combine_region_tables(output, prefixes, samples):
    "the regions' tables as the run's, OUTPUT.spliceevents.*"
    combine_temp_files_by_suffix(output, prefixes, ['.spliceevents.bed'])
    for suffix, columns in (('.spliceevents.counts.tsv', EVENT_ROW_COLUMNS), ('.spliceevents.PSIjunc.tsv', EVENT_ROW_COLUMNS),
                            ('.spliceevents.PSItot.tsv', EVENT_ROW_COLUMNS), ('.spliceevents.totals.tsv', EVENT_TOTALS_COLUMNS)):
        write_with_header(output, prefixes, suffix, columns + samples)


def read_manifest(path):
    "a manifest's samples, (name, bam file), a line each but comments"
    allsamples = []
    with open(path) as fh:
        for line in fh:
            if line[0] != '#':
                sample, bamfile = line.rstrip().split('\t')
                allsamples.append((sample, bamfile))
    return allsamples


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
        cuts = {_region_key(region): intergenic_cuts(region, gtf_data.transcripts, PARTITION_EXONS) for region in report_regions}
        plan = CorrectionPlan([], {}, {}, [], gtf_data)
        for region, owned, report, nreads in partition_regions(report_regions, bamfiles, cuts):
            plan.regions.append(region)
            plan.reports[_region_key(region)] = [report]
            plan.read_starts[_region_key(region)] = owned
            plan.weights.append(nreads)
        return plan
    if args.region_bed:
        regions, reports = correction_regions(report_regions, bamfiles)
        # only the annotation of the regions reads are corrected in
        return CorrectionPlan(regions, reports, {}, None, load_gtf(args.annot, regions))
    return CorrectionPlan(report_regions, {_region_key(region): [region] for region in report_regions}, {}, None,
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
               read_starts=plan.read_starts)
    runner.run(_run_region, args=args, allsamples=allsamples, region_reports=plan.reports,
               annot_afe_ss=annot_afe_ss, annot_ale_ss=annot_ale_ss, read_starts=plan.read_starts)

    combine_region_tables(args.output, [p.file_prefix for p in runner], [sample for sample, _ in allsamples])
    if not args.keep_intermediate:
        shutil.rmtree(temp_dir)
