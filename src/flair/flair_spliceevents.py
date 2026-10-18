import os
from bisect import bisect_left
import os.path as osp
from flair import FlairInputDataError
import pipettor
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


def zero_counts(allsamples):
    "no reads in any sample"
    return {s: 0 for s in allsamples}


def add_reads(counts, key, allsamples, sample, numreads):
    "add a sample's reads to counts[key], its reads per sample"
    if key not in counts:
        counts[key] = zero_counts(allsamples)
    counts[key][sample] += numreads


def add_counts_at(counts, key, other, allsamples):
    "add other, reads per sample, to counts[key]"
    counts[key] = add_counts_to_dict(counts.get(key, zero_counts(allsamples)), other)


def sum_counts(all_counts, allsamples):
    "the sum of reads per sample"
    total = zero_counts(allsamples)
    for counts in all_counts:
        total = add_counts_to_dict(total, counts)
    return total


def supported(counts, allsamples, min_reads):
    "do some sample's reads reach min_reads"
    return any(counts[s] >= min_reads for s in allsamples)


def sample_psis(counts, totals, allsamples, min_total):
    "the PSIs, counts over totals, of the samples with at least min_total reads in totals"
    return [counts[s] / totals[s] for s in allsamples if totals[s] >= min_total]


def junc_between(ss1, ss2):
    "the junction between two splice sites"
    return (min(ss1, ss2), max(ss1, ss2))


def some_psi_between(counts, totals, allsamples, junc_frac_of_event):
    "is counts, of totals, between junc_frac_of_event and 1 minus it in some sample"
    return any(junc_frac_of_event <= counts[s] / totals[s] <= 1 - junc_frac_of_event for s in allsamples if totals[s] > 0)


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


def _event_bed(chrom, blocks, name, strand, color):
    "an event's BED record of blocks, (start, end) pairs, scored when the event is written"
    start, end = blocks[0][0], blocks[-1][1]
    return Bed(chrom, start, end, name, score=0, strand=strand, thickStart=start, thickEnd=end, itemRgb=color,
               blocks=[BedBlock(s, e) for s, e in blocks])


def junction_bed(chrom, junc, exons, name, strand, color):
    """an event's BED record of a junction, as 10 bases each side of it, and the exons
    included between its ends"""
    return _event_bed(chrom, [(junc[0] - 10, junc[0])] + list(exons) + [(junc[1], junc[1] + 10)], name, strand, color)


def site_bed(chrom, start, end, name, strand, color):
    "an event's BED record of one block"
    return _event_bed(chrom, [(start, end)], name, strand, color)


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

def extract_35ss_info(juncs, strand, allsamples, sample, numreads, ss5to3, ss3to5, alljuncs):
    for j in juncs:
        ss5, ss3 = (j[0], j[1]) if strand == '+' else (j[1], j[0])
        add_reads(ss5to3.setdefault(ss5, {}), ss3, allsamples, sample, numreads)
        add_reads(ss3to5.setdefault(ss3, {}), ss5, allsamples, sample, numreads)
        add_reads(alljuncs, j, allsamples, sample, numreads)
    return ss5to3, ss3to5, alljuncs


def extract_exon_usage_info(juncs, allsamples, sample, numreads, allblocks, exonjpairs):
    # just for middle exons
    for i in range(len(juncs) - 1):
        add_reads(allblocks, (juncs[i][1], juncs[i + 1][0]), allsamples, sample, numreads)
        add_reads(exonjpairs, (juncs[i], juncs[i + 1]), allsamples, sample, numreads)
    return allblocks, exonjpairs


def extract_end_coverage_info(juncs, readinfo, allsamples, sample, allblocks):
    for r in readinfo:
        # getting all block coverage to measure intron retention
        add_reads(allblocks, (r.start, juncs[0][0]), allsamples, sample, 1)
        add_reads(allblocks, (juncs[-1][1], r.end), allsamples, sample, 1)
    return allblocks

def group_ends(ends, allsamples, window):
    """read ends, (position, sample), clustered where they lie within window of the
    next, each cluster's reads counted at its median"""
    end_to_counts = {}
    curr_group = []
    for end in sorted(ends):
        if curr_group and end[0] - curr_group[-1][0] > window:
            _add_end_group(end_to_counts, curr_group, allsamples)
            curr_group = []
        curr_group.append(end)
    if curr_group:
        _add_end_group(end_to_counts, curr_group, allsamples)
    return end_to_counts


def _add_end_group(end_to_counts, group, allsamples):
    counts = end_to_counts.setdefault(int(median([e[0] for e in group])), zero_counts(allsamples))
    for _, sample in group:
        counts[sample] += 1


def extract_terminal(juncs, strand, terminal, allsamples, sample, numreads, terminal_exons):
    "a chain's first (afe) or last (ale) exon, as its splice site and junction"
    if (terminal == 'afe') == (strand == '+'):
        junc, site = juncs[0], juncs[0][0]
    else:
        junc, site = juncs[-1], juncs[-1][1]
    add_reads(terminal_exons.setdefault(site, {}), junc, allsamples, sample, numreads)


def extract_tandem_splicing_info(sjc, alljuncs, outer_junc_to_exons, allsamples, sample, numreads):
    flat_ss = [ss for junc in sjc for ss in junc]
    for j in alljuncs:
        if j not in sjc and j[0] in flat_ss and j[1] in flat_ss:
            inner_to_counts = outer_junc_to_exons.setdefault(j, {})
            inner_juncs = []
            outer_left, outer_right = False, False
            for j2 in sjc:
                if j[0] <= j2[0] and j2[1] <= j[1]:
                    inner_juncs.append(j2)
                if j[0] == j2[0]:
                    outer_left = True
                if j[1] == j2[1]:
                    outer_right = True
            # this checks for weird junction patterns - requires outer junction to be used as an outer junction
            if outer_left and outer_right:
                add_reads(inner_to_counts, tuple(inner_juncs), allsamples, sample, numreads)
    return outer_junc_to_exons


def extract_terminal_info(allsamples, chains_by_key, key, strand, terminal):
    """the junctions, first (afe) or last (ale) exons and read starts (or ends) of the
    chains of key, a 5' (or 3') gene: (position, sample) per read.  Its fragments,
    under no junctions, have only their starts (or ends)"""
    alljuncs, terminal_exons, ends = {}, {}, []
    take_start = (terminal == 'afe') == (strand == '+')
    for sample, chains in zip(allsamples, chains_by_key):
        for juncs, readinfo in chains.get(key, {}).items():
            numreads = len(readinfo)
            if juncs:
                for j in juncs:
                    add_reads(alljuncs, j, allsamples, sample, numreads)
                extract_terminal(juncs, strand, terminal, allsamples, sample, numreads, terminal_exons)
            ends.extend(((r.start if take_start else r.end), sample) for r in readinfo)
    return alljuncs, terminal_exons, ends


def extract_splicing_info(allsamples, chains_by_key, key, strand):
    """the splice sites, junctions, exons and blocks of the chains of key, a cluster,
    for the events of splicing"""
    ss5to3, ss3to5, alljuncs, exonjpairs, allblocks = {}, {}, {}, {}, {}
    outer_junc_to_exons = {}
    for sample, chains in zip(allsamples, chains_by_key):
        for juncs, readinfo in chains.get(key, {}).items():
            numreads = len(readinfo)
            extract_35ss_info(juncs, strand, allsamples, sample, numreads, ss5to3, ss3to5, alljuncs)
            extract_exon_usage_info(juncs, allsamples, sample, numreads, allblocks, exonjpairs)
            extract_end_coverage_info(juncs, readinfo, allsamples, sample, allblocks)

    # the exons included between the ends of a skipping junction, once all the
    # junctions are known: a chain counted before the skipping junction was seen, as
    # in a sample listed before the samples with skipping, was left out
    for sample, chains in zip(allsamples, chains_by_key):
        for juncs, readinfo in chains.get(key, {}).items():
            extract_tandem_splicing_info(juncs, alljuncs, outer_junc_to_exons, allsamples, sample, len(readinfo))

    return ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, outer_junc_to_exons


def write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, gene, mycolor, min_read_support, interval_to_reads, junc_frac_of_event, event_support, outer_junc_to_exons):  # noqa: C901 - FIXME: reduce complexity
    exon_to_outer_juncs = {}
    exon_to_all_inc = {}
    exon_to_all_exc = {}
    event_to_info = {}
    for prevjunc, nextjunc in exonjpairs:
        outerjunc = (prevjunc[0], nextjunc[1])
        exon = (prevjunc[1], nextjunc[0])
        if exon not in exon_to_outer_juncs:
            exon_to_outer_juncs[exon] = set()
            exon_to_all_exc[exon] = zero_counts(allsamples)
        exon_to_outer_juncs[exon].add(outerjunc)
        add_counts_at(exon_to_all_inc, exon, exonjpairs[(prevjunc, nextjunc)], allsamples)
        if outerjunc in alljuncs:
            add_counts_at(exon_to_all_exc, exon, alljuncs[outerjunc], allsamples)

    seen_junc_combos = set()
    esjuncs = {}
    for exon in exon_to_outer_juncs:
        inc, exc = exon_to_all_inc[exon], exon_to_all_exc[exon]
        tot_counts = add_counts_to_dict(inc, exc)
        if supported(inc, allsamples, min_read_support) and supported(exc, allsamples, min_read_support) \
                and psi_varies(sample_psis(inc, tot_counts, allsamples, event_support), junc_frac_of_event):
            outer_juncs = exon_to_outer_juncs[exon]
            ename = event_name('es', thischrom, exon[0], exon[1], strand)

            # the outer junctions whose inclusion is a large enough part of the event
            goodouterjuncs = set()
            innerjuncs = set()
            for outerjunc in outer_juncs:
                ejpair = exonjpairs[((outerjunc[0], exon[0]), (exon[1], outerjunc[1]))]
                if supported(ejpair, allsamples, min_read_support) \
                        and some_psi_between(ejpair, tot_counts, allsamples, junc_frac_of_event):
                    goodouterjuncs.add(outerjunc)
                    innerjuncs.add((outerjunc[0], exon[0]))
                    innerjuncs.add((exon[1], outerjunc[1]))

            beds = []
            for outerjunc in goodouterjuncs:
                beds.append(junction_bed(thischrom, outerjunc, [exon], row_name(ename, 'inc'), strand, mycolor))
                inc_counts = exonjpairs[((outerjunc[0], exon[0]), (exon[1], outerjunc[1]))]
                add_counts_at(esjuncs.setdefault((outerjunc[0], exon[0]), {}), outerjunc, inc_counts, allsamples)
                add_counts_at(esjuncs.setdefault((exon[1], outerjunc[1]), {}), outerjunc, inc_counts, allsamples)

            event = event_to_info[ename] = SplicingEvent(ename, 'es', gene, thischrom, strand, tot_counts, allsamples)
            event.introns = {j for outerjunc in outer_juncs
                             for j in ((outerjunc[0], exon[0]), (exon[1], outerjunc[1]))}
            event.introns.update(outerjunc for outerjunc in outer_juncs if outerjunc in alljuncs)
            event.events['inc'] = SplicingEventJunction('inc', beds, inc, innerjuncs, set(), goodouterjuncs, {exon, })
            event.events['exc'] = SplicingEventJunction('exc', [], exc, set(), innerjuncs, goodouterjuncs, {exon, })
            ref_juncs = goodouterjuncs if len(goodouterjuncs) > 0 else outer_juncs
            event.totoverlap = get_overlapping_reads((max(x[0] for x in ref_juncs), min(x[1] for x in ref_juncs)), interval_to_reads, allsamples)
            seen_junc_combos.add(frozenset(goodouterjuncs | innerjuncs))

    # exons skipped together: the chains between the ends of a skipping junction
    for outerjunc, innerjuncs_to_counts in outer_junc_to_exons.items():
        exc = alljuncs[outerjunc]
        if not supported(exc, allsamples, min_read_support):
            continue
        tot_counts = sum_counts([exc, *innerjuncs_to_counts.values()], allsamples)
        if not psi_varies(sample_psis(exc, tot_counts, allsamples, event_support), junc_frac_of_event):
            continue
        goodB = set()
        all_exons = set()
        my_juncs = set()
        othercounts = zero_counts(allsamples)
        for innerjuncs, counts in innerjuncs_to_counts.items():
            if supported(counts, allsamples, min_read_support) \
                    and psi_varies(sample_psis(counts, tot_counts, allsamples, event_support), junc_frac_of_event):
                goodB.add(innerjuncs)
                all_exons.update(internal_exons(innerjuncs))
                my_juncs.update(innerjuncs)
            else:
                othercounts = add_counts_to_dict(othercounts, counts)

        # only require one exon relative to outer junction
        if len(goodB) > 0 and frozenset(my_juncs | {outerjunc, }) not in seen_junc_combos:
            ename = event_name('ces', thischrom, outerjunc[0], outerjunc[1], strand)
            event = event_to_info[ename] = SplicingEvent(ename, 'ces', gene, thischrom, strand, tot_counts, allsamples)
            event.introns = {outerjunc}.union(*innerjuncs_to_counts)
            bed = junction_bed(thischrom, outerjunc, [], row_name(ename, 'exc'), strand, mycolor)
            event.events['exc'] = SplicingEventJunction('exc', [bed], exc, set(), my_juncs, {outerjunc, }, all_exons)

            for innerjuncs in goodB:
                my_exons = internal_exons(innerjuncs)
                # the exons included
                jname = coord_ranges(my_exons)
                bed = junction_bed(thischrom, outerjunc, my_exons, row_name(ename, jname), strand, mycolor)
                event.events[jname] = SplicingEventJunction(jname, [bed], innerjuncs_to_counts[innerjuncs],
                                                            set(innerjuncs), my_juncs - set(innerjuncs), {outerjunc, }, set(my_exons))
                for j in innerjuncs:
                    add_counts_at(esjuncs.setdefault(j, {}), outerjunc, innerjuncs_to_counts[innerjuncs], allsamples)

            event.keep_other(othercounts, min_read_support)
            event.totoverlap = get_overlapping_reads(outerjunc, interval_to_reads, allsamples)

    return event_to_info, esjuncs


class ReadCoverage:
    """The reads of a strand of a region, per sample, by their spans rounded to 10
    bases, for the reads covering an event's locus, whatever their gene"""
    def __init__(self, allsamples):
        self.allsamples = allsamples
        self.spans = {}
        self.index = None

    def add(self, start, end, sample, numreads):
        add_reads(self.spans, (round(start, -1), round(end, -1)), self.allsamples, sample, numreads)

    def covering(self, ref):
        "reads per sample whose spans contain ref, (start, end)"
        if self.index is None:
            self.index = IntervalIndex()
            for span, counts in self.spans.items():
                self.index.add(span[0], span[1], (span, counts))
            self.index.build()
        # a span containing ref overlaps the range between its ends, which can be in either
        # order, as a terminal exon event's (its introns' last start, first end)
        low, high = min(ref), max(ref)
        return sum_counts([counts for span, counts in self.index.overlap(low, max(high, low + 1))
                           if span[0] < ref[0] and ref[1] < span[1]], self.allsamples)


def get_overlapping_reads(ref_junc, coverage, allsamples):
    "the reads per sample covering ref_junc, of a ReadCoverage"
    return coverage.covering(ref_junc)


def add_counts_to_dict(a, b):
    return {s: a[s] + b[s] for s in a.keys()}

def near_annotated_site(site, annot_sites):
    "is site within TERMINAL_JUNCTION_DIST_FROM_ANNOT of one of annot_sites, a sorted list"
    i = bisect_left(annot_sites, site - TERMINAL_JUNCTION_DIST_FROM_ANNOT)
    return i < len(annot_sites) and annot_sites[i] <= site + TERMINAL_JUNCTION_DIST_FROM_ANNOT


def is_terminal_site(junc_counts, alljuncs, allsamples, min_read_support, annotated):
    """Is a splice site of first (or last) junctions really a first (or last) exon's:
    pooled over the samples, at least min_read_support reads have its junctions as
    their first (or last), and few of its junctions' reads have an exon before (or
    after) it, more being allowed of an annotated one"""
    terminal = sum(counts[s] for counts in junc_counts.values() for s in allsamples)
    total = sum(alljuncs[junc][s] for junc in junc_counts for s in allsamples)
    max_fraction = MAX_NONTERMINAL_FRACTION_ANNOTATED if annotated else MAX_NONTERMINAL_FRACTION
    return terminal >= min_read_support and total - terminal <= max_fraction * total


def process_terminal_exons(termExon, eventtype, thischrom, strand, gene, allsamples, mycolor, min_read_support, junc_frac_of_event, event_support, interval_to_reads, alljuncs, annot_terminal_ss):  # noqa: C901 - FIXME: reduce complexity
    """Alternative first (afe) or last (ale) exons of a gene, from termExon, its reads'
    first or last junctions by their splice site.  Only the sites really of a first
    or last exon (is_terminal_site) are its alternatives, and only their reads are
    the event's: a truncated read's first junction, after an internal exon, is
    neither"""
    event_to_info = {}
    junc_combos = set()
    annot_ss = annot_terminal_ss.get(strand, [])
    ss_to_counts = {termSS: sum_counts(junc_counts.values(), allsamples) for termSS, junc_counts in termExon.items()
                    if is_terminal_site(junc_counts, alljuncs, allsamples, min_read_support,
                                        near_annotated_site(termSS, annot_ss))}
    if len(ss_to_counts) <= 1:
        return event_to_info, junc_combos
    tot_counts = sum_counts(ss_to_counts.values(), allsamples)
    if not supported(tot_counts, allsamples, event_support):
        return event_to_info, junc_combos

    goodB, goodBstrict = [], []
    my_juncs = set()
    othercounts = zero_counts(allsamples)
    for termSS in ss_to_counts:
        counts = ss_to_counts[termSS]
        if supported(counts, allsamples, min_read_support) \
                and psi_varies(sample_psis(counts, tot_counts, allsamples, event_support), junc_frac_of_event):
            goodB.append(termSS)
            if any(counts[s] / tot_counts[s] >= MIN_TERMINAL_SS_FRACTION for s in allsamples if tot_counts[s] > 0):
                goodBstrict.append(termSS)
            for termJunc, junc_counts in termExon[termSS].items():
                if supported(junc_counts, allsamples, min_read_support) \
                        and some_psi_between(junc_counts, tot_counts, allsamples, junc_frac_of_event):
                    my_juncs.add(termJunc)
        else:
            othercounts = add_counts_to_dict(othercounts, counts)

    # check that the junctions involved are distant enough from each other
    if len(goodB) > 1 and len(goodBstrict) > 1 and max(goodBstrict) - min(goodBstrict) >= MIN_TERMINAL_JUNCTION_SEPARATION:
        # each first or last exon by its introns: those of the event's, or failing that,
        # all its own
        ss_juncs = {ssB: ({x for x in my_juncs if ssB in x} or set(termExon[ssB])) for ssB in goodB}
        ename = event_name(eventtype, thischrom, *span([j for juncs in ss_juncs.values() for j in juncs]), strand)
        event = event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)
        for ssB in goodB:
            jname = coord_ranges(ss_juncs[ssB])
            # 100 bases of the exon, on its side of the splice site
            if (eventtype == 'afe' and strand == '+') or (eventtype == 'ale' and strand == '-'):
                bs, be = ssB - 100, ssB
            else:
                bs, be = ssB, ssB + 100
            inc_juncs = {x for x in my_juncs if ssB in x}
            event.events[jname] = SplicingEventJunction(jname, [site_bed(thischrom, bs, be, row_name(ename, jname), strand, mycolor)],
                                                        ss_to_counts[ssB], inc_juncs, my_juncs - inc_juncs)
        event.keep_other(othercounts, min_read_support)
        junc_combos.add(frozenset(my_juncs))
        ref_junc = (max(x[0] for x in my_juncs), min(x[1] for x in my_juncs))
        event.totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, allsamples)

    # TODO: add initial check for all reads with junction are equal to all terminal exon reads!
    # TODO: identify good juncs for each terminal SS, save junc combo
    # Good juncs must be significant proportion of all reads for event!
    return event_to_info, junc_combos


def process_junction_events(ssAtoB, esjuncs, eventtype, thischrom, strand, gene, allsamples, colordict, min_read_support, interval_to_reads, junc_frac_of_event, event_support, afe_ale_junc_comb):  # noqa: C901 - FIXME: reduce complexity
    """Splice sites B alternative to each other from one site A (alt5 or alt3, by
    eventtype), and the sites A sharing the same alternative sites B (alt3ss or alt5ss)"""
    event_to_info = {}
    ssB_groups_to_ssA = {}
    for ssA, ssB_counts in ssAtoB.items():
        if len(ssB_counts) <= 1:  # don't report when no alternative splicing at all
            continue
        goodB, goodBstrict = [], []
        my_juncs = set()
        othercounts = zero_counts(allsamples)
        tot_counts = sum_counts(ssB_counts.values(), allsamples)
        for ssB, counts in ssB_counts.items():
            junc = junc_between(ssA, ssB)
            if psi_varies(sample_psis(counts, tot_counts, allsamples, event_support), junc_frac_of_event):
                goodB.append(ssB)
                my_juncs.add(junc)
                if junc not in esjuncs:
                    goodBstrict.append(ssB)
                else:
                    # a junction of exon skipping counts only for its reads apart from the skipping
                    other_juncs = [junc_between(ssA, ssB2) for ssB2 in ssB_counts if ssB2 != ssB]
                    escounts = sum_counts([esjuncs[junc][j] for j in other_juncs if j in esjuncs[junc]], allsamples)
                    if any(counts[s] > 0 and (counts[s] - escounts[s]) / counts[s] > 0.1 and counts[s] - escounts[s] > min_read_support
                           for s in allsamples):
                        goodBstrict.append(ssB)
            else:
                othercounts = add_counts_to_dict(othercounts, counts)

        if frozenset(my_juncs) in afe_ale_junc_comb:
            continue
        if len(goodB) > 1:
            ssB_groups_to_ssA.setdefault(frozenset(goodB), {})[ssA] = tot_counts
        if supported(tot_counts, allsamples, event_support) and len(goodBstrict) > 1:  # only if alt splicing, not including junctions involved in exon skipping
            # from the shared site to the furthest alternative one, so unique to the shared site
            ename = event_name(eventtype, thischrom, *span([junc_between(ssA, ssB) for ssB in goodB]), strand)
            event = event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)
            event.introns = {junc_between(ssA, ssB) for ssB in ssB_counts}
            for ssB in goodB:
                junc = junc_between(ssA, ssB)
                jname = coord_ranges([junc])
                bed = site_bed(thischrom, junc[0], junc[1], row_name(ename, jname), strand, colordict[eventtype])
                event.events[jname] = SplicingEventJunction(jname, [bed], ssB_counts[ssB], {junc, }, my_juncs - {junc, })
            event.keep_other(othercounts, min_read_support)
            # all goodB lie on the same side of ssA
            ref_junc = (ssA, min(goodB)) if ssA < min(goodB) else (max(goodB), ssA)
            event.totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, allsamples)

    eventtype = 'alt3' if eventtype == 'alt5' else 'alt5'
    for ssBgroup, ssA_counts in ssB_groups_to_ssA.items():
        if len(ssA_counts) <= 1:
            continue
        tot_counts = sum_counts(ssA_counts.values(), allsamples)
        if not supported(tot_counts, allsamples, event_support):
            continue
        good = []
        my_juncs = set()
        othercounts = zero_counts(allsamples)
        for ssA, counts in ssA_counts.items():
            if psi_varies(sample_psis(counts, tot_counts, allsamples, event_support), junc_frac_of_event):
                good.append(ssA)
                my_juncs.update(junc_between(ssA, ssB) for ssB in ssBgroup)
            else:
                othercounts = add_counts_to_dict(othercounts, counts)
        if len(good) > 1:
            # unique: each site A is in one group, and the outer bounds end at the group's
            # outermost site A
            ename = event_name(f'{eventtype}ss', thischrom, *span(my_juncs), strand)
            event = event_to_info[ename] = SplicingEvent(ename, f'{eventtype}ss', gene, thischrom, strand, tot_counts, allsamples)
            event.introns = {junc_between(ssA, ssB) for ssA in ssA_counts for ssB in ssBgroup}
            for ssA in good:
                if (eventtype == 'alt5' and strand == '+') or (eventtype == 'alt3' and strand == '-'):
                    bs, be = ssA - 10, ssA
                else:
                    bs, be = ssA, ssA + 10
                inc_juncs = {junc_between(ssA, ssB) for ssB in ssBgroup}
                jname = coord_ranges(inc_juncs)
                bed = site_bed(thischrom, bs, be, row_name(ename, jname), strand, colordict[eventtype])
                event.events[jname] = SplicingEventJunction(jname, [bed], ssA_counts[ssA], inc_juncs, my_juncs - inc_juncs)
            event.keep_other(othercounts, min_read_support)
            ref_junc = (max(good), min(ssBgroup)) if median(good) < median(ssBgroup) else (max(ssBgroup), min(good))
            event.totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, allsamples)

    return event_to_info


class SplicingEvent:
    """Represents metadata about a splicing event"""
    def __init__(self, name, eventtype, gene, chrom, strand, totjunc, allsamples):
        self.name = name
        self.chrom = chrom
        self.eventtype = eventtype
        self.gene = gene
        self.strand = strand
        self.totjunc = totjunc
        self.events = {}
        self.totoverlap = None
        self.other = {s: 0 for s in allsamples}
        # all the introns of its alternatives, those of other too, for its gene
        self.introns = set()
        # the other annotated genes sharing an intron or splice site with it
        self.additional_genes = []

    def keep_other(self, othercounts, min_reads):
        "count the reads of the event's other junctions or sites, if some sample has min_reads of them"
        if any(othercounts[s] >= min_reads for s in othercounts):
            self.other = othercounts

class SplicingEventJunction:
    """represents metadata about a specific junction in a splicing event"""
    def __init__(self, name, beds, samplecounts, inc_junc, exc_junc, outer_junc=None, inc_exon=None):
        self.name = name
        self.beds = beds
        self.inc_juncs = inc_junc
        self.exc_juncs = exc_junc
        self.samplecounts = samplecounts
        self.inc_exon = set()
        self.outer_juncs = set()
        if outer_junc is not None:
            self.outer_juncs = outer_junc
        if inc_exon is not None:
            self.inc_exon = inc_exon


def write_intron_retention(alljuncs, allsamples, allblocks, thischrom, strand, gene, mycolor, min_read_support, interval_to_reads, junc_frac_of_event, event_support):
    event_to_info = {}
    for j, splicedcounts in alljuncs.items():
        # this nested loop could get brutal
        retainedcounts = sum_counts([counts for b, counts in allblocks.items() if b[0] < j[0] and j[1] < b[1]], allsamples)
        tot_counts = add_counts_to_dict(splicedcounts, retainedcounts)
        if supported(splicedcounts, allsamples, min_read_support) and supported(retainedcounts, allsamples, min_read_support) \
                and psi_varies(sample_psis(retainedcounts, tot_counts, allsamples, event_support), junc_frac_of_event):
            ename = event_name('ir', thischrom, j[0], j[1], strand)
            event = event_to_info[ename] = SplicingEvent(ename, 'ir', gene, thischrom, strand, tot_counts, allsamples)
            event.introns = {j}
            bed = junction_bed(thischrom, j, [], row_name(ename, 'spliced'), strand, mycolor)
            event.events['spliced'] = SplicingEventJunction('spliced', [bed], splicedcounts, {j, }, set())
            bed = site_bed(thischrom, j[0], j[1], row_name(ename, 'retained'), strand, mycolor)
            event.events['retained'] = SplicingEventJunction('retained', [bed], retainedcounts, set(), {j, })
            event.totoverlap = get_overlapping_reads(j, interval_to_reads, allsamples)
    return event_to_info


def write_ends(grouped_ends, allsamples, thischrom, strand, gene, eventtype, mycolor, support, junc_frac_of_event, event_support):
    event_to_info = {}
    tot_counts = sum_counts(grouped_ends.values(), allsamples)
    if not supported(tot_counts, allsamples, event_support):
        return event_to_info
    good_ends = set()
    othercounts = zero_counts(allsamples)
    for e, counts in grouped_ends.items():
        if supported(counts, allsamples, support) \
                and psi_varies(sample_psis(counts, tot_counts, allsamples, event_support), junc_frac_of_event):
            good_ends.add(e)
        else:
            othercounts = add_counts_to_dict(othercounts, counts)
    if len(good_ends) > 1:
        # the 1-based position of each end's base: the reads' starts are 0-based, on +
        # their tss and on - their tts
        from_starts = (eventtype == 'tss') == (strand == '+')
        positions = {e: e + 1 if from_starts else e for e in good_ends}
        ename = event_name(eventtype, thischrom, min(positions.values()) - 1, max(positions.values()), strand)
        event = event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)
        for e, pos in positions.items():
            jname = str(pos)
            bed = site_bed(thischrom, pos - 1, pos, row_name(ename, jname), strand, mycolor)
            event.events[jname] = SplicingEventJunction(jname, [bed], grouped_ends[e], {pos, }, set(positions.values()) - {pos, })
        event.keep_other(othercounts, support)
    return event_to_info


def write_counts_psi(info, ncounts, junctot, fulltot, allsamples, outcounts, outpsijunc, outpsitot, event_support):
    outline = info + [str(ncounts[s]) for s in allsamples]
    outcounts.write('\t'.join(outline) + '\n')
    juncpsi = [ncounts[s] / junctot[s] if junctot[s] >= event_support else 'NA' for s in allsamples]
    outline = info + [str(round(ncounts[s] / junctot[s], 4)) if junctot[s] >= event_support else '' for s in allsamples]
    outpsijunc.write('\t'.join(outline) + '\n')
    outline = info + [str(round(ncounts[s] / fulltot[s], 4)) if fulltot[s] >= event_support else '' for s in allsamples]
    outpsitot.write('\t'.join(outline) + '\n')
    return juncpsi

def get_junc_string(chrom, juncs):
    if len(juncs) == 0:
        return ''
    else:
        juncs = sorted(list(juncs))
        if isinstance(juncs[0], int):
            return ','.join([f'{chrom}:{x}' for x in juncs])
        else:
            return ','.join([f'{chrom}:{x[0] + 1}-{x[1]}' for x in juncs])


def write_events(event_to_info, allsamples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot,  # noqa: C901 - FIXME: reduce complexity
                 outtotals, event_support, annotated):
    """Write the events expressed enough: each row's counts and two PSIs, over the
    event's reads (PSIjunc) and over the reads covering its locus (PSItot), with the
    row's introns that aren't annotated, its BED lines scored by its median PSIjunc,
    and the event's two totals, which spliceevents_outliers divides by"""
    for ename in event_to_info:
        event = event_to_info[ename]
        if event.totoverlap is None:
            event.totoverlap = event.totjunc

        if any([event_frac_of_tot <= event.totjunc[s] / event.totoverlap[s] if event.totoverlap[s] > 0 else False for s in allsamples]):  # check that junction/event is expressed enough relative to all reads covering locus
            # do an extra check that one junction that is not exclusion is at least junc_frac_of_event of total locus
            one_event_pass = False
            for jname in event.events:
                jinfo = event.events[jname]
                inc_frac_list = sample_psis(jinfo.samplecounts, event.totoverlap, allsamples, 1)

                if (event.eventtype not in {'es', 'ces'} or jname != 'exc') \
                        and psi_varies(inc_frac_list, junc_frac_of_event):
                    one_event_pass = True

            if one_event_pass:
                genes = [event.gene, annotated.gene_names(event.gene.split(',')),
                         ','.join(event.additional_genes), annotated.gene_names(event.additional_genes)]
                for jname in event.events:
                    jinfo = event.events[jname]

                    novel = {j for j in row_introns(event, jname) if not annotated.is_annotated(j, event.strand)}
                    outinfo = [ename, jname, event.eventtype, *genes, get_junc_string(event.chrom, jinfo.inc_juncs),
                               get_junc_string(event.chrom, jinfo.exc_juncs), get_junc_string(event.chrom, jinfo.outer_juncs),
                               get_junc_string(event.chrom, jinfo.inc_exon), get_junc_string(event.chrom, novel)]
                    juncpsi = write_counts_psi(outinfo, jinfo.samplecounts, event.totjunc, event.totoverlap, allsamples, outcounts, outpsijunc, outpsitot, event_support)

                    # every sample can be below event_support, leaving no PSI to
                    # take a median of and nothing to score the BED with
                    vals = [x for x in juncpsi if x != 'NA']
                    if len(vals) > 0:
                        med = median(vals)
                        for bed in jinfo.beds:
                            bed.score = round(med * 100)
                            bed.write(outbed)
                if any([event.other[s] > 0 for s in allsamples]):
                    outinfo = [ename, 'other', event.eventtype, *genes, '', '', '', '', '']
                    write_counts_psi(outinfo, event.other, event.totjunc, event.totoverlap, allsamples, outcounts, outpsijunc, outpsitot, event_support)
                for total, counts in (('event_reads', event.totjunc), ('locus_reads', event.totoverlap)):
                    outtotals.write('\t'.join([ename, event.eventtype, *genes, total] + [str(counts[s]) for s in allsamples]) + '\n')


ETYPE_TO_COLOR = {'skipped_exons': '66,105,245', 'retained_introns': '144,66,245',
                  'alt3': '245,215,66', 'alt5': '43,184,39', 'tss': '255,0,0', 'tts': '0,0,255'}


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


def group_chains(chrom, allsamples, sample_chains, novel_gene):
    """Each sample's reads grouped three ways, each a list per sample of
    key -> junctions -> [ReadEnds]: by cluster of junction chains linked by shared
    splice sites, on one strand, for the events of splicing, which don't depend on
    genes; by each chain's 5'-most gene, for the first exons and starts; and by its
    3'-most gene, for the last exons and ends.  A gene's 5' fragments are with its
    starts and its 3' fragments with its ends, under no junctions, and nowhere else.
    Also each key's strand, each chain's genes, by (strand, junctions), a novel gene's
    if of none annotated, and each strand's ReadCoverage, of its spliced reads"""
    chain_genes = {}
    chain_isoforms = {}
    coverage = {'+': ReadCoverage(allsamples), '-': ReadCoverage(allsamples)}
    for sample, chains in zip(allsamples, sample_chains):
        for chain in chains:
            if chain.fragment:
                continue
            key = (chain.strand, chain.juncs)
            chain_genes[key] = chain.genes or (novel_gene[key],)
            if key not in chain_isoforms:
                chain_isoforms[key] = Isoform(chrom, chain.strand, tuple(Junc(*j) for j in chain.juncs))
            chain_isoforms[key].reads.extend(chain.ends)
            for r in chain.ends:
                coverage[chain.strand].add(r.start, r.end, sample, 1)
    for isoform in chain_isoforms.values():
        isoform.start, isoform.end = int(median(isoform.starts)), int(median(isoform.ends))
    cluster_of, key_strand = _splicing_clusters(chain_isoforms)
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
    row = event.events[alternative]
    if event.eventtype in ('es', 'ces') and alternative == 'exc':
        return set(row.outer_juncs)
    return {j for j in row.inc_juncs if isinstance(j, tuple)}


def event_gene_scores(introns, strand, annotated, junc_reads):
    """the genes sharing an intron or splice site with an event, on its strand, and their
    scores: 2 for an intron the gene has annotated, or failing that 1 for each of its
    splice sites the gene has, each intron weighted by its share of the reads of all
    the introns, junc_reads, so a rare alternative counts little"""
    total = sum(junc_reads.get(junc, 0) for junc in introns)
    scores = Counter()
    for junc in introns:
        weight = junc_reads.get(junc, 0) / total if total > 0 else 1
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


def process_region_events(temp_prefix, thischrom, allsamples, sample_chains, novel_gene, junc_support, event_frac_of_tot, junc_frac_of_event, event_support,  # noqa: C901 - FIXME: reduce complexity
                          annot_afe_ss, annot_ale_ss, report_regions, gene_juncs, annotated=None):
    """A region's events, from each sample's ReadChains, written to its tables.  The
    events of splicing are found per cluster of junction chains linked by shared
    splice sites, without genes; the alternative first exons and starts per each
    chain's 5'-most gene, and the last exons and ends per its 3'-most gene, so a
    readthrough's reads count with the gene each end lies in.  An event's locus
    reads are all the strand's reads covering it.  A splicing event's gene is the one
    annotated with most of its introns, those of all its alternatives, and their
    splice sites, weighted by their reads (annotated_event_genes), or, with none, that
    of most of its reads
    (read_event_genes); several only on a tie.
    Its first and last exons and ends are of the gene they were grouped by"""
    annotated = annotated or AnnotatedSplicing()
    views, key_strand, chain_genes, coverage = group_chains(thischrom, allsamples, sample_chains, novel_gene)
    events = []  # (event, the span of its gene's junctions, for the regions)
    afe_ale_junc_comb = {'+': set(), '-': set()}
    # first exons and starts by 5' gene, last exons and ends by 3' gene
    for view, terminal, annot_terminal_ss, end_type in (('five', 'afe', annot_afe_ss, 'tss'),
                                                        ('three', 'ale', annot_ale_ss, 'tts')):
        for gene in _keys(views[view]):
            strand = key_strand[gene]
            cov = coverage[strand]
            alljuncs, terminal_exons, ends = extract_terminal_info(allsamples, views[view], gene, strand, terminal)
            gene_span = gene_junction_span(gene_juncs.get(gene, ()))
            info, junc_comb = process_terminal_exons(terminal_exons, terminal, thischrom, strand, gene, allsamples,
                                                     ETYPE_TO_COLOR[end_type], junc_support, junc_frac_of_event, event_support, cov,
                                                     alljuncs, annot_terminal_ss)
            afe_ale_junc_comb[strand] |= junc_comb
            info |= write_ends(group_ends(ends, allsamples, 100), allsamples, thischrom, strand, gene, end_type, ETYPE_TO_COLOR[end_type],
                               junc_support, junc_frac_of_event, event_support)
            events.extend((event, gene_span) for event in info.values())
    # the events of splicing, by cluster
    for cluster in sorted(_keys(views['cluster']), key=lambda c: int(c[1:])):
        strand = key_strand[cluster]
        label = ''  # each event's gene, once it is found
        cov = coverage[strand]
        ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, outer_junc_to_exons = extract_splicing_info(allsamples, views['cluster'], cluster, strand)
        es_info, esjuncs = write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, label,
                                               ETYPE_TO_COLOR['skipped_exons'], junc_support, cov, junc_frac_of_event, event_support, outer_junc_to_exons)
        a3_info = process_junction_events(ss5to3, esjuncs, 'alt3', thischrom, strand, label, allsamples, ETYPE_TO_COLOR,
                                          junc_support, cov, junc_frac_of_event, event_support, afe_ale_junc_comb[strand])
        a5_info = process_junction_events(ss3to5, esjuncs, 'alt5', thischrom, strand, label, allsamples, ETYPE_TO_COLOR,
                                          junc_support, cov, junc_frac_of_event, event_support, afe_ale_junc_comb[strand])
        ir_info = write_intron_retention(alljuncs, allsamples, allblocks, thischrom, strand, label,
                                         ETYPE_TO_COLOR['retained_introns'], junc_support, cov, junc_frac_of_event, event_support)
        cluster_chains = [chains.get(cluster, {}) for chains in views['cluster']]
        junc_reads = {junc: sum(counts.values()) for junc, counts in alljuncs.items()}
        for event in (es_info | a3_info | a5_info | ir_info).values():
            scores = event_gene_scores(event.introns, strand, annotated, junc_reads)
            genes = annotated_event_genes(scores, annotated) or read_event_genes(event.introns, cluster_chains, chain_genes, strand)
            event.gene = ','.join(genes)
            # the other genes sharing a splice site, best first
            event.additional_genes = sorted(set(scores) - set(genes), key=lambda gene: (-scores[gene], gene))
            events.append((event, (0, 0)))

    with open(temp_prefix + '.spliceevents.bed', 'w') as outbed, open(temp_prefix + '.spliceevents.counts.tsv', 'w') as outcounts, \
            open(temp_prefix + '.spliceevents.PSIjunc.tsv', 'w') as outpsijunc, open(temp_prefix + '.spliceevents.PSItot.tsv', 'w') as outpsitot, \
            open(temp_prefix + '.spliceevents.totals.tsv', 'w') as outtotals:
        # only the events inside the regions asked for, whose reads are all there
        reported = {}
        for event, gene_span in events:
            if event_in_regions(event, report_regions, gene_span):
                if event.name in reported:
                    raise AssertionError(f'two events named {event.name}, of {reported[event.name].gene} and {event.gene}')
                reported[event.name] = event
        write_events(reported, allsamples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot,
                     outtotals, event_support, annotated)


def correct_sample_reads(args, region, work_dir, bamfile_name, genome, annots, junction_corrector):
    """The sample's reads in region, corrected as flair transcriptome corrects them,
    grouped by junction chain: aligned to the annotated transcripts, then their
    junctions corrected, the unspliced reads under no junctions.  The alignment's
    intermediates are written in work_dir, removed after"""
    os.makedirs(work_dir, exist_ok=True)
    temp_prefix = osp.join(work_dir, 'reads')
    pipettor.run([('samtools', 'view', '-h', bamfile_name, f'{region.name}:{region.start}-{region.end}'),
                  ('samtools', 'fasta', '-')],
                 stdout=temp_prefix + '.reads.fasta')
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
                                   check_motifs=not args.trust_junctions)
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


def get_juncs_single_sample(args, region, path, bamfile_name, genome, junction_corrector, annots, terminal_exons):
    """correct a sample's reads in region, give their junction chains annotated genes,
    and write them to path, with the unspliced reads that are fragments of one gene's
    first or last exons, terminal_exons; the chains of no annotated gene get novel
    genes in event calling, over all samples"""
    sj_to_ends = correct_sample_reads(args, region, path + '.correction', bamfile_name, genome, annots, junction_corrector)
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
    novel = {}
    for chains in sample_chains:
        for chain in chains:
            if not chain.genes:
                key = (chain.strand, chain.juncs)
                if key not in novel:
                    novel[key] = Isoform(chrom, chain.strand, tuple(Junc(*j) for j in chain.juncs))
                novel[key].reads.extend(chain.ends)
    for isoform in novel.values():
        isoform.start, isoform.end = int(median(isoform.starts)), int(median(isoform.ends))
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
    """chrom -> the first and last exons, TerminalExons, of the spliced transcripts, as
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
            site = exon.end if (exon.kind == 'first') == (exon.strand == '+') else exon.start
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
            if (exon.kind == 'last') == (exon.strand == '+'):  # the splice site is the exon's start
                at_site = exon.start - FRAGMENT_SPLICE_SITE_WINDOW <= start < exon.end
            else:
                at_site = exon.start < end <= exon.end + FRAGMENT_SPLICE_SITE_WINDOW
            if at_site:
                found.add(exon.kind)
        if len(found) != 1:
            return None
        return found.pop(), gene, gene_strand


def _correct_region_sample(*, partition, gtf_data, junction_corrector, args, item, terminal_exons):
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
            get_juncs_single_sample(args, region, path, bamfile, genome, junction_corrector, annots, region_terminal)


def _run_region(*, partition, gtf_data, junction_corrector, args, allsamples, region_reports, annot_afe_ss, annot_ale_ss):  # noqa: C901 - FIXME: reduce complexity
    annots = annot_data_from_gtf(gtf_data, partition.region)
    if len(annots.transcripts) > 0:  # check if any annotated transcripts in region
        # each sample's corrected reads, from _correct_region_sample
        sample_chains = [read_read_chains(read_chains_path(partition, sample)) for sample, _ in allsamples]
        novel_gene = novel_genes_of_chains(partition.region.name, sample_chains)
        process_region_events(
            partition.file_prefix, partition.region.name, [x[0] for x in allsamples],
            sample_chains, novel_gene, args.feature_support,
            args.event_frac_of_tot, args.junc_frac_of_event, args.event_support,
            annot_afe_ss.get(partition.region.name, {}), annot_ale_ss.get(partition.region.name, {}),
            region_reports[_region_key(partition.region)], gene_junctions(sample_chains, novel_gene, annots),
            AnnotatedSplicing(annots))

        if not args.keep_intermediate:
            for sample, _ in allsamples:
                os.remove(read_chains_path(partition, sample))
    else:
        # no annotated genes, so no events: empty tables, as every partition's are combined
        for suffix in REGION_OUTPUT_SUFFIXES:
            open(partition.file_prefix + suffix, 'w').close()


# the tables each region writes, combined into the run's
REGION_OUTPUT_SUFFIXES = ('.spliceevents.bed', '.spliceevents.counts.tsv', '.spliceevents.PSIjunc.tsv',
                          '.spliceevents.PSItot.tsv', '.spliceevents.totals.tsv')


def _region_key(region):
    return (region.name, region.start, region.end)


def correction_regions(report_regions, bamfiles):
    """The regions to correct reads in, for the regions events are reported in: each
    report region extended over the reads overlapping it in any sample, so a read
    running past its end is corrected with the annotation and junctions it reaches,
    as anywhere else.  Extended regions that overlap are merged.  Returns the
    regions, and a map of each one's key to the report regions in it."""
    extended = []
    for region in report_regions:
        start, end = region.start, region.end
        for bamfile in bamfiles:
            with pysam.AlignmentFile(bamfile, 'rb') as bam:
                for read in bam.fetch(region.name, region.start, region.end):
                    if read.is_mapped and not read.is_secondary:
                        start, end = min(start, read.reference_start), max(end, read.reference_end)
        extended.append((SeqRange(region.name, start, end), region))
    merged = combine_regions([ext for ext, _ in extended])
    reports = {_region_key(region): [] for region in merged}
    for ext, region in extended:
        container = next(m for m in merged if m.name == ext.name and m.start <= ext.start and ext.end <= m.end)
        reports[_region_key(container)].append(region)
    return merged, reports


# events whose reads are all of a gene's: its starts and ends, first and last exons
GENE_WIDE_EVENTS = {'tss', 'tts', 'afe', 'ale'}


def _event_positions(event):
    "the splice sites and ends of an event's junctions, exons and ends"
    positions = []
    for jinfo in event.events.values():
        for feature in (*jinfo.inc_juncs, *jinfo.exc_juncs, *jinfo.outer_juncs, *jinfo.inc_exon):
            positions.extend(feature if isinstance(feature, tuple) else (feature,))
    return positions


def gene_junction_span(juncs):
    """The span of a gene's junctions, its annotated ones, or a novel gene's of its
    reads, from the first one's start to the last one's end"""
    if not juncs:
        return (0, 0)
    return (min(j[0] for j in juncs), max(j[1] for j in juncs))


def event_in_regions(event, regions, gene_span):
    """Is an event inside one of the regions events are reported in: all its
    junctions and exons, so every read reaching it overlaps the region and its
    counts are complete.  An event of all its gene's reads (GENE_WIDE_EVENTS) only
    when all the gene's annotated junctions are, gene_span (gene_junction_span),
    as a user choosing a region around a gene would think of it: its reads are
    spliced, so nearly all have a junction in the region and were fetched.  The
    gene's first and last exons may run past the region"""
    if event.eventtype in GENE_WIDE_EVENTS:
        low, high = gene_span
    else:
        positions = _event_positions(event)
        if not positions:
            return False
        low, high = min(positions), max(positions)
    return any(region.name == event.chrom and region.start <= low and high <= region.end for region in regions)


def combine_regions(regions, buffersize=0):
    regions.sort()
    new_regions = []
    lastchrom, laststart, lastend = -1, -1, -1
    for range in regions:
        c, s, e = range.name, range.start, range.end
        if c != lastchrom or s > lastend + buffersize:
            if lastchrom != -1:
                new_regions.append(SeqRange(lastchrom, laststart, lastend))
            lastchrom, laststart, lastend = c, s, e
        else:
            lastend = max((lastend, e))
    if lastchrom != -1:
        new_regions.append(SeqRange(lastchrom, laststart, lastend))
    return new_regions


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


def spliceevents(args):  # noqa: C901 - FIXME: reduce complexity
    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)

    # FIXME: make this a common function somewhere
    tempDir = make_temp_dir(args.output)
    logging.info('temp directory: %s', tempDir)

    allsamples = []
    for line in open(args.manifest):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            sample, bamfile = line
            allsamples.append((sample, bamfile))

    if args.region_bed:
        report_regions = [SeqRange(bed.chrom, bed.chromStart, bed.chromEnd) for bed in BedReader(args.region_bed)]
        all_regions, region_reports = correction_regions(combine_regions(report_regions),
                                                         [bamfile for _, bamfile in allsamples])
    else:
        all_regions = [SeqRange(chrom, 0, genome.get_reference_length(chrom))
                       for chrom in genome.references]
        region_reports = {_region_key(region): [region] for region in all_regions}
    logging.info(f'input regions: {len(all_regions)}')

    # with --region_bed, only the annotation of the regions reads are corrected in, the
    # transcripts each region's partition would take, rather than the whole GTF's
    annot_regions = all_regions if args.region_bed else None
    logging.info('loading annotation GTF')
    annot_gtf_data = gtf_data_parser(args.annot, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES,
                                     regions=annot_regions)
    logging.info('building intron support database')
    junction_corrector = junction_corrector_factory(args.ss_window, args.junction_support,
                                                    annot_gtf_data=annot_gtf_data,
                                                    intron_beds=args.junction_bed,
                                                    star_sj_tabs=args.junction_tab)

    # first and last exons, which unspliced fragments lie in, and, of --annot_basic,
    # which first and last exons may also be internal ones
    annot_afe_ss, annot_ale_ss = {}, {}
    if args.annot_basic:
        logging.info('loading first and last exons of --annot_basic')
        terminal_exons = annotated_terminal_exons(gtf_data_parser(args.annot_basic, attrs=GtfAttrsSet.FLAIR,
                                                                  include_features=TRANSCRIPT_EXON_FEATURES,
                                                                  regions=annot_regions).transcripts)
        annot_afe_ss, annot_ale_ss = annotated_terminal_splice_sites(terminal_exons)
    else:
        terminal_exons = annotated_terminal_exons(annot_gtf_data.transcripts)

    logging.info(f'correcting reads of {len(allsamples)} samples in {len(all_regions)} regions, then finding events, with {args.threads} threads')

    runner = PartitionRunner(all_regions, tempDir, gtf_data=annot_gtf_data, junction_corrector=junction_corrector,
                             threads=args.threads)
    runner.run(_correct_region_sample, items=allsamples, args=args, terminal_exons=terminal_exons)
    runner.run(_run_region, args=args, allsamples=allsamples, region_reports=region_reports,
               annot_afe_ss=annot_afe_ss, annot_ale_ss=annot_ale_ss)

    prefixes = [p.file_prefix for p in runner]
    combine_temp_files_by_suffix(args.output, prefixes, ['.spliceevents.bed'])
    samples = [x[0] for x in allsamples]
    write_with_header(args.output, prefixes, '.spliceevents.counts.tsv', EVENT_ROW_COLUMNS + samples)
    write_with_header(args.output, prefixes, '.spliceevents.PSIjunc.tsv', EVENT_ROW_COLUMNS + samples)
    write_with_header(args.output, prefixes, '.spliceevents.PSItot.tsv', EVENT_ROW_COLUMNS + samples)
    write_with_header(args.output, prefixes, '.spliceevents.totals.tsv', EVENT_TOTALS_COLUMNS + samples)

    if not args.keep_intermediate:
        shutil.rmtree(tempDir)

    genome.close()
