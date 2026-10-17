import os
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
from flair.isoform_data import ReadRec
from flair.io_utils import make_temp_dir
from flair.read_processing import generate_genomic_alignment_read_to_clipping_file
from flair.read_correction import filter_correct_group_reads
from flair.annotation_alignment import identify_good_match_to_annot
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import annot_data_from_gtf
from flair.pycbio.hgdata.bed import Bed, BedBlock, BedReader


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
                             'first and last splice sites mark alternative first and last exons.  A '
                             'comprehensive annotation includes truncated transcripts, whose ends are not '
                             'real first and last exons')
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
    parser.add_argument('--output_read_ends', action='store_true',
                        help='write OUTPUT.spliceevents.readends.bed, the corrected ends of the reads')
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
MIN_TERMINAL_JUNCTION_SEPARATION = 100
MIN_TERMINAL_SS_FRACTION = 0.1

def get_annot_gene_hits(gene_to_exons, exons):
    gene_hits = []
    for annotgene in gene_to_exons:
        annotexons = gene_to_exons[annotgene]
        if min((annotexons[-1][1], exons[-1][1])) > max(
                (annotexons[0][0], exons[0][0])):  # there is overlap in the genes
            coveredpos = set()
            for s, e in exons:
                for ast, ae in annotexons:
                    for p in range(max((ast, s)), min((ae, e))):
                        coveredpos.add(p)
            if len(coveredpos) > sum([x[1] - x[0] for x in exons]) * 0.5:
                gene_hits.append([len(coveredpos), annotgene])
    return gene_hits

def get_juncs_to_gene(juncs, isoinfo, sjc_to_gene, junc_to_gene, gene_to_exons, gene_to_juncs):
    thisgene = None
    if juncs in sjc_to_gene:
        thisgene = sjc_to_gene[juncs]
    else:
        gene_hits = {}
        for j in juncs:
            if j in junc_to_gene:
                gene_id = junc_to_gene[j]
                if gene_id not in gene_hits:
                    gene_hits[gene_id] = [0, -1 * len(gene_to_juncs[gene_id])]
                gene_hits[gene_id][0] += 1
        if gene_hits:
            sortedgenes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)
            thisgene = sortedgenes[0][0]
        else:
            # look for exon overlap.  isoinfo is an Isoform, which is not iterable;
            # its reads are what was wanted.  The median is how isoform_data turns a
            # read cluster into one isoform, and neither extreme is outlier proof
            mystart = int(median(isoinfo.starts))
            myend = int(median(isoinfo.ends))
            exons = [(mystart, juncs[0][0])] + [(juncs[i][1], juncs[i + 1][0]) for i in range(len(juncs) - 1)] + [
                (juncs[-1][1], myend)]
            # strand = isoinfo[0][2]  # not a super robust strand picking, assumes well stranded reads
            gene_hits = get_annot_gene_hits(gene_to_exons, exons)
            if len(gene_hits) > 0:
                gene_hits.sort(reverse=True)
                thisgene = gene_hits[0][1]
    return thisgene

def group_juncs_by_annot_gene(sjtoends, sjc_to_gene, junc_to_gene, gene_to_exons, gene_to_juncs):
    genetojuncs, nogenejuncs, sereads = {}, {}, []
    no_gene_reads, se_reads_cnt = 0, 0
    for (chrom, juncs) in sjtoends:
        if len(juncs) > 0:  # skip single-exon reads
            thisgene = get_juncs_to_gene(juncs, sjtoends[(chrom, juncs)], sjc_to_gene, junc_to_gene, gene_to_exons, gene_to_juncs)
            if thisgene:
                if thisgene not in genetojuncs:
                    genetojuncs[thisgene] = {}
                genetojuncs[thisgene][juncs] = sjtoends[(chrom, juncs)]
            else:
                no_gene_reads += sjtoends[(chrom, juncs)].num_reads
                nogenejuncs[juncs] = sjtoends[(chrom, juncs)]
        else:
            se_reads_cnt += sjtoends[(chrom, juncs)].num_reads
            sereads.extend(sjtoends[(chrom, juncs)].reads)
    logging.debug(f"gene assignment: {no_gene_reads} spliced reads not assigned to a gene, {se_reads_cnt} single-exon reads discarded")
    return genetojuncs, nogenejuncs, sereads


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


def extract_end_coverage_info(juncs, readinfo, allsamples, sample, thischrom, gene, strand, allblocks, outends):
    for r in readinfo:
        # getting all block coverage to measure intron retention
        add_reads(allblocks, (r.start, juncs[0][0]), allsamples, sample, 1)
        add_reads(allblocks, (juncs[-1][1], r.end), allsamples, sample, 1)
        if outends:
            Bed(thischrom, r.start, r.end, name=gene + '|' + r.name,
                score=0, strand=strand).write(outends)
    return allblocks

def determine_juncs_are_subset(juncs, alljuncs):
    is_subset = False
    unique_seq_bound = [[], []]
    for otherjuncs in alljuncs:
        if juncs != otherjuncs and len(juncs) < len(otherjuncs):
            if str(juncs)[1:-1].rstrip(',') in str(otherjuncs):
                is_subset = True
                other_internal_exons = [(otherjuncs[x][1], otherjuncs[x + 1][0]) for x in range(len(otherjuncs) - 1)]
                for other_exon in other_internal_exons:
                    if juncs[0][0] == other_exon[1]:
                        unique_seq_bound[0].append(other_exon[1] - other_exon[0])
                    if juncs[-1][1] == other_exon[0]:
                        unique_seq_bound[1].append(other_exon[1] - other_exon[0])
    if is_subset:
        unique_seq_bound[0] = max(unique_seq_bound[0]) if len(unique_seq_bound[0]) > 0 else None
        unique_seq_bound[1] = max(unique_seq_bound[1]) if len(unique_seq_bound[1]) > 0 else None
        return unique_seq_bound
    else:
        return [None, None]

def get_start_end_counts(readinfo, subset_info, sample, first_junc, last_junc, t_starts_ends, t_first_last_sj, allsamples):
    for r in readinfo:
        # if subset_info[0] == None or first_junc-start > subset_info[0]:
        t_starts_ends[0].append((r.start, sample))
        add_reads(t_first_last_sj[0], first_junc, allsamples, sample, 1)
        # if subset_info[1] == None or end-last_junc > subset_info[1]:
        t_starts_ends[1].append((r.end, sample))
        add_reads(t_first_last_sj[1], last_junc, allsamples, sample, 1)
    return t_starts_ends, t_first_last_sj


def group_ends(t_starts_ends, allsamples, window):
    grouped_ends = [[], []]
    for i in range(2):
        end_to_counts = {}
        last_end = 0
        curr_group = []
        # all_groups = []
        for end in sorted(t_starts_ends[i]):
            if end[0] - last_end > window:
                if len(curr_group) > 0:
                    # all_groups.append(curr_group)
                    g = int(median([x[0] for x in curr_group]))
                    end_to_counts[g] = {s: 0 for s in allsamples}
                    for e in curr_group:
                        end_to_counts[g][e[1]] += 1
                curr_group = []
            curr_group.append(end)
            last_end = end[0]
        # all_groups.append(curr_group)
        # end_to_counts[median(curr_group)] = len(curr_group)
        g = int(median([x[0] for x in curr_group]))
        end_to_counts[g] = {s: 0 for s in allsamples}
        for e in curr_group:
            end_to_counts[g][e[1]] += 1

        grouped_ends[i] = end_to_counts
    return grouped_ends


def extract_afe_ale(juncs, strand, allsamples, sample, numreads, alljuncs, afe, ale):
    if strand == '+':
        afe_junc, ale_junc = juncs[0], juncs[-1]
        afe_ss, ale_ss = juncs[0][0], juncs[-1][1]
    else:
        ale_junc, afe_junc = juncs[0], juncs[-1]
        ale_ss, afe_ss = juncs[0][0], juncs[-1][1]

    add_reads(afe.setdefault(afe_ss, {}), afe_junc, allsamples, sample, numreads)
    add_reads(ale.setdefault(ale_ss, {}), ale_junc, allsamples, sample, numreads)

    return afe, ale


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


def extract_splicing_info(allsamples, all_genes_to_juncs, gene, strand, thischrom, outends):
    ss5to3, ss3to5, alljuncs, exonjpairs, allblocks = {}, {}, {}, {}, {}
    afe, ale = {}, {}
    t_starts_ends, t_first_last_sj = [[], []], [{}, {}]
    interval_to_reads = {}
    outer_junc_to_exons = {}
    # get valid read ends by excluding subset junction chains (check for whether end is within exon)
    # only check relative to other reads, not reference?
    # get list of valid read ends after excluding subset junction chains
    # cluster within window, take median/mode of window (use code from transcriptome)?

    # process each sample for each gene - only ever storing data for one gene at a time
    for s in range(len(allsamples)):
        sample = allsamples[s]
        # TODO: change to genetojuncs being all juncs from all samples!!
        genetojuncs = all_genes_to_juncs[s]
        if gene in genetojuncs:
            for juncs in genetojuncs[gene]:
                readinfo = genetojuncs[gene][juncs]
                numreads = len(readinfo)

                # get intervals for calculating coverage
                for r in readinfo:
                    add_reads(interval_to_reads, (round(r.start, -1), round(r.end, -1)), allsamples, sample, 1)

                ss5to3, ss3to5, alljuncs = extract_35ss_info(juncs, strand, allsamples, sample, numreads,
                                                             ss5to3, ss3to5, alljuncs)
                afe, ale = extract_afe_ale(juncs, strand, allsamples, sample, numreads, alljuncs, afe, ale)
                allblocks, exonjpairs = extract_exon_usage_info(juncs, allsamples, sample, numreads,
                                                                allblocks, exonjpairs)
                allblocks = extract_end_coverage_info(juncs, readinfo, allsamples, sample, thischrom, gene, strand,
                                                      allblocks, outends)
                subset_info = determine_juncs_are_subset(juncs, genetojuncs[gene])
                t_starts_ends, t_first_last_sj = get_start_end_counts(readinfo, subset_info, sample, juncs[0][0], juncs[-1][1],
                                                                      t_starts_ends, t_first_last_sj, allsamples)

    # the exons included between the ends of a skipping junction, once all the gene's
    # junctions are known: a chain counted before the skipping junction was seen, as
    # in a sample listed before the samples with skipping, was left out
    for s, sample in enumerate(allsamples):
        for juncs, readinfo in all_genes_to_juncs[s].get(gene, {}).items():
            outer_junc_to_exons = extract_tandem_splicing_info(juncs, alljuncs, outer_junc_to_exons, allsamples, sample, len(readinfo))

    return ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj, interval_to_reads, afe, ale, outer_junc_to_exons

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
            ename = f'es-of-{thischrom}:{exon[0]}-{exon[1]}'

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
                beds.append(junction_bed(thischrom, outerjunc, [exon], f'inc_{ename}', strand, mycolor))
                inc_counts = exonjpairs[((outerjunc[0], exon[0]), (exon[1], outerjunc[1]))]
                add_counts_at(esjuncs.setdefault((outerjunc[0], exon[0]), {}), outerjunc, inc_counts, allsamples)
                add_counts_at(esjuncs.setdefault((exon[1], outerjunc[1]), {}), outerjunc, inc_counts, allsamples)

            event = event_to_info[ename] = SplicingEvent(ename, 'es', gene, thischrom, strand, tot_counts, allsamples)
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
            ename = f'ces-relTo-{thischrom}:{outerjunc[0]}-{outerjunc[1]}({strand})-{gene}'
            event = event_to_info[ename] = SplicingEvent(ename, 'ces', gene, thischrom, strand, tot_counts, allsamples)
            bed = junction_bed(thischrom, outerjunc, [], f'exc_{ename}', strand, mycolor)
            event.events['exc'] = SplicingEventJunction('exc', [bed], exc, set(), my_juncs, {outerjunc, }, all_exons)

            for innerjuncs in goodB:
                my_exons = internal_exons(innerjuncs)
                jname = 'inc-of-' + ','.join([f'{thischrom}:{x[0]}-{x[1]}' for x in my_exons])
                bed = junction_bed(thischrom, outerjunc, my_exons, f'{jname}_{ename}', strand, mycolor)
                event.events[jname] = SplicingEventJunction(jname, [bed], innerjuncs_to_counts[innerjuncs],
                                                            set(innerjuncs), my_juncs - set(innerjuncs), {outerjunc, }, set(my_exons))
                for j in innerjuncs:
                    add_counts_at(esjuncs.setdefault(j, {}), outerjunc, innerjuncs_to_counts[innerjuncs], allsamples)

            event.keep_other(othercounts, min_read_support)
            event.totoverlap = get_overlapping_reads(outerjunc, interval_to_reads, allsamples)

    return event_to_info, esjuncs


def get_overlapping_reads(ref_junc, interval_to_reads, allsamples):
    tot_counts = {s: 0 for s in allsamples}
    for i in interval_to_reads:
        if i[0] < ref_junc[0] and ref_junc[1] < i[1]:
            tot_counts = add_counts_to_dict(tot_counts, interval_to_reads[i])
    return tot_counts


def add_counts_to_dict(a, b):
    return {s: a[s] + b[s] for s in a.keys()}

def process_terminal_exons(termExon, eventtype, thischrom, strand, gene, allsamples, mycolor, min_read_support, junc_frac_of_event, event_support, interval_to_reads, alljuncs, annot_terminal_ss):  # noqa: C901 - FIXME: reduce complexity
    event_to_info = {}
    junc_combos = set()
    if len(termExon) <= 1:
        return event_to_info, junc_combos

    # the splice sites really of a terminal exon: annotated as one, or every read of each
    # of their junctions starting or ending there
    annot_ss = annot_terminal_ss.get(gene, ())
    ss_to_counts = {}
    actual_terminal = []
    othercounts = zero_counts(allsamples)
    for termSS, junc_counts in termExon.items():
        ss_to_counts[termSS] = sum_counts(junc_counts.values(), allsamples)
        if any(abs(termSS - ss) <= TERMINAL_JUNCTION_DIST_FROM_ANNOT for ss in annot_ss) \
                or all(any(counts[s] == alljuncs[termJunc][s] and counts[s] >= min_read_support for s in allsamples)
                       for termJunc, counts in junc_counts.items()):
            actual_terminal.append(termSS)
        else:
            othercounts = add_counts_to_dict(othercounts, ss_to_counts[termSS])
    tot_counts = sum_counts(ss_to_counts.values(), allsamples)
    if not supported(tot_counts, allsamples, event_support):
        return event_to_info, junc_combos

    goodB, goodBstrict = [], []
    my_juncs = set()
    for termSS in actual_terminal:
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
        ename = f'{eventtype}-({strand})-{gene}'
        event = event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)
        for ssB in goodB:
            jname = f'{thischrom}:{ssB}'
            # 100 bases of the exon, on its side of the splice site
            if (eventtype == 'AFE' and strand == '+') or (eventtype == 'ALE' and strand == '-'):
                bs, be = ssB - 100, ssB
            else:
                bs, be = ssB, ssB + 100
            inc_juncs = {x for x in my_juncs if ssB in x}
            event.events[jname] = SplicingEventJunction(jname, [site_bed(thischrom, bs, be, f'{jname}_{ename}', strand, mycolor)],
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
            ename = f'{eventtype}-relTo-{thischrom}:{ssA}({strand})-{gene}'
            event = event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)
            for ssB in goodB:
                junc = junc_between(ssA, ssB)
                jname = f'{thischrom}:{junc[0]}-{junc[1]}'
                bed = site_bed(thischrom, junc[0], junc[1], f'{jname}_{ename}', strand, colordict[eventtype])
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
            oj_string = [f'{thischrom}:{x}' for x in sorted(list(ssBgroup))]
            ename = f'{eventtype}ss-relTo-{",".join(oj_string)}({strand})-{gene}'
            event = event_to_info[ename] = SplicingEvent(ename, f'{eventtype}ss', gene, thischrom, strand, tot_counts, allsamples)
            for ssA in good:
                if (eventtype == 'alt5' and strand == '+') or (eventtype == 'alt3' and strand == '-'):
                    bs, be = ssA - 10, ssA
                else:
                    bs, be = ssA, ssA + 10
                jname = f'{thischrom}:{ssA}'
                inc_juncs = {junc_between(ssA, ssB) for ssB in ssBgroup}
                bed = site_bed(thischrom, bs, be, f'{jname}_{ename}', strand, colordict[eventtype])
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
            ename = f'ir-of-{thischrom}:{j[0]}-{j[1]}({strand})-{gene}'
            event = event_to_info[ename] = SplicingEvent(ename, 'ir', gene, thischrom, strand, tot_counts, allsamples)
            bed = junction_bed(thischrom, j, [], f'spliced_{ename}', strand, mycolor)
            event.events['spliced'] = SplicingEventJunction('spliced', [bed], splicedcounts, {j, }, set())
            bed = site_bed(thischrom, j[0], j[1], f'retained_{ename}', strand, mycolor)
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
        ename = f'{eventtype}-({strand})-{gene}'
        event = event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)
        for e in good_ends:
            jname = f'{thischrom}:{e}'
            bed = site_bed(thischrom, e - 1, e, f'{jname}_{ename}', strand, mycolor)
            event.events[jname] = SplicingEventJunction(jname, [bed], grouped_ends[e], {e, }, good_ends - {e, })
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
            return ','.join([f'{chrom}:{x[0]}-{x[1]}' for x in juncs])


def write_events(event_to_info, allsamples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot,  # noqa: C901 - FIXME: reduce complexity
                 outtotals, event_support):
    """Write the events expressed enough: each row's counts and two PSIs, over the
    event's reads (PSIjunc) and over the reads covering its locus (PSItot), its BED
    lines scored by its median PSIjunc, and the event's two totals, which
    spliceevents_outliers divides by"""
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
                for jname in event.events:
                    jinfo = event.events[jname]

                    outinfo = [f'{jname}_{ename}', ename, event.eventtype, event.gene, get_junc_string(event.chrom, jinfo.inc_juncs),
                               get_junc_string(event.chrom, jinfo.exc_juncs), get_junc_string(event.chrom, jinfo.outer_juncs), get_junc_string(event.chrom, jinfo.inc_exon)]
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
                    outinfo = [f'other_{ename}', ename, event.eventtype, event.gene, '', '', '', '']
                    write_counts_psi(outinfo, event.other, event.totjunc, event.totoverlap, allsamples, outcounts, outpsijunc, outpsitot, event_support)
                for total, counts in (('event_reads', event.totjunc), ('locus_reads', event.totoverlap)):
                    outtotals.write('\t'.join([ename, event.eventtype, event.gene, total] + [str(counts[s]) for s in allsamples]) + '\n')


def process_gene_to_events(temp_prefix, thischrom, allsamples, all_genes_to_juncs, genetostrand, junc_support, output_read_ends, event_frac_of_tot, junc_frac_of_event, event_support, annot_afe_ss, annot_ale_ss,
                           report_regions, gene_to_annot_juncs):
    etypetocolor = {'skipped_exons': '66,105,245', 'retained_introns': '144,66,245',
                    'alt3': '245,215,66', 'alt5': '43,184,39', 'tss': '255,0,0', 'tts': '0,0,255'}
    # outends is just for writing out read ends for Harrison - he can group them more intelligently
    outends = None
    if output_read_ends:
        outends = open(temp_prefix + '.spliceevents.readends.bed', 'w')

    allgenes = set.union(*[set(all_genes_to_juncs[s].keys()) for s in range(len(all_genes_to_juncs))])
    with open(temp_prefix + '.spliceevents.bed', 'w') as outbed, open(temp_prefix + '.spliceevents.counts.tsv', 'w') as outcounts, \
            open(temp_prefix + '.spliceevents.PSIjunc.tsv', 'w') as outpsijunc, open(temp_prefix + '.spliceevents.PSItot.tsv', 'w') as outpsitot, \
            open(temp_prefix + '.spliceevents.totals.tsv', 'w') as outtotals:

        # sorted, as a set of names is ordered differently in each run
        for gene in sorted(allgenes):

            # get gene strand - this is not optimal, def need to revamp how getting annot info
            strand = genetostrand[gene]
            ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj, interval_to_reads, afe, ale, outer_junc_to_exons = extract_splicing_info(allsamples, all_genes_to_juncs,
                                                                                                                                                                      gene, strand, thischrom, outends)

            # exon skipping
            # TODO combine exon skipping at multiple junctions (so all junctions that skip exon are combined)
            es_info, esjuncs = write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, gene,
                                                   etypetocolor['skipped_exons'], junc_support, interval_to_reads, junc_frac_of_event, event_support, outer_junc_to_exons)
            afe_info, afe_junc_comb = process_terminal_exons(afe, 'AFE', thischrom, strand, gene, allsamples, etypetocolor['tss'], junc_support, junc_frac_of_event, event_support, interval_to_reads, alljuncs, annot_afe_ss)
            ale_info, ale_junc_comb = process_terminal_exons(ale, 'ALE', thischrom, strand, gene, allsamples, etypetocolor['tts'], junc_support, junc_frac_of_event, event_support, interval_to_reads, alljuncs, annot_ale_ss)
            afe_ale_junc_comb = afe_junc_comb | ale_junc_comb
            a3_info = process_junction_events(ss5to3, esjuncs, 'alt3', thischrom, strand, gene, allsamples, etypetocolor,
                                              junc_support, interval_to_reads, junc_frac_of_event, event_support, afe_ale_junc_comb)
            a5_info = process_junction_events(ss3to5, esjuncs, 'alt5', thischrom, strand, gene, allsamples, etypetocolor,
                                              junc_support, interval_to_reads, junc_frac_of_event, event_support, afe_ale_junc_comb)
            # intron retention
            ir_info = write_intron_retention(alljuncs, allsamples, allblocks, thischrom, strand, gene,
                                             etypetocolor['retained_introns'], junc_support, interval_to_reads, junc_frac_of_event, event_support)

            grouped_ends = group_ends(t_starts_ends, allsamples, 100)
            if strand == '-':
                grouped_ends = grouped_ends[::-1]
                t_first_last_sj = t_first_last_sj[::-1]
            tss_info = write_ends(grouped_ends[0], allsamples, thischrom, strand, gene, 'tss', etypetocolor['tss'], junc_support, junc_frac_of_event, event_support)
            tts_info = write_ends(grouped_ends[1], allsamples, thischrom, strand, gene, 'tts', etypetocolor['tts'], junc_support, junc_frac_of_event, event_support)
            # write_ends(t_first_last_sj[0], allsamples, thischrom, strand, gene, 'tss', etypetocolor['tss'], out, outbed, support)
            # write_ends(t_first_last_sj[1], allsamples, thischrom, strand, gene, 'tts', etypetocolor['tts'], out, outbed, support)
            all_info = (((a3_info | a5_info) | (es_info | ir_info)) | (tts_info | tss_info)) | (afe_info | ale_info)
            # only the events inside the regions asked for, whose reads are all there
            gene_span = gene_junction_span(gene_to_annot_juncs.get(gene, ()))
            all_info = {ename: event for ename, event in all_info.items()
                        if event_in_regions(event, report_regions, gene_span)}

            write_events(all_info, allsamples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot,
                         outtotals, event_support)

    if output_read_ends:
        outends.close()


def correct_sample_reads(args, region, work_dir, bamfile_name, genome, annots, junction_corrector):
    """The sample's spliced reads in region, corrected as flair transcriptome corrects
    them, grouped by junction chain: aligned to the annotated transcripts, then their
    junctions corrected, single-exon reads, often fragments, left out.  The
    alignment's intermediates are written in work_dir, removed after"""
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
                                   keep_single_exon=False, trust_strand=args.trust_strand,
                                   check_motifs=not args.trust_junctions)
    shutil.rmtree(work_dir)
    return sj_to_ends


def get_juncs_single_sample(args, region, temp_prefix, sample, bamfile_name, genome, junction_corrector, annots):
    temp_prefix = temp_prefix + '_' + sample
    sj_to_ends = correct_sample_reads(args, region, temp_prefix + '.correction', bamfile_name, genome, annots, junction_corrector)
    genetojuncs, nogenejuncs, sereads = group_juncs_by_annot_gene(sj_to_ends, annots.sjc_to_gene, annots.junc_to_gene_id, annots.gene_to_exons, annots.gene_to_annot_juncs)

    with open(temp_prefix + '_gene_to_juncs.txt', 'w') as out:
        for gene in genetojuncs:
            for juncs in genetojuncs[gene]:
                juncstring = ','.join(['.'.join([str(y) for y in x]) for x in juncs])
                for read_info in genetojuncs[gene][juncs].reads:
                    outline = [gene, juncstring, str(read_info.start), str(read_info.end), genetojuncs[gene][juncs].strand, read_info.name]
                    out.write('\t'.join(outline) + '\n')


def annotated_terminal_splice_sites(gtf_file, regions=None):
    """gene -> the splice sites of the first exons, and of the last exons, of its spliced
    transcripts in gtf_file, such as GENCODE basic, which mark its alternative first
    and last exons; given regions, of the transcripts overlapping them"""
    annot_afe_ss, annot_ale_ss = {}, {}
    for trans in gtf_data_parser(gtf_file, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES,
                                 regions=regions).transcripts:
        exons = sorted((exon.start, exon.end) for exon in trans.exons)
        if len(exons) > 1:
            afe_ss, ale_ss = exons[0][1], exons[-1][0]
            if trans.strand == '-':
                afe_ss, ale_ss = ale_ss, afe_ss
            annot_afe_ss.setdefault(trans.gene_id, set()).add(afe_ss)
            annot_ale_ss.setdefault(trans.gene_id, set()).add(ale_ss)
    return annot_afe_ss, annot_ale_ss


def _correct_region_sample(*, partition, gtf_data, junction_corrector, args, item):
    """Correct one sample's reads in a region, for _run_region.  A task of its own,
    so the samples of a region are corrected in parallel"""
    sample, bamfile = item
    if os.path.exists(partition.file_prefix + '_' + sample + '_gene_to_juncs.txt'):
        return
    annots = annot_data_from_gtf(gtf_data, partition.region)
    if len(annots.transcripts) > 0:
        with pysam.FastaFile(args.genome) as genome:
            get_juncs_single_sample(args, partition.region, partition.file_prefix, sample, bamfile, genome, junction_corrector, annots)


def _run_region(*, partition, gtf_data, junction_corrector, args, allsamples, region_reports, annot_afe_ss, annot_ale_ss):  # noqa: C901 - FIXME: reduce complexity
    annots = annot_data_from_gtf(gtf_data, partition.region)
    if len(annots.transcripts) > 0:  # check if any annotated transcripts in region
        # each sample's corrected reads, from _correct_region_sample
        all_genes_to_juncs = []
        for sample, bamfile in allsamples:
            gene_to_juncs = {}
            for line in open(partition.file_prefix + '_' + sample + '_gene_to_juncs.txt'):
                line = line.rstrip().split('\t')
                gene, juncstring, start, end, strand, readname = line
                juncs = [x.split('.') for x in juncstring.split(',')]
                juncs = tuple([(int(x[0]), int(x[1])) for x in juncs])
                if gene not in gene_to_juncs:
                    gene_to_juncs[gene] = {}
                if juncs not in gene_to_juncs[gene]:
                    gene_to_juncs[gene][juncs] = []
                gene_to_juncs[gene][juncs].append(ReadRec(None, strand, (), int(start), int(end), readname))
            all_genes_to_juncs.append(gene_to_juncs)

        process_gene_to_events(
            partition.file_prefix, partition.region.name, [x[0] for x in allsamples],
            all_genes_to_juncs, annots.gene_to_strand, args.feature_support, args.output_read_ends,
            args.event_frac_of_tot, args.junc_frac_of_event, args.event_support,
            annot_afe_ss, annot_ale_ss,
            region_reports[_region_key(partition.region)], annots.gene_to_annot_juncs)

        if not args.keep_intermediate:
            for sample, bamfile in allsamples:
                pipettor.run([('rm', partition.file_prefix + '_' + sample + '_gene_to_juncs.txt')])
    else:
        # no annotated genes, so no events: empty tables, as every partition's are combined
        suffixes = REGION_OUTPUT_SUFFIXES + (('.spliceevents.readends.bed',) if args.output_read_ends else ())
        for suffix in suffixes:
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
GENE_WIDE_EVENTS = {'tss', 'tts', 'AFE', 'ALE'}


def _event_positions(event):
    "the splice sites and ends of an event's junctions, exons and ends"
    positions = []
    for jinfo in event.events.values():
        for feature in (*jinfo.inc_juncs, *jinfo.exc_juncs, *jinfo.outer_juncs, *jinfo.inc_exon):
            positions.extend(feature if isinstance(feature, tuple) else (feature,))
    return positions


def gene_junction_span(annot_juncs):
    """The span of a gene's annotated junctions, from the first one's start to the
    last one's end"""
    if not annot_juncs:
        return (0, 0)
    return (min(j.start for j in annot_juncs), max(j.end for j in annot_juncs))


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
EVENT_ROW_COLUMNS = ['eventname', 'event', 'eventtype', 'gene', 'junctions_included', 'junctions_excluded',
                     'outer_junctions', 'exons']
EVENT_TOTALS_COLUMNS = ['event', 'eventtype', 'gene', 'total']


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

    annot_afe_ss, annot_ale_ss = {}, {}
    if args.annot_basic:
        logging.info('loading first and last exons of --annot_basic')
        annot_afe_ss, annot_ale_ss = annotated_terminal_splice_sites(args.annot_basic, annot_regions)

    logging.info(f'correcting reads of {len(allsamples)} samples in {len(all_regions)} regions, then finding events, with {args.threads} threads')

    runner = PartitionRunner(all_regions, tempDir, gtf_data=annot_gtf_data, junction_corrector=junction_corrector,
                             threads=args.threads)
    runner.run(_correct_region_sample, items=allsamples, args=args)
    runner.run(_run_region, args=args, allsamples=allsamples, region_reports=region_reports,
               annot_afe_ss=annot_afe_ss, annot_ale_ss=annot_ale_ss)

    prefixes = [p.file_prefix for p in runner]
    combine_temp_files_by_suffix(args.output, prefixes, ['.spliceevents.bed'])
    if args.output_read_ends:
        combine_temp_files_by_suffix(args.output, prefixes, ['.spliceevents.readends.bed'])
    samples = [x[0] for x in allsamples]
    write_with_header(args.output, prefixes, '.spliceevents.counts.tsv', EVENT_ROW_COLUMNS + samples)
    write_with_header(args.output, prefixes, '.spliceevents.PSIjunc.tsv', EVENT_ROW_COLUMNS + samples)
    write_with_header(args.output, prefixes, '.spliceevents.PSItot.tsv', EVENT_ROW_COLUMNS + samples)
    write_with_header(args.output, prefixes, '.spliceevents.totals.tsv', EVENT_TOTALS_COLUMNS + samples)

    if not args.keep_intermediate:
        shutil.rmtree(tempDir)

    genome.close()
