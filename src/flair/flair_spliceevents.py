import os
import os.path as osp
import sys
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
from flair.pycbio.hgdata.bed import Bed, BedReader

pkgdir = osp.dirname(osp.realpath(__file__))
gtf_to_bed_prog = osp.join(pkgdir, "gtf_to_bed.py")

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
        if strand == '+':
            ss5, ss3 = j[0], j[1]
        else:
            ss5, ss3 = j[1], j[0]
        if ss5 not in ss5to3:
            ss5to3[ss5] = {}
        if ss3 not in ss3to5:
            ss3to5[ss3] = {}
        if ss3 not in ss5to3[ss5]:
            ss5to3[ss5][ss3] = {s: 0 for s in allsamples}
        if ss5 not in ss3to5[ss3]:
            ss3to5[ss3][ss5] = {s: 0 for s in allsamples}
        ss5to3[ss5][ss3][sample] += numreads
        ss3to5[ss3][ss5][sample] += numreads
        if j not in alljuncs:
            alljuncs[j] = {s: 0 for s in allsamples}
        alljuncs[j][sample] += numreads
    return ss5to3, ss3to5, alljuncs


def extract_exon_usage_info(juncs, allsamples, sample, numreads, allblocks, exonjpairs):
    # just for middle exons
    for i in range(len(juncs) - 1):
        exon = (juncs[i][1], juncs[i + 1][0])
        prevjunc, nextjunc = juncs[i], juncs[i + 1]
        if exon not in allblocks:
            allblocks[exon] = {s: 0 for s in allsamples}
        allblocks[exon][sample] += numreads
        if (prevjunc, nextjunc) not in exonjpairs:
            exonjpairs[(prevjunc, nextjunc)] = {s: 0 for s in allsamples}
        exonjpairs[(prevjunc, nextjunc)][sample] += numreads
    return allblocks, exonjpairs


def extract_end_coverage_info(juncs, readinfo, allsamples, sample, thischrom, gene, strand, allblocks, outends):
    for r in readinfo:
        # getting all block coverage to measure intron retention
        firstexon = (r.start, juncs[0][0])
        lastexon = (juncs[-1][1], r.end)
        if firstexon not in allblocks:
            allblocks[firstexon] = {s: 0 for s in allsamples}
        if lastexon not in allblocks:
            allblocks[lastexon] = {s: 0 for s in allsamples}
        allblocks[firstexon][sample] += 1
        allblocks[lastexon][sample] += 1
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
        if first_junc not in t_first_last_sj[0]:
            t_first_last_sj[0][first_junc] = {s: 0 for s in allsamples}
        t_first_last_sj[0][first_junc][sample] += 1
        # if subset_info[1] == None or end-last_junc > subset_info[1]:
        t_starts_ends[1].append((r.end, sample))
        if last_junc not in t_first_last_sj[1]:
            t_first_last_sj[1][last_junc] = {s: 0 for s in allsamples}
        t_first_last_sj[1][last_junc][sample] += 1
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

    if afe_ss not in afe:
        afe[afe_ss] = {}
    if afe_junc not in afe[afe_ss]:
        afe[afe_ss][afe_junc] = {s: 0 for s in allsamples}
    afe[afe_ss][afe_junc][sample] += numreads

    if ale_ss not in ale:
        ale[ale_ss] = {}
    if ale_junc not in ale[ale_ss]:
        ale[ale_ss][ale_junc] = {s: 0 for s in allsamples}
    ale[ale_ss][ale_junc][sample] += numreads

    return afe, ale


def extract_tandem_splicing_info(sjc, alljuncs, outer_junc_to_exons, allsamples, sample, numreads):
    flat_ss = [ss for junc in sjc for ss in junc]
    for j in alljuncs:
        if j not in sjc and j[0] in flat_ss and j[1] in flat_ss:
            if j not in outer_junc_to_exons:
                outer_junc_to_exons[j] = {}
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
                inner_juncs = tuple(inner_juncs)
                if inner_juncs not in outer_junc_to_exons[j]:
                    outer_junc_to_exons[j][inner_juncs] = {s: 0 for s in allsamples}
                outer_junc_to_exons[j][inner_juncs][sample] += numreads
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
                    interval = (round(r.start, -1), round(r.end, -1))
                    if interval not in interval_to_reads:
                        interval_to_reads[interval] = {s: 0 for s in allsamples}
                    interval_to_reads[interval][sample] += 1

                ss5to3, ss3to5, alljuncs = extract_35ss_info(juncs, strand, allsamples, sample, numreads,
                                                             ss5to3, ss3to5, alljuncs)
                afe, ale = extract_afe_ale(juncs, strand, allsamples, sample, numreads, alljuncs, afe, ale)
                outer_junc_to_exons = extract_tandem_splicing_info(juncs, alljuncs, outer_junc_to_exons, allsamples, sample, numreads)
                allblocks, exonjpairs = extract_exon_usage_info(juncs, allsamples, sample, numreads,
                                                                allblocks, exonjpairs)
                allblocks = extract_end_coverage_info(juncs, readinfo, allsamples, sample, thischrom, gene, strand,
                                                      allblocks, outends)
                subset_info = determine_juncs_are_subset(juncs, genetojuncs[gene])
                t_starts_ends, t_first_last_sj = get_start_end_counts(readinfo, subset_info, sample, juncs[0][0], juncs[-1][1],
                                                                      t_starts_ends, t_first_last_sj, allsamples)

    return ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj, interval_to_reads, afe, ale, outer_junc_to_exons

def write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, gene, mycolor, min_read_support, interval_to_reads, junc_frac_of_event, event_support, outer_junc_to_exons):  # noqa: C901 - FIXME: reduce complexity
    # outer_junc_to_exons = {}
    # outer_junc_to_all_inc = {}
    exon_to_outer_juncs = {}
    exon_to_all_inc = {}
    exon_to_all_exc = {}
    event_to_info = {}
    for prevjunc, nextjunc in exonjpairs:
        outerjunc = (prevjunc[0], nextjunc[1])
        exon = (prevjunc[1], nextjunc[0])

        if exon not in exon_to_outer_juncs:
            exon_to_outer_juncs[exon] = set()
            exon_to_all_inc[exon] = {s: 0 for s in allsamples}
            exon_to_all_exc[exon] = {s: 0 for s in allsamples}
        exon_to_outer_juncs[exon].add(outerjunc)
        exon_to_all_inc[exon] = add_counts_to_dict(exon_to_all_inc[exon], exonjpairs[(prevjunc, nextjunc)])
        if outerjunc in alljuncs:
            # if outerjunc not in outer_junc_to_exons:
            #     outer_junc_to_exons[outerjunc] = set()
            #     outer_junc_to_all_inc[outerjunc] = {s: 0 for s in allsamples}
            exon_to_all_exc[exon] = add_counts_to_dict(exon_to_all_exc[exon], alljuncs[outerjunc])
            # outer_junc_to_exons[outerjunc].add(exon)
            # outer_junc_to_all_inc[outerjunc] = add_counts_to_dict(outer_junc_to_all_inc[outerjunc], exonjpairs[(prevjunc, nextjunc)])

    seen_junc_combos = set()
    esjuncs = {}
    for exon in exon_to_outer_juncs:
        tot_counts = exon_to_all_inc[exon]
        tot_counts = add_counts_to_dict(tot_counts, exon_to_all_exc[exon])

        inc_frac_list = [exon_to_all_inc[exon][s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]
        if len(inc_frac_list) >= 1 \
                and any([exon_to_all_inc[exon][s] >= min_read_support for s in allsamples]) \
                and any([exon_to_all_exc[exon][s] >= min_read_support for s in allsamples]) \
                and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                     or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
            outer_juncs = exon_to_outer_juncs[exon]

            bedlines = []
            goodouterjuncs = set()
            innerjuncs = set()
            for outerjunc in outer_juncs:
                ejpair = exonjpairs[((outerjunc[0], exon[0]), (exon[1], outerjunc[1]))]
                if any([ejpair[s] >= min_read_support for s in allsamples]) \
                        and any([junc_frac_of_event <= ejpair[s] / (exon_to_all_inc[exon][s] + exon_to_all_exc[exon][s]) <= 1 - junc_frac_of_event
                                 if (exon_to_all_inc[exon][s] + exon_to_all_exc[exon][s]) > 0 else False for s in allsamples]):
                    goodouterjuncs.add(outerjunc)
                    innerjuncs.add((outerjunc[0], exon[0]))
                    innerjuncs.add((exon[1], outerjunc[1]))

            ename = f'es-of-{thischrom}:{exon[0]}-{exon[1]}'

            for outerjunc in goodouterjuncs:
                # FIXME: use BED class
                bedlines.append([thischrom, outerjunc[0] - 10, outerjunc[1] + 10, f'inc_{ename}', 0, strand,
                                 outerjunc[0] - 10, outerjunc[1] + 10, mycolor, 3,
                                 f'10,{exon[1] - exon[0]},10',
                                 f'0,{10 + exon[0] - outerjunc[0]},{10 + outerjunc[1] - outerjunc[0]}'])
                if (outerjunc[0], exon[0]) not in esjuncs:
                    esjuncs[(outerjunc[0], exon[0])] = {}
                if outerjunc not in esjuncs[(outerjunc[0], exon[0])]:
                    esjuncs[(outerjunc[0], exon[0])][outerjunc] = {s: 0 for s in allsamples}
                esjuncs[(outerjunc[0], exon[0])][outerjunc] = add_counts_to_dict(esjuncs[(outerjunc[0], exon[0])][outerjunc], exonjpairs[((outerjunc[0], exon[0]), (exon[1], outerjunc[1]))])
                if (exon[1], outerjunc[1]) not in esjuncs:
                    esjuncs[(exon[1], outerjunc[1])] = {}
                if outerjunc not in esjuncs[(exon[1], outerjunc[1])]:
                    esjuncs[(exon[1], outerjunc[1])][outerjunc] = {s: 0 for s in allsamples}
                esjuncs[(exon[1], outerjunc[1])][outerjunc] = add_counts_to_dict(esjuncs[(exon[1], outerjunc[1])][outerjunc], exonjpairs[((outerjunc[0], exon[0]), (exon[1], outerjunc[1]))])

            event_to_info[ename] = SplicingEvent(ename, 'es', gene, thischrom, strand, tot_counts, allsamples)

            event_to_info[ename].events['inc'] = SplicingEventJunction('inc', bedlines, exon_to_all_inc[exon], innerjuncs, set(), goodouterjuncs, {exon, })

            event_to_info[ename].events['exc'] = SplicingEventJunction('exc', [], exon_to_all_exc[exon], set(), innerjuncs, goodouterjuncs, {exon, })

            if len(goodouterjuncs) > 0:
                event_to_info[ename].totoverlap = get_overlapping_reads((max([x[0] for x in goodouterjuncs]), min([x[1] for x in goodouterjuncs])), interval_to_reads, allsamples)
            else:
                event_to_info[ename].totoverlap = get_overlapping_reads((max([x[0] for x in outer_juncs]), min([x[1] for x in outer_juncs])), interval_to_reads, allsamples)

            finaljuncs = goodouterjuncs | innerjuncs
            seen_junc_combos.add(frozenset(finaljuncs))

    for outerjunc in outer_junc_to_exons:
        innerjuncs_to_counts = outer_junc_to_exons[outerjunc]
        if any([alljuncs[outerjunc][s] >= min_read_support for s in allsamples]):
            goodB = set()
            all_exons = set()
            my_juncs = set()  # {outerjunc,}
            othercounts = {s: 0 for s in allsamples}
            tot_counts = alljuncs[outerjunc]
            for innerjuncs in innerjuncs_to_counts:
                tot_counts = add_counts_to_dict(tot_counts, innerjuncs_to_counts[innerjuncs])
            inc_frac_list = [alljuncs[outerjunc][s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]
            if len(inc_frac_list) > 0 \
                and any([alljuncs[outerjunc][s] >= min_read_support for s in allsamples]) \
                and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                     or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                for innerjuncs in innerjuncs_to_counts:
                    inc_frac_list = [innerjuncs_to_counts[innerjuncs][s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]
                    if any([innerjuncs_to_counts[innerjuncs][s] >= min_read_support for s in allsamples]) \
                        and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                             or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                        goodB.add(innerjuncs)
                        my_exons = [(innerjuncs[i][1], innerjuncs[i + 1][0]) for i in range(len(innerjuncs) - 1)]
                        all_exons.update(set(my_exons))
                        for j in innerjuncs:
                            my_juncs.add(j)
                    else:
                        othercounts = add_counts_to_dict(othercounts, innerjuncs_to_counts[innerjuncs])

                if len(goodB) > 0 \
                        and frozenset(my_juncs | {outerjunc, }) not in seen_junc_combos:  # only require one exon relative to outer junction
                    ename = f'ces-relTo-{thischrom}:{outerjunc[0]}-{outerjunc[1]}({strand})-{gene}'
                    event_to_info[ename] = SplicingEvent(ename, 'ces', gene, thischrom, strand, tot_counts, allsamples)

                    # FIXME: use BED class
                    bedline = [thischrom, outerjunc[0] - 10, outerjunc[1] + 10, f'exc_{ename}', 0, strand, outerjunc[0] - 10, outerjunc[1] + 10,
                               mycolor, 2, '10,10', f'0,{10 + outerjunc[1] - outerjunc[0]}']
                    event_to_info[ename].events['exc'] = SplicingEventJunction('exc', [bedline], alljuncs[outerjunc], set(), my_juncs, {outerjunc, }, all_exons)

                    for innerjuncs in goodB:
                        my_exons = [(innerjuncs[i][1], innerjuncs[i + 1][0]) for i in range(len(innerjuncs) - 1)]
                        exonstring = ','.join([f'{thischrom}:{x[0]}-{x[1]}' for x in my_exons])
                        jname = f'inc-of-{exonstring}'
                        esizes = ','.join([str(x[1] - x[0]) for x in my_exons])
                        estarts = ','.join([str(10 + x[0] - outerjunc[0]) for x in my_exons])
                        # FIXME: use BED class
                        bedline = [thischrom, outerjunc[0] - 10, outerjunc[1] + 10, f'{jname}_{ename}', 0, strand, outerjunc[0] - 10, outerjunc[1] + 10,
                                   mycolor, 3, f'10,{esizes},10', f'0,{estarts},{10 + outerjunc[1] - outerjunc[0]}']
                        event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], innerjuncs_to_counts[innerjuncs],
                                                                                   set(innerjuncs), my_juncs - set(innerjuncs), {outerjunc, }, set(my_exons))

                        for j in innerjuncs:
                            if j not in esjuncs:
                                esjuncs[j] = {}
                            if outerjunc not in esjuncs[j]:
                                esjuncs[j][outerjunc] = {s: 0 for s in allsamples}
                            esjuncs[j][outerjunc] = add_counts_to_dict(esjuncs[j][outerjunc], innerjuncs_to_counts[innerjuncs])

                    if any([othercounts[s] >= min_read_support for s in allsamples]):
                        event_to_info[ename].other = othercounts
                    event_to_info[ename].totoverlap = get_overlapping_reads(outerjunc, interval_to_reads, allsamples)

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
    if len(termExon) > 1:
        tot_counts = {s: 0 for s in allsamples}
        othercounts = {s: 0 for s in allsamples}
        ss_to_counts = {}
        ss_to_junc_qual = {}
        actual_terminal = []
        for termSS in termExon:
            ss_to_counts[termSS] = {s: 0 for s in allsamples}
            ss_to_junc_qual[termSS] = []
            is_annot_terminal = False
            if gene in annot_terminal_ss:
                for ss in annot_terminal_ss[gene]:
                    if ss - 20 <= termSS <= ss + 20:
                        is_annot_terminal = True
                        break
            for termJunc in termExon[termSS]:
                ss_to_counts[termSS] = add_counts_to_dict(ss_to_counts[termSS], termExon[termSS][termJunc])
                # here we determine whether this splice site is really a terminal exon
                # either all reads terminate here, or it's annotated as a terminal junction
                if is_annot_terminal:
                    ss_to_junc_qual[termSS].append(True)
                else:
                    ss_to_junc_qual[termSS].append(any([termExon[termSS][termJunc][s] == alljuncs[termJunc][s] and termExon[termSS][termJunc][s] >= min_read_support for s in allsamples]))

        for termSS in termExon:
            if all(ss_to_junc_qual[termSS]):
                actual_terminal.append(termSS)
            else:
                othercounts = add_counts_to_dict(othercounts, ss_to_counts[termSS])
            tot_counts = add_counts_to_dict(tot_counts, ss_to_counts[termSS])

        if any([tot_counts[s] >= event_support for s in allsamples]):
            goodB = []
            goodBstrict = []
            my_juncs = set()
            my_juncs = set()
            for termSS in actual_terminal:
                inc_frac_list = [ss_to_counts[termSS][s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]

                if any([ss_to_counts[termSS][s] >= min_read_support for s in allsamples]) \
                    and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                         or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                    goodB.append(termSS)
                    if any([ss_to_counts[termSS][s] / tot_counts[s] >= MIN_TERMINAL_SS_FRACTION if tot_counts[s] > 0 else False for s in allsamples]):
                        goodBstrict.append(termSS)
                    for termJunc in termExon[termSS]:
                        if any([termExon[termSS][termJunc][s] >= min_read_support for s in allsamples]) \
                                and any([junc_frac_of_event <= termExon[termSS][termJunc][s] / tot_counts[s] <= 1 - junc_frac_of_event if tot_counts[s] > 0 else False for s in allsamples]):
                            my_juncs.add(termJunc)
                else:
                    othercounts = add_counts_to_dict(othercounts, ss_to_counts[termSS])

            # check that the junctions involved are distant enough from each other
            if len(goodB) > 1 and len(goodBstrict) > 1 and max(goodBstrict) - min(goodBstrict) >= MIN_TERMINAL_JUNCTION_SEPARATION:
                ename = f'{eventtype}-({strand})-{gene}'
                event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)

                for ssB in goodB:
                    jname = f'{thischrom}:{ssB}'
                    if (eventtype == 'AFE' and strand == '+') or (eventtype == 'ALE' and strand == '-'):
                        bs, be = ssB - 100, ssB
                    else:
                        bs, be = ssB, ssB + 100
                    # FIXME: use BED class
                    bedline = [thischrom, bs, be, f'{jname}_{ename}', 0, strand, bs, be, mycolor, 1, be - bs, 0]
                    inc_juncs = {x for x in my_juncs if ssB in x}
                    event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], ss_to_counts[ssB], inc_juncs, my_juncs - inc_juncs)

                if any([othercounts[s] >= min_read_support for s in allsamples]):
                    event_to_info[ename].other = othercounts
                junc_combos.add(frozenset(my_juncs))
                ref_junc = (max([x[0] for x in my_juncs]), min([x[1] for x in my_juncs]))
                event_to_info[ename].totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, allsamples)

            # TODO: add initial check for all reads with junction are equal to all terminal exon reads!
            # TODO: identify good juncs for each terminal SS, save junc combo
            # Good juncs must be significant proportion of all reads for event!
    return event_to_info, junc_combos


def process_junction_events(ssAtoB, esjuncs, eventtype, thischrom, strand, gene, allsamples, colordict, min_read_support, interval_to_reads, junc_frac_of_event, event_support, afe_ale_junc_comb):  # noqa: C901 - FIXME: reduce complexity
    event_to_info = {}
    ssB_groups_to_ssA = {}
    for ssA in ssAtoB:
        if len(ssAtoB[ssA]) > 1:  # don't report when no alternative splicing at all
            goodB, goodBstrict = [], []
            junccomb = set()
            my_juncs = set()
            othercounts = {s: 0 for s in allsamples}
            tot_counts = {s: 0 for s in allsamples}
            for ssB in ssAtoB[ssA]:
                tot_counts = add_counts_to_dict(tot_counts, ssAtoB[ssA][ssB])

            for ssB in ssAtoB[ssA]:

                junc = (min((ssA, ssB)), max((ssA, ssB)))
                inc_frac_list = [ssAtoB[ssA][ssB][s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]

                if len(inc_frac_list) >= 1 \
                    and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                         or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):

                    goodB.append(ssB)
                    my_juncs.add(junc)
                    junccomb.add(junc)
                    if junc not in esjuncs:
                        goodBstrict.append(ssB)
                    else:
                        escounts = {s: 0 for s in allsamples}
                        for ssB2 in ssAtoB[ssA]:
                            if ssB2 != ssB:
                                junc2 = (min((ssA, ssB2)), max((ssA, ssB2)))
                                if junc2 in esjuncs[junc]:
                                    escounts = add_counts_to_dict(escounts, esjuncs[junc][junc2])
                        deltafrac = [(ssAtoB[ssA][ssB][s] - escounts[s]) / ssAtoB[ssA][ssB][s] if ssAtoB[ssA][ssB][s] > 0 else 0 for s in allsamples]
                        delta = [ssAtoB[ssA][ssB][s] - escounts[s] for s in allsamples]
                        if any([deltafrac[x] > 0.1 and delta[x] > min_read_support for x in range(len(delta))]):
                            goodBstrict.append(ssB)
                else:
                    othercounts = add_counts_to_dict(othercounts, ssAtoB[ssA][ssB])

            if frozenset(junccomb) not in afe_ale_junc_comb:
                if len(goodB) > 1:
                    ssBgroup = frozenset(goodB)
                    if ssBgroup not in ssB_groups_to_ssA:
                        ssB_groups_to_ssA[ssBgroup] = {}
                    ssB_groups_to_ssA[ssBgroup][ssA] = tot_counts

                if any([tot_counts[s] >= event_support for s in allsamples]) \
                        and len(goodBstrict) > 1:  # only if alt splicing, not including junctions involved in exon skipping
                    ename = f'{eventtype}-relTo-{thischrom}:{ssA}({strand})-{gene}'
                    event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)
                    for ssB in goodB:
                        junc = (min((ssA, ssB)), max((ssA, ssB)))
                        jname = f'{thischrom}:{junc[0]}-{junc[1]}'
                        bs, be = min((ssA, ssB)), max((ssA, ssB))
                        # FIXME: use BED class
                        bedline = [thischrom, bs, be, f'{jname}_{ename}', 0, strand, bs, be, colordict[eventtype], 1, be - bs, 0]
                        event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], ssAtoB[ssA][ssB], {junc, }, my_juncs - {junc, })

                    if any([othercounts[s] >= min_read_support for s in allsamples]):
                        event_to_info[ename].other = othercounts
                    # against min(goodB), not the ssB the loop above left behind.  All
                    # goodB lie on the same side of ssA, so this is the same junction
                    # the leaked value gave, said properly
                    ref_junc = (ssA, min(goodB)) if ssA < min(goodB) else (max(goodB), ssA)
                    event_to_info[ename].totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, allsamples)
    eventtype = 'alt3' if eventtype == 'alt5' else 'alt5'
    for ssBgroup in ssB_groups_to_ssA:
        if len(ssB_groups_to_ssA[ssBgroup]) > 1:
            tot_counts = {s: 0 for s in allsamples}
            othercounts = {s: 0 for s in allsamples}
            for ssA in ssB_groups_to_ssA[ssBgroup]:
                tot_counts = add_counts_to_dict(tot_counts, ssB_groups_to_ssA[ssBgroup][ssA])
            if any([tot_counts[s] >= event_support for s in allsamples]):
                good = []
                my_juncs = set()
                for ssA in ssB_groups_to_ssA[ssBgroup]:
                    inc_frac_list = [ssB_groups_to_ssA[ssBgroup][ssA][s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]
                    if (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                            or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                        good.append(ssA)
                        for ssB in ssBgroup:
                            my_juncs.add((min((ssA, ssB)), max((ssA, ssB))))
                    else:
                        othercounts = add_counts_to_dict(othercounts, ssB_groups_to_ssA[ssBgroup][ssA])
                if len(good) > 1:
                    oj_string = [f'{thischrom}:{x}' for x in sorted(list(ssBgroup))]
                    ename = f'{eventtype}ss-relTo-{",".join(oj_string)}({strand})-{gene}'
                    event_to_info[ename] = SplicingEvent(ename, f'{eventtype}ss', gene, thischrom, strand, tot_counts, allsamples)
                    for ssA in good:
                        if (eventtype == 'alt5' and strand == '+') or (eventtype == 'alt3' and strand == '-'):
                            bs, be = ssA - 10, ssA
                        else:
                            bs, be = ssA, ssA + 10
                        jname = f'{thischrom}:{ssA}'
                        inc_juncs = set()
                        for ssB in ssBgroup:
                            inc_juncs.add((min((ssA, ssB)), max((ssA, ssB))))

                        # FIXME: use BED class
                        bedline = [thischrom, bs, be, f'{jname}_{ename}', 0, strand, bs, be, colordict[eventtype], 1, be - bs, 0]
                        event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], ssB_groups_to_ssA[ssBgroup][ssA], inc_juncs, my_juncs - inc_juncs)
                    if any([othercounts[s] >= min_read_support for s in allsamples]):
                        event_to_info[ename].other = othercounts
                    ref_junc = (max(good), min(ssBgroup)) if median(good) < median(ssBgroup) else (max(ssBgroup), min(good))
                    event_to_info[ename].totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, allsamples)

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

class SplicingEventJunction:
    """represents metadata about a specific junction in a splicing event"""
    def __init__(self, name, bedlines, samplecounts, inc_junc, exc_junc, outer_junc=None, inc_exon=None):
        # FIXME: use BED class
        self.name = name
        self.bedlines = bedlines
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
    for j in alljuncs:
        retainedcounts = {s: 0 for s in allsamples}
        for b in allblocks:  # this nested for loop could get brutal
            if b[0] < j[0] and j[1] < b[1]:
                for s in allsamples:
                    retainedcounts[s] += allblocks[b][s]
        splicedcounts = alljuncs[j]
        # if sum(allretainedcounts.values()) > 0:  # only report if junction retained at all
        tot_counts = splicedcounts
        tot_counts = add_counts_to_dict(tot_counts, retainedcounts)

        inc_frac_list = [retainedcounts[s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]
        if len(inc_frac_list) >= 1 \
            and any([splicedcounts[s] >= min_read_support for s in allsamples]) \
            and any([retainedcounts[s] >= min_read_support for s in allsamples]) \
            and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                 or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
            ename = f'ir-of-{thischrom}:{j[0]}-{j[1]}({strand})-{gene}'
            event_to_info[ename] = SplicingEvent(ename, 'ir', gene, thischrom, strand, tot_counts, allsamples)

            # FIXME: use BED class
            bedline = [thischrom, j[0] - 10, j[1] + 10, f'spliced_{ename}', 0, strand, j[0] - 10, j[1] + 10, mycolor, 2, '10,10', f'0,{10 + j[1] - j[0]}']
            event_to_info[ename].events['spliced'] = SplicingEventJunction('spliced', [bedline], splicedcounts, {j, }, {})

            # FIXME: use BED class
            bedline = [thischrom, j[0], j[1], f'retained_{ename}', 0, strand, j[0], j[1], mycolor, 1, j[1] - j[0], 0]
            event_to_info[ename].events['retained'] = SplicingEventJunction('retained', [bedline], retainedcounts, {}, {j, })

            event_to_info[ename].totoverlap = get_overlapping_reads(j, interval_to_reads, allsamples)
    return event_to_info

def write_ends(grouped_ends, allsamples, thischrom, strand, gene, eventtype, mycolor, support, junc_frac_of_event, event_support):
    good_ends = set()
    event_to_info = {}
    othercounts = {s: 0 for s in allsamples}
    tot_counts = {s: 0 for s in allsamples}
    for e in grouped_ends:
        tot_counts = add_counts_to_dict(tot_counts, grouped_ends[e])
    if any([tot_counts[s] >= event_support for s in allsamples]):
        for e in grouped_ends:
            inc_frac_list = [grouped_ends[e][s] / tot_counts[s] for s in allsamples if tot_counts[s] >= event_support]

            if any([grouped_ends[e][s] >= support for s in allsamples]) \
                and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                     or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                good_ends.add(e)
            else:
                othercounts = add_counts_to_dict(othercounts, grouped_ends[e])
        if len(good_ends) > 1:
            ename = f'{eventtype}-({strand})-{gene}'
            event_to_info[ename] = SplicingEvent(ename, eventtype, gene, thischrom, strand, tot_counts, allsamples)

            for e in good_ends:
                jname = f'{thischrom}:{e}'
                # FIXME: use BED class
                bedline = [thischrom, e - 1, e, f'{jname}_{ename}', 0, strand, e - 1, e, mycolor, 1, 1, 0]
                event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], grouped_ends[e], {e, }, good_ends - {e, })
            if any([othercounts[s] >= support for s in allsamples]):
                event_to_info[ename].other = othercounts
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
                inc_frac_list = [jinfo.samplecounts[s] / event.totoverlap[s] for s in allsamples if event.totoverlap[s] > 0]

                if (event.eventtype not in {'es', 'ces'} or jname != 'exc') \
                    and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                         or (len(allsamples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
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
                        # FIXME: use BED class
                        for line in jinfo.bedlines:
                            line[4] = round(med * 100)
                            outbed.write('\t'.join([str(x) for x in line]) + '\n')
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

        for gene in allgenes:

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


def process_bed_line(bed_rec):
    transcript = bed_rec.name
    gene = bed_rec.name.split('_')[-1]
    exons = [(blk.start, blk.end) for blk in bed_rec.blocks]
    junctions = [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)]
    strand = bed_rec.strand
    return gene, transcript, exons, junctions, strand


def _run_region(*, partition, gtf_data, junction_corrector, args, allsamples, region_reports):  # noqa: C901 - FIXME: reduce complexity
    # FIXMEL what are these files being created, just do in memory
    region_bed = partition.output_path('region.bed')
    out = open(region_bed, 'w')
    out.write('\t'.join([partition.region.name, str(partition.region.start), str(partition.region.end)]) + '\n')
    out.close()

    annots = annot_data_from_gtf(gtf_data, partition.region)
    if len(annots.transcripts) > 0:  # check if any annotated transcripts in region
        if args.annot_basic:
            region_annot_basic = partition.file_prefix + '.annotation.basic.bed'
            pipettor.run([('bedtools', 'intersect', '-wa', '-a', args.annot_basic, '-b', region_bed)], stdout=region_annot_basic)

        annot_afe_ss, annot_ale_ss = {}, {}
        if args.annot_basic:
            for bed_rec in BedReader(region_annot_basic, fixScores=True):
                gene, transcript, exons, junctions, strand = process_bed_line(bed_rec)
                if gene not in annot_afe_ss:
                    annot_afe_ss[gene] = set()
                    annot_ale_ss[gene] = set()
                if len(junctions) > 0:
                    afe_ss, ale_ss = junctions[0][0], junctions[-1][1]
                    if strand == '-':
                        afe_ss, ale_ss = ale_ss, afe_ss
                    annot_afe_ss[gene].add(afe_ss)
                    annot_ale_ss[gene].add(ale_ss)

        # align reads to annot [transcripts +- 1000bp], filter to only good aligns, convert to genomic coords
        # get reads from that bam file (no need to correct), save junctions - am actually doing correct, could probably remove that
        # add junctions from other reads that did not match ref transcriptome well after correction
        # load splice junctions for chrom

        with pysam.FastaFile(args.genome) as genome:
            for sample, bamfile in allsamples:
                if not os.path.exists(partition.file_prefix + '_' + sample + '_gene_to_juncs.txt'):
                    get_juncs_single_sample(args, partition.region, partition.file_prefix, sample, bamfile, genome, junction_corrector, annots)

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
        report_regions = []
        for line in open(args.region_bed):
            # FIXME: use BED class
            line = line.rstrip().split('\t')
            report_regions.append(SeqRange(line[0], int(line[1]), int(line[2])))
        all_regions, region_reports = correction_regions(combine_regions(report_regions),
                                                         [bamfile for _, bamfile in allsamples])
    else:
        all_regions = [SeqRange(chrom, 0, genome.get_reference_length(chrom))
                       for chrom in genome.references]
        region_reports = {_region_key(region): [region] for region in all_regions}
    logging.info(f'input regions: {len(all_regions)}')

    logging.info('loading annotation GTF')
    annot_gtf_data = gtf_data_parser(args.annot, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)
    logging.info('building intron support database')
    junction_corrector = junction_corrector_factory(args.ss_window, args.junction_support,
                                                    annot_gtf_data=annot_gtf_data,
                                                    intron_beds=args.junction_bed,
                                                    star_sj_tabs=args.junction_tab)

    # FIXME: we appear to do this just to do bed intersect; just do in memory.
    if args.annot_basic:
        annot_basic_bed = tempDir + '/annotation.basic.bed'
        if not os.path.exists(annot_basic_bed):
            pipettor.run([(sys.executable, gtf_to_bed_prog, args.annot_basic, annot_basic_bed, '--include_gene')])
        args.annot_basic = annot_basic_bed

    logging.info(f'running regions with {args.threads} threads')

    runner = PartitionRunner(all_regions, tempDir, gtf_data=annot_gtf_data, junction_corrector=junction_corrector,
                             threads=args.threads)
    runner.run(_run_region, args=args, allsamples=allsamples, region_reports=region_reports)

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
