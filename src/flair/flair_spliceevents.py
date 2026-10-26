#! /usr/bin/env python3

import argparse
import os
import os.path as osp
import sys
import pipettor
import shutil
import pysam
import logging
import scipy.stats as sps
from flair.partition_runner import PartitionRunner, combine_temp_files_by_suffix
from flair.gene_juncs_tsv import GeneJuncsReader, GeneJuncsWriter
from flair.read_ends_tsv import ReadEndsReader
from flair import SeqRange, resolve_deprecated_option
from flair import thread_share
from statistics import median
from flair.junction_correct import junction_corrector_factory
from flair.isoform_data import ReadRec
from flair.io_utils import check_genome_fasta, make_temp_dir, sort_bed_file
from flair.read_processing import get_sequence_from_bed, generate_genomic_alignment_read_to_clipping_file
from flair.read_correction import filter_correct_group_reads
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import annot_data_from_gtf
from flair.gtf_to_bed import GENE_ID_EXTRA_COL
from flair.splicing_events_tsv import SplicingEventsWriter
from flair.splicing_outliers_tsv import SplicingOutliersWriter
from flair.pycbio.hgdata.bed import Bed, BedBlock, BedReader
from flair.count_sam_transcripts import run_count_sam_transcripts
from flair.pycbio.sys import fileOps

pkgdir = osp.dirname(osp.realpath(__file__))
gtf_to_bed_prog = osp.join(pkgdir, "gtf_to_bed.py")

# joined files whose parts each carry a header
EVENT_SUFFIXES = ('.diffsplice.counts.tsv', '.diffsplice.PSIjunc.tsv',
                  '.diffsplice.PSItot.tsv')
OUTLIER_SUFFIXES = ('.diffsplice.outliers.tsv', '.diffsplice.outliers.filtered.tsv')

# Bases of the flanking exons included either side of an event in the BED that
# spliceevents writes for viewing, <output>.diffsplice.bed.
#
# An event is defined by its junctions, so its own extent is the intron: drawn
# as it stands, the two sides of an event would be a bare gap with nothing to
# anchor them, and the inclusion and exclusion rows would be hard to tell apart.
# Each row therefore starts FLANK_STUB bases before the first junction and ends
# FLANK_STUB after the last, and carries those two stubs as its outer blocks, so
# inclusion shows as stub, exon, stub and exclusion as stub, stub with the
# junction as the gap between them.
#
# It is display padding and nothing else: the value does not reach any count,
# PSI or test, and the coordinates the event is named and reported by are the
# junctions themselves, not these bounds.
FLANK_STUB = 10

def get_args():
    parser = argparse.ArgumentParser(description='identifies counts of different splicing events directly from a '
                                                 'bam file of long-rna-seq alignments')
    parser.add_argument('-m', '--manifest', required=True,
                        help='manifest tsv with sample name and path to bam file aligned to genome')
    parser.add_argument('-g', '--genome', required=True, type=str,
                        help='FASTA of reference genome, can be minimap2 indexed')
    parser.add_argument('-o', '--output', default='flair',
                        help='output file name base for FLAIR isoforms (default: flair)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='the number of threads to use (4)')
    # FIXME: this is different than all other modules that use --gtf.
    parser.add_argument('-f', '--annot', required=True, default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms '
                             'to annotated isoforms and adjusting TSS/TESs')
    # FIXME: what is the difference than the above, same description
    parser.add_argument('--annot_basic', default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms '
                             'to annotated isoforms and adjusting TSS/TESs')
    parser.add_argument('--junction_star', help='splice junctions in STAR SJ.out.tab format, as STAR writes '
                                                'when it aligns reads')
    parser.add_argument('--junction_tab', help='deprecated name for --junction_star')
    parser.add_argument('--junction_bed', help='splice junctions as BED6 to BED9 with the number of supporting '
                                               'reads in the score column, as intron-prospector writes from a BAM '
                                               'of either short or long reads.  An alternative format to '
                                               '--junction_star, not a different kind of evidence')
    parser.add_argument('--region_bed',
                        help='bed file with regions to parallelize by; if not specified, all chromosomes are used')
    parser.add_argument('--junction_support', type=int, default=1,
                        help='minimum number of supporting reads to keep a junction, from either --junction_star '
                             'or --junction_bed; for bed the score column holds the count (default: %(default)s)')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')
    parser.add_argument('--junc_support', type=int, default=2,
                        help='''minimum number of supporting reads for a specific junction''')
    parser.add_argument('--event_support', type=int, default=20,
                        help='''minimum number of supporting reads for an event for PSI to be calculated''')
    parser.add_argument('--event_frac_of_tot', type=float, default=0.1,
                        help='''minimum (total reads in splicing event)/(total coverage over locus) to call event''')
    parser.add_argument('--junc_frac_of_event', type=float, default=0.02,
                        help='''minimum (reads with junction/total reads in splicing event) Set to 0 for max recall''')
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.
                        Intermediate files include: promoter-supported reads file,
                        read assignments to firstpass isoforms''')
    parser.add_argument('--check_outliers', default=False, action='store_true',
                        help='''whether to run a statistical analysis to identify and filter outliers from the input samples''')
    parser.add_argument('--keep_sup', default=False, action='store_true',
                        help='''specify if you want to keep supplementary alignments to define isoforms''')
    parser.add_argument('--output_read_ends', default=False, action='store_true',
                        help='''specify if you want to output corrected read ends bed file''')
    parser.add_argument('--noaligntoannot', default=False, action='store_true',
                        help='''specify if you don't want
                        an initial alignment to the annotated sequences and only want transcript
                        detection from the genomic alignment.
                        Will be slightly faster but less accurate if the annotation is good''')

    args = parser.parse_args()
    resolve_deprecated_option(args, deprecated="junction_tab", current="junction_star")

    args.trust_ends = False
    args.remove_internal_priming = False
    args.quality = 0  # should only be used for genomic alignment

    if not os.path.exists(args.genome):
        parser.error(f'Genome file path does not exist: {args.genome}')
    return args


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

def get_juncs_to_gene(juncs, isoinfo, sjc_to_gene, junc_to_gene, gene_to_exons, gene_to_annot_juncs):
    gene = None
    if juncs in sjc_to_gene:
        gene = sjc_to_gene[juncs]
    else:
        gene_hits = {}
        for j in juncs:
            if j in junc_to_gene:
                gene_id = junc_to_gene[j]
                if gene_id not in gene_hits:
                    gene_hits[gene_id] = [0, -1 * len(gene_to_annot_juncs[gene_id])]
                gene_hits[gene_id][0] += 1
        if gene_hits:
            sortedgenes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)
            gene = sortedgenes[0][0]
        else:
            # look for exon overlap.  isoinfo is an Isoform, which is not iterable;
            # its reads are what was wanted.  The median is how isoform_data turns a
            # read cluster into one isoform, and neither extreme is outlier proof
            mystart = int(median(isoinfo.starts))
            myend = int(median(isoinfo.ends))
            exons = ([(mystart, juncs[0].start)]
                     + [(juncs[i].end, juncs[i + 1].start) for i in range(len(juncs) - 1)]
                     + [(juncs[-1].end, myend)])
            # strand = isoinfo[0][2]  # not a super robust strand picking, assumes well stranded reads
            gene_hits = get_annot_gene_hits(gene_to_exons, exons)
            if len(gene_hits) > 0:
                gene_hits.sort(reverse=True)
                gene = gene_hits[0][1]
    return gene

def group_juncs_by_annot_gene(sjtoends, sjc_to_gene, junc_to_gene, gene_to_exons, gene_to_annot_juncs):
    gene_to_juncs, nogenejuncs, sereads = {}, {}, []
    no_gene_reads, se_reads_cnt = 0, 0
    for (chrom, juncs) in sjtoends:
        if len(juncs) > 0:  # skip single-exon reads
            gene = get_juncs_to_gene(juncs, sjtoends[(chrom, juncs)], sjc_to_gene, junc_to_gene, gene_to_exons, gene_to_annot_juncs)
            if gene:
                if gene not in gene_to_juncs:
                    gene_to_juncs[gene] = {}
                gene_to_juncs[gene][juncs] = sjtoends[(chrom, juncs)]
            else:
                no_gene_reads += sjtoends[(chrom, juncs)].num_reads
                nogenejuncs[juncs] = sjtoends[(chrom, juncs)]
        else:
            se_reads_cnt += sjtoends[(chrom, juncs)].num_reads
            sereads.extend(sjtoends[(chrom, juncs)].reads)
    logging.debug(f"gene assignment: {no_gene_reads} spliced reads not assigned to a gene, {se_reads_cnt} single-exon reads discarded")
    return gene_to_juncs, nogenejuncs, sereads


def extract_35ss_info(juncs, strand, all_samples, sample, numreads, ss5to3, ss3to5, alljuncs):
    for j in juncs:
        if strand == '+':
            ss5, ss3 = j.start, j.end
        else:
            ss5, ss3 = j.end, j.start
        if ss5 not in ss5to3:
            ss5to3[ss5] = {}
        if ss3 not in ss3to5:
            ss3to5[ss3] = {}
        if ss3 not in ss5to3[ss5]:
            ss5to3[ss5][ss3] = {s: 0 for s in all_samples}
        if ss5 not in ss3to5[ss3]:
            ss3to5[ss3][ss5] = {s: 0 for s in all_samples}
        ss5to3[ss5][ss3][sample] += numreads
        ss3to5[ss3][ss5][sample] += numreads
        if j not in alljuncs:
            alljuncs[j] = {s: 0 for s in all_samples}
        alljuncs[j][sample] += numreads
    return ss5to3, ss3to5, alljuncs


def extract_exon_usage_info(juncs, all_samples, sample, numreads, allblocks, exonjpairs):
    # just for middle exons
    for i in range(len(juncs) - 1):
        exon = (juncs[i][1], juncs[i + 1][0])
        prevjunc, nextjunc = juncs[i], juncs[i + 1]
        if exon not in allblocks:
            allblocks[exon] = {s: 0 for s in all_samples}
        allblocks[exon][sample] += numreads
        if (prevjunc, nextjunc) not in exonjpairs:
            exonjpairs[(prevjunc, nextjunc)] = {s: 0 for s in all_samples}
        exonjpairs[(prevjunc, nextjunc)][sample] += numreads
    return allblocks, exonjpairs


def extract_end_coverage_info(juncs, readinfo, all_samples, sample, chrom, gene, strand, allblocks, outends):
    for r in readinfo:
        # getting all block coverage to measure intron retention
        firstexon = (r.start, juncs[0].start)
        lastexon = (juncs[-1].end, r.end)
        if firstexon not in allblocks:
            allblocks[firstexon] = {s: 0 for s in all_samples}
        if lastexon not in allblocks:
            allblocks[lastexon] = {s: 0 for s in all_samples}
        allblocks[firstexon][sample] += 1
        allblocks[lastexon][sample] += 1
        if outends:
            Bed(chrom, r.start, r.end, name=gene + '|' + r.name,
                score=0, strand=strand).write(outends)
    return allblocks

def determine_juncs_are_subset(juncs, alljuncs):
    is_subset = False
    unique_seq_bound = [[], []]
    for otherjuncs in alljuncs:
        if juncs != otherjuncs and len(juncs) < len(otherjuncs):
            if str(juncs)[1:-1].rstrip(',') in str(otherjuncs):
                is_subset = True
                other_internal_exons = [(otherjuncs[x].end, otherjuncs[x + 1].start) for x in range(len(otherjuncs) - 1)]
                for other_exon in other_internal_exons:
                    if juncs[0].start == other_exon[1]:
                        unique_seq_bound[0].append(other_exon[1] - other_exon[0])
                    if juncs[-1].end == other_exon[0]:
                        unique_seq_bound[1].append(other_exon[1] - other_exon[0])
    if is_subset:
        unique_seq_bound[0] = max(unique_seq_bound[0]) if len(unique_seq_bound[0]) > 0 else None
        unique_seq_bound[1] = max(unique_seq_bound[1]) if len(unique_seq_bound[1]) > 0 else None
        return unique_seq_bound
    else:
        return [None, None]

def get_start_end_counts(readinfo, subset_info, sample, first_junc, last_junc, t_starts_ends, t_first_last_sj, all_samples):
    for r in readinfo:
        # if subset_info[0] == None or first_junc-start > subset_info[0]:
        t_starts_ends[0].append((r.start, sample))
        if first_junc not in t_first_last_sj[0]:
            t_first_last_sj[0][first_junc] = {s: 0 for s in all_samples}
        t_first_last_sj[0][first_junc][sample] += 1
        # if subset_info[1] == None or end-last_junc > subset_info[1]:
        t_starts_ends[1].append((r.end, sample))
        if last_junc not in t_first_last_sj[1]:
            t_first_last_sj[1][last_junc] = {s: 0 for s in all_samples}
        t_first_last_sj[1][last_junc][sample] += 1
    return t_starts_ends, t_first_last_sj


def group_ends(t_starts_ends, all_samples, window):
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
                    end_to_counts[g] = {s: 0 for s in all_samples}
                    for e in curr_group:
                        end_to_counts[g][e[1]] += 1
                curr_group = []
            curr_group.append(end)
            last_end = end[0]
        # all_groups.append(curr_group)
        # end_to_counts[median(curr_group)] = len(curr_group)
        g = int(median([x[0] for x in curr_group]))
        end_to_counts[g] = {s: 0 for s in all_samples}
        for e in curr_group:
            end_to_counts[g][e[1]] += 1

        grouped_ends[i] = end_to_counts
    return grouped_ends


def extract_afe_ale(juncs, strand, all_samples, sample, numreads, alljuncs, afe, ale):
    if strand == '+':
        afe_junc, ale_junc = juncs[0], juncs[-1]
        afe_ss, ale_ss = juncs[0][0], juncs[-1][1]
    else:
        ale_junc, afe_junc = juncs[0], juncs[-1]
        ale_ss, afe_ss = juncs[0][0], juncs[-1][1]

    if afe_ss not in afe:
        afe[afe_ss] = {}
    if afe_junc not in afe[afe_ss]:
        afe[afe_ss][afe_junc] = {s: 0 for s in all_samples}
    afe[afe_ss][afe_junc][sample] += numreads

    if ale_ss not in ale:
        ale[ale_ss] = {}
    if ale_junc not in ale[ale_ss]:
        ale[ale_ss][ale_junc] = {s: 0 for s in all_samples}
    ale[ale_ss][ale_junc][sample] += numreads

    return afe, ale


def extract_tandem_splicing_info(sjc, alljuncs, outer_junc_to_exons, all_samples, sample, numreads):
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
                    outer_junc_to_exons[j][inner_juncs] = {s: 0 for s in all_samples}
                outer_junc_to_exons[j][inner_juncs][sample] += numreads
    return outer_junc_to_exons


def extract_splicing_info(all_samples, all_genes_to_juncs, gene, strand, chrom, outends):
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
    for s in range(len(all_samples)):
        sample = all_samples[s]
        # TODO: change to gene_to_juncs being all juncs from all samples!!
        gene_to_juncs = all_genes_to_juncs[s]
        if gene in gene_to_juncs:
            for juncs in gene_to_juncs[gene]:
                readinfo = gene_to_juncs[gene][juncs]
                numreads = len(readinfo)

                # get intervals for calculating coverage
                for r in readinfo:
                    interval = (round(r.start, -1), round(r.end, -1))
                    if interval not in interval_to_reads:
                        interval_to_reads[interval] = {s: 0 for s in all_samples}
                    interval_to_reads[interval][sample] += 1

                ss5to3, ss3to5, alljuncs = extract_35ss_info(juncs, strand, all_samples, sample, numreads,
                                                             ss5to3, ss3to5, alljuncs)
                afe, ale = extract_afe_ale(juncs, strand, all_samples, sample, numreads, alljuncs, afe, ale)
                outer_junc_to_exons = extract_tandem_splicing_info(juncs, alljuncs, outer_junc_to_exons, all_samples, sample, numreads)
                allblocks, exonjpairs = extract_exon_usage_info(juncs, all_samples, sample, numreads,
                                                                allblocks, exonjpairs)
                allblocks = extract_end_coverage_info(juncs, readinfo, all_samples, sample, chrom, gene, strand,
                                                      allblocks, outends)
                subset_info = determine_juncs_are_subset(juncs, gene_to_juncs[gene])
                t_starts_ends, t_first_last_sj = get_start_end_counts(readinfo, subset_info, sample, juncs[0][0], juncs[-1][1],
                                                                      t_starts_ends, t_first_last_sj, all_samples)

    return ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj, interval_to_reads, afe, ale, outer_junc_to_exons

def write_exon_skipping(exonjpairs, alljuncs, all_samples, chrom, strand, gene, mycolor, min_read_support, interval_to_reads, junc_frac_of_event, event_support, outer_junc_to_exons):  # noqa: C901 - FIXME: reduce complexity
    # outer_junc_to_exons = {}
    # outer_junc_to_all_inc = {}
    exon_to_outer_juncs = {}
    exon_to_all_inc = {}
    exon_to_all_exc = {}
    event_to_info = {}
    for prevjunc, nextjunc in exonjpairs:
        outer_junc = (prevjunc[0], nextjunc[1])
        exon = (prevjunc[1], nextjunc[0])

        if exon not in exon_to_outer_juncs:
            exon_to_outer_juncs[exon] = set()
            exon_to_all_inc[exon] = {s: 0 for s in all_samples}
            exon_to_all_exc[exon] = {s: 0 for s in all_samples}
        exon_to_outer_juncs[exon].add(outer_junc)
        exon_to_all_inc[exon] = add_counts_to_dict(exon_to_all_inc[exon], exonjpairs[(prevjunc, nextjunc)])
        if outer_junc in alljuncs:
            # if outer_junc not in outer_junc_to_exons:
            #     outer_junc_to_exons[outer_junc] = set()
            #     outer_junc_to_all_inc[outer_junc] = {s: 0 for s in all_samples}
            exon_to_all_exc[exon] = add_counts_to_dict(exon_to_all_exc[exon], alljuncs[outer_junc])
            # outer_junc_to_exons[outer_junc].add(exon)
            # outer_junc_to_all_inc[outer_junc] = add_counts_to_dict(outer_junc_to_all_inc[outer_junc], exonjpairs[(prevjunc, nextjunc)])

    seen_junc_combos = set()
    esjuncs = {}
    for exon in exon_to_outer_juncs:
        tot_counts = exon_to_all_inc[exon]
        tot_counts = add_counts_to_dict(tot_counts, exon_to_all_exc[exon])

        inc_frac_list = [exon_to_all_inc[exon][s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]
        if len(inc_frac_list) >= 1 \
                and any([exon_to_all_inc[exon][s] >= min_read_support for s in all_samples]) \
                and any([exon_to_all_exc[exon][s] >= min_read_support for s in all_samples]) \
                and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                     or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
            outer_juncs = exon_to_outer_juncs[exon]

            bedlines = []
            good_outer_juncs = set()
            inner_juncs = set()
            for outer_junc in outer_juncs:
                ejpair = exonjpairs[((outer_junc[0], exon[0]), (exon[1], outer_junc[1]))]
                if any([ejpair[s] >= min_read_support for s in all_samples]) \
                        and any([junc_frac_of_event <= ejpair[s] / (exon_to_all_inc[exon][s] + exon_to_all_exc[exon][s]) <= 1 - junc_frac_of_event
                                 if (exon_to_all_inc[exon][s] + exon_to_all_exc[exon][s]) > 0 else False for s in all_samples]):
                    good_outer_juncs.add(outer_junc)
                    inner_juncs.add((outer_junc[0], exon[0]))
                    inner_juncs.add((exon[1], outer_junc[1]))

            ename = f'es-of-{chrom}:{exon[0]}-{exon[1]}'

            for outer_junc in good_outer_juncs:
                # the skipped exon between the two flanking stubs
                bed_start, bed_end = outer_junc[0] - FLANK_STUB, outer_junc[1] + FLANK_STUB
                blocks = [BedBlock(bed_start, outer_junc[0]),
                          BedBlock(exon[0], exon[1]),
                          BedBlock(outer_junc[1], bed_end)]
                bedlines.append(Bed(chrom, bed_start, bed_end, name=f'inc_{ename}',
                                    score=0, strand=strand, thickStart=bed_start,
                                    thickEnd=bed_end, itemRgb=mycolor, blocks=blocks))
                if (outer_junc[0], exon[0]) not in esjuncs:
                    esjuncs[(outer_junc[0], exon[0])] = {}
                if outer_junc not in esjuncs[(outer_junc[0], exon[0])]:
                    esjuncs[(outer_junc[0], exon[0])][outer_junc] = {s: 0 for s in all_samples}
                esjuncs[(outer_junc[0], exon[0])][outer_junc] = add_counts_to_dict(esjuncs[(outer_junc[0], exon[0])][outer_junc], exonjpairs[((outer_junc[0], exon[0]), (exon[1], outer_junc[1]))])
                if (exon[1], outer_junc[1]) not in esjuncs:
                    esjuncs[(exon[1], outer_junc[1])] = {}
                if outer_junc not in esjuncs[(exon[1], outer_junc[1])]:
                    esjuncs[(exon[1], outer_junc[1])][outer_junc] = {s: 0 for s in all_samples}
                esjuncs[(exon[1], outer_junc[1])][outer_junc] = add_counts_to_dict(esjuncs[(exon[1], outer_junc[1])][outer_junc], exonjpairs[((outer_junc[0], exon[0]), (exon[1], outer_junc[1]))])

            event_to_info[ename] = SplicingEvent(ename, 'es', gene, chrom, strand, tot_counts, all_samples)

            event_to_info[ename].events['inc'] = SplicingEventJunction('inc', bedlines, exon_to_all_inc[exon], inner_juncs, set(), good_outer_juncs, {exon, })

            event_to_info[ename].events['exc'] = SplicingEventJunction('exc', [], exon_to_all_exc[exon], set(), inner_juncs, good_outer_juncs, {exon, })

            if len(good_outer_juncs) > 0:
                event_to_info[ename].totoverlap = get_overlapping_reads((max([x[0] for x in good_outer_juncs]), min([x[1] for x in good_outer_juncs])), interval_to_reads, all_samples)
            else:
                event_to_info[ename].totoverlap = get_overlapping_reads((max([x[0] for x in outer_juncs]), min([x[1] for x in outer_juncs])), interval_to_reads, all_samples)

            finaljuncs = good_outer_juncs | inner_juncs
            seen_junc_combos.add(frozenset(finaljuncs))

    for outer_junc in outer_junc_to_exons:
        inner_juncs_to_counts = outer_junc_to_exons[outer_junc]
        if any([alljuncs[outer_junc][s] >= min_read_support for s in all_samples]):
            good_b = set()
            all_exons = set()
            my_juncs = set()  # {outer_junc,}
            othercounts = {s: 0 for s in all_samples}
            tot_counts = alljuncs[outer_junc]
            for inner_juncs in inner_juncs_to_counts:
                tot_counts = add_counts_to_dict(tot_counts, inner_juncs_to_counts[inner_juncs])
            inc_frac_list = [alljuncs[outer_junc][s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]
            if len(inc_frac_list) > 0 \
                and any([alljuncs[outer_junc][s] >= min_read_support for s in all_samples]) \
                and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                     or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                for inner_juncs in inner_juncs_to_counts:
                    inc_frac_list = [inner_juncs_to_counts[inner_juncs][s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]
                    if any([inner_juncs_to_counts[inner_juncs][s] >= min_read_support for s in all_samples]) \
                        and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                             or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                        good_b.add(inner_juncs)
                        my_exons = [(inner_juncs[i][1], inner_juncs[i + 1][0]) for i in range(len(inner_juncs) - 1)]
                        all_exons.update(set(my_exons))
                        for j in inner_juncs:
                            my_juncs.add(j)
                    else:
                        othercounts = add_counts_to_dict(othercounts, inner_juncs_to_counts[inner_juncs])

                if len(good_b) > 0 \
                        and frozenset(my_juncs | {outer_junc, }) not in seen_junc_combos:  # only require one exon relative to outer junction
                    ename = f'ces-relTo-{chrom}:{outer_junc[0]}-{outer_junc[1]}({strand})-{gene}'
                    event_to_info[ename] = SplicingEvent(ename, 'ces', gene, chrom, strand, tot_counts, all_samples)

                    bed_start, bed_end = outer_junc[0] - FLANK_STUB, outer_junc[1] + FLANK_STUB
                    bedline = Bed(chrom, bed_start, bed_end, name=f'exc_{ename}', score=0,
                                  strand=strand, thickStart=bed_start, thickEnd=bed_end,
                                  itemRgb=mycolor,
                                  blocks=[BedBlock(bed_start, outer_junc[0]),
                                          BedBlock(outer_junc[1], bed_end)])
                    event_to_info[ename].events['exc'] = SplicingEventJunction('exc', [bedline], alljuncs[outer_junc], set(), my_juncs, {outer_junc, }, all_exons)

                    for inner_juncs in good_b:
                        my_exons = [(inner_juncs[i][1], inner_juncs[i + 1][0]) for i in range(len(inner_juncs) - 1)]
                        exonstring = ','.join([f'{chrom}:{x[0]}-{x[1]}' for x in my_exons])
                        jname = f'inc-of-{exonstring}'
                        bed_start, bed_end = outer_junc[0] - FLANK_STUB, outer_junc[1] + FLANK_STUB
                        blocks = ([BedBlock(bed_start, outer_junc[0])]
                                  + [BedBlock(x[0], x[1]) for x in my_exons]
                                  + [BedBlock(outer_junc[1], bed_end)])
                        bedline = Bed(chrom, bed_start, bed_end, name=f'{jname}_{ename}',
                                      score=0, strand=strand, thickStart=bed_start,
                                      thickEnd=bed_end, itemRgb=mycolor, blocks=blocks)
                        event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], inner_juncs_to_counts[inner_juncs],
                                                                                   set(inner_juncs), my_juncs - set(inner_juncs), {outer_junc, }, set(my_exons))

                        for j in inner_juncs:
                            if j not in esjuncs:
                                esjuncs[j] = {}
                            if outer_junc not in esjuncs[j]:
                                esjuncs[j][outer_junc] = {s: 0 for s in all_samples}
                            esjuncs[j][outer_junc] = add_counts_to_dict(esjuncs[j][outer_junc], inner_juncs_to_counts[inner_juncs])

                    if any([othercounts[s] >= min_read_support for s in all_samples]):
                        event_to_info[ename].other = othercounts
                    event_to_info[ename].totoverlap = get_overlapping_reads(outer_junc, interval_to_reads, all_samples)

    return event_to_info, esjuncs


def get_overlapping_reads(ref_junc, interval_to_reads, all_samples):
    tot_counts = {s: 0 for s in all_samples}
    for i in interval_to_reads:
        if i[0] < ref_junc[0] and ref_junc[1] < i[1]:
            tot_counts = add_counts_to_dict(tot_counts, interval_to_reads[i])
    return tot_counts


def add_counts_to_dict(a, b):
    return {s: a[s] + b[s] for s in a.keys()}

def process_terminal_exons(term_exon, eventtype, chrom, strand, gene, all_samples, mycolor, min_read_support, junc_frac_of_event, event_support, interval_to_reads, alljuncs, annot_terminal_ss):  # noqa: C901 - FIXME: reduce complexity
    event_to_info = {}
    junc_combos = set()
    if len(term_exon) > 1:
        tot_counts = {s: 0 for s in all_samples}
        othercounts = {s: 0 for s in all_samples}
        ss_to_counts = {}
        ss_to_junc_qual = {}
        actual_terminal = []
        for term_ss in term_exon:
            ss_to_counts[term_ss] = {s: 0 for s in all_samples}
            ss_to_junc_qual[term_ss] = []
            is_annot_terminal = False
            if gene in annot_terminal_ss:
                for ss in annot_terminal_ss[gene]:
                    if ss - 20 <= term_ss <= ss + 20:
                        is_annot_terminal = True
                        break
            for term_junc in term_exon[term_ss]:
                ss_to_counts[term_ss] = add_counts_to_dict(ss_to_counts[term_ss], term_exon[term_ss][term_junc])
                # here we determine whether this splice site is really a terminal exon
                # either all reads terminate here, or it's annotated as a terminal junction
                if is_annot_terminal:
                    ss_to_junc_qual[term_ss].append(True)
                else:
                    ss_to_junc_qual[term_ss].append(any([term_exon[term_ss][term_junc][s] == alljuncs[term_junc][s] and term_exon[term_ss][term_junc][s] >= min_read_support for s in all_samples]))

        for term_ss in term_exon:
            if all(ss_to_junc_qual[term_ss]):
                actual_terminal.append(term_ss)
            else:
                othercounts = add_counts_to_dict(othercounts, ss_to_counts[term_ss])
            tot_counts = add_counts_to_dict(tot_counts, ss_to_counts[term_ss])

        if any([tot_counts[s] >= event_support for s in all_samples]):
            good_b = []
            good_b_strict = []
            my_juncs = set()
            my_juncs = set()
            for term_ss in actual_terminal:
                inc_frac_list = [ss_to_counts[term_ss][s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]

                if any([ss_to_counts[term_ss][s] >= min_read_support for s in all_samples]) \
                    and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                         or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                    good_b.append(term_ss)
                    if any([ss_to_counts[term_ss][s] / tot_counts[s] >= MIN_TERMINAL_SS_FRACTION if tot_counts[s] > 0 else False for s in all_samples]):
                        good_b_strict.append(term_ss)
                    for term_junc in term_exon[term_ss]:
                        if any([term_exon[term_ss][term_junc][s] >= min_read_support for s in all_samples]) \
                                and any([junc_frac_of_event <= term_exon[term_ss][term_junc][s] / tot_counts[s] <= 1 - junc_frac_of_event if tot_counts[s] > 0 else False for s in all_samples]):
                            my_juncs.add(term_junc)
                else:
                    othercounts = add_counts_to_dict(othercounts, ss_to_counts[term_ss])

            # check that the junctions involved are distant enough from each other
            if len(good_b) > 1 and len(good_b_strict) > 1 and max(good_b_strict) - min(good_b_strict) >= MIN_TERMINAL_JUNCTION_SEPARATION:
                ename = f'{eventtype}-({strand})-{gene}'
                event_to_info[ename] = SplicingEvent(ename, eventtype, gene, chrom, strand, tot_counts, all_samples)

                for ss_b in good_b:
                    jname = f'{chrom}:{ss_b}'
                    if (eventtype == 'AFE' and strand == '+') or (eventtype == 'ALE' and strand == '-'):
                        bs, be = ss_b - 100, ss_b
                    else:
                        bs, be = ss_b, ss_b + 100
                    bedline = Bed(chrom, bs, be, name=f'{jname}_{ename}', score=0,
                                  strand=strand, thickStart=bs, thickEnd=be,
                                  itemRgb=mycolor, blocks=[BedBlock(bs, be)])
                    inc_juncs = {x for x in my_juncs if ss_b in x}
                    event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], ss_to_counts[ss_b], inc_juncs, my_juncs - inc_juncs)

                if any([othercounts[s] >= min_read_support for s in all_samples]):
                    event_to_info[ename].other = othercounts
                junc_combos.add(frozenset(my_juncs))
                ref_junc = (max([x[0] for x in my_juncs]), min([x[1] for x in my_juncs]))
                event_to_info[ename].totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, all_samples)

            # TODO: add initial check for all reads with junction are equal to all terminal exon reads!
            # TODO: identify good juncs for each terminal SS, save junc combo
            # Good juncs must be significant proportion of all reads for event!
    return event_to_info, junc_combos


def process_junction_events(ss_a_to_b, esjuncs, eventtype, chrom, strand, gene, all_samples, colordict, min_read_support, interval_to_reads, junc_frac_of_event, event_support, afe_ale_junc_comb):  # noqa: C901 - FIXME: reduce complexity
    event_to_info = {}
    ss_b_groups_to_ss_a = {}
    for ss_a in ss_a_to_b:
        if len(ss_a_to_b[ss_a]) > 1:  # don't report when no alternative splicing at all
            good_b, good_b_strict = [], []
            junccomb = set()
            my_juncs = set()
            othercounts = {s: 0 for s in all_samples}
            tot_counts = {s: 0 for s in all_samples}
            for ss_b in ss_a_to_b[ss_a]:
                tot_counts = add_counts_to_dict(tot_counts, ss_a_to_b[ss_a][ss_b])

            for ss_b in ss_a_to_b[ss_a]:

                junc = (min((ss_a, ss_b)), max((ss_a, ss_b)))
                inc_frac_list = [ss_a_to_b[ss_a][ss_b][s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]

                if len(inc_frac_list) >= 1 \
                    and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                         or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):

                    good_b.append(ss_b)
                    my_juncs.add(junc)
                    junccomb.add(junc)
                    if junc not in esjuncs:
                        good_b_strict.append(ss_b)
                    else:
                        escounts = {s: 0 for s in all_samples}
                        for ss_b2 in ss_a_to_b[ss_a]:
                            if ss_b2 != ss_b:
                                junc2 = (min((ss_a, ss_b2)), max((ss_a, ss_b2)))
                                if junc2 in esjuncs[junc]:
                                    escounts = add_counts_to_dict(escounts, esjuncs[junc][junc2])
                        deltafrac = [(ss_a_to_b[ss_a][ss_b][s] - escounts[s]) / ss_a_to_b[ss_a][ss_b][s] if ss_a_to_b[ss_a][ss_b][s] > 0 else 0 for s in all_samples]
                        delta = [ss_a_to_b[ss_a][ss_b][s] - escounts[s] for s in all_samples]
                        if any([deltafrac[x] > 0.1 and delta[x] > min_read_support for x in range(len(delta))]):
                            good_b_strict.append(ss_b)
                else:
                    othercounts = add_counts_to_dict(othercounts, ss_a_to_b[ss_a][ss_b])

            if frozenset(junccomb) not in afe_ale_junc_comb:
                if len(good_b) > 1:
                    ss_b_group = frozenset(good_b)
                    if ss_b_group not in ss_b_groups_to_ss_a:
                        ss_b_groups_to_ss_a[ss_b_group] = {}
                    ss_b_groups_to_ss_a[ss_b_group][ss_a] = tot_counts

                if any([tot_counts[s] >= event_support for s in all_samples]) \
                        and len(good_b_strict) > 1:  # only if alt splicing, not including junctions involved in exon skipping
                    ename = f'{eventtype}-relTo-{chrom}:{ss_a}({strand})-{gene}'
                    event_to_info[ename] = SplicingEvent(ename, eventtype, gene, chrom, strand, tot_counts, all_samples)
                    for ss_b in good_b:
                        junc = (min((ss_a, ss_b)), max((ss_a, ss_b)))
                        jname = f'{chrom}:{junc[0]}-{junc[1]}'
                        bs, be = min((ss_a, ss_b)), max((ss_a, ss_b))
                        bedline = Bed(chrom, bs, be, name=f'{jname}_{ename}', score=0,
                                      strand=strand, thickStart=bs, thickEnd=be,
                                      itemRgb=colordict[eventtype], blocks=[BedBlock(bs, be)])
                        event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], ss_a_to_b[ss_a][ss_b], {junc, }, my_juncs - {junc, })

                    if any([othercounts[s] >= min_read_support for s in all_samples]):
                        event_to_info[ename].other = othercounts
                    # against min(good_b), not the ss_b the loop above left behind.  All
                    # good_b lie on the same side of ss_a, so this is the same junction
                    # the leaked value gave, said properly
                    ref_junc = (ss_a, min(good_b)) if ss_a < min(good_b) else (max(good_b), ss_a)
                    event_to_info[ename].totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, all_samples)
    eventtype = 'alt3' if eventtype == 'alt5' else 'alt5'
    for ss_b_group in ss_b_groups_to_ss_a:
        if len(ss_b_groups_to_ss_a[ss_b_group]) > 1:
            tot_counts = {s: 0 for s in all_samples}
            othercounts = {s: 0 for s in all_samples}
            for ss_a in ss_b_groups_to_ss_a[ss_b_group]:
                tot_counts = add_counts_to_dict(tot_counts, ss_b_groups_to_ss_a[ss_b_group][ss_a])
            if any([tot_counts[s] >= event_support for s in all_samples]):
                good = []
                my_juncs = set()
                for ss_a in ss_b_groups_to_ss_a[ss_b_group]:
                    inc_frac_list = [ss_b_groups_to_ss_a[ss_b_group][ss_a][s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]
                    if (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                            or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                        good.append(ss_a)
                        for ss_b in ss_b_group:
                            my_juncs.add((min((ss_a, ss_b)), max((ss_a, ss_b))))
                    else:
                        othercounts = add_counts_to_dict(othercounts, ss_b_groups_to_ss_a[ss_b_group][ss_a])
                if len(good) > 1:
                    oj_string = [f'{chrom}:{x}' for x in sorted(list(ss_b_group))]
                    ename = f'{eventtype}ss-relTo-{",".join(oj_string)}({strand})-{gene}'
                    event_to_info[ename] = SplicingEvent(ename, f'{eventtype}ss', gene, chrom, strand, tot_counts, all_samples)
                    for ss_a in good:
                        if (eventtype == 'alt5' and strand == '+') or (eventtype == 'alt3' and strand == '-'):
                            bs, be = ss_a - 10, ss_a
                        else:
                            bs, be = ss_a, ss_a + 10
                        jname = f'{chrom}:{ss_a}'
                        inc_juncs = set()
                        for ss_b in ss_b_group:
                            inc_juncs.add((min((ss_a, ss_b)), max((ss_a, ss_b))))

                        bedline = Bed(chrom, bs, be, name=f'{jname}_{ename}', score=0,
                                      strand=strand, thickStart=bs, thickEnd=be,
                                      itemRgb=colordict[eventtype], blocks=[BedBlock(bs, be)])
                        event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bedline], ss_b_groups_to_ss_a[ss_b_group][ss_a], inc_juncs, my_juncs - inc_juncs)
                    if any([othercounts[s] >= min_read_support for s in all_samples]):
                        event_to_info[ename].other = othercounts
                    ref_junc = (max(good), min(ss_b_group)) if median(good) < median(ss_b_group) else (max(ss_b_group), min(good))
                    event_to_info[ename].totoverlap = get_overlapping_reads(ref_junc, interval_to_reads, all_samples)

    return event_to_info

class SplicingEvent:
    """Represents metadata about a splicing event"""
    def __init__(self, name, eventtype, gene, chrom, strand, totjunc, all_samples):
        self.name = name
        self.chrom = chrom
        self.eventtype = eventtype
        self.gene = gene
        self.strand = strand
        self.totjunc = totjunc
        self.events = {}
        self.totoverlap = None
        self.other = {s: 0 for s in all_samples}

class SplicingEventJunction:
    """represents metadata about a specific junction in a splicing event"""
    def __init__(self, name, bedlines, samplecounts, inc_junc, exc_junc, outer_junc=None, inc_exon=None):
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


def write_intron_retention(alljuncs, all_samples, allblocks, chrom, strand, gene, mycolor, min_read_support, interval_to_reads, junc_frac_of_event, event_support):
    event_to_info = {}
    for j in alljuncs:
        retainedcounts = {s: 0 for s in all_samples}
        for b in allblocks:  # this nested for loop could get brutal
            if b[0] < j[0] and j[1] < b[1]:
                for s in all_samples:
                    retainedcounts[s] += allblocks[b][s]
        splicedcounts = alljuncs[j]
        # if sum(allretainedcounts.values()) > 0:  # only report if junction retained at all
        tot_counts = splicedcounts
        tot_counts = add_counts_to_dict(tot_counts, retainedcounts)

        inc_frac_list = [retainedcounts[s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]
        if len(inc_frac_list) >= 1 \
            and any([splicedcounts[s] >= min_read_support for s in all_samples]) \
            and any([retainedcounts[s] >= min_read_support for s in all_samples]) \
            and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                 or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
            ename = f'ir-of-{chrom}:{j[0]}-{j[1]}({strand})-{gene}'
            event_to_info[ename] = SplicingEvent(ename, 'ir', gene, chrom, strand, tot_counts, all_samples)

            bed_start, bed_end = j[0] - FLANK_STUB, j[1] + FLANK_STUB
            bedline = Bed(chrom, bed_start, bed_end, name=f'spliced_{ename}', score=0,
                          strand=strand, thickStart=bed_start, thickEnd=bed_end,
                          itemRgb=mycolor,
                          blocks=[BedBlock(bed_start, j[0]), BedBlock(j[1], bed_end)])
            event_to_info[ename].events['spliced'] = SplicingEventJunction('spliced', [bedline], splicedcounts, {j, }, {})

            bedline = Bed(chrom, j[0], j[1], name=f'retained_{ename}', score=0,
                          strand=strand, thickStart=j[0], thickEnd=j[1],
                          itemRgb=mycolor, blocks=[BedBlock(j[0], j[1])])
            event_to_info[ename].events['retained'] = SplicingEventJunction('retained', [bedline], retainedcounts, {}, {j, })

            event_to_info[ename].totoverlap = get_overlapping_reads(j, interval_to_reads, all_samples)
    return event_to_info

def write_ends(grouped_ends, all_samples, chrom, strand, gene, eventtype, mycolor, support, junc_frac_of_event, event_support):
    good_ends = set()
    event_to_info = {}
    othercounts = {s: 0 for s in all_samples}
    tot_counts = {s: 0 for s in all_samples}
    for e in grouped_ends:
        tot_counts = add_counts_to_dict(tot_counts, grouped_ends[e])
    if any([tot_counts[s] >= event_support for s in all_samples]):
        for e in grouped_ends:
            inc_frac_list = [grouped_ends[e][s] / tot_counts[s] for s in all_samples if tot_counts[s] >= event_support]

            if any([grouped_ends[e][s] >= support for s in all_samples]) \
                and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                     or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                good_ends.add(e)
            else:
                othercounts = add_counts_to_dict(othercounts, grouped_ends[e])
        if len(good_ends) > 1:
            ename = f'{eventtype}-({strand})-{gene}'
            event_to_info[ename] = SplicingEvent(ename, eventtype, gene, chrom, strand, tot_counts, all_samples)

            for e in good_ends:
                jname = f'{chrom}:{e}'
                bed = Bed(chrom, e - 1, e, name=f'{jname}_{ename}', score=0, strand=strand,
                          thickStart=e - 1, thickEnd=e, itemRgb=mycolor,
                          blocks=[BedBlock(e - 1, e)])
                event_to_info[ename].events[jname] = SplicingEventJunction(jname, [bed], grouped_ends[e], {e, }, good_ends - {e, })
            if any([othercounts[s] >= support for s in all_samples]):
                event_to_info[ename].other = othercounts
    return event_to_info


def psi_values(ncounts, totals, all_samples, event_support):
    "PSI of each sample, empty where the sample has too few reads to take one from"
    return [str(round(ncounts[s] / totals[s], 4)) if totals[s] >= event_support else ''
            for s in all_samples]

def write_counts_psi(info, ncounts, junctot, fulltot, all_samples, outcounts, outpsijunc, outpsitot, event_support):
    outcounts.writeEvent(info, [str(ncounts[s]) for s in all_samples])
    juncpsi = [ncounts[s] / junctot[s] if junctot[s] >= event_support else 'NA' for s in all_samples]
    outpsijunc.writeEvent(info, psi_values(ncounts, junctot, all_samples, event_support))
    outpsitot.writeEvent(info, psi_values(ncounts, fulltot, all_samples, event_support))
    return juncpsi

def get_junc_coords(chrom, juncs):
    "the coordinates of a set of junctions or exons, as the event columns hold them"
    juncs = sorted(juncs)
    if len(juncs) == 0:
        return []
    elif isinstance(juncs[0], int):
        return [f'{chrom}:{x}' for x in juncs]
    else:
        return [f'{chrom}:{x[0]}-{x[1]}' for x in juncs]


def get_psi_and_filter(event_to_info, all_samples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot, event_support, outoutlier, outolfilt):  # noqa: C901 - FIXME: reduce complexity
    sig_events = []
    etype_to_sig = {}
    for ename in event_to_info:
        event = event_to_info[ename]
        if event.totoverlap is None:
            event.totoverlap = event.totjunc

        if any([event_frac_of_tot <= event.totjunc[s] / event.totoverlap[s] if event.totoverlap[s] > 0 else False for s in all_samples]):  # check that junction/event is expressed enough relative to all reads covering locus
            # do an extra check that one junction that is not exclusion is at least junc_frac_of_event of total locus
            one_event_pass = False
            for jname in event.events:
                jinfo = event.events[jname]
                inc_frac_list = [jinfo.samplecounts[s] / event.totoverlap[s] for s in all_samples if event.totoverlap[s] > 0]

                if (event.eventtype not in {'es', 'ces'} or jname != 'exc') \
                    and (max(inc_frac_list) - min(inc_frac_list) >= junc_frac_of_event
                         or (len(all_samples) == 1 and junc_frac_of_event <= inc_frac_list[0] <= 1 - junc_frac_of_event)):
                    one_event_pass = True

            if one_event_pass:
                for jname in event.events:
                    jinfo = event.events[jname]

                    outinfo = [f'{jname}_{ename}', event.eventtype, event.gene, get_junc_coords(event.chrom, jinfo.inc_juncs),
                               get_junc_coords(event.chrom, jinfo.exc_juncs), get_junc_coords(event.chrom, jinfo.outer_juncs), get_junc_coords(event.chrom, jinfo.inc_exon)]
                    juncpsi = write_counts_psi(outinfo, jinfo.samplecounts, event.totjunc, event.totoverlap, all_samples, outcounts, outpsijunc, outpsitot, event_support)

                    # every sample can be below event_support, leaving no PSI to
                    # take a median of and nothing to score the BED with
                    vals_for_outlier = [x for x in juncpsi if x != 'NA']
                    if len(vals_for_outlier) > 0:
                        med = median(vals_for_outlier)
                        for bed in jinfo.bedlines:
                            bed.score = round(med * 100)
                            bed.write(outbed)

                        if (outoutlier is not None) and (len(vals_for_outlier) >= 5):
                            dev = sps.iqr(vals_for_outlier) / 2
                            for s in all_samples:
                                count = jinfo.samplecounts[s]
                                tot = event.totjunc[s]
                                psi = count / tot if tot > event_support else 'NA'
                                if psi != 'NA' and abs(psi - med) >= 0.1:  # 10% PSI
                                    delta = abs(psi - med)
                                    read_counts = [count, tot]
                                    # event_name, gene, median, dev, psi, numsamplesused, count/junctot, absdeltapsi, numberofdev
                                    if (dev == 0 and psi != med) or delta > dev * 3:
                                        outline = outinfo[:3] + [s, round(med, 6), round(dev, 6), round(psi, 6), len(vals_for_outlier), read_counts, round(delta, 6), 100000]
                                        if dev != 0:
                                            outline[-1] = round(delta / dev, 6)
                                        sig_events.append(outline)
                                        if event.eventtype not in etype_to_sig:
                                            etype_to_sig[event.eventtype] = {}
                                        all_these_juncs = (frozenset(jinfo.inc_juncs), frozenset(jinfo.exc_juncs), frozenset(jinfo.outer_juncs))
                                        if all_these_juncs not in etype_to_sig[event.eventtype]:
                                            etype_to_sig[event.eventtype][all_these_juncs] = {}
                                        etype_to_sig[event.eventtype][all_these_juncs][s] = outline
                if any([event.other[s] > 0 for s in all_samples]):
                    outinfo = [f'other_{ename}', event.eventtype, event.gene, [], [], [], []]
                    write_counts_psi(outinfo, event.other, event.totjunc, event.totoverlap, all_samples, outcounts, outpsijunc, outpsitot, event_support)
    if outoutlier is not None:
        sig_events.sort(reverse=True, key=lambda x: x[::-1])
        for line in sig_events:
            outoutlier.writeRow(line)

        seen_junctions, condensed_sig = {}, []
        seen_es_juncs = {}
        TEMP_OUTLIER_DEV_THRESHOLD = 6

        for et in ['es', 'AFE', 'ALE', 'ir', 'alt5ss', 'alt3ss']:  # alt5ss and alt3ss overlaps with ALE and AFE should have already been filtered out when generating events
            if et in etype_to_sig:
                for juncs in etype_to_sig[et]:
                    sample_to_vals = etype_to_sig[et][juncs]
                    if et != 'ir':
                        inc_juncs, exc_juncs, outer_juncs = juncs
                        for j1 in inc_juncs:
                            for j2 in outer_juncs:
                                k = frozenset((j1, j2))
                                if k not in seen_junctions:
                                    seen_junctions[k] = {}
                                seen_junctions[k][juncs] = sample_to_vals
                        if et == 'es':
                            all_juncs = (inc_juncs | exc_juncs) | outer_juncs
                            seen_es_juncs[all_juncs] = sample_to_vals
                    for s in sample_to_vals:
                        condensed_sig.append(sample_to_vals[s])

        for et in ['ces', 'alt5', 'alt3']:
            if et in etype_to_sig:
                for juncs in etype_to_sig[et]:
                    inc_juncs, exc_juncs, outer_juncs = juncs

                    sample_to_vals = etype_to_sig[et][juncs]
                    good_s = set()
                    has_overlap = False
                    if et == 'ces':  # only for ces, check for subet of full es juncs
                        if len(inc_juncs) == 0:
                            all_juncs = exc_juncs | outer_juncs  # exclusion
                        else:
                            all_juncs = inc_juncs | outer_juncs  # include only the specific exon that is significatn
                        for j2 in seen_es_juncs:
                            if all_juncs - j2 == frozenset():  # is subset of other
                                has_overlap = True
                                for s in sample_to_vals:
                                    if s not in seen_es_juncs[j2]:
                                        good_s.add(s)
                                    else:
                                        event_dev = sample_to_vals[s][-1]
                                        other_dev = seen_es_juncs[j2][s][-1]
                                        if other_dev / event_dev < 0.8:  # deviation is 20% better for this event
                                            good_s.add(s)
                    else:
                        all_juncs = inc_juncs | exc_juncs
                        for jkey in seen_junctions:
                            if all_juncs - jkey == frozenset():  # is subset of other
                                has_overlap = True
                                for j2 in seen_junctions[jkey]:
                                    for s in sample_to_vals:
                                        if s not in seen_junctions[jkey][j2]:
                                            good_s.add(s)
                                        else:
                                            event_dev = sample_to_vals[s][-1]
                                            other_dev = seen_junctions[jkey][j2][s][-1]
                                            if other_dev / event_dev < 0.8:  # deviation is 20% better for this event
                                                good_s.add(s)

                    # adding good values to output
                    if not has_overlap:
                        good_s = sample_to_vals.keys()
                        # adding junctions to reference
                        if all_juncs not in seen_junctions:
                            seen_junctions[all_juncs] = {}
                        seen_junctions[all_juncs][juncs] = sample_to_vals
                    for s in good_s:
                        condensed_sig.append(sample_to_vals[s])
        condensed_sig.sort(reverse=True, key=lambda x: x[::-1])
        for line in condensed_sig:
            if line[-1] > TEMP_OUTLIER_DEV_THRESHOLD:
                outolfilt.writeRow(line)


def process_gene_to_events(temp_prefix, chrom, all_samples, all_genes_to_juncs, gene_to_strand, junc_support, output_read_ends, event_frac_of_tot, junc_frac_of_event, event_support, annot_afe_ss, annot_ale_ss, check_outliers):
    etype_to_color = {'skipped_exons': '66,105,245', 'retained_introns': '144,66,245',
                      'alt3': '245,215,66', 'alt5': '43,184,39', 'tss': '255,0,0', 'tts': '0,0,255'}
    # outends is just for writing out read ends for Harrison - he can group them more intelligently
    outends = None
    if output_read_ends:
        outends = open(temp_prefix + '.diffsplice.readends.bed', 'w')

    allgenes = set.union(*[set(all_genes_to_juncs[s].keys()) for s in range(len(all_genes_to_juncs))])
    outoutlier, outolfilt = None, None
    if check_outliers:
        outoutlier = SplicingOutliersWriter(temp_prefix + '.diffsplice.outliers.tsv')
        outolfilt = SplicingOutliersWriter(temp_prefix + '.diffsplice.outliers.filtered.tsv')
    with open(temp_prefix + '.diffsplice.bed', 'w') as outbed, \
            SplicingEventsWriter(temp_prefix + '.diffsplice.counts.tsv', all_samples) as outcounts, \
            SplicingEventsWriter(temp_prefix + '.diffsplice.PSIjunc.tsv', all_samples) as outpsijunc, \
            SplicingEventsWriter(temp_prefix + '.diffsplice.PSItot.tsv', all_samples) as outpsitot:

        for gene in allgenes:

            # get gene strand - this is not optimal, def need to revamp how getting annot info
            strand = gene_to_strand[gene]
            ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj, interval_to_reads, afe, ale, outer_junc_to_exons = extract_splicing_info(all_samples, all_genes_to_juncs,
                                                                                                                                                                      gene, strand, chrom, outends)

            # exon skipping
            # TODO combine exon skipping at multiple junctions (so all junctions that skip exon are combined)
            es_info, esjuncs = write_exon_skipping(exonjpairs, alljuncs, all_samples, chrom, strand, gene,
                                                   etype_to_color['skipped_exons'], junc_support, interval_to_reads, junc_frac_of_event, event_support, outer_junc_to_exons)
            afe_info, afe_junc_comb = process_terminal_exons(afe, 'AFE', chrom, strand, gene, all_samples, etype_to_color['tss'], junc_support, junc_frac_of_event, event_support, interval_to_reads, alljuncs, annot_afe_ss)
            ale_info, ale_junc_comb = process_terminal_exons(ale, 'ALE', chrom, strand, gene, all_samples, etype_to_color['tts'], junc_support, junc_frac_of_event, event_support, interval_to_reads, alljuncs, annot_ale_ss)
            afe_ale_junc_comb = afe_junc_comb | ale_junc_comb
            a3_info = process_junction_events(ss5to3, esjuncs, 'alt3', chrom, strand, gene, all_samples, etype_to_color,
                                              junc_support, interval_to_reads, junc_frac_of_event, event_support, afe_ale_junc_comb)
            a5_info = process_junction_events(ss3to5, esjuncs, 'alt5', chrom, strand, gene, all_samples, etype_to_color,
                                              junc_support, interval_to_reads, junc_frac_of_event, event_support, afe_ale_junc_comb)
            # intron retention
            ir_info = write_intron_retention(alljuncs, all_samples, allblocks, chrom, strand, gene,
                                             etype_to_color['retained_introns'], junc_support, interval_to_reads, junc_frac_of_event, event_support)

            grouped_ends = group_ends(t_starts_ends, all_samples, 100)
            if strand == '-':
                grouped_ends = grouped_ends[::-1]
                t_first_last_sj = t_first_last_sj[::-1]
            tss_info = write_ends(grouped_ends[0], all_samples, chrom, strand, gene, 'tss', etype_to_color['tss'], junc_support, junc_frac_of_event, event_support)
            tts_info = write_ends(grouped_ends[1], all_samples, chrom, strand, gene, 'tts', etype_to_color['tts'], junc_support, junc_frac_of_event, event_support)
            # write_ends(t_first_last_sj[0], all_samples, chrom, strand, gene, 'tss', etype_to_color['tss'], out, outbed, support)
            # write_ends(t_first_last_sj[1], all_samples, chrom, strand, gene, 'tts', etype_to_color['tts'], out, outbed, support)
            all_info = (((a3_info | a5_info) | (es_info | ir_info)) | (tts_info | tss_info)) | (afe_info | ale_info)

            get_psi_and_filter(all_info, all_samples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot, event_support, outoutlier, outolfilt)

    if check_outliers:
        outoutlier.close()
        outolfilt.close()
    if output_read_ends:
        outends.close()


def generate_good_match_to_annot(args, temp_prefix, region, bamfile_name, region_annot, region_annot_fa, clipping_file):
    if not args.noaligntoannot:
        pipettor.run([('samtools', 'view', '-h', bamfile_name, f'{region.name}:{region.start}-{region.end}'),
                      ('samtools', 'fasta', '-')],
                     stdout=temp_prefix + '.reads.fasta')
        mm2_cmd = ['minimap2', '-a', '-N', '4', '--MD', region_annot_fa, temp_prefix + '.reads.fasta']
        # 1 thread because already multithreaded here
        # minimap2 borrows threads left idle in the thread budget while it runs
        with thread_share.program_threads() as mm2_threads:
            run_count_sam_transcripts(
                mm2_cmd=mm2_cmd[:1] + ['-t', str(mm2_threads)] + mm2_cmd[1:],
                output=temp_prefix + '.matchannot.counts.tsv',
                threads=1,
                check_splice=True,
                isoforms=region_annot,
                trimmedreads=clipping_file,
                allow_UTR_indels=True,
                soft_clipping_buffer=10,
                output_endpos=temp_prefix + '.readtoends.tsv')
        return temp_prefix + '.readtoends.tsv'
    else:
        return None


def read_annot_alignments(good_annot_aligns):
    """Read-to-transcript assignments from the ends file.  With --noaligntoannot
    there is no file and no assignments."""
    read_to_transcript = {}
    if good_annot_aligns is None:
        return read_to_transcript
    for row in ReadEndsReader(good_annot_aligns):
        # an empty start_sj_index is a single exon transcript, with no junction to
        # measure the read's ends from
        if row.start_sj_index is not None:
            read_to_transcript[row.read] = (row.transcript, row.start_sj_index,
                                            row.start_sj_dist, row.end_sj_index,
                                            row.end_sj_dist)
    return read_to_transcript


def remove_region_temp_files(temp_prefix, aligned_to_annot):
    """Drop the per-region intermediates.  The annotation-alignment ones only exist
    when that step ran."""
    temp_files = [f'{temp_prefix}.reads.genomicclipping.tsv']
    if aligned_to_annot:
        temp_files += [f'{temp_prefix}.matchannot.counts.tsv',
                       f'{temp_prefix}.readtoends.tsv',
                       f'{temp_prefix}.reads.fasta']
    pipettor.run([tuple(['rm'] + temp_files)])


def get_juncs_single_sample(args, region, temp_prefix, sample, bamfile_name, region_annot, region_annot_fa, junction_corrector, annots):
    temp_prefix = temp_prefix + '_' + sample

    bam_file = pysam.AlignmentFile(bamfile_name, 'rb')
    num_reads, clipping_file = generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, region.name, region.start, region.end)
    bam_file.close()

    # FIXME: what format is this? reading into memory, why write file?
    good_annot_aligns = generate_good_match_to_annot(args, temp_prefix, region, bamfile_name, region_annot, region_annot_fa, clipping_file)

    read_to_transcript = read_annot_alignments(good_annot_aligns)

    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')
    sj_to_ends = {}
    filter_correct_group_reads(bam_file=bamfile, region=region,
                               read_to_annot_transcript=read_to_transcript,
                               annots=annots, junction_corrector=junction_corrector,
                               genome=None,
                               quality=args.quality, keep_sup=args.keep_sup,
                               sj_to_ends=sj_to_ends,
                               allow_outside_range=True,
                               keep_single_exon=False)
    bamfile.close()

    gene_to_juncs, nogenejuncs, sereads = group_juncs_by_annot_gene(sj_to_ends, annots.sjc_to_gene, annots.junc_to_gene_id, annots.gene_to_exons, annots.gene_to_annot_juncs)

    c = 0
    with GeneJuncsWriter(temp_prefix + '_gene_to_juncs.tsv') as out:
        for gene in gene_to_juncs:
            for juncs in gene_to_juncs[gene]:
                for read_info in gene_to_juncs[gene][juncs].reads:
                    c += 1
                    out.writeRow((gene, juncs, read_info.start, read_info.end,
                                  gene_to_juncs[gene][juncs].strand, read_info.name))
    remove_region_temp_files(temp_prefix, good_annot_aligns is not None)


def process_bed_line(bed_rec):
    transcript = bed_rec.name
    gene = bed_rec.extraCols[GENE_ID_EXTRA_COL]
    exons = [(blk.start, blk.end) for blk in bed_rec.blocks]
    junctions = [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)]
    strand = bed_rec.strand
    return gene, transcript, exons, junctions, strand


def _run_region(*, partition, gtf_data, junction_corrector, args, all_samples):  # noqa: C901 - FIXME: reduce complexity
    # FIXMEL what are these files being created, just do in memory
    region_bed = partition.output_path('region.bed')
    with open(region_bed, 'w') as out:
        Bed(partition.region.name, partition.region.start, partition.region.end).write(out)

    region_annot = partition.file_prefix + '.annotation.bed'
    pipettor.run([('bedtools', 'intersect', '-wa', '-a', args.annot, '-b', region_bed)], stdout=region_annot)

    if os.path.getsize(region_annot) > 0:  # check if any annotated transcripts in region
        if args.annot_basic:
            region_annot_basic = partition.file_prefix + '.annotation.basic.bed'
            pipettor.run([('bedtools', 'intersect', '-wa', '-a', args.annot_basic, '-b', region_bed)], stdout=region_annot_basic)

        region_annot_fa = None
        if not args.noaligntoannot:
            region_annot_fa = partition.file_prefix + '.annotation.fa'
            get_sequence_from_bed(args.genome, region_annot, region_annot_fa)

        region_juncs = None
        if args.junction_bed:
            region_juncs = partition.file_prefix + '.juncbed.bed'
            pipettor.run([('bedtools', 'intersect', '-wa', '-a', args.junction_bed, '-b', region_bed)], stdout=region_juncs)

        annots = annot_data_from_gtf(gtf_data, partition.region)

        annot_afe_ss, annot_ale_ss = {}, {}
        if args.annot_basic:
            for bed_rec in BedReader(region_annot_basic, numStdCols=12, fixScores=True):
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

        for sample, bamfile in all_samples:
            if not os.path.exists(partition.file_prefix + '_' + sample + '_gene_to_juncs.tsv'):
                get_juncs_single_sample(args, partition.region, partition.file_prefix, sample, bamfile, region_annot, region_annot_fa, junction_corrector, annots)

        all_genes_to_juncs = []
        for sample, bamfile in all_samples:
            gene_to_juncs = {}
            for row in GeneJuncsReader(partition.file_prefix + '_' + sample + '_gene_to_juncs.tsv'):
                juncs = row.juncs
                if row.gene not in gene_to_juncs:
                    gene_to_juncs[row.gene] = {}
                if juncs not in gene_to_juncs[row.gene]:
                    gene_to_juncs[row.gene][juncs] = []
                gene_to_juncs[row.gene][juncs].append(
                    ReadRec(None, row.strand, (), row.start, row.end, row.read))
            all_genes_to_juncs.append(gene_to_juncs)

        process_gene_to_events(
            partition.file_prefix, partition.region.name, [x[0] for x in all_samples],
            all_genes_to_juncs, annots.gene_to_strand, args.junc_support, args.output_read_ends,
            args.event_frac_of_tot, args.junc_frac_of_event, args.event_support,
            annot_afe_ss, annot_ale_ss, args.check_outliers)

        if not args.keep_intermediate:
            for sample, bamfile in all_samples:
                pipettor.run([('rm', partition.file_prefix + '_' + sample + '_gene_to_juncs.tsv')])


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


def main():  # noqa: C901 - FIXME: reduce complexity
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    check_genome_fasta(args.genome)
    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)

    # FIXME: make this a common function somewhere
    temp_dir = make_temp_dir(args.output)
    logging.info('temp directory: %s', temp_dir)

    if args.region_bed:
        all_regions = [SeqRange(bed.chrom, bed.chromStart, bed.chromEnd)
                       for bed in BedReader(args.region_bed)]
        all_regions = combine_regions(all_regions)
    else:
        all_regions = [SeqRange(chrom, 0, genome.get_reference_length(chrom))
                       for chrom in genome.references]
    logging.info(f'input regions: {len(all_regions)}')

    all_samples = []
    for line in fileOps.iterLines(args.manifest):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            sample, bamfile = line
            all_samples.append((sample, bamfile))

    if args.annot:
        logging.info('loading annotation GTF')
        annot_gtf_data = gtf_data_parser(args.annot, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)
    logging.info('building intron support database')
    junction_corrector = junction_corrector_factory(args.ss_window, args.junction_support,
                                                    annot_gtf_data=annot_gtf_data,
                                                    intron_beds=args.junction_bed,
                                                    star_sj_tabs=args.junction_star)

    # FIXME: we appear to do this just to do bed intersect; just do in memory.
    annot_bed = temp_dir + '/annotation.bed'
    if not os.path.exists(annot_bed):
        pipettor.run([(sys.executable, gtf_to_bed_prog, args.annot, annot_bed, '--include_gene')])
    args.annot = annot_bed
    if args.annot_basic:
        annot_basic_bed = temp_dir + '/annotation.basic.bed'
        if not os.path.exists(annot_basic_bed):
            pipettor.run([(sys.executable, gtf_to_bed_prog, args.annot_basic, annot_basic_bed, '--include_gene')])
        args.annot_basic = annot_basic_bed

    logging.info(f'running regions with {args.threads} threads')

    runner = PartitionRunner(all_regions, temp_dir, gtf_data=annot_gtf_data, junction_corrector=junction_corrector,
                             threads=args.threads)
    runner.run(_run_region, args=args, all_samples=all_samples)

    combine_temp_files_by_suffix(args.output, [p.file_prefix for p in runner],
                                 ['.diffsplice.bed'] + list(EVENT_SUFFIXES),
                                 headers=EVENT_SUFFIXES)
    sort_bed_file(args.output + '.diffsplice.bed')
    if args.output_read_ends:
        combine_temp_files_by_suffix(args.output, [p.file_prefix for p in runner], ['.diffsplice.readends.bed'])
        sort_bed_file(args.output + '.diffsplice.readends.bed')
    if args.check_outliers:
        combine_temp_files_by_suffix(args.output, [p.file_prefix for p in runner],
                                     list(OUTLIER_SUFFIXES), headers=OUTLIER_SUFFIXES)

    if not args.keep_intermediate:
        shutil.rmtree(temp_dir)

    genome.close()


if __name__ == "__main__":
    main()
