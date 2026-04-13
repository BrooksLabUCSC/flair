#! /usr/bin/env python3

import argparse
import os
import pipettor
import shutil
import pysam
import logging
import scipy.stats as sps
from flair.partition_runner import PartitionRunner
from flair import SeqRange
from statistics import median
from flair.intron_support import IntronSupport
from flair.junction_correct import JunctionCorrector
from flair.isoform_data import (Junc, ReadRec)
from flair.read_processing import (add_corrected_read_to_groups, get_sequence_from_bed)
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import annot_data_from_gtf

# FIXME: this is temp, need to move into a
import flair.flair_transcriptome as ft

def get_args():
    parser = argparse.ArgumentParser(description='identifies counts of different splicing events directly from a '
                                                 'bam file of long-rna-seq alignments')
    parser.add_argument('-m', '--manifest', required=True,
                        help='manifest tsv with sample name and path to bam file aligned to genome')
    parser.add_argument('-g', '--genome', required=True, type=str,
                        help='FastA of reference genome, can be minimap2 indexed')
    parser.add_argument('-o', '--output', default='flair',
                        help='output file name base for FLAIR isoforms (default: flair)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (4)')
    # FIXME: this is different than all other modules that use --gtf.
    parser.add_argument('-f', '--annot', required=True, default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms '
                             'to annotated isoforms and adjusting TSS/TESs')
    # FIXME: what is the difference than the above, same description
    parser.add_argument('--annot_basic', default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms '
                             'to annotated isoforms and adjusting TSS/TESs')

    # FIXME:
    # parser.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
    #                                            'Use this option if you aligned your short-reads with STAR, '
    #                                            'STAR will automatically output this file')
    parser.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from short-read alignment with junctions_from_sam)')
    parser.add_argument('--region_bed',
                        help='bed file with regions to parallelize by; if not specified, all chromosomes are used')
    parser.add_argument('--junction_support', type=int, default=1,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support.')
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

    # parser.add_argument('--parallel_mode', default='auto:1GB',
    #                     help='''parallelization mode. Default: "auto:1GB" This indicates an automatic threshold where
    #                          if the file is less than 1GB, parallelization is done by chromosome, but if it's larger,
    #                          parallelization is done by region of non-overlapping reads. Other modes: bychrom, byregion,
    #                          auto:xGB - for setting the auto threshold, it must be in units of GB.''')

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
            # look for exon overlap
            mystart = max(x.start for x in isoinfo)
            myend = min(x.end for x in isoinfo)
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
            outends.write('\t'.join([thischrom, str(r.start), str(r.end),
                                     gene + '|' + r.name, '.', strand]) + '\n')
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


def extract_splicing_info(allsamples, allgenetojuncs, gene, strand, thischrom, outends):
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
        genetojuncs = allgenetojuncs[s]
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
                    ref_junc = (ssA, min(goodB)) if ssA < ssB else (max(goodB), ssA)
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


def get_psi_and_filter(event_to_info, allsamples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot, event_support, outoutlier, outolfilt):  # noqa: C901 - FIXME: reduce complexity
    sig_events = []
    etype_to_sig = {}
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

                    outinfo = [f'{jname}_{ename}', event.eventtype, event.gene, get_junc_string(event.chrom, jinfo.inc_juncs),
                               get_junc_string(event.chrom, jinfo.exc_juncs), get_junc_string(event.chrom, jinfo.outer_juncs), get_junc_string(event.chrom, jinfo.inc_exon)]
                    juncpsi = write_counts_psi(outinfo, jinfo.samplecounts, event.totjunc, event.totoverlap, allsamples, outcounts, outpsijunc, outpsitot, event_support)

                    if outoutlier is not None:
                        vals_for_outlier = [x for x in juncpsi if x != 'NA']
                        med = median(vals_for_outlier)
                        # FIXME: use BED class
                        for line in jinfo.bedlines:
                            line[4] = round(med * 100)
                            outbed.write('\t'.join([str(x) for x in line]) + '\n')

                        if len(vals_for_outlier) >= 5:
                            dev = sps.iqr(vals_for_outlier) / 2
                            for s in allsamples:
                                thiscount = jinfo.samplecounts[s]
                                thistot = event.totjunc[s]
                                thispsi = thiscount / thistot if thistot > event_support else 'NA'
                                if thispsi != 'NA' and abs(thispsi - med) >= 0.1:  # 10% PSI
                                    delta = abs(thispsi - med)
                                    countsstr = f'{thiscount};{thistot}'
                                    # event_name, gene, median, dev, psi, numsamplesused, thiscount/junctot, absdeltapsi, numberofdev
                                    if (dev == 0 and thispsi != med) or delta > dev * 3:
                                        outline = outinfo[:3] + [s, round(med, 6), round(dev, 6), round(thispsi, 6), len(vals_for_outlier), countsstr, round(delta, 6), 100000]
                                        if dev != 0:
                                            outline[-1] = round(delta / dev, 6)
                                        sig_events.append(outline)
                                        if event.eventtype not in etype_to_sig:
                                            etype_to_sig[event.eventtype] = {}
                                        all_these_juncs = (frozenset(jinfo.inc_juncs), frozenset(jinfo.exc_juncs), frozenset(jinfo.outer_juncs))
                                        if all_these_juncs not in etype_to_sig[event.eventtype]:
                                            etype_to_sig[event.eventtype][all_these_juncs] = {}
                                        etype_to_sig[event.eventtype][all_these_juncs][s] = outline
                if any([event.other[s] > 0 for s in allsamples]):
                    outinfo = [f'other_{ename}', event.eventtype, event.gene, '', '', '', '']
                    write_counts_psi(outinfo, event.other, event.totjunc, event.totoverlap, allsamples, outcounts, outpsijunc, outpsitot, event_support)
    if outoutlier is not None:
        sig_events.sort(reverse=True, key=lambda x: x[::-1])
        for line in sig_events:
            outoutlier.write('\t'.join([str(x) for x in line]) + '\n')

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
                                        this_dev = sample_to_vals[s][-1]
                                        other_dev = seen_es_juncs[j2][s][-1]
                                        if other_dev / this_dev < 0.8:  # deviation is 20% better for this event
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
                                            this_dev = sample_to_vals[s][-1]
                                            other_dev = seen_junctions[jkey][j2][s][-1]
                                            if other_dev / this_dev < 0.8:  # deviation is 20% better for this event
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
                outolfilt.write('\t'.join([str(x) for x in line]) + '\n')


def process_gene_to_events(temp_prefix, thischrom, allsamples, allgenetojuncs, genetostrand, junc_support, output_read_ends, event_frac_of_tot, junc_frac_of_event, event_support, annot_afe_ss, annot_ale_ss, check_outliers):
    etypetocolor = {'skipped_exons': '66,105,245', 'retained_introns': '144,66,245',
                    'alt3': '245,215,66', 'alt5': '43,184,39', 'tss': '255,0,0', 'tts': '0,0,255'}
    # outends is just for writing out read ends for Harrison - he can group them more intelligently
    outends = None
    if output_read_ends:
        outends = open(temp_prefix + '.diffsplice.readends.bed', 'w')

    allgenes = set.union(*[set(allgenetojuncs[s].keys()) for s in range(len(allgenetojuncs))])
    outoutlier, outolfilt = None, None
    if check_outliers:
        outoutlier = open(temp_prefix + '.diffsplice.outliers.tsv', 'w')
        outolfilt = open(temp_prefix + '.diffsplice.outliers.filtered.tsv', 'w')
    with open(temp_prefix + '.diffsplice.bed', 'w') as outbed, open(temp_prefix + '.diffsplice.counts.tsv', 'w') as outcounts, \
            open(temp_prefix + '.diffsplice.PSIjunc.tsv', 'w') as outpsijunc, open(temp_prefix + '.diffsplice.PSItot.tsv', 'w') as outpsitot:

        for gene in allgenes:

            # get gene strand - this is not optimal, def need to revamp how getting annot info
            strand = genetostrand[gene]
            ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj, interval_to_reads, afe, ale, outer_junc_to_exons = extract_splicing_info(allsamples, allgenetojuncs,
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

            get_psi_and_filter(all_info, allsamples, event_frac_of_tot, junc_frac_of_event, outbed, outcounts, outpsijunc, outpsitot, event_support, outoutlier, outolfilt)

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
        mm2_cmd = ('minimap2', '-a', '-N', '4', '--MD',
                   region_annot_fa, temp_prefix + '.reads.fasta')
        flairpath = '/'.join(os.path.realpath(__file__).split('/')[:-1])
        # count_cmd = ('python3', flairpath + '/filter_transcriptome_align.py',
        count_cmd = ('python3', flairpath + '/count_sam_transcripts.py',
                     '--sam', '-',
                     '-o', temp_prefix + '.matchannot.counts.tsv',
                     '-t', 1,  # feeding 1 thread in because this is already multithreaded here
                     '--quality', 0,
                     #  '--generate_map', temp_prefix + '.matchannot.read.map.txt',
                     '--check_splice',
                     '-i', region_annot,
                     '--trimmedreads', clipping_file,
                     '--allow_UTR_indels',
                     '--soft_clipping_buffer', 10,
                     '--output_endpos', temp_prefix + '.readtoends.txt',
                     )

        pipettor.run([mm2_cmd, count_cmd])
        return temp_prefix + '.readtoends.txt'  # temp_prefix + '.matchannot.read.map.txt'
    else:
        return None


def get_juncs_single_sample(listofargs):  # noqa: C901 - FIXME: reduce complexity
    args, region, temp_prefix, sample, bamfile_name, region_annot, region_annot_fa, region_juncs, annots = listofargs

    # FIXME: convert to using PartitionRunner
    intron_support = IntronSupport()
    if region_juncs is not None:
        intron_support.load_introns_bed(region_juncs)
    intron_support.load_annot_bed(region_annot)
    junction_corrector = JunctionCorrector(intron_support, args.ss_window, args.junction_support)

    genome = pysam.FastaFile(args.genome)

    temp_prefix = temp_prefix + '_' + sample

    # print(region.name, region.start, region.end, sample, 'getting annot match')
    bam_file = pysam.AlignmentFile(bamfile_name, 'rb')
    num_reads, clipping_file = ft.generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, region)
    bam_file.close()

    goodannotaligns = generate_good_match_to_annot(args, temp_prefix, region, bamfile_name, region_annot, region_annot_fa, clipping_file)

    read_to_transcript = {}
    for line in open(goodannotaligns):
        line = line.rstrip().split('\t')
        read, transcript = line[:2]
        startindex, startdist, endindex, enddist = [int(x) for x in line[2:]]
        read_to_transcript[read] = (transcript, startindex, startdist, endindex, enddist)

    # print(region.name, region.start, region.end, sample, 'correcting reads')

    sj_to_ends = {}
    bamfile = pysam.AlignmentFile(bamfile_name, 'rb')
    dropped_secondary_sup, dropped_quality, dropped_correction = 0, 0, 0
    annot_match, corrected_cnt = 0, 0
    for read in bamfile.fetch(region.name, region.start, region.end):
        if not read.is_secondary and (not read.is_supplementary or args.keep_sup):
            readrec = ReadRec.from_read(read)
            corrected = False
            if read.query_name in read_to_transcript:
                transcript, startindex, startdist, endindex, enddist = read_to_transcript[read.query_name]
                juncs = annots.transcript_to_sjc[transcript]
                if len(juncs) > 0:
                    newstart = juncs[startindex][0] - startdist
                    newend = juncs[endindex][1] + enddist
                    juncs = tuple([Junc(x[0], x[1]) for x in juncs[startindex:endindex + 1]])
                    strand = annots.gene_to_strand[transcript.split('_')[-1]]
                    readrec.correct_from_annotation(newstart, newend, strand, juncs)
                    corrected = True
                    annot_match += 1
                else:
                    annot_match += 1
            elif read.mapping_quality >= args.quality:
                corrected = junction_corrector.correct_readrec(readrec)
                if corrected:
                    corrected_cnt += 1
                else:
                    dropped_correction += 1
                    logging.debug(f"read dropped: junction correction failed: {read.query_name}")
            else:
                dropped_quality += 1
                logging.debug(f"read dropped: low quality ({read.mapping_quality} < {args.quality}): {read.query_name}")
            if corrected:
                add_corrected_read_to_groups(readrec, sj_to_ends)
        else:
            dropped_secondary_sup += 1
            logging.debug(f"read dropped: secondary or supplementary: {read.query_name}")
    bamfile.close()
    logging.debug(f"spliceevents region {region.name}:{region.start}-{region.end}: "
                  f"annot_match={annot_match}, corrected={corrected_cnt}, "
                  f"dropped_quality={dropped_quality}, dropped_correction={dropped_correction}, "
                  f"dropped_secondary_sup={dropped_secondary_sup}")

    genetojuncs, nogenejuncs, sereads = group_juncs_by_annot_gene(sj_to_ends, annots.sjc_to_gene, annots.junc_to_gene_id, annots.gene_to_exons, annots.gene_to_annot_juncs)

    c = 0
    with open(temp_prefix + '_gene_to_juncs.txt', 'w') as out:
        for gene in genetojuncs:
            for juncs in genetojuncs[gene]:
                juncstring = ','.join(['.'.join([str(y) for y in x]) for x in juncs])
                for read_info in genetojuncs[gene][juncs].reads:
                    c += 1
                    outline = [gene, juncstring, str(read_info.start), str(read_info.end), genetojuncs[gene][juncs].strand, read_info.name]
                    out.write('\t'.join(outline) + '\n')
    # print(sample, c, 'final multi-junction reads')
    pipettor.run([('rm', f'{temp_prefix}.matchannot.counts.tsv', f'{temp_prefix}.readtoends.txt', f'{temp_prefix}.reads.fasta', f'{temp_prefix}.reads.genomicclipping.txt')])
    genome.close()

def process_bed_line(line):
    # FIXME: use BED class
    line = line.rstrip().split('\t')
    transcript = line[3]
    gene = line[3].split('_')[-1]
    start, esizes, estarts = int(line[1]), [int(x) for x in line[-2].rstrip(',').split(',')], \
        [int(x) for x in line[-1].rstrip(',').split(',')]
    exons = [(start + estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
    junctions = [(start + estarts[i] + esizes[i], start + estarts[i + 1]) for i in range(len(esizes) - 1)]
    strand = line[5]
    return gene, transcript, exons, junctions, strand


def _run_region(*, partition, gtf_data, intron_support, args, allsamples):  # noqa: C901 - FIXME: reduce complexity
    region_bed = partition.output_path('region.bed')
    out = open(region_bed, 'w')
    out.write('\t'.join([partition.region.name, str(partition.region.start), str(partition.region.end)]) + '\n')
    out.close()

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
            for line in open(region_annot_basic):
                gene, transcript, exons, junctions, strand = process_bed_line(line)
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

        for sample, bamfile in allsamples:
            if not os.path.exists(partition.file_prefix + '_' + sample + '_gene_to_juncs.txt'):
                get_juncs_single_sample([args, partition.region, partition.file_prefix, sample, bamfile, region_annot, region_annot_fa, region_juncs, annots])

        # p = multiprocessing.Pool(4)
        # childErrs = set()
        # c = 1
        # for i in p.imap(get_juncs_single_sample, chunkcmds):
        #     # logging.info(f'\r{region.name} {region.start} {region.end} done running sample chunk {c} of {len(chunkcmds)}')
        #     childErrs.add(i)
        #     c += 1
        # p.close()
        # p.join()
        # if len(childErrs) > 1:
        #     raise ValueError(childErrs)

        allgenetojuncs = []
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
            allgenetojuncs.append(gene_to_juncs)

        process_gene_to_events(
            partition.file_prefix, partition.region.name, [x[0] for x in allsamples],
            allgenetojuncs, annots.gene_to_strand, args.junc_support, args.output_read_ends,
            args.event_frac_of_tot, args.junc_frac_of_event, args.event_support,
            annot_afe_ss, annot_ale_ss, args.check_outliers)

        if not args.keep_intermediate:
            for sample, bamfile in allsamples:
                pipettor.run([('rm', partition.file_prefix + '_' + sample + '_gene_to_juncs.txt')])


# def preprocess_single_sample(listofargs):
#     args, region, temp_prefix, sample, bamfile_name = listofargs

#     temp_prefix = temp_prefix + '_' + sample

#     bamfile = pysam.AlignmentFile(bamfile_name, 'rb')
#     starts, ends = set(), set()
#     for read in bamfile.fetch(region.name, region.start, region.end):
#         if not read.is_secondary and (not read.is_supplementary or args.keep_sup) and read.mapping_quality >= args.quality:
#             starts.add(read.reference_start)
#             ends.add(read.reference_end)
#     bamfile.close()

#     with open(temp_prefix + '_ends.txt', 'w') as out:
#         if len(starts) > 0:
#             out.write(f'{min(starts)}\t{max(ends)}')
#         else:
#             out.write(f'\t')


# def preprocess_by_region(listofargs):
#     args, region, temp_prefix, allsamples = listofargs
#     chunkcmds = []
#     for sample, bamfile in allsamples:
#         chunkcmds.append([args, region, temp_prefix, sample, bamfile])

#     p = multiprocessing.Pool(4)
#     childErrs = set()
#     for i in p.imap(preprocess_single_sample, chunkcmds):
#         childErrs.add(i)
#     p.close()
#     p.join()
#     if len(childErrs) > 1:
#         raise ValueError(childErrs)


#     starts, ends = set(), set()
#     for sample, bamfile in allsamples:
#         for line in open(temp_prefix + '_' + sample + '_ends.txt'):
#             s, e = line.split('\t')
#             if s != '':
#                 starts.add(int(s))
#                 ends.add(int(e))
#         pipettor.run([('rm', temp_prefix + '_' + sample + '_ends.txt')])

#     with open(temp_prefix + '_ends.txt', 'w') as out:
#         if len(starts) > 0:
#             out.write(f'{min(starts)}\t{max(ends)}')
#         else:
#             out.write(f'\t')


# class NoDaemonProcess(multiprocessing.Process):
#     # make 'daemon' attribute always return False
#     @property
#     def daemon(self):
#         return False

#     @daemon.setter
#     def daemon(self, val):
#         pass

# # We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# # because the latter is only a wrapper function, not a proper class.
# class NoDaemonProcessPool(multiprocessing.pool.Pool):
#     def Process(self, *args, **kwds):
#         proc = super(NoDaemonProcessPool, self).Process(*args, **kwds)
#         proc.__class__ = NoDaemonProcess

#         return proc

def combine_regions(regions, buffersize=0):
    regions.sort()
    new_regions = []
    lastchrom, laststart, lastend = -1, -1, -1
    for range in regions:
        c, s, e = range.name, range.start, range.end
        if c != lastchrom or s > lastend + buffersize:
            if lastchrom != -1:
                new_regions.append(ft.SeqRange(lastchrom, laststart, lastend))
            lastchrom, laststart, lastend = c, s, e
        else:
            lastend = max((lastend, e))
    if lastchrom != -1:
        new_regions.append(ft.SeqRange(lastchrom, laststart, lastend))
    return new_regions


def main():  # noqa: C901 - FIXME: reduce complexity
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)
    logging.info('making temp dir')

    tempDir = ft.make_temp_dir(args.output)
    print('temp directory:', tempDir)

    if args.region_bed:
        all_regions = []
        for line in open(args.region_bed):
            # FIXME: use BED class
            line = line.rstrip().split('\t')
            all_regions.append(SeqRange(line[0], int(line[1]), int(line[2])))
        all_regions = combine_regions(all_regions)
    else:
        all_regions = [SeqRange(chrom, 0, genome.get_reference_length(chrom))
                       for chrom in genome.references]
    logging.info(f'input regions: {len(all_regions)}')

    allsamples = []
    for line in open(args.manifest):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            sample, bamfile = line
            allsamples.append((sample, bamfile))

    out = open(tempDir + '0000.header.diffsplice.counts.tsv', 'w')
    out.write('\t'.join(['featureID'] + [x[0] for x in allsamples]) + '\n')
    out.close()

    logging.info('pre-processing annotation')
    logging.info('loading annotation GTF')
    annot_gtf_data = gtf_data_parser(args.annot, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)

    annot_bed = tempDir + '/annotation.bed'
    if not os.path.exists(annot_bed):
        pipettor.run([('gtf_to_bed', args.annot, annot_bed, '--include_gene')])
    args.annot = annot_bed
    if args.annot_basic:
        annot_basic_bed = tempDir + '/annotation.basic.bed'
        if not os.path.exists(annot_basic_bed):
            pipettor.run([('gtf_to_bed', args.annot_basic, annot_basic_bed, '--include_gene')])
        args.annot_basic = annot_basic_bed

    logging.info(f'running regions with {args.threads} threads')

    runner = PartitionRunner(all_regions, tempDir, gtf_data=annot_gtf_data, threads=args.threads)
    runner.run(_run_region, args=args, allsamples=allsamples)

    ft.combine_temp_files_by_suffix(args.output, [p.file_prefix for p in runner],
                                    ['.diffsplice.bed', '.diffsplice.counts.tsv', '.diffsplice.PSIjunc.tsv', '.diffsplice.PSItot.tsv'])
    if args.output_read_ends:
        ft.combine_temp_files_by_suffix(args.output, [p.file_prefix for p in runner], ['.diffsplice.readends.bed'])
    if args.check_outliers:
        ft.combine_temp_files_by_suffix(args.output, [p.file_prefix for p in runner], ['.diffsplice.outliers.tsv', '.diffsplice.outliers.filtered.tsv'])

    counts_header = ['eventname', 'eventtype', 'gene', 'junctions_included', 'junctions_excluded', 'outer_junctions', 'exons']
    for suffix in ['.diffsplice.counts', '.diffsplice.PSIjunc', '.diffsplice.PSItot']:
        with open(args.output + suffix + '.new.tsv', 'w') as out:
            out.write('\t'.join(counts_header + [x[0] for x in allsamples]) + '\n')
            for line in open(args.output + suffix + '.tsv'):
                out.write(line)
        # FIXME: use os.rename
        pipettor.run([('mv', args.output + suffix + '.new.tsv', args.output + suffix + '.tsv')])

    if args.check_outliers:
        outlier_header = ['eventname', 'eventtype', 'gene', 'sample', 'medianPSI', 'dev(IQR/2)', 'sample_val', 'tot_not_NA_samples', 'event_reads;total_locus_reads', 'delta_PSI_to_med', 'dev_from_med']
        for suffix in ['.diffsplice.outliers', '.diffsplice.outliers.filtered']:
            with open(args.output + suffix + '.new.tsv', 'w') as out:
                out.write('\t'.join(outlier_header) + '\n')
                for line in open(args.output + suffix + '.tsv'):
                    out.write(line)
            pipettor.run([('mv', args.output + suffix + '.new.tsv', args.output + suffix + '.tsv')])

    if not args.keep_intermediate:
        shutil.rmtree(tempDir)

    genome.close()


if __name__ == "__main__":
    main()
