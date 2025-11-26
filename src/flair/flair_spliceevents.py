#! /usr/bin/env python3

import argparse
import os
import pipettor
import shutil
import pysam
import logging
from flair.ssPrep import buildIntervalTree
import multiprocessing as mp
from flair import FlairInputDataError
import flair.flair_transcriptome as ft
from statistics import median


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
    parser.add_argument('-f', '--gtf', default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms '
                             'to annotated isoforms and adjusting TSS/TESs')

    mutexc = parser.add_mutually_exclusive_group(required=False)
    mutexc.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    mutexc.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from short-read alignment with junctions_from_sam)')
    parser.add_argument('--junction_support', type=int, default=1,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support.')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')

    parser.add_argument('-s', '--support', type=int, default=2,
                        help='''minimum number of supporting reads for an isoform (3)''')

    parser.add_argument('--parallelmode', default='auto:1GB',
                        help='''parallelization mode. Default: "auto:1GB" This indicates an automatic threshold where
                             if the file is less than 1GB, parallelization is done by chromosome, but if it's larger,
                             parallelization is done by region of non-overlapping reads. Other modes: bychrom, byregion,
                             auto:xGB - for setting the auto threshold, it must be in units of GB.''')

    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.
                Intermediate files include: promoter-supported reads file,
                read assignments to firstpass isoforms''')
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

    check_file_paths(args)
    args = ft.add_preset_args(args)
    args.quality = 0 #should only be used for genomic alignment
    return args


def check_file_paths(args):
    if not os.path.exists(args.manifest):
        raise FlairInputDataError(f'Manifest file path does not exist: {args.manifest}')
    if not os.path.exists(args.genome):
        raise FlairInputDataError(f'Genome file path does not exist: {args.genome}')
    if not (args.parallelmode in {'bychrom', 'byregion'}
            or (args.parallelmode[:5] == 'auto:'
                and ((args.parallelmode[-2:] == 'GB' and args.parallelmode[5:-2].replace(".", "").isnumeric())
                     or args.parallelmode[5:].replace(".", "").isnumeric()))):
        raise FlairInputDataError(
            f'parallelmode {args.parallelmode} is not in an allowed format. See docs for allowed formats')


def get_annot_tstart_tend(tinfo):
    tstart, tend, strand = tinfo[0]
    if tstart is None:
        tstart = min([x[0] for x in tinfo[1]])
        tend = max([x[1] for x in tinfo[1]])
    return tstart, tend, strand

def get_annot_gene_hits(annotgenes, exons):
    gene_hits = []
    for annotgene in annotgenes:
        annotexons = sorted(list(annotgenes[annotgene]))
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

def get_juncs_to_gene(juncs, isoinfo, juncstotranscript, junctogene, genetoannotjuncs, allsplicedexons):
    thisgene, thistranscript = None, None
    if juncs in juncstotranscript:
        thistranscript, thisgene = juncstotranscript[juncs]
    else:
        thisgene = None
        gene_hits = ft.get_genes_with_shared_juncs(juncs, junctogene, genetoannotjuncs)
        if gene_hits:
            sortedgenes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)
            thisgene = sortedgenes[0][0]
        else:
            # look for exon overlap
            mystart = max(x[0] for x in isoinfo)
            myend = min(x[1] for x in isoinfo)
            exons = [(mystart, juncs[0][0])] + [(juncs[i][1], juncs[i + 1][0]) for i in range(len(juncs) - 1)] + [
                (juncs[-1][1], myend)]
            strand = isoinfo[0][2]  # not a super robust strand picking, assumes well stranded reads
            gene_hits = get_annot_gene_hits(allsplicedexons[strand], exons)
            if len(gene_hits) > 0:
                gene_hits.sort(reverse=True)
                thisgene = gene_hits[0][1]
    return thisgene, thistranscript

def group_juncs_by_annot_gene(sjtoends, juncstotranscript, junctogene, genetoannotjuncs, allsplicedexons):
    genetojuncs, nogenejuncs, sereads = {}, {}, []
    for juncs in sjtoends:  # assuming unspliced reads already removed
        if len(juncs) > 0:  # remove unspliced reads
            thisgene, thistranscript = get_juncs_to_gene(juncs, sjtoends[juncs], juncstotranscript,
                                                         junctogene, genetoannotjuncs, allsplicedexons)

            if thisgene:
                if thisgene not in genetojuncs:
                    genetojuncs[thisgene] = {}
                genetojuncs[thisgene][juncs] = sjtoends[juncs]
            else:
                nogenejuncs[juncs] = sjtoends[juncs]
        else:
            sereads.extend(sjtoends[juncs])
    return genetojuncs, nogenejuncs, sereads


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


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
    for start, end, rstrand, name in readinfo:
        # getting all block coverage to measure intron retention
        firstexon = (start, juncs[0][0])
        lastexon = (juncs[-1][1], end)
        if firstexon not in allblocks:
            allblocks[firstexon] = {s: 0 for s in allsamples}
        if lastexon not in allblocks:
            allblocks[lastexon] = {s: 0 for s in allsamples}
        allblocks[firstexon][sample] += 1
        allblocks[lastexon][sample] += 1
        if outends:
            outends.write('\t'.join([thischrom, str(start), str(end),
                                     gene + '|' + name, '.', strand]) + '\n')
    return allblocks

def determine_juncs_are_subset(juncs, alljuncs):
    is_subset = False
    unique_seq_bound = [[], []]
    for otherjuncs in alljuncs:
        if juncs != otherjuncs and len(juncs) < len(otherjuncs):
            # print(str(juncs)[1:-1].rstrip(','), str(otherjuncs))
            if str(juncs)[1:-1].rstrip(',') in str(otherjuncs):
                is_subset = True
                other_internal_exons = [(otherjuncs[x][1], otherjuncs[x+1][0]) for x in range(len(otherjuncs)-1)]
                # if juncs[0][0] == otherjuncs[0][0]:
                #     terminal_exon_is_subset[0] = 1
                # if juncs[-1][1] == otherjuncs[-1][1]:
                #     terminal_exon_is_subset[1] = 1
                for other_exon in other_internal_exons:
                    if juncs[0][0] == other_exon[1]:
                        unique_seq_bound[0].append(other_exon[1] - other_exon[0])
                    if juncs[-1][1] == other_exon[0]:
                        unique_seq_bound[1].append(other_exon[1] - other_exon[0])
    if is_subset:
        unique_seq_bound[0] = max(unique_seq_bound[0]) if len(unique_seq_bound[0])>0 else None
        unique_seq_bound[1] = max(unique_seq_bound[1]) if len(unique_seq_bound[1])>0 else None
        return unique_seq_bound
    else:
        return [None, None]

def get_start_end_counts(readinfo, subset_info, sample, first_junc, last_junc, t_starts_ends, t_first_last_sj, allsamples):
    for start, end, rstrand, name in readinfo:
        # print(end, subset_info)
        if subset_info[0] == None or first_junc-start > subset_info[0]:
            t_starts_ends[0].append((start, sample))
            if first_junc not in t_first_last_sj[0]:
                t_first_last_sj[0][first_junc] = {s:0 for s in allsamples}
            t_first_last_sj[0][first_junc][sample] += 1
        if subset_info[1] == None or end-last_junc > subset_info[1]:
            t_starts_ends[1].append((end, sample))
            if last_junc not in t_first_last_sj[1]:
                t_first_last_sj[1][last_junc] = {s:0 for s in allsamples}
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
            last_end=end[0]
        # all_groups.append(curr_group)
        # end_to_counts[median(curr_group)] = len(curr_group)
        g = int(median([x[0] for x in curr_group]))
        end_to_counts[g] = {s: 0 for s in allsamples}
        for e in curr_group:
            end_to_counts[g][e[1]] += 1
        grouped_ends[i] = end_to_counts
    return grouped_ends



def extract_splicing_info(allsamples, allgenetojuncs, gene, strand, thischrom, outends):
    ss5to3, ss3to5, alljuncs, exonjpairs, allblocks = {}, {}, {}, {}, {}
    t_starts_ends, t_first_last_sj = [[],[]], [{},{}]
    ##get valid read ends by excluding subset junction chains (check for whether end is within exon)
        ##only check relative to other reads, not reference?
    ##get list of valid read ends after excluding subset junction chains
    ##cluster within window, take median/mode of window (use code from transcriptome)?

    # process each sample for each gene - only ever storing data for one gene at a time
    for s in range(len(allsamples)):
        sample = allsamples[s]
        genetojuncs = allgenetojuncs[s]
        if gene in genetojuncs:
            for juncs in genetojuncs[gene]:
                readinfo = genetojuncs[gene][juncs]
                numreads = len(readinfo)

                ss5to3, ss3to5, alljuncs = extract_35ss_info(juncs, strand, allsamples, sample, numreads,
                                                             ss5to3, ss3to5, alljuncs)

                allblocks, exonjpairs = extract_exon_usage_info(juncs, allsamples, sample, numreads,
                                                                allblocks, exonjpairs)
                allblocks = extract_end_coverage_info(juncs, readinfo, allsamples, sample, thischrom, gene, strand,
                                                      allblocks, outends)
                subset_info = determine_juncs_are_subset(juncs, genetojuncs[gene])
                t_starts_ends, t_first_last_sj = get_start_end_counts(readinfo, subset_info, sample, juncs[0][0], juncs[-1][1],
                                                                      t_starts_ends, t_first_last_sj, allsamples)


    return ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj


def write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, gene, mycolor, out, outbed, min_read_support):
    esjuncs = set()
    for prevjunc, nextjunc in exonjpairs:
        outerjunc = (prevjunc[0], nextjunc[1])
        if outerjunc in alljuncs:  # if there's any reads that skip this exon
            increads = [exonjpairs[(prevjunc, nextjunc)][s] for s in allsamples]
            excreads = [alljuncs[outerjunc][s] for s in allsamples]
            if any([x >= min_read_support for x in increads]) and any([x >= min_read_support for x in excreads]):

                esname = f'es-of-{thischrom}:{prevjunc[1]}-{nextjunc[0]}' \
                         f'({strand})-at-{outerjunc[0]}-{outerjunc[1]}-{gene}'
                incline = [thischrom, prevjunc[1], nextjunc[0], f'inc_{esname}', str(min(sum(increads), 1000)), strand,
                           prevjunc[1], nextjunc[0], mycolor]
                outbed.write('\t'.join([str(x) for x in incline]) + '\n')
                incline = [f'inc_{esname}'] + increads
                out.write('\t'.join([str(x) for x in incline]) + '\n')
                excline = [thischrom, outerjunc[0], outerjunc[1], f'exc_{esname}', str(min(sum(excreads), 1000)), strand,
                           outerjunc[0], outerjunc[1], mycolor]
                outbed.write('\t'.join([str(x) for x in excline]) + '\n')
                excline = [f'exc_{esname}'] + excreads
                out.write('\t'.join([str(x) for x in excline]) + '\n')

                esjuncs.add(prevjunc)
                esjuncs.add(nextjunc)
    return esjuncs


def process_junction_events(ssAtoB, esjuncs, type, thischrom, strand, gene, allsamples, mycolor, out, outbed, min_read_support):
    for ssA in ssAtoB:
        if len(ssAtoB[ssA]) > 1:  # don't report when no alternative splicing at all
            goodB = []
            for ssB in ssAtoB[ssA]:
                junc = (min((ssA, ssB)), max((ssA, ssB)))
                # TODO should I instead be trying to subtract the exon skipping read support from the alt 5/3 read support?
                sscounts = [ssAtoB[ssA][ssB][s] for s in allsamples]
                if junc not in esjuncs and any([x >= min_read_support for x in sscounts]):
                    goodB.append(ssB)
            if len(goodB) > 1:  # only if alt splicing, not including junctions involved in exon skipping
                ssAname = f'{type}-relTo-{thischrom}:{ssA}({strand})-{gene}'
                for ssB in goodB:
                    sscounts = [ssAtoB[ssA][ssB][s] for s in allsamples]
                    ssBname = f'{ssB}_{ssAname}'
                    outline = [thischrom, min((ssA, ssB)), max((ssA, ssB)), ssBname, str(min(sum(sscounts), 1000)), strand,
                               min((ssA, ssB)), max((ssA, ssB)), mycolor]
                    outbed.write('\t'.join([str(x) for x in outline]) + '\n')
                    outline = [ssBname] + sscounts
                    out.write('\t'.join([str(x) for x in outline]) + '\n')


def write_intron_retention(alljuncs, allsamples, allblocks, thischrom, strand, gene, mycolor, out, outbed, min_read_support):
    for j in alljuncs:
        allretainedcounts = {s: 0 for s in allsamples}
        for b in allblocks:  # this nested for loop could get brutal
            if b[0] < j[0] and j[1] < b[1]:
                for s in allsamples:
                    allretainedcounts[s] += allblocks[b][s]
        splicedcounts = [alljuncs[j][s] for s in allsamples]
        retainedcounts = [allretainedcounts[s] for s in allsamples]
        # if sum(allretainedcounts.values()) > 0:  # only report if junction retained at all
        if any([x >= min_read_support for x in splicedcounts]) and any([x >= min_read_support for x in retainedcounts]):
            irname = f'ir-of-{thischrom}:{j[0]}-{j[1]}({strand})-{gene}'
            incline = [thischrom, j[0], j[1], f'spliced_{irname}', str(min(sum(splicedcounts), 1000)), strand,
                       j[0], j[1], mycolor]

            outbed.write('\t'.join([str(x) for x in incline]) + '\n')
            incline = [f'spliced_{irname}'] + splicedcounts
            out.write('\t'.join([str(x) for x in incline]) + '\n')
            excline = [thischrom, j[0], j[1], f'retained_{irname}', str(min(sum(retainedcounts), 1000)), strand,
                       j[0], j[1], mycolor]

            outbed.write('\t'.join([str(x) for x in excline]) + '\n')
            excline = [f'retained_{irname}'] + retainedcounts
            out.write('\t'.join([str(x) for x in excline]) + '\n')

def write_ends(grouped_ends, allsamples, thischrom, strand, gene, eventname, mycolor, out, outbed, support):
    good_ends = []
    for e in grouped_ends:
        outcounts = [grouped_ends[e][s] for s in allsamples]
        if any([x>=support for x in outcounts]):
            good_ends.append(e)
    if len(good_ends) > 1:
        for e in good_ends:
            outcounts = [grouped_ends[e][s] for s in allsamples]
            name = f'{e}_{eventname}-{gene}'
            bedline = [thischrom,e-1, e, name, str(min(sum(outcounts), 1000)), strand,
                           e-1, e, mycolor]
            outbed.write('\t'.join([str(x) for x in bedline]) + '\n')
            countsline = [name] + outcounts
            out.write('\t'.join([str(x) for x in countsline]) + '\n')



def process_gene_to_events(tempprefix, thischrom, allsamples, allgenetojuncs, allannottranscripts, genetostrand, support, output_read_ends):
    etypetocolor = {'skipped_exons': '66,105,245', 'retained_introns': '144,66,245',
                    'alt3': '245,215,66', 'alt5': '43,184,39', 'tss':'255,0,0', 'tts':'0,0,255'}
    out = open(tempprefix + '.diffsplice.counts.tsv', 'w')
    outbed = open(tempprefix + '.diffsplice.bed', 'w')
    ##outends is just for writing out read ends for Harrison - he can group them more intelligently
    outends = None
    if output_read_ends:
        outends = open(tempprefix + '.diffsplice.readends.bed', 'w')

    allgenes = set.union(*[set(allgenetojuncs[s].keys()) for s in range(len(allgenetojuncs))])
    for gene in allgenes:
        # get gene strand - this is not optimal, def need to revamp how getting annot info
        strand = genetostrand[gene]
        ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, t_starts_ends, t_first_last_sj = extract_splicing_info(allsamples, allgenetojuncs,
                                                                                gene, strand, thischrom, outends)

        # exon skipping
        # TODO combine exon skipping at multiple junctions (so all junctions that skip exon are combined)
        esjuncs = write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, gene,
                                      etypetocolor['skipped_exons'], out, outbed, support)
        process_junction_events(ss5to3, esjuncs, 'alt3', thischrom, strand, gene, allsamples, etypetocolor['alt3'],
                                out, outbed, support)
        process_junction_events(ss3to5, esjuncs, 'alt5', thischrom, strand, gene, allsamples, etypetocolor['alt5'],
                                out, outbed, support)
        # intron retention
        write_intron_retention(alljuncs, allsamples, allblocks, thischrom, strand, gene,
                               etypetocolor['retained_introns'], out, outbed, support)

        grouped_ends = group_ends(t_starts_ends, allsamples, 100)
        if strand == '-':
            grouped_ends = grouped_ends[::-1]
            t_first_last_sj = t_first_last_sj[::-1]
        write_ends(grouped_ends[0], allsamples, thischrom, strand, gene, 'tss', etypetocolor['tss'], out, outbed, support)
        write_ends(grouped_ends[1], allsamples, thischrom, strand, gene, 'tts', etypetocolor['tts'], out, outbed, support)
        # write_ends(t_first_last_sj[0], allsamples, thischrom, strand, gene, 'tss', etypetocolor['tss'], out, outbed, support)
        # write_ends(t_first_last_sj[1], allsamples, thischrom, strand, gene, 'tts', etypetocolor['tts'], out, outbed, support)

    out.close()
    outbed.close()
    if output_read_ends:
        outends.close()


def generate_good_match_to_annot(args, tempprefix, thischrom, annottranscripttoexons, alltranscripts, genome, bamfile, junctogene):
    if not args.noaligntoannot and len(alltranscripts) > 0:
        # FIXME Instead of randomly addding to transcript ends, standardize ends based on gene/splice junction?

        _ = ft.generate_transcriptome_reference(tempprefix, alltranscripts, annottranscripttoexons,
                                                thischrom, genome, junctogene, add_length_at_ends=500, normalize_ends=True)

        clippingfile = tempprefix + '.reads.genomicclipping.txt'
        mm2_cmd = ('minimap2', '-a', '-N', '4', '--MD',
                   tempprefix + '.annotated_transcripts.fa', tempprefix + '.reads.fasta')
        flairpath = '/'.join(os.path.realpath(__file__).split('/')[:-1])

        count_cmd = ('python3', flairpath + '/filter_transcriptome_align.py',
                     '--sam', '-',
                     '-o', tempprefix + '.matchannot.counts.tsv',
                     '-t', 1,  # feeding 1 thread in because this is already multithreaded here
                     '--quality', 0,
                     '--output_bam', tempprefix + '.matchannot.bam',
                     '--stringent',
                     '--check_splice',
                     '-i', tempprefix + '.annotated_transcripts.bed',
                     '--trimmedreads', clippingfile)

        liftovercmd = ('python3', os.path.dirname(os.path.abspath(__file__)) + '/lift_transcriptome_aligned_bam_to_genome.py',
                       tempprefix + '.annotated_transcripts.bed',
                       tempprefix + '.matchannot.bam',
                       bamfile,
                       tempprefix + '.matchannot')

        removefilescmd = ('rm',
                          tempprefix + '.matchannot.counts.tsv',
                          tempprefix + '.matchannot.bam',
                          tempprefix + '.matchannot.genomelift.bam')
        pipettor.run([mm2_cmd, count_cmd])
        pipettor.run([liftovercmd])
        pipettor.run([removefilescmd])
    else:
        out = pysam.AlignmentFile(tempprefix + '.matchannot.genomelift.sorted.bam', 'wb', template=bamfile)
        out.close()
    return tempprefix + '.matchannot.genomelift.sorted.bam'


def runcollapsebychrom(listofargs):
    args, tempprefix, splicesiteannot_chrom, juncstotranscript, junctogene, allannotse, allsplicedexons, \
        genetoannotjuncs, genetostrand, annottranscripttoexons, allannottranscripts = listofargs
    # first extract reads for chrom as fasta
    tempsplit = tempprefix.split('/')[-1].split('-')
    rchrom, rstart, rend = '-'.join(tempsplit[:-2]), tempsplit[-2], tempsplit[-1]

    # align reads to annot [transcripts +- 1000bp], filter to only good aligns, convert to genomic coords
    # get reads from that bam file (no need to correct), save junctions - am actually doing correct, could probably remove that
    # add junctions from other reads that did not match ref transcriptome well after correction

    # load splice junctions for chrom
    intervalTree, junctionBoundaryDict = buildIntervalTree(splicesiteannot_chrom, args.ss_window, rchrom, False)
    genome = pysam.FastaFile(args.genome)

    allgenetojuncs, allsamples = [], []
    for line in open(args.manifest):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            sample, bamfile = line

            pipettor.run([('samtools', 'view', '-h', bamfile, rchrom + ':' + rstart + '-' + rend),
                          ('samtools', 'fasta', '-')],
                         stdout=open(tempprefix + '.reads.fasta', 'w'))

            samfile = pysam.AlignmentFile(bamfile, 'rb')
            ft.generate_genomic_clipping_reference(tempprefix, samfile, rchrom, rstart, rend)
            samfile.close()

            goodannotaligns = generate_good_match_to_annot(args, tempprefix, rchrom, annottranscripttoexons,
                                                           allannottranscripts, genome, bamfile, junctogene)

            samfile = pysam.AlignmentFile(goodannotaligns, 'rb')
            sjtoends, goodaligntoannot = ft.filter_correct_group_reads(args, tempprefix, rchrom, rstart, rend, samfile,
                                                                    set(), intervalTree, junctionBoundaryDict,
                                                                    generate_fasta=False, sj_to_ends={},
                                                                    return_used_reads=True, allow_secondary=True)
            samfile.close()
            samfile = pysam.AlignmentFile(bamfile, 'rb')
            sjtoends = ft.filter_correct_group_reads(args, tempprefix, rchrom, rstart, rend, samfile, goodaligntoannot,
                                                  intervalTree, junctionBoundaryDict,
                                                  generate_fasta=False, sj_to_ends=sjtoends, return_used_reads=False)
            samfile.close()
            genetojuncs, nogenejuncs, sereads = group_juncs_by_annot_gene(sjtoends, juncstotranscript, junctogene,
                                                                          genetoannotjuncs, allsplicedexons)
            allsamples.append(sample)
            allgenetojuncs.append(genetojuncs)

    process_gene_to_events(tempprefix, rchrom, allsamples, allgenetojuncs, allannottranscripts, genetostrand, args.support, args.output_read_ends)


def collapsefrombam():
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)
    allchrom = genome.references
    logging.info('making temp dir')
    tempDir = ft.make_correct_temp_dir()

    allregions = []
    for chrom in allchrom:
        chromsize = genome.get_reference_length(chrom)
        allregions.append((chrom, 0, chromsize))
    logging.info(f'Running by chroms: {len(allregions)}')
    logging.info('Generating splice site database')
    knownchromosomes, annotationFiles = ft.generate_known_SS_database(args, tempDir)

    regionstoannotdata = {}
    if args.gtf:
        logging.info('Extracting annotation from GTF')
        regionstoannotdata = ft.get_annot_info(args.gtf, allregions)
    logging.info('splitting by chunk')

    chunkcmds = []
    tempprefixes = []
    allsamples = []
    for line in open(args.manifest):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            allsamples.append(line[0])
    out = open(tempDir + '0000.header.diffsplice.counts.tsv', 'w')
    out.write('\t'.join(['featureID'] + allsamples) + '\n')
    out.close()

    for rchrom, rstart, rend in allregions:
        if rchrom in knownchromosomes:
            juncstotranscript, junctogene, allannotse, allsplicedexons, genetoannotjuncs, \
                genetostrand, annottranscripttoexons, allannottranscripts = {}, {}, [], {'+':{},'-':{}}, {}, {}, {}, []
            if args.gtf:
                juncstotranscript, junctogene, allannotse, allsplicedexons, genetoannotjuncs, \
                    genetostrand, annottranscripttoexons, allannottranscripts = \
                    regionstoannotdata[(rchrom, rstart, rend)].return_data()

            splicesiteannot_chrom = annotationFiles[rchrom]
            tempprefix = tempDir + '-'.join([rchrom, str(rstart), str(rend)])
            chunkcmds.append([args, tempprefix, splicesiteannot_chrom, juncstotranscript, junctogene,
                              allannotse, allsplicedexons, genetoannotjuncs, genetostrand,
                              annottranscripttoexons, allannottranscripts])
            tempprefixes.append(tempprefix)
    mp.set_start_method('fork')
    p = mp.Pool(args.threads)
    childErrs = set()
    c = 1
    for i in p.imap(runcollapsebychrom, chunkcmds):
        logging.info(f'\rdone running chunk {c} of {len(chunkcmds)}')
        childErrs.add(i)
        c += 1
    p.close()
    p.join()
    if len(childErrs) > 1:
        raise ValueError(childErrs)

    ft.combine_temp_files_by_suffix(args, tempprefixes,
                                ['.diffsplice.bed', '.diffsplice.counts.tsv'])
    if args.output_read_ends:
        ft.combine_temp_files_by_suffix(args, tempprefixes, ['.diffsplice.readends.bed'])

    with open(args.output + '.diffsplice.counts.new.tsv', 'w') as out:
        out.write('\t'.join(['ids'] + allsamples) + '\n')
        for line in open(args.output + '.diffsplice.counts.tsv'):
            out.write(line)
    pipettor.run([('mv', args.output + '.diffsplice.counts.new.tsv', args.output + '.diffsplice.counts.tsv')])


    if not args.keep_intermediate:
        shutil.rmtree(tempDir)

    genome.close()


if __name__ == "__main__":
    collapsefrombam()
