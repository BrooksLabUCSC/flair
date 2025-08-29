#! /usr/bin/env python3

import sys
import argparse
import os
import pipettor
import glob
import uuid
import shutil
import pysam
import logging
from flair.flair_align import inferMM2JuncStrand, intronChainToestarts
from flair.ssUtils import addOtherJuncs, gtfToSSBed
from flair.ssPrep import buildIntervalTree, ssCorrect
import multiprocessing as mp
import time
from collections import Counter
from flair import FlairInputDataError
from flair.flair_transcriptome import *


# export PATH="/private/groups/brookslab/cafelton/git-flair/flair/bin:/private/groups/brookslab/cafelton/git-flair/flair/src/flair:$PATH"

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--manifest',
                        help='manifest tsv with sample name and path to bam file aligned to genome')
    parser.add_argument('-g', '--genome', type=str,
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

    parser.add_argument('-s', '--support', type=float, default=3.0,
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
    parser.add_argument('--noaligntoannot', default=False, action='store_true',
                        help='''specify if you don't want 
                            an initial alignment to the annotated sequences and only want transcript 
                            detection from the genomic alignment.
                             Will be slightly faster but less accurate if the annotation is good''')

    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        parser.error('No arguments passed. Please provide bam file, genome, and gtf')

    args, unknown = parser.parse_known_args()
    if unknown:
        logging.info(f'unrecognized arguments: {" ".join(unknown)}')

    check_file_paths(args)
    args = add_preset_args(args)
    return args


def check_file_paths(args):
    if not args.manifest:
        raise FlairInputDataError(f'Please include the --manifest option')
    if not args.genome:
        raise FlairInputDataError(f'Please include the --genome option\n')
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
    if tstart == None:
        tstart = min([x[0] for x in tinfo[1]])
        tend = max([x[1] for x in tinfo[1]])
    return tstart, tend, strand



# def filtercorrectgroupreads(args, tempprefix, rchrom, rstart, rend, samfile, goodaligntoannot, intervalTree,
#                             junctionBoundaryDict):
#     sjtoends = {}
#     c = 0
#     for read in samfile.fetch(rchrom, int(rstart), int(rend)):
#         if not read.is_secondary and (not read.is_supplementary or args.keep_sup):
#             c += 1
#             if read.mapping_quality >= args.quality:
#                 juncstrand = inferMM2JuncStrand(read)
#                 bedread = BedRead()
#                 bedread.generate_from_cigar(read.reference_start, read.is_reverse, read.cigartuples,
#                                             read.query_name,
#                                             read.reference_name, read.mapping_quality, juncstrand)
#                 correctedread = correctsingleread(bedread, intervalTree, junctionBoundaryDict)
#                 if correctedread:
#                     junckey = tuple(sorted(correctedread.juncs))
#                     if len(junckey) > 0: ##FILTER OUT UNSPLICED READS
#                         if junckey not in sjtoends:
#                             sjtoends[junckey] = []
#                         sjtoends[junckey].append((correctedread.start, correctedread.end, correctedread.strand, correctedread.name))
#     return sjtoends




def group_juncs_by_annot_gene(sjtoends, juncstotranscript, junctogene, genetoannotjuncs):
    genetojuncs, nogenejuncs, sereads = {}, {}, []
    for juncs in sjtoends: ##assuming unspliced reads already removed
        if len(juncs) > 0: ##remove unspliced reads
            if juncs in juncstotranscript:
                thistranscript, thisgene = juncstotranscript[juncs]
            else:
                thisgene = None
                gene_hits = getsplicedisogenehits(juncs, junctogene, genetoannotjuncs)
                if gene_hits:
                    sortedgenes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)
                    thisgene = sortedgenes[0][0]
            if thisgene:
                if thisgene not in genetojuncs: genetojuncs[thisgene] = {}
                genetojuncs[thisgene][juncs] = sjtoends[juncs]
            else:
                nogenejuncs[juncs] = sjtoends[juncs]
        else:
            sereads.extend(sjtoends[juncs])
    return genetojuncs, nogenejuncs, sereads


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def process_gene_to_events(tempprefix, thischrom, allsamples, allgenetojuncs, allannottranscripts, genetostrand):
    etypetocolor = {'skipped_exons':'66,105,245', 'retained_introns':'144,66,245', 'alt3':'245,215,66', 'alt5':'43,184,39'}
    out = open(tempprefix + '.diffsplice.counts.tsv', 'w')
    outbed = open(tempprefix + '.diffsplice.bed', 'w')
    outends = open(tempprefix + '.diffsplice.readends.bed', 'w')

    allgenes = set.union(*[set(allgenetojuncs[s].keys()) for s in range(len(allgenetojuncs))])
    for gene in allgenes:
        ##get gene strand - this is not optimal, def need to revamp how getting annot info
        strand = genetostrand[gene]
        ss5to3, ss3to5, alljuncs, exonjpairs, allblocks, startpos, endpos = {}, {}, {}, {}, {}, {}, {}
        ##process each sample for each gene - only ever storing data for one gene at a time
        for s in range(len(allsamples)):
            sample = allsamples[s]
            genetojuncs = allgenetojuncs[s]
            if gene in genetojuncs:
                for juncs in genetojuncs[gene]:
                    readinfo = genetojuncs[gene][juncs]
                    numreads = len(readinfo)
                    for j in juncs:
                        if strand == '+': ss5, ss3 = j[0], j[1]
                        else: ss5, ss3 = j[1], j[0]
                        if ss5 not in ss5to3: ss5to3[ss5] = {}
                        if ss3 not in ss3to5: ss3to5[ss3] = {}
                        if ss3 not in ss5to3[ss5]: ss5to3[ss5][ss3] = {s:0 for s in allsamples}
                        if ss5 not in ss3to5[ss3]: ss3to5[ss3][ss5] = {s:0 for s in allsamples}
                        ss5to3[ss5][ss3][sample] += numreads
                        ss3to5[ss3][ss5][sample] += numreads
                        if j not in alljuncs: alljuncs[j] = {s:0 for s in allsamples}
                        alljuncs[j][sample] += numreads

                    ##just for middle exons
                    for i in range(len(juncs) - 1):
                        exon = (juncs[i][1], juncs[i + 1][0])
                        prevjunc, nextjunc = juncs[i], juncs[i+1]
                        if exon not in allblocks: allblocks[exon] = {s:0 for s in allsamples}
                        allblocks[exon][sample] += numreads
                        if (prevjunc, nextjunc) not in exonjpairs: exonjpairs[(prevjunc, nextjunc)] = {s:0 for s in allsamples}
                        exonjpairs[(prevjunc, nextjunc)][sample] += numreads


                    for start, end, rstrand, name in readinfo:
                        ##getting all block coverage to measure intron retention
                        firstexon = (start, juncs[0][0])
                        lastexon = (juncs[-1][1], end)
                        if firstexon not in allblocks: allblocks[firstexon] = {s:0 for s in allsamples}
                        if lastexon not in allblocks: allblocks[lastexon] = {s:0 for s in allsamples}
                        allblocks[firstexon][sample] += 1
                        allblocks[lastexon][sample] += 1
                        outends.write('\t'.join([thischrom, str(start), str(end), gene + '|' + name, '.', strand]) + '\n')
                        # ##for alt TSS, TTS
                        # rstart, rend = round(start, -2), round(end, -2)
                        # if strand == '-': rstart, rend = rend, rstart
                        # if rstart not in startpos: startpos[rstart] = {s:0 for s in allsamples}
                        # if rend not in endpos: endpos[rend] = {s:0 for s in allsamples}
                        # startpos[rstart][sample] += 1
                        # endpos[rend][sample] += 1

        # if len(startpos) > 1:
        #     name = f'TSS-{gene}-{thischrom}({strand})'
        #     for pos in startpos:
        #         outline = [thischrom, pos-1, pos, f'{pos}_{name}',
        #                    '.', strand, pos-1, pos, 'orange']
        #         outbed.write('\t'.join([str(x) for x in outline]) + '\n')
        #         outline = [f'{pos}_{name}'] + [startpos[pos][s] for s in allsamples]
        #         out.write('\t'.join([str(x) for x in outline]) + '\n')
        # if len(endpos) > 1:
        #     name = f'TTS-{gene}-{thischrom}({strand})'
        #     for pos in endpos:
        #         outline = [thischrom, pos-1, pos, f'{pos}_{name}',
        #                    '.', strand, pos-1, pos, 'pink']
        #         outbed.write('\t'.join([str(x) for x in outline]) + '\n')
        #         outline = [f'{pos}_{name}'] + [endpos[pos][s] for s in allsamples]
        #         out.write('\t'.join([str(x) for x in outline]) + '\n')

        ##exon skipping
        ## TODO combine exon skipping at multiple junctions (so all junctions that skip exon are combined)
        esjuncs = set()
        for prevjunc, nextjunc in exonjpairs:
            outerjunc = (prevjunc[0], nextjunc[1])
            if outerjunc in alljuncs:  ##if there's any reads that skip this exon
                increads = [exonjpairs[(prevjunc, nextjunc)][s] for s in allsamples]
                excreads = [alljuncs[outerjunc][s] for s in allsamples]
                esname = f'es-of-{thischrom}:{prevjunc[1]}-{nextjunc[0]}({strand})-at-{outerjunc[0]}-{outerjunc[1]}-{gene}'
                incline = [thischrom, prevjunc[1], nextjunc[0], f'inc_{esname}', '.', strand,
                           prevjunc[1], nextjunc[0], etypetocolor['skipped_exons']]
                outbed.write('\t'.join([str(x) for x in incline]) + '\n')
                incline = [f'inc_{esname}'] + increads
                out.write('\t'.join([str(x) for x in incline]) + '\n')
                excline = [thischrom, outerjunc[0], outerjunc[1], f'exc_{esname}', '.', strand,
                           outerjunc[0], outerjunc[1], etypetocolor['skipped_exons']]
                outbed.write('\t'.join([str(x) for x in excline]) + '\n')
                excline = [f'exc_{esname}'] + excreads
                out.write('\t'.join([str(x) for x in excline]) + '\n')

                esjuncs.add(prevjunc)
                esjuncs.add(nextjunc)

        for ss5 in ss5to3:
            if len(ss5to3[ss5]) > 1: ##don't report when no alternative splicing at all
                good3 = []
                for ss3 in ss5to3[ss5]:
                    junc = (min((ss3, ss5)), max((ss3, ss5)))
                    if junc not in esjuncs: ## TODO should I instead be trying to subtract the exon skipping read support from the alt 5/3 read support?
                        good3.append(ss3)
                if len(good3) > 1: ##only if alt splicing, not including junctions involved in exon skipping
                    ss5name = f'alt3-relTo-{thischrom}:{ss5}({strand})-{gene}'
                    for ss3 in good3:
                        sscounts = [ss5to3[ss5][ss3][s] for s in allsamples]
                        ss3name = f'{ss3}_{ss5name}'
                        outline = [thischrom, min((ss3, ss5)), max((ss3, ss5)), ss3name, '.', strand, min((ss3, ss5)), max((ss3, ss5)), etypetocolor['alt3']]
                        outbed.write('\t'.join([str(x) for x in outline]) + '\n')
                        outline = [ss3name] + sscounts
                        out.write('\t'.join([str(x) for x in outline]) + '\n')
        for ss3 in ss3to5:
            if len(ss3to5[ss3]) > 1: ##don't report when no alternative splicing at all
                good5 = []
                for ss5 in ss3to5[ss3]:
                    junc = (min((ss3, ss5)), max((ss3, ss5)))
                    if junc not in esjuncs:
                        good5.append(ss5)
                if len(good5) > 1:  ##only if alt splicing, not including junctions involved in exon skipping
                    ss3name = f'alt5-relTo-{thischrom}:{ss3}({strand})-{gene}'
                    for ss5 in good5:
                        sscounts = [ss3to5[ss3][ss5][s] for s in allsamples]
                        ss5name = f'{ss5}_{ss3name}'
                        outline = [thischrom, min((ss3, ss5)), max((ss3, ss5)), ss5name, '.', strand, min((ss3, ss5)), max((ss3, ss5)), etypetocolor['alt5']]
                        outbed.write('\t'.join([str(x) for x in outline]) + '\n')
                        outline = [ss5name] + sscounts
                        out.write('\t'.join([str(x) for x in outline]) + '\n')

        ##intron retention
        for j in alljuncs:
            allretainedcounts = {s:0 for s in allsamples}
            for b in allblocks: ##this nested for loop could get brutal
                if b[0] < j[0] and j[1] < b[1]:
                    for s in allsamples:
                        allretainedcounts[s] += allblocks[b][s]
            if sum(allretainedcounts.values()) > 0: ##only report if junction retained at all
                splicedcounts = [alljuncs[j][s] for s in allsamples]
                irname = f'ir-of-{thischrom}:{j[0]}-{j[1]}({strand})-{gene}'
                incline = [thischrom, j[0], j[1], f'spliced_{irname}', '.', strand,
                           j[0], j[1], etypetocolor['retained_introns']]
                outbed.write('\t'.join([str(x) for x in incline]) + '\n')
                incline = [f'spliced_{irname}'] + splicedcounts
                out.write('\t'.join([str(x) for x in incline]) + '\n')
                excline = [thischrom, j[0], j[1], f'retained_{irname}', '.', strand,
                           j[0], j[1], etypetocolor['retained_introns']]
                outbed.write('\t'.join([str(x) for x in excline]) + '\n')
                excline = [f'retained_{irname}'] + [allretainedcounts[s] for s in allsamples]
                out.write('\t'.join([str(x) for x in excline]) + '\n')
    out.close()
    outbed.close()
    outends.close()


def generate_good_match_to_annot(args, tempprefix, thischrom, annottranscripttoexons, alltranscripts, genome, bamfile):
    if not args.noaligntoannot and len(alltranscripts) > 0:
        ##FIXME Instead of randomly addding to transcript ends, standardize ends based on gene/splice junction?

        transcripttostrand = generate_transcriptome_reference(tempprefix, alltranscripts, annottranscripttoexons, thischrom, genome, addseqatends=500, normalizeends=True)
        clippingfile = tempprefix + '.reads.genomicclipping.txt'
        mm2_cmd = ('minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD', tempprefix + '.annotated_transcripts.fa', tempprefix + '.reads.fasta')
        flairpath = '/'.join(os.path.realpath(__file__).split('/')[:-1])
        count_cmd = ('python3', flairpath + '/filter_transcriptome_align.py', '--sam', '-', '-o', tempprefix + '.matchannot.counts.tsv', '-t', 1,  ###feeding 1 thread in because this is already multithreaded here
                     '--quality', 0, '--output_bam', tempprefix + '.matchannot.bam', '--stringent', '--check_splice', '-i', tempprefix + '.annotated_transcripts.bed', '--trimmedreads', clippingfile)
        liftovercmd = ('python3', '/private/groups/brookslab/cafelton/fusions-code/lift_transcriptome_aligned_bam_to_genome.py', tempprefix + '.annotated_transcripts.bed',
                       tempprefix + '.matchannot.bam', bamfile, tempprefix + '.matchannot')
        removefilescmd = ('rm', tempprefix + '.matchannot.counts.tsv', tempprefix + '.matchannot.bam', tempprefix + '.matchannot.genomelift.bam')
        pipettor.run([mm2_cmd, count_cmd])
        pipettor.run([liftovercmd])
        pipettor.run([removefilescmd])
    else:
        out = pysam.AlignmentFile(tempprefix + '.matchannot.genomelift.sorted.bam', 'wb', template=bamfile)
        out.close()
    return tempprefix + '.matchannot.genomelift.sorted.bam'







def runcollapsebychrom(listofargs):
    args, tempprefix, splicesiteannot_chrom, juncstotranscript, junctogene, allannotse, genetoannotjuncs, genetostrand, annottranscripttoexons, allannottranscripts = listofargs
    # first extract reads for chrom as fasta
    tempsplit = tempprefix.split('/')[-1].split('-')
    rchrom, rstart, rend = '-'.join(tempsplit[:-2]), tempsplit[-2], tempsplit[-1]


    ##align reads to annot [transcripts +- 1000bp], filter to only good aligns, convert to genomic coords
    ##get reads from that bam file (no need to correct), save junctions - am actually doing correct, could probably remove that
    ##add junctions from other reads that did not match ref transcriptome well after correction

    # load splice junctions for chrom
    intervalTree, junctionBoundaryDict = buildIntervalTree(splicesiteannot_chrom, args.ss_window, rchrom, False)
    genome = pysam.FastaFile(args.genome)

    allgenetojuncs, allsamples = [], []
    for line in open(args.manifest):
        line = line.rstrip().split('\t')
        sample, bamfile = line

        pipettor.run([('samtools', 'view', '-h', bamfile, rchrom + ':' + rstart + '-' + rend),
                      ('samtools', 'fasta', '-')],
                     stdout=open(tempprefix + '.reads.fasta', 'w'))

        samfile = pysam.AlignmentFile(bamfile, 'rb')
        generate_genomic_clipping_reference(tempprefix, samfile, rchrom, rstart, rend)
        samfile.close()

        goodannotaligns = generate_good_match_to_annot(args, tempprefix, rchrom, annottranscripttoexons,
                                                       allannottranscripts, genome, bamfile)

        samfile = pysam.AlignmentFile(goodannotaligns, 'rb')
        sjtoends, goodaligntoannot = filtercorrectgroupreads(args, tempprefix, rchrom, rstart, rend, samfile, set(), intervalTree,
                                           junctionBoundaryDict, generatefasta=False, sjtoends={}, returnusedreads=True, allowsecondary=True)
        samfile.close()
        samfile = pysam.AlignmentFile(bamfile, 'rb')
        sjtoends = filtercorrectgroupreads(args, tempprefix, rchrom, rstart, rend, samfile, goodaligntoannot, intervalTree,
                                           junctionBoundaryDict, generatefasta=False, sjtoends=sjtoends, returnusedreads=False)
        samfile.close()
        genetojuncs, nogenejuncs, sereads = group_juncs_by_annot_gene(sjtoends, juncstotranscript, junctogene, genetoannotjuncs)
        allsamples.append(sample)
        allgenetojuncs.append(genetojuncs)

    process_gene_to_events(tempprefix, rchrom, allsamples, allgenetojuncs, allannottranscripts, genetostrand)




def collapsefrombam():
    logging.basicConfig(level=logging.INFO)
    args = get_args()
    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)
    allchrom = genome.references
    logging.info('making temp dir')
    tempDir = makecorrecttempdir()

    allregions = []
    for chrom in allchrom:
        chromsize = genome.get_reference_length(chrom)
        allregions.append((chrom, 0, chromsize))
    logging.info(f'Running by chroms: {len(allregions)}')
    logging.info('Generating splice site database')
    knownchromosomes, annotationFiles = generateKnownSSDatabase(args, tempDir)

    regionstoannotdata = {}
    if args.gtf:
        logging.info('Extracting annotation from GTF')
        regionstoannotdata = getannotinfo(args.gtf, allregions)
    logging.info('splitting by chunk')

    chunkcmds = []
    tempprefixes = []
    allsamples = []
    for line in open(args.manifest):
        line = line.rstrip().split('\t')
        allsamples.append(line[0])
    out = open(tempDir + '0000.header.diffsplice.counts.tsv', 'w')
    out.write('\t'.join(['featureID'] + allsamples) + '\n')
    out.close()

    for rchrom, rstart, rend in allregions:
        if rchrom in knownchromosomes:
            juncstotranscript, junctogene, allannotse, genetoannotjuncs, annottranscripttoexons, allannottranscripts = {}, {}, [], {}, {}, []
            if args.gtf:
                juncstotranscript, junctogene, allannotse, genetoannotjuncs, genetostrand, annottranscripttoexons, allannottranscripts = \
                                 regionstoannotdata[(rchrom, rstart, rend)].returndata()

            splicesiteannot_chrom = annotationFiles[rchrom]
            tempprefix = tempDir + '-'.join([rchrom, str(rstart), str(rend)])
            chunkcmds.append([args, tempprefix, splicesiteannot_chrom, juncstotranscript,
                              junctogene, allannotse, genetoannotjuncs, genetostrand, annottranscripttoexons,
                              allannottranscripts])
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
        raise Error(childErrs)

    combinetempfilesbysuffix(args, tempprefixes, ['.diffsplice.bed', '.diffsplice.counts.tsv', '.diffsplice.readends.bed'])

    if not args.keep_intermediate:
        shutil.rmtree(tempDir)


    genome.close()


if __name__ == "__main__":
    collapsefrombam()
