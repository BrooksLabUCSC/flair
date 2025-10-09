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

    args = parser.parse_args()

    check_file_paths(args)
    args = ft.add_preset_args(args)
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
        outends.write('\t'.join([thischrom, str(start), str(end),
                                 gene + '|' + name, '.', strand]) + '\n')
    return allblocks


def extract_splicing_info(allsamples, allgenetojuncs, gene, strand, thischrom, outends):
    ss5to3, ss3to5, alljuncs, exonjpairs, allblocks = {}, {}, {}, {}, {}
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

    return ss5to3, ss3to5, alljuncs, exonjpairs, allblocks


def write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, gene, mycolor, out, outbed):
    esjuncs = set()
    for prevjunc, nextjunc in exonjpairs:
        outerjunc = (prevjunc[0], nextjunc[1])
        if outerjunc in alljuncs:  # if there's any reads that skip this exon
            increads = [exonjpairs[(prevjunc, nextjunc)][s] for s in allsamples]
            excreads = [alljuncs[outerjunc][s] for s in allsamples]
            esname = f'es-of-{thischrom}:{prevjunc[1]}-{nextjunc[0]}' \
                     f'({strand})-at-{outerjunc[0]}-{outerjunc[1]}-{gene}'
            incline = [thischrom, prevjunc[1], nextjunc[0], f'inc_{esname}', '.', strand,
                       prevjunc[1], nextjunc[0], mycolor]
            outbed.write('\t'.join([str(x) for x in incline]) + '\n')
            incline = [f'inc_{esname}'] + increads
            out.write('\t'.join([str(x) for x in incline]) + '\n')
            excline = [thischrom, outerjunc[0], outerjunc[1], f'exc_{esname}', '.', strand,
                       outerjunc[0], outerjunc[1], mycolor]
            outbed.write('\t'.join([str(x) for x in excline]) + '\n')
            excline = [f'exc_{esname}'] + excreads
            out.write('\t'.join([str(x) for x in excline]) + '\n')

            esjuncs.add(prevjunc)
            esjuncs.add(nextjunc)
    return esjuncs


def process_junction_events(ssAtoB, esjuncs, type, thischrom, strand, gene, allsamples, mycolor, out, outbed):
    for ssA in ssAtoB:
        if len(ssAtoB[ssA]) > 1:  # don't report when no alternative splicing at all
            goodB = []
            for ssB in ssAtoB[ssA]:
                junc = (min((ssA, ssB)), max((ssA, ssB)))
                # TODO should I instead be trying to subtract the exon skipping read support from the alt 5/3 read support?
                if junc not in esjuncs:
                    goodB.append(ssB)
            if len(goodB) > 1:  # only if alt splicing, not including junctions involved in exon skipping
                ssAname = f'{type}-relTo-{thischrom}:{ssA}({strand})-{gene}'
                for ssB in goodB:
                    sscounts = [ssAtoB[ssA][ssB][s] for s in allsamples]
                    ssBname = f'{ssB}_{ssAname}'
                    outline = [thischrom, min((ssA, ssB)), max((ssA, ssB)), ssBname, '.', strand,
                               min((ssA, ssB)), max((ssA, ssB)), mycolor]
                    outbed.write('\t'.join([str(x) for x in outline]) + '\n')
                    outline = [ssBname] + sscounts
                    out.write('\t'.join([str(x) for x in outline]) + '\n')


def write_intron_retention(alljuncs, allsamples, allblocks, thischrom, strand, gene, mycolor, out, outbed):
    for j in alljuncs:
        allretainedcounts = {s: 0 for s in allsamples}
        for b in allblocks:  # this nested for loop could get brutal
            if b[0] < j[0] and j[1] < b[1]:
                for s in allsamples:
                    allretainedcounts[s] += allblocks[b][s]
        if sum(allretainedcounts.values()) > 0:  # only report if junction retained at all
            splicedcounts = [alljuncs[j][s] for s in allsamples]
            irname = f'ir-of-{thischrom}:{j[0]}-{j[1]}({strand})-{gene}'
            incline = [thischrom, j[0], j[1], f'spliced_{irname}', '.', strand,
                       j[0], j[1], mycolor]

            outbed.write('\t'.join([str(x) for x in incline]) + '\n')
            incline = [f'spliced_{irname}'] + splicedcounts
            out.write('\t'.join([str(x) for x in incline]) + '\n')
            excline = [thischrom, j[0], j[1], f'retained_{irname}', '.', strand,
                       j[0], j[1], mycolor]

            outbed.write('\t'.join([str(x) for x in excline]) + '\n')
            excline = [f'retained_{irname}'] + [allretainedcounts[s] for s in allsamples]
            out.write('\t'.join([str(x) for x in excline]) + '\n')


def process_gene_to_events(tempprefix, thischrom, allsamples, allgenetojuncs, allannottranscripts, genetostrand):
    etypetocolor = {'skipped_exons': '66,105,245', 'retained_introns': '144,66,245',
                    'alt3': '245,215,66', 'alt5': '43,184,39'}
    out = open(tempprefix + '.diffsplice.counts.tsv', 'w')
    outbed = open(tempprefix + '.diffsplice.bed', 'w')
    outends = open(tempprefix + '.diffsplice.readends.bed', 'w')

    allgenes = set.union(*[set(allgenetojuncs[s].keys()) for s in range(len(allgenetojuncs))])
    for gene in allgenes:
        # get gene strand - this is not optimal, def need to revamp how getting annot info
        strand = genetostrand[gene]
        ss5to3, ss3to5, alljuncs, exonjpairs, allblocks = extract_splicing_info(allsamples, allgenetojuncs,
                                                                                gene, strand, thischrom, outends)

        # exon skipping
        # TODO combine exon skipping at multiple junctions (so all junctions that skip exon are combined)
        esjuncs = write_exon_skipping(exonjpairs, alljuncs, allsamples, thischrom, strand, gene,
                                      etypetocolor['skipped_exons'], out, outbed)
        process_junction_events(ss5to3, esjuncs, 'alt3', thischrom, strand, gene, allsamples, etypetocolor['alt3'],
                                out, outbed)
        process_junction_events(ss3to5, esjuncs, 'alt5', thischrom, strand, gene, allsamples, etypetocolor['alt5'],
                                out, outbed)
        # intron retention
        write_intron_retention(alljuncs, allsamples, allblocks, thischrom, strand, gene,
                               etypetocolor['retained_introns'], out, outbed)

    out.close()
    outbed.close()
    outends.close()


def generate_good_match_to_annot(args, tempprefix, thischrom, annottranscripttoexons, alltranscripts, genome, bamfile):
    if not args.noaligntoannot and len(alltranscripts) > 0:
        # FIXME Instead of randomly addding to transcript ends, standardize ends based on gene/splice junction?

        _ = ft.generate_transcriptome_reference(tempprefix, alltranscripts, annottranscripttoexons,
                                                thischrom, genome, addseqatends=500, normalizeends=True)

        clippingfile = tempprefix + '.reads.genomicclipping.txt'
        mm2_cmd = ('minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD',
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

        liftovercmd = ('python3',
                       '/private/groups/brookslab/cafelton/fusions-code/lift_transcriptome_aligned_bam_to_genome.py',
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
        line = line.rstrip().split('\t')
        sample, bamfile = line

        pipettor.run([('samtools', 'view', '-h', bamfile, rchrom + ':' + rstart + '-' + rend),
                      ('samtools', 'fasta', '-')],
                     stdout=open(tempprefix + '.reads.fasta', 'w'))

        samfile = pysam.AlignmentFile(bamfile, 'rb')
        ft.generate_genomic_clipping_reference(tempprefix, samfile, rchrom, rstart, rend)
        samfile.close()

        goodannotaligns = generate_good_match_to_annot(args, tempprefix, rchrom, annottranscripttoexons,
                                                       allannottranscripts, genome, bamfile)

        samfile = pysam.AlignmentFile(goodannotaligns, 'rb')
        sjtoends, goodaligntoannot = ft.filter_correct_group_reads(args, tempprefix, rchrom, rstart, rend, samfile,
                                                                set(), intervalTree, junctionBoundaryDict,
                                                                generatefasta=False, sjtoends={},
                                                                returnusedreads=True, allowsecondary=True)
        samfile.close()
        samfile = pysam.AlignmentFile(bamfile, 'rb')
        sjtoends = ft.filter_correct_group_reads(args, tempprefix, rchrom, rstart, rend, samfile, goodaligntoannot,
                                              intervalTree, junctionBoundaryDict,
                                              generatefasta=False, sjtoends=sjtoends, returnusedreads=False)
        samfile.close()
        genetojuncs, nogenejuncs, sereads = group_juncs_by_annot_gene(sjtoends, juncstotranscript, junctogene,
                                                                      genetoannotjuncs, allsplicedexons)
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
        line = line.rstrip().split('\t')
        allsamples.append(line[0])
    out = open(tempDir + '0000.header.diffsplice.counts.tsv', 'w')
    out.write('\t'.join(['featureID'] + allsamples) + '\n')
    out.close()

    for rchrom, rstart, rend in allregions:
        if rchrom in knownchromosomes:
            juncstotranscript, junctogene, allannotse, allsplicedexons, genetoannotjuncs, \
                genetostrand, annottranscripttoexons, allannottranscripts = {}, {}, [], {}, {}, {}, {}, []
            if args.gtf:
                juncstotranscript, junctogene, allannotse, allsplicedexons, genetoannotjuncs, \
                    genetostrand, annottranscripttoexons, allannottranscripts = \
                    regionstoannotdata[(rchrom, rstart, rend)].returndata()

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
                                ['.diffsplice.bed', '.diffsplice.counts.tsv', '.diffsplice.readends.bed'])

    if not args.keep_intermediate:
        shutil.rmtree(tempDir)

    genome.close()


if __name__ == "__main__":
    collapsefrombam()
