#!/usr/bin/env python3

import sys
import argparse
from bisect import bisect_left
import pysam

# import pybedtools as pbt

###python3 /private/groups/brookslab/cafelton/git-flair/flair/src/flair/remove_internal_priming.py -r /private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa -a /private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf -b /private/groups/brookslab/tswon/flair_runs/flair3_runs/Berger_Patient_Data/BergerPatientData110124/flnc-1_flair/flair.align.bam -o flnc-1.intprimingremoved.bam

def parseargs():
    parser = argparse.ArgumentParser(description='''for removing internal priming reads from bam file
		usage=python3 remove_internal_priming.py -b bamfile -o outputfile -a annotgtf''')
    parser.add_argument('-o', '--output', help='output file name (also .bam)')
    parser.add_argument('-b', '--bam',
                        help='input bam file, can be aligned to genome or transcriptome')
    parser.add_argument('-a', '--annot',
                        help='gtf file of annotated transcripts, not required, only use if genomic alignment')
    parser.add_argument('--annotbed',
                        help='bed file of annotated transcripts, not required, only use if transcriptomic alignment')
    parser.add_argument('-r', '--reference',
                        help='fa file of reference genome or transcriptome')
    parser.add_argument('--intprimingthreshold', type=int, default=12,
                        help='number of bases that are at leas 75%% As required to call read as internal priming')
    parser.add_argument('--intprimingfracAs', type=float, default=0.6,
                        help='number of bases that are at leas 75%% As required to call read as internal priming')
    args = parser.parse_args()
    return args


def checkIsNearAnnotEnd(read3endpos, annotends):
    pos1 = bisect_left(annotends, read3endpos)
    if pos1 == len(annotends): return abs(annotends[pos1 - 1] - read3endpos) <= 200
    disttoend = min(abs(annotends[pos1 - 1] - read3endpos), abs(annotends[pos1] - read3endpos))
    return disttoend <= 200


def getannotends(annotfile):
    if annotfile:
        annottranscriptends = {}
        lastexon = None
        for line in open(annotfile):
            if line[0] != '#':
                line = line.split('\t')
                if line[2] == 'transcript' and 'basic' in line[8]:
                    chr, start, stop, strand = line[0], int(line[3]), int(line[4]), line[6]
                    if chr not in annottranscriptends: annottranscriptends[chr] = set()
                    if strand == '+':
                        annottranscriptends[chr].add(stop)
                    else:
                        annottranscriptends[chr].add(start)
                    if lastexon: ###assuming that gtf files have exons in order
                        if lastexon[2] == '+':
                            for i in range(lastexon[0]+200, lastexon[1], 200):
                                annottranscriptends[chr].add(i)
                        else:
                            for i in range(lastexon[0], lastexon[1]-200, 200):
                                annottranscriptends[chr].add(i)
                elif line[2] == 'exon' and 'basic' in line[8]:
                    lastexon = (int(line[3]), int(line[4]), line[6])
        for chr in annottranscriptends:
            annottranscriptends[chr] = sorted(list(annottranscriptends[chr]))
        return annottranscriptends
    else: return None

###add annotation-reliant check for transcript end, implement binary search
def checkInternalPriming(read3endpos, thischr, genome, reqfreq, threshold):  # , printoutput):#, threshold):
    genomeseqnearend = genome.fetch(thischr, max(read3endpos - 30, 0), min(read3endpos + 30, genome.get_reference_length(thischr))).upper()
    maxlen, maxfreq = 0, 0
    if len(genomeseqnearend) > threshold*2:
        halfseqlen = int(len(genomeseqnearend)/2)
        for i in list(range(-1 * halfseqlen, -1*threshold)) + list(range(threshold, halfseqlen)):
            thisseq = genomeseqnearend[min(i + halfseqlen, halfseqlen): max(i + halfseqlen, halfseqlen)]
            thiscount = max(thisseq.count('A'), thisseq.count('T'))
            thisfreq = thiscount / len(thisseq) if len(thisseq) > 0 else 0
            if thisfreq >= reqfreq and len(thisseq) > maxlen: maxlen, maxfreq = len(thisseq), thisfreq
            # if 55205880 < read3endpos < 55205919: print(len(genomeseqnearend), i, min(i + halfseqlen, halfseqlen), max(i + halfseqlen, halfseqlen))
            # if 55205880 < read3endpos < 55205919: print(thisseq, thiscount, len(thisseq), thisfreq, thisfreq >= reqfreq, len(thisseq) > maxlen, maxlen, maxfreq)
        # if 32034696 < read3endpos < 32034835 and maxlen < 10:
        # if 55205880 < read3endpos < 55205919 and maxlen < threshold:
        #     print(genomeseqnearend, maxlen, maxfreq)
    return maxlen  # >= threshold

def removeinternalpriming(refname, refstart, refend, isrev, genome, annottranscriptends, annotexons, threshold, fracAs):
    read3endpos = refend if not isrev else refstart
    if not annottranscriptends:  # if aligned to transcriptome, check distance to transcript end
        # reflen = genome.get_reference_length(refname)
        if annotexons and refname in annotexons:
            theseexons = annotexons[refname]
            # print('near end', refname, read3endpos > sum(theseexons) - theseexons[-1], read3endpos, sum(theseexons) - theseexons[-1])
            if len(theseexons) > 1 and read3endpos > sum(theseexons) - theseexons[-1]:
                return True
            elif len(theseexons) == 1 and read3endpos >= theseexons[0]-200:
                return True
    intprimlen = checkInternalPriming(read3endpos, refname, genome, fracAs, threshold)  # , c == 1000)
    if intprimlen < threshold:  ##read doesn't have stretch of As beyond threshold
        # print('passes threshold', intprimlen)
        return True
    # if 32034696 < read3endpos < 32034835 and intprimlen < 10: break
    # if 55205880 < read3endpos < 55205919 and intprimlen < threshold: break
    elif annottranscriptends and refname in annottranscriptends:  ##if annot provided + read has strech of As, check if end near annot end
        isnearannotend = checkIsNearAnnotEnd(read3endpos, annottranscriptends[refname])
        if isnearannotend:
            return True
    return False

def processbamfile(bamfile, outfile, annottranscriptends, annotexons, reference, threshold, fracAs):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    genome = pysam.FastaFile(reference)
    outbam = pysam.AlignmentFile(outfile, 'wb', template=samfile)
    for read in samfile:
        if read.is_mapped:
            isgood = removeinternalpriming(read.reference_name, read.reference_start, read.reference_end, read.is_reverse, genome, annottranscriptends, annotexons, threshold, fracAs)
            if isgood:
                outbam.write(read)
    samfile.close()
    outbam.close()
    pysam.index(outfile)


if __name__ == '__main__':
    args = parseargs()
    annottranscriptends = getannotends(args.annot)
    processbamfile(args.bam, args.output, annottranscriptends, None, args.reference, args.intprimingthreshold, args.intprimingfracAs)
