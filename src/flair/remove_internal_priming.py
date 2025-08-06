#!/usr/bin/env python3

import sys
import argparse
from bisect import bisect_left
import pysam


def checkIsNearAnnotEnd(read3endpos, annotends):
    """
    Implements binary search for nearest transcript end, return true if pos is <=200bp from nearest end
    """
    pos1 = bisect_left(annotends, read3endpos)
    if pos1 == len(annotends): return abs(annotends[pos1 - 1] - read3endpos) <= 200
    disttoend = min(abs(annotends[pos1 - 1] - read3endpos), abs(annotends[pos1] - read3endpos))
    return disttoend <= 200


def getannotends(annotfile):
    """
    Takes in a gtf annotation file, returns a dictionary od chromsome to sorted list of transcript end positions
    """
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
def checkInternalPriming(read3endpos, thischr, genome, reqfreq, threshold):
    """
    Checks the genomic sequence adjacent to the read end position for a stretch of As with
    a frequency >= reqfreq and a length >= threshold
    """
    genomeseqnearend = genome.fetch(thischr, max(read3endpos - 30, 0), min(read3endpos + 30, genome.get_reference_length(thischr))).upper()
    maxlen, maxfreq = 0, 0
    if len(genomeseqnearend) > threshold*2:
        halfseqlen = int(len(genomeseqnearend)/2)
        for i in list(range(-1 * halfseqlen, -1*threshold)) + list(range(threshold, halfseqlen)):
            thisseq = genomeseqnearend[min(i + halfseqlen, halfseqlen): max(i + halfseqlen, halfseqlen)]
            thiscount = max(thisseq.count('A'), thisseq.count('T'))
            thisfreq = thiscount / len(thisseq) if len(thisseq) > 0 else 0
            if thisfreq >= reqfreq and len(thisseq) > maxlen: maxlen, maxfreq = len(thisseq), thisfreq
    return maxlen >= threshold

def removeinternalpriming(refname, refstart, refend, isrev, genome, annottranscriptends, annotexons, threshold, fracAs):
    """
    Given info from a an aligned bam read, check whether it has internal priming
    refname, refstart, refend, isrev - all info about read alignment
    genome: pysam.FastaFile object
    annottranscriptends: [genomic alignment only] dictionary of chrom to sorted list of transcript end pos
    annotexons: [transcriptomic alignment only] dictionary of transcript name to list of exon lengths
    threshold: max length of stretch of As before something is internal priming
    fracAs: minimum frequency of As in sequence to qualify as polyA (6/8 bp=A -> threshold=0.75)
    """
    read3endpos = refend if not isrev else refstart
    # if aligned to transcriptome, check distance to transcript end
    # if read end is close enough to transcript end, return True (no internal priming)
    if not annottranscriptends:
        if annotexons and refname in annotexons:
            theseexons = annotexons[refname]
            # multi exon transcript
            if len(theseexons) > 1 and read3endpos > sum(theseexons) - theseexons[-1]:
                return True
            # single exon transcript
            elif len(theseexons) == 1 and read3endpos >= theseexons[0]-200:
                return True
    # if read doesn't have stretch of As beyond threshold, doesn't have internal priming
    if not checkInternalPriming(read3endpos, refname, genome, fracAs, threshold):
        return True
    # if annot provided + read has strech of As, check if end near annot end
    elif annottranscriptends and refname in annottranscriptends:
        isnearannotend = checkIsNearAnnotEnd(read3endpos, annottranscriptends[refname])
        if isnearannotend:
            return True
    return False
