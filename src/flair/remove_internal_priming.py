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
                        help='gtf file of annotated transcripts, not required')
    parser.add_argument('-r', '--reference',
                        help='fa file of reference genome or transcriptome')
    parser.add_argument('-t', '--threshold', type=int, default=12,
                        help='number of bases that are at leas 75% As required to call read as internal priming')
    parser.add_argument('-f', '--fracAs', type=float, default=0.6,
                        help='number of bases that are at leas 75% As required to call read as internal priming')
    args = parser.parse_args()
    return args


def checkIsNearAnnotEnd(read3endpos, annotends):
    pos1 = bisect_left(annotends, read3endpos)
    if pos1 == len(annotends): return abs(annotends[pos1 - 1] - read3endpos) < 200
    disttoend = min(abs(annotends[pos1 - 1] - read3endpos), abs(annotends[pos1] - read3endpos))
    return disttoend < 200


def getannotends(annotfile):
    if annotfile:
        annottranscriptends = {}
        for line in open(annotfile):
            if line[0] != '#':
                line = line.split('\t', 7)
                if line[2] == 'transcript':
                    chr, start, stop, strand = line[0], int(line[3]), int(line[4]), line[6]
                    if chr not in annottranscriptends: annottranscriptends[chr] = set()
                    if strand == '+':
                        annottranscriptends[chr].add(stop)
                    else:
                        annottranscriptends[chr].add(start)
        for chr in annottranscriptends:
            annottranscriptends[chr] = sorted(list(annottranscriptends[chr]))
        return annottranscriptends
    else: return None


###add annotation-reliant check for transcript end, implement binary search
def checkInternalPriming(read3endpos, thischr, genome, reqfreq, threshold):  # , printoutput):#, threshold):
    genomeseqnearend = genome.fetch(thischr, read3endpos - 30, read3endpos + 30).upper()
    maxlen, maxfreq = 0, 0
    if len(genomeseqnearend) > threshold*2:
        halfseqlen = int(len(genomeseqnearend)/2)
        for i in list(range(-1 * halfseqlen, -1*threshold)) + list(range(threshold, halfseqlen)):
            otheredge = halfseqlen-threshold if i < 1 else halfseqlen + threshold
            thisseq = genomeseqnearend[min(i + halfseqlen, otheredge): max(i + halfseqlen, otheredge)]
            thiscount = max(thisseq.count('A'), thisseq.count('T'))
            thisfreq = thiscount / len(thisseq) if len(thisseq) > 0 else 0
            if thisfreq >= reqfreq and len(thisseq) > maxlen: maxlen, maxfreq = len(thisseq), thisfreq
        # if 32034696 < read3endpos < 32034835 and maxlen < 10:
        #     for i in list(range(-30, 0)) + list(range(1, 30)):
        #         thisseq = genomeseqnearend[min(i + 30, 30): max(i + 30, 30)]
        #         thiscount = max(thisseq.count('A'), thisseq.count('T'))
        #         thisfreq = thiscount / len(thisseq)
        #         print(thisseq, thiscount, thisfreq)
        #     print(genomeseqnearend, maxlen, maxfreq)
    return maxlen  # >= threshold

def removeinternalpriming(bamfile, outfile, annottranscriptends, reference, threshold, fracAs):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    genome = pysam.FastaFile(reference)
    outbam = pysam.AlignmentFile(outfile, 'wb', template=samfile)
    for read in samfile:
        if read.is_mapped:
            read3endpos = read.reference_end if not read.is_reverse else read.reference_start
            thischr = read.reference_name
            if not annottranscriptends: #if aligned to transcriptome, check distance to transcript end
                reflen = genome.get_reference_length(thischr)
                if reflen - read3endpos < 200:
                    outbam.write(read)
                    continue
            intprimlen = checkInternalPriming(read3endpos, thischr, genome, fracAs, threshold)  # , c == 1000)
            if intprimlen < threshold: ##read doesn't have stretch of As beyond threshold
                outbam.write(read)
            # if 32034696 < read3endpos < 32034835 and intprimlen < 10: break
            elif annottranscriptends and thischr in annottranscriptends: ##if annot provided + read has strech of As, check if end near annot end
                isnearannotend = checkIsNearAnnotEnd(read3endpos, annottranscriptends[thischr])
                if isnearannotend:
                    outbam.write(read)
    samfile.close()
    outbam.close()
    pysam.index(outfile)


if __name__ == '__main__':
    args = parseargs()
    annottranscriptends = getannotends(args.annot)
    removeinternalpriming(args.bam, args.output, annottranscriptends, args.reference, args.threshold, args.fracAs)
