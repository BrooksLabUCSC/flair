#!/usr/bin/env python3

import sys
import pysam
import re
import logging
import os

basetocomp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
def revcomp(seq):
    seq = list(seq.upper())[::-1]
    return ''.join([basetocomp[x] for x in seq])


def addOtherJuncs(juncs, filetype, bedJuncs, minsup, chromosomes, printErrFname, known, verbose, printErr):
    # guess what kind of bedFile
    if os.path.getsize(bedJuncs) == 0:
        # raise Exception("Empty junctions BED file, not supported")
        logging.warn('WARNING: orthogonal junctions bed file is EMPTY')
        return juncs, chromosomes, False

    if filetype == 'bed':
        # normal bed
        strandCol = 5
        starOffset = 0
        scorecol=4
    else: #filetype == 'tab'
        strandCol = 3
        starOffset = 1
        scorecol = 6

    tempJuncs = list()
    addedFlag = False
    with open(bedJuncs,'r') as bedLines:
        for line in bedLines:
            cols = line.rstrip().split()
            if len(cols) == 12:
                raise Exception("Bed12 not currently supported for short-read junctions. Please convert to bed6 or bed9.")

            chrom, c1, c2, strand, score = cols[0], int(cols[1])-starOffset, int(cols[2]), cols[strandCol], int(cols[scorecol])

            if score >= minsup:
                if chrom not in juncs:
                    juncs[chrom] = dict()

                if c2-c1 < 5:
                    continue

                if filetype == 'tab':
                    if strand == "1": strand = "+"
                    elif strand == "2": strand = "-"
                    else: continue

                chromosomes.add(chrom)
                key = (c1, c2, strand)
                if key in juncs[chrom]:
                    juncs[chrom][key] = "both"
                    continue
                tempJuncs.append((chrom,c1,c2,strand))
                addedFlag = True
    if addedFlag == False:
        return juncs, chromosomes, addedFlag

    for chrom,c1,c2,strand in tempJuncs:
        key = (c1, c2, strand)
        known1, known2 = known.get((chrom, c1), None), known.get((chrom, c2), None)
        if known1 is not None:
            if known1 != strand:
                continue
        if known2 is not None:
            if known2 != strand:
                continue

        if key not in juncs[chrom]:
            juncs[chrom][key] = "sr"

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** GTF Juncs + other juncs now total %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)

    return juncs, chromosomes, addedFlag


def gtfToSSBed(gtffile, knownSS, printErr, printErrFname, verbose):
    ''' Convenience function, reformats GTF to bed'''

    # First: get all exons per transcript.
    exons = dict()
    chromosomes = set()
    with open(gtffile,'r') as lines:
        for l in lines:
            if l[0] == "#": # skip header lines
                continue

            cols = l.split("\t")

            if "exon" == cols[2]:

                # -1 for 1 to 0 based conversion
                chrom, c1, c2, strand = cols[0], int(cols[3])-1, int(cols[4]), cols[6]
                chromosomes.add(chrom)
                try:
                    txn = re.search('transcript_id "([^\"]+)"', l).group(1)
                except Exception as ex:
                    raise Exception("** ERROR expect transcript_id in GTF format, cannot read %s" % gtffile) from ex

                key = (chrom, txn, strand)

                if key not in exons:
                    exons[key] = list()
                exons[key].append(c1)
                exons[key].append(c2)

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Read GTF. Got %s transcripts" % len(list(exons.keys())), file=fo)
            print("** Getting introns...Read GTF", file=fo)

    # Second: get junction and splice sites from transcript exons.
    txnList = list(exons.keys())
    juncs = dict()

    for exonInfo in txnList:
        chrom, txn, strand = exonInfo

        if chrom not in juncs:
            juncs[chrom] = dict()

        coords = list(exons[exonInfo])

        # assume lowest and highest as TSS and TES, and remove them
        coords.sort()
        coords = coords[1:-1]

        # Coords is list of exons, so a list less than 2 is a single exon gene.
        if len(coords) < 2: continue

        for pos in range(0,len(coords)-1,2):
            c1 = coords[pos]
            c2 = coords[pos+1]

            if abs(c2 - c1) <= 5:
                continue

            juncs[chrom][(c1,c2,strand)] = "gtf"
            knownSS[(chrom, c1)] = strand
            knownSS[(chrom, c2)] = strand

    if printErr:
        with open(printErrFname,'a+') as fo:
            print("** Created %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)
    return juncs, chromosomes, knownSS
