#!/usr/bin/env python3

import sys
import csv
import os
from flair import FlairInputDataError

def pull_starts(bedfile, outfilename, nvrna=False, reverse=False):
    bedfh = open(bedfile)

    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom, name, start, end, strand, numblocks = line[0], line[3], int(line[1]), int(line[2]), line[5], line[9]

            if not nvrna and numblocks == 1:  # single exon gene, add both sides for cdna since strand is hard
                writer.writerow([chrom, start, start, name])
                writer.writerow([chrom, end, end, name])
                continue
            elif not reverse and '+' in strand or reverse and '-' in strand:
                tss = start
            elif not reverse and '-' in strand or reverse and '+' in strand:
                tss = end
            else:  # ambiguous strand, write both
                writer.writerow([chrom, start, start, name])
                tss = end
            writer.writerow([chrom, tss, tss, name])

if __name__ == "__main__":
    try:
        bedfile = sys.argv[1]
        outfilename = sys.argv[2]
        nvrna = reverse = False
        if len(sys.argv) > 3:
            nvrna = 'nvrna' in sys.argv[3]  # specify if stranded protocol
            reverse = 'reverse' in sys.argv[3]
    except:
        raise FlairInputDataError('pull_starts.py bed outfilename [nvrna]')

    pull_starts(bedfile, outfilename, nvrna=nvrna, reverse=reverse)
