#!/usr/bin/env python3
import sys
import csv
import os
from flair import FlairInputDataError

try:
    fasta = open(sys.argv[1])
    outfilename = sys.argv[2]
    if len(sys.argv) > 3:
        outfilename2 = sys.argv[3]
    else:
        outfilename2 = ''
except:
    raise FlairInputDataError('usage: fasta_seq_lengths fasta outfilename [outfilename2]\n')

length_frequencies = {}
with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    seqlen = 0
    name = None
    for line in fasta:
        line = line.rstrip()
        if line.startswith('>'):
            if seqlen:
                writer.writerow([name, seqlen])
                if seqlen not in length_frequencies:
                    length_frequencies[seqlen] = 0
                length_frequencies[seqlen] += 1
            name = line[1:]
            seqlen = 0
            continue
        seqlen += len(line.rstrip())
    writer.writerow([name, seqlen])

if outfilename2:
    alllengths = sorted(length_frequencies.keys())
    with open(outfilename2, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        for length in alllengths:
            writer.writerow([length, length_frequencies[length]])
