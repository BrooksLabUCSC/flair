#!/usr/bin/env python3
import sys
import csv
import os

def select_from_bed(inputbed, queryfile, outputfile):
    keepnames = set()

    for line in open(inputbed):  # bed of promoter-supported reads
        line = line.rstrip().split('\t')
        keepnames.add(line[3])

    with open(outputfile, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        for line in open(queryfile):
            line = line.rstrip().split('\t')
            if line[3] in keepnames:
                writer.writerow(line)
