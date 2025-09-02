#!/usr/bin/env python3

import csv
import os
import argparse

# TODO: update readthedocs with flags
def main():
    parser = argparse.ArgumentParser(description='script for collapse-range')
    required = parser.add_argument_group('required named arguments')
    required.add_argument('--counts_file', type=str, required=True,
            help='counts file from count_sam_transcripts')
    required.add_argument('--bed', type=str, required=True,
            help='bed file to filter by counts')
    parser.add_argument('--min_reads', type=float, default=0,
            help='minimum number of supporting reads')
    required.add_argument('--out_bed', type=str, required=True,
            help='output file name, bed')
    parser.add_argument('--generate_map',  dest='generate_map', default=False,
            help='''name of read-isoform mapping file to be filtered''')
    parser.add_argument('--append_counts', action='store_true',
            help='''specify this argument to append the supporting
            read counts as a new column in the BED file''')
    args = parser.parse_args()

    match_counts(counts_file=args.counts_file, output_file=args.out_bed, bed=args.bed, append_counts=args.append_counts,
          min_reads=args.min_reads, isoform_file=args.generate_map)

def match_counts(counts_file, output_file, bed, append_counts=False, min_reads=0, isoform_file=False):
    counts = {}
    for line in open(counts_file):
        line = line.rstrip().split('\t')
        counts[line[0]] = float(line[1])

    if isoform_file:
        map_file = {}
        for line in open(isoform_file):
            line = line.rstrip().split('\t')
            map_file[line[0]] = line[1]
        out_map = open(isoform_file, 'wt')  # writing to a file of the same name TODO fix this abomination
        writer_map = csv.writer(out_map, delimiter='\t', lineterminator=os.linesep)

    out_bedfh = open(output_file, 'wt')
    writer = csv.writer(out_bedfh, delimiter='\t', lineterminator=os.linesep)
    for line in open(bed):
        line = line.rstrip().split('\t')
        name = line[3]

        if name in counts:
            count = counts[name]
        else:
            count = 0
        if count >= min_reads:
            if append_counts:
                writer.writerow(line + [count])
            else:
                writer.writerow(line)
            if isoform_file:
                writer_map.writerow([name, map_file[name]])

if __name__ == "__main__":
    main()
