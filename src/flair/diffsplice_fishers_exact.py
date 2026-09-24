#!/usr/bin/env python3

import argparse
import csv
import os
from flair import FlairInputDataError
from flair.pycbio.sys import cli

def build_parser():
    desc = "Fisher's exact test of two samples' inclusion and exclusion counts for each splicing event"
    parser = argparse.ArgumentParser(prog='diffsplice_fishers_exact', description=desc)
    parser.add_argument('events_quant_tsv', help='event inclusion/exclusion counts from flair diffsplice')
    parser.add_argument('colname1', help='column name of the first sample to compare')
    parser.add_argument('colname2', help='column name of the second sample to compare')
    parser.add_argument('fishers_tsv', help='output TSV of the per-event test results')
    return parser

def diffsplice_fishers_exact(events_quant_tsv, colname1, colname2, fishers_tsv):
    # imported here rather than at module scope; scipy takes 0.7s to load and this
    # is the only use of it
    import scipy.stats as sps
    events_quant = open(events_quant_tsv)
    outfilename = fishers_tsv
    header = events_quant.readline().rstrip().split('\t')

    if colname1 in header:
        col1 = header.index(colname1)
    else:
        raise FlairInputDataError('Could not find {} in {}\n'.format(colname1, ' '.join(header)))

    if colname2 in header:
        col2 = header.index(colname2)
    else:
        raise FlairInputDataError('Could not find {} in {}\n'.format(colname2, ' '.join(header)))

    events = {}
    for line in events_quant:
        line = line.rstrip().split('\t')
        feature = line[0][line[0].find('_') + 1:]
        if feature not in events:
            events[feature] = {}
            events[feature]['entries'] = []
            events[feature]['counts'] = []
        events[feature]['entries'] += [line]
        events[feature]['counts'] += [[float(line[col1]), float(line[col2])]]

    features_sorted = sorted(events.keys())
    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        writer.writerow(header + [colname1 + '-' + colname2 + '_pval'])
        for feature in features_sorted:
            for line in events[feature]['entries']:
                writer.writerow(line + [sps.fisher_exact(events[feature]['counts'])[1]])


def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        diffsplice_fishers_exact(args.events_quant_tsv, args.colname1, args.colname2,
                                 args.fishers_tsv)


if __name__ == "__main__":
    main()
