#!/usr/bin/env python3

import argparse
from flair import FlairInputDataError
from flair.diffsplice_fishers_tsv import diffsplice_fishers_writer
from flair.event_quant_tsv import EventQuantReader
from flair.pycbio.sys import cli

def build_parser():
    desc = "Fisher's exact test of two samples' inclusion and exclusion counts for each splicing event"
    parser = argparse.ArgumentParser(prog='diffsplice_fishers_exact', description=desc)
    parser.add_argument('events_quant_tsv', help='event inclusion/exclusion counts from flair diffsplice')
    parser.add_argument('colname1', help='column name of the first sample to compare')
    parser.add_argument('colname2', help='column name of the second sample to compare')
    parser.add_argument('fishers_tsv', help='output TSV of the per-event test results')
    return parser

def sample_column_index(reader, colname, events_quant_tsv):
    "index of a named sample column within a row"
    if colname not in reader.sample_columns:
        raise FlairInputDataError(
            f"{events_quant_tsv} has no sample column named {colname}; it has: "
            f"{' '.join(reader.sample_columns)}")
    return reader.columns.index(colname)

def feature_of(row):
    "the feature_id without its inclusion or exclusion prefix"
    return row.feature_id[row.feature_id.find('_') + 1:]

def group_events(reader, col1, col2):
    """The rows of each event, with the two by two table the test needs.  Events are
    keyed by the feature the two sides share, which orders the output by locus."""
    events = {}
    for row in reader:
        feature = feature_of(row)
        if feature not in events:
            events[feature] = {'entries': [], 'counts': []}
        events[feature]['entries'].append(row)
        events[feature]['counts'].append([float(row[col1]), float(row[col2])])
    return events

def diffsplice_fishers_exact(events_quant_tsv, colname1, colname2, fishers_tsv):
    # imported here rather than at module scope; scipy takes 0.7s to load and this
    # is the only use of it
    import scipy.stats as sps
    with EventQuantReader(events_quant_tsv) as reader:
        col1 = sample_column_index(reader, colname1, events_quant_tsv)
        col2 = sample_column_index(reader, colname2, events_quant_tsv)
        columns = list(reader.columns)
        events = group_events(reader, col1, col2)

    with diffsplice_fishers_writer(fishers_tsv, columns, colname1, colname2) as writer:
        for feature in sorted(events.keys()):
            pval = sps.fisher_exact(events[feature]['counts'])[1]
            for row in events[feature]['entries']:
                writer.writeRow(list(row) + [pval])

def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        diffsplice_fishers_exact(args.events_quant_tsv, args.colname1, args.colname2,
                                 args.fishers_tsv)


if __name__ == "__main__":
    main()
