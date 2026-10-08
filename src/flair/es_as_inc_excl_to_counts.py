#!/usr/bin/env python3
"""Attach counts to the skipped exon events es_as.py calls, writing the
es events.quant.tsv."""

import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np  # noqa: E402 - openblas setting must be before numpy import
from flair.pycbio.sys import cli  # noqa: E402


def read_counts(counts_matrix_tsv):
    "the counts of each isoform, and the sample column names"
    counts = {}
    with open(counts_matrix_tsv) as fh:
        header = next(fh).rstrip().split()
        for line in fh:
            cols = line.rstrip().split()
            counts[cols[0]] = np.asarray(cols[1:], dtype=np.float32)
    return counts, header[1:]


def event_counts(counts, isos, num_samples):
    "the counts of an event, summed over the isoforms supporting it"
    return np.sum(np.asarray([counts.get(iso, np.zeros(num_samples))
                              for iso in isos.split(',')]), axis=0)


def write_event(exon, counts, inc_isos, exc_isos, num_samples):
    for label, isos in (('inclusion', inc_isos), ('exclusion', exc_isos)):
        values = event_counts(counts, isos, num_samples)
        print(f'{label}_{exon}', exon, '\t'.join(str(v) for v in values), isos, sep='\t')


def write_events(es_events_tsv, counts, sample_names):
    print('\t'.join(['feature_id', 'coordinate'] + sample_names + ['isoform_ids']))
    with open(es_events_tsv) as fh:
        for line in fh:
            # split on tabs, not on whitespace: an empty isoform list is an empty
            # field, which a whitespace split drops, leaving a short row
            exon, _, _, _, inc_isos, exc_isos = line.rstrip('\n').split('\t')
            # an exon needs an isoform on each side to be a skipping event
            if inc_isos and exc_isos:
                write_event(exon, counts, inc_isos, exc_isos, len(sample_names))


def main():
    counts_matrix_tsv, es_events_tsv = sys.argv[1], sys.argv[2]
    with cli.ErrorHandler():
        counts, sample_names = read_counts(counts_matrix_tsv)
        write_events(es_events_tsv, counts, sample_names)


if __name__ == '__main__':
    main()
