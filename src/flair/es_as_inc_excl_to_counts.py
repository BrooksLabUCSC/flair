#!/usr/bin/env python3
"""Attach counts to the skipped exon events es_as.py calls, writing the
es events.quant.tsv."""

import os
import sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np  # noqa: E402 - openblas setting must be before numpy import
from flair.counts_matrix_tsv import read_sample_columns, read_counts_rows  # noqa: E402
from flair.es_events_tsv import EsEventsReader  # noqa: E402
from flair.event_quant_tsv import EventQuantWriter  # noqa: E402
from flair.pycbio.sys import cli  # noqa: E402

def side_counts(data, n_samps, isos):
    "the counts of one side of an event, summed over its isoforms"
    vals = np.asarray([data.get(iso, np.zeros(n_samps)) for iso in isos])
    return np.sum(vals, axis=0)

def write_es_events(counts_matrix_tsv, es_events_tsv):
    sample_names = read_sample_columns(counts_matrix_tsv)
    n_samps = len(sample_names)
    data = {row.isoform_id: np.asarray(row.counts, dtype=np.float32)
            for row in read_counts_rows(counts_matrix_tsv)}

    with EventQuantWriter(None, sample_names, outFh=sys.stdout) as writer:
        for row in EsEventsReader(es_events_tsv):
            # an exon needs an isoform on each side to be a skipping event
            if row.inclusion_isos and row.exclusion_isos:
                writer.writeSide('inclusion', row.exon, row.exon,
                                 side_counts(data, n_samps, row.inclusion_isos),
                                 row.inclusion_isos)
                writer.writeSide('exclusion', row.exon, row.exon,
                                 side_counts(data, n_samps, row.exclusion_isos),
                                 row.exclusion_isos)

def main():
    with cli.ErrorHandler():
        write_es_events(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
