#!/usr/bin/env python3

import argparse
from flair import FlairInputDataError
from flair.counts_matrix_tsv import ID_COLUMNS, CountsRow
from flair.cpm_tsv import cpm_writer
from flair.pycbio.sys import cli
from flair.pycbio.sys import fileOps
from flair.cli_help import FlairHelpFormatter

def parse_args():
    parser = argparse.ArgumentParser(prog='counts_to_cpm',
                                     formatter_class=FlairHelpFormatter,
                                     description='write a counts matrix as counts per million')
    parser.add_argument('counts_matrix', help='the counts matrix of flair quantify')
    parser.add_argument('cpm_matrix', help='the CPM matrix to write')
    return parser.parse_args()


def counts_to_cpm(counts_matrix, outfilename, empty_ok=False):
    """Write a counts matrix as counts per million: each count over the sample's
    counted reads, times a million.  A sample with no counts in any row has no CPM:
    an error, or with empty_ok, 0 in every row.  Returns the names of such samples."""
    header = counts_matrix.readline().rstrip().split('\t')
    num_id_cols = len(ID_COLUMNS)
    num_samples = len(header[num_id_cols:])
    matrix_data = [header]
    totals = [0] * num_samples
    for line in counts_matrix:
        line = line.rstrip().split('\t')
        ids, counts = line[:num_id_cols], [float(count) for count in line[num_id_cols:]]
        for n in range(num_samples):
            totals[n] += counts[n]
        matrix_data += [ids + counts]

    totals = [total / 1e6 for total in totals]
    empty = [header[n + num_id_cols] for n, total in enumerate(totals) if total == 0]
    if len(empty) > 0 and not empty_ok:
        raise FlairInputDataError(f"{outfilename}: samples {', '.join(empty)} have no counts in any row, "
                                  "so CPM cannot be computed; drop them from the counts matrix")

    sample_columns = matrix_data[0][num_id_cols:]
    with cpm_writer(outfilename, sample_columns) as writer:
        for line in matrix_data[1:]:
            cpm = [round(line[n + num_id_cols] / totals[n], 5) if totals[n] > 0 else 0
                   for n in range(num_samples)]
            writer.writeCountsRow(CountsRow(line[0], line[1], cpm))
    return empty


def main():
    args = parse_args()
    with cli.ErrorHandler():
        with fileOps.opengz(args.counts_matrix) as counts_matrix:
            counts_to_cpm(counts_matrix, args.cpm_matrix)


if __name__ == '__main__':
    main()
