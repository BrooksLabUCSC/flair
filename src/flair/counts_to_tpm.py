#!/usr/bin/env python3

import argparse
from flair import FlairInputDataError
from flair.counts_matrix_tsv import ID_COLUMNS, CountsRow
from flair.seq_lengths_tsv import SeqLengthsReader
from flair.tpm_tsv import tpm_writer
from flair.pycbio.sys import cli
from flair.pycbio.sys import fileOps
from flair.cli_help import FlairHelpFormatter

def parse_args():
    parser = argparse.ArgumentParser(prog='counts_to_tpm',
                                     formatter_class=FlairHelpFormatter,
                                     description='write a counts matrix as transcripts per million')
    parser.add_argument('counts_matrix', help='the counts matrix of flair quantify')
    parser.add_argument('tpm_matrix', help='the TPM matrix to write')
    parser.add_argument('isoform_sizes', nargs='?',
                        help='each isoform\'s length, from fasta_seq_lengths; without it the counts '
                             'are not length normalized, giving reads per million')
    return parser.parse_args()

def counts_to_tpm(counts_matrix, outfilename, sizefile=None, empty_ok=False):
    """Write a counts matrix as transcripts per million.  A sample with no counts in
    any row has no TPM: an error, or with empty_ok, 0 in every row.  Returns the names
    of such samples."""
    sizes = {}
    if sizefile:
        sizes = {row.name: row.length for row in SeqLengthsReader(sizefile)}

    header = counts_matrix.readline().rstrip().split('\t')
    num_id_cols = len(ID_COLUMNS)
    num_samples = len(header[num_id_cols:])
    matrix_data = [header]
    all_rpk = [0] * num_samples
    for line in counts_matrix:
        line = line.rstrip().split('\t')
        ids, counts = line[:num_id_cols], line[num_id_cols:]
        isoform_id = ids[-1]
        if sizes:
            rpk = [float(count) / sizes[isoform_id] for count in counts]
        else:
            rpk = [float(count) for count in counts]
        for n in range(num_samples):
            all_rpk[n] += rpk[n]
        matrix_data += [ids + rpk]

    all_rpk = [rpk / 1e6 for rpk in all_rpk]
    empty = [header[n + num_id_cols] for n, rpk in enumerate(all_rpk) if rpk == 0]
    if len(empty) > 0 and not empty_ok:
        raise FlairInputDataError(f"{outfilename}: samples {', '.join(empty)} have no counts in any row, "
                                  "so TPM cannot be computed; drop them from the counts matrix")

    sample_columns = matrix_data[0][num_id_cols:]
    with tpm_writer(outfilename, sample_columns) as writer:
        for line in matrix_data[1:]:
            tpm = [round(line[n + num_id_cols] / all_rpk[n], 5) if all_rpk[n] > 0 else 0
                   for n in range(num_samples)]
            writer.writeCountsRow(CountsRow(line[0], line[1], tpm))
    return empty


def main():
    args = parse_args()
    with cli.ErrorHandler():
        with fileOps.opengz(args.counts_matrix) as counts_matrix:
            counts_to_tpm(counts_matrix, args.tpm_matrix, args.isoform_sizes)


if __name__ == '__main__':
    main()
