#!/usr/bin/env python3

import sys
from flair import FlairInputDataError
from flair.counts_matrix_tsv import ID_COLUMNS, CountsRow
from flair.seq_lengths_tsv import SeqLengthsReader
from flair.tpm_tsv import tpm_writer
from flair.pycbio.sys import cli
from flair.pycbio.sys import fileOps

def parse_input():
    try:
        counts_matrix = fileOps.opengz(sys.argv[1])
        outfilename = sys.argv[2]
        sizefile = sys.argv[3] if len(sys.argv) > 3 else None
    except Exception:
        raise FlairInputDataError('usage: counts_to_tpm.py counts_matrix.tsv count_matrix.tpm.tsv [iso.sizes]\n'
                                  'convenience script for obtaining a file of isoform sizes: bin/fasta_seq_lengths.py\n'
                                  'if no isoform size file is provided, no length normalization will be done (just reads per million)\n')
    return counts_matrix, outfilename, sizefile

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
    counts_matrix, outfilename, sizefile = parse_input()
    with cli.ErrorHandler():
        counts_to_tpm(counts_matrix, outfilename, sizefile)


if __name__ == '__main__':
    main()
