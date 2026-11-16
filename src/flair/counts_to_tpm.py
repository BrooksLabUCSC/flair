#!/usr/bin/env python3

import sys
import csv
import os
from flair import FlairInputDataError


def read_isoform_sizes(sizefile):
    "the length of every isoform, from the file fasta_seq_lengths writes"
    sizes = {}
    with open(sizefile) as fh:
        for line in fh:
            name, length = line.rstrip().split('\t')
            sizes[name] = float(length)
    return sizes


def _isoform_size(isoform_sizes, isoform_id):
    try:
        return isoform_sizes[isoform_id]
    except KeyError as ex:
        raise FlairInputDataError(f'no length for isoform {isoform_id}; a length is needed for every isoform of '
                                  f'the counts matrix.  From fasta_seq_lengths the names are the whole FASTA '
                                  f'header line, so that FASTA must be the one the counts matrix names, with no '
                                  f'description after the id') from ex


def counts_to_tpm(counts_matrix, outfilename, isoform_sizes):
    """Write a counts matrix as transcripts per million, isoform_sizes giving the
    length in bases of each isoform of the matrix."""
    header = counts_matrix.readline().rstrip().split('\t')
    num_samples = len(header[1:])
    matrix_data = [header]
    all_rpk = [0] * num_samples
    for line in counts_matrix:
        line = line.rstrip().split('\t')
        isoform_id = line[0]
        size = _isoform_size(isoform_sizes, isoform_id)
        rpk = [float(count) / size for count in line[1:]]
        for n in range(num_samples):
            all_rpk[n] += rpk[n]
        matrix_data += [[isoform_id] + rpk]

    all_rpk = [rpk / 1e6 for rpk in all_rpk]

    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        writer.writerow(matrix_data[0])
        for line in matrix_data[1:]:
            for n in range(num_samples):
                line[n + 1] = line[n + 1] / all_rpk[n]
            writer.writerow(line)


def main():
    if len(sys.argv) != 4:
        raise FlairInputDataError('usage: counts_to_tpm.py counts_matrix.tsv count_matrix.tpm.tsv iso.sizes\n'
                                  'convenience script for obtaining a file of isoform sizes: fasta_seq_lengths\n')
    counts_matrix, outfilename, sizefile = sys.argv[1:4]
    with open(counts_matrix) as fh:
        counts_to_tpm(fh, outfilename, read_isoform_sizes(sizefile))


if __name__ == '__main__':
    main()
