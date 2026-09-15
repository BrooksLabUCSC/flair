#!/usr/bin/env python3

import sys
import csv
import os
from flair import FlairInputDataError

def parse_input():
    try:
        counts_matrix = open(sys.argv[1])
        outfilename = sys.argv[2]
        if len(sys.argv) > 3:
            sizefile = open(sys.argv[3])
        else:
            sizefile = None
    except Exception:
        raise FlairInputDataError('usage: counts_to_tpm.py counts_matrix.tsv count_matrix.tpm.tsv [iso.sizes]\n'
                                  'convenience script for obtaining a file of isoform sizes: bin/fasta_seq_lengths.py\n'
                                  'if no isoform size file is provided, no length normalization will be done (just reads per million)\n')
    return counts_matrix, outfilename, sizefile

def counts_to_tpm(counts_matrix, outfilename, sizefile=None):
    sizes = {}
    if sizefile:
        for line in sizefile:
            line = line.rstrip().split('\t')
            sizes[line[0]] = float(line[1])

    header = counts_matrix.readline().rstrip().split('\t')
    num_samples = len(header[1:])
    matrix_data = [header]
    all_rpk = [0] * num_samples
    for line in counts_matrix:
        line = line.rstrip().split('\t')
        isoform_id = line[0]
        if sizes:
            rpk = [float(count) / sizes[isoform_id] for count in line[1:]]
        else:
            rpk = [float(count) for count in line[1:]]
        for n in range(num_samples):
            all_rpk[n] += rpk[n]
        matrix_data += [[isoform_id] + rpk]

    all_rpk = [rpk / 1e6 for rpk in all_rpk]
    empty = [n for n, rpk in enumerate(all_rpk) if rpk == 0]
    if len(empty) > 0:
        raise FlairInputDataError(f"{outfilename}: samples in columns {empty} have no counts in any row, "
                                  "so TPM cannot be computed; drop them from the counts matrix")

    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        writer.writerow(matrix_data[0])
        for line in matrix_data[1:]:
            for n in range(num_samples):
                line[n + 1] = round(line[n + 1] / all_rpk[n], 5)
            writer.writerow(line)


if __name__ == '__main__':
    counts_matrix, outfilename, sizefile = parse_input()
    counts_to_tpm(counts_matrix, outfilename, sizefile)
