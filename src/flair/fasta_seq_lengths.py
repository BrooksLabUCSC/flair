#!/usr/bin/env python3
import sys
import csv
import os
from flair import FlairInputDataError

# FIXME: use argparse

def _write_seq_length(writer, length_frequencies, name, seqlen):
    "one sequence's length row, and its contribution to the histogram"
    if name is not None:
        writer.writerow([name, seqlen])
        length_frequencies[seqlen] = length_frequencies.get(seqlen, 0) + 1


def main():
    try:
        fasta = open(sys.argv[1])
        outfilename = sys.argv[2]
        if len(sys.argv) > 3:
            outfilename2 = sys.argv[3]
        else:
            outfilename2 = ''
    except Exception:
        raise FlairInputDataError('usage: fasta_seq_lengths fasta outfilename [outfilename2]\n')

    length_frequencies = {}
    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        seqlen = 0
        name = None
        for line in fasta:
            line = line.rstrip()
            if line.startswith('>'):
                _write_seq_length(writer, length_frequencies, name, seqlen)
                name = line[1:]
                seqlen = 0
            else:
                seqlen += len(line)
        # the last sequence: it used to be written without being counted, and an
        # empty file wrote a row of None 0
        _write_seq_length(writer, length_frequencies, name, seqlen)

    if outfilename2:
        alllengths = sorted(length_frequencies.keys())
        with open(outfilename2, 'wt') as outfile:
            writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
            for length in alllengths:
                writer.writerow([length, length_frequencies[length]])


if __name__ == "__main__":
    main()
