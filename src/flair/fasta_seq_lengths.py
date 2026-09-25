#!/usr/bin/env python3
import argparse
import csv
import os
from flair.pycbio.sys import cli

def _write_seq_length(writer, length_frequencies, name, seqlen):
    "one sequence's length row, and its contribution to the histogram"
    if name is not None:
        writer.writerow([name, seqlen])
        length_frequencies[seqlen] = length_frequencies.get(seqlen, 0) + 1


def build_parser():
    desc = "Write the length of every sequence in a FASTA, and optionally a histogram of those lengths"
    parser = argparse.ArgumentParser(prog='fasta_seq_lengths', description=desc)
    parser.add_argument('fasta', help='FASTA of sequences to measure')
    parser.add_argument('lengths_tsv', help='output TSV of sequence name and length')
    parser.add_argument('histogram_tsv', nargs='?', default='',
                        help='optional output TSV of length and the number of sequences of that length')
    return parser

def fasta_seq_lengths(fasta_file, lengths_tsv, histogram_tsv):
    fasta = open(fasta_file)
    outfilename = lengths_tsv
    outfilename2 = histogram_tsv
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


def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        fasta_seq_lengths(args.fasta, args.lengths_tsv, args.histogram_tsv)


if __name__ == "__main__":
    main()
