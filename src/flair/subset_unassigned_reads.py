#!/usr/bin/env python3

import sys, os, csv
s = float(sys.argv[3])
assigned_names = set()
for line in open(sys.argv[1]):  # map
    iso, reads = line.rstrip().split('\t')
    reads = reads.split(',')
    # if len(reads) < s:
    #   continue
    for r in reads:
        assigned_names.add(r)

headers_keep = set()
isbed = sys.argv[2][-3:].lower() != 'psl'
with open(sys.argv[4], 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    for line in open(sys.argv[2]):  # bed
        line = line.rstrip().split('\t')
        if isbed:
            name = line[3][:line[3].rfind(';')]
        else:
            name = line[9][:line[9].rfind(';')]

        if name not in assigned_names:
            writer.writerow(line)
            headers_keep.add(name)

import mappy as mm

headers_used = set()
for fle in sys.argv[5:]:
    for read in mm.fastx_read(fle):
        header, seq, qual = read
        if header in headers_keep:
            print('>'+header)
            print(seq)
            headers_used.add(header)

diff = len(headers_keep -  headers_used)
if diff > 0:
    sys.stderr.write('{} names do not match any names in fastq file(s)'.format(diff))
    sys.stderr.write('e.g. {} in bed but not in fastq\n'.format(list(headers_keep - headers_used)[0]))
    sys.exit(1)
