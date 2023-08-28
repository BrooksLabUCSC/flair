#!/usr/bin/env python3

import sys
import os
import csv
import mappy as mm

def main():
    readmap = sys.argv[1]
    query = sys.argv[2]
    support = float(sys.argv[3])
    output = sys.argv[4]
    isbed = sys.argv[2][-3:].lower() != 'psl'
    fle = sys.argv[5:]
    subset_unassigned_reads(readmap=readmap, query=query, support=support, output=output, fastx=fle, outfa='/dev/stdout', isbed=isbed)

def subset_unassigned_reads(readmap, query, support, output, fastx, outfa, isbed=True):
    assigned_names = set()
    for line in open(readmap):  # map
        iso, reads = line.rstrip().split('\t')
        reads = reads.split(',')
        # if len(reads) < support:
        #   continue
        for r in reads:
            assigned_names.add(r)

    headers_keep = set()
    
    with open(output, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        for line in open(query):  # bed
            line = line.rstrip().split('\t')
            if isbed:
                name = line[3]
            else:
                name = line[9]

            if name not in assigned_names:
                writer.writerow(line)
                headers_keep.add(name)

    headers_used = set()
    with open(outfa, 'w') as outf:
        for fle in fastx:
            for read in mm.fastx_read(fle):
                header, seq, qual = read
                if header in headers_keep:
                    print('>'+header, file=outf)
                    print(seq, file=outf)
                    headers_used.add(header)
    outf.close()

    diff = len(headers_keep - headers_used)
    if diff > 0:
        sys.stderr.write('{} names do not match any names in fastq file(s)'.format(diff))
        sys.stderr.write('e.g. {} in bed but not in fastq\n'.format(list(headers_keep - headers_used)[0]))
        sys.exit(1)

if __name__ == "__main__":
    main()
