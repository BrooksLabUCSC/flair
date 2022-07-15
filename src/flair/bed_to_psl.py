#!/usr/bin/env python3
import sys, csv, os

try:
    chromsizesfile = open(sys.argv[1])
    bed = open(sys.argv[2])
    outfilename = sys.argv[3]
except:
    sys.stderr.write('usage: bed_to_psl.py chromsizes bedfile pslfile\n')
    sys.exit(1)

chromsizes = {}
for line in chromsizesfile:
    line = line.rstrip().split('\t')
    chromsizes[line[0]] = line[1]

with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    for line in bed:
        line = line.rstrip().split('\t')
        chrom, start, end, name, score, strand = line[:6]
        blocknum, blocksizes, relblockstarts = line[9:12]
        sizes = [int(n) for n in blocksizes.split(',')[:-1]]
        starts = [int(n) for n in relblockstarts.split(',')[:-1]]
        blockstarts = ','.join([str(int(start)+relstart) for relstart in starts]) + ','
        
        qblockstarts, qbs = '', 0
        for s in sizes:
            qblockstarts += str(qbs) + ','
            qbs += s
        qstart, qend = starts[0], qbs

        tNumInsert, tBaseInsert = 0, 0
        for i in range(len(starts)-1):
            gap_size = starts[i+1] - (starts[i]+sizes[i])
            if gap_size:
                tNumInsert += 1
                tBaseInsert += gap_size

        pslline = [sum(sizes)] + [0]*5 + [tNumInsert, tBaseInsert]
        pslline += [strand, name, qend-qstart, qstart, qend]
        if chrom not in chromsizes:
            pslline += [chrom, 0, start, end]
        else:
            pslline += [chrom, chromsizes[chrom], start, end]
        pslline += [blocknum, blocksizes, qblockstarts, blockstarts]
        pslline += line[12:]
        writer.writerow(pslline)