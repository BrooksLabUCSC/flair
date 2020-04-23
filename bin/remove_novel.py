#!/usr/bin/env python3
import sys, csv

try:
    gp = open(sys.argv[1])
    psl = open(sys.argv[2])
    outpsl = sys.argv[3]
except:
    sys.stderr.write('script.py junctions.gp psl nonovels.psl \n')
    sys.exit(1)

annotated = {}
for line in gp:
    line = line.rstrip().split('\t')
    chrom = line[1]
    if chrom not in annotated:
        annotated[chrom] = set()
    blockstarts = [int(n) for n in line[8].split(',')[:-1]][1:]
    blockends = [int(n) for n in line[9].split(',')[:-1]][:-1]
    for start, end in zip(blockstarts, blockends):
        annotated[chrom].add((end, start))

lastchrom, lastjunc = '', ''
notfound = set()
with open(outpsl, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    for line in psl:
        line = line.rstrip().split('\t')
        strand, chrom = line[8], line[13]
        if chrom not in annotated:
            continue
        starts = [int(n) for n in line[20].split(',')[:-1]]
        sizes = [int(n) for n in line[18].split(',')[:-1]]      
        validjuncs = True
        for b in range(len(starts)-1):
            junction = (starts[b]+sizes[b], starts[b+1])  # accounting for lots of indexing errors
            junction2 = (junction[0]-1, junction[1])
            junction3 = (junction[0]+1, junction[1])
            if junction not in annotated[chrom] and junction2 not in annotated[chrom] and junction3 not in annotated[chrom]:
                validjuncs = False
                lastchrom = chrom
                lastjunc = junction3
                notfound.add((junction, chrom, strand))
                break
        if validjuncs:
            writer.writerow(line)