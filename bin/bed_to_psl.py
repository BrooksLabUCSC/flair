import sys, csv

try:
    chromsizesfile = open(sys.argv[1])
    bed = open(sys.argv[2])
    outfilename = sys.argv[3]
except:
    print('usage: script.py bedfile chromsizes pslfile')
    sys.exit()

def get_junctions_psl(starts, sizes):
    junctions = set()
    if len(starts) != 1:
        for b in range(len(starts)-1):
            junctions.add((starts[b]+sizes[b], starts[b+1]))
        return junctions

chromsizes = {}
for line in chromsizesfile:
    line = line.rstrip().split('\t')
    chromsizes[line[0]] = line[1]

with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    for line in bed:
        line = line.rstrip().split('\t')
        chrom, start, end, name, score, strand = line[:6]
        blocknum, blocksizes, relblockstarts = line[9:]
        sizes = [int(n) for n in blocksizes.split(',')[:-1]]
        starts = [int(n) for n in relblockstarts.split(',')[:-1]]
        blockstarts = ','.join([str(int(start)+relstart) for relstart in starts]) + ','
        pslline = [0] * 8
        pslline += [strand, name, sum(sizes), 0, sum(sizes)]
        pslline += [chrom, chromsizes[chrom], start, end]
        pslline += [blocknum, blocksizes, relblockstarts, blockstarts]
        writer.writerow(pslline)