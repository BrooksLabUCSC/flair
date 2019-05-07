import sys, csv

try:
    chromsizesfile = open(sys.argv[1])
    bed = open(sys.argv[2])
    outfilename = sys.argv[3]
except:
    print('usage: script.py chromsizes bedfile pslfile')
    sys.exit()

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
        qstart, qend = starts[0], starts[-1]+sizes[-1]
        blockstarts = ','.join([str(int(start)+relstart) for relstart in starts]) + ','
        qblockstarts = ''
        qbs = 0
        for s in sizes:
            qblockstarts += str(qbs) + ','
            qbs += s
        pslline = [0] * 8
        pslline += [strand, name, qend-qstart, qstart, qend]
        pslline += [chrom, chromsizes[chrom], start, end]
        pslline += [blocknum, blocksizes, qblockstarts, blockstarts]
        writer.writerow(pslline)