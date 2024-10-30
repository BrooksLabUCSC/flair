
import sys


isoformsbed = sys.argv[1]
readmapfile = sys.argv[2]
file = sys.argv[3]
breakpointfile = sys.argv[4]
outname = sys.argv[5]


isoreadsup = {}
freadsfinal = set()
for line in open(readmapfile):
    line = line.split('\t')
    reads = line[1].split(',')
    isoreadsup[line[0]] = len(reads)
    if len(reads) > 1:
        freadsfinal.update(reads)
# print('done reading read map')








last = False
readsfile = open(file)
temp = file.split('.')
temp[-2] = 'fusionreads'
freads = open('.'.join(temp), 'w')

if file.split('.')[-1] == 'fasta' or file.split('.')[-1] == 'fa':
    for line in readsfile:
        if line[0] == '>':
            readname = line[1:].rstrip().split()[0]
            if readname in freadsfinal: last = True
            else: last = False
        if last: freads.write(line)
else:
    linecount = 0
    for line in readsfile:
        if linecount % 4 == 0:
            readname = line[1:].rstrip().split()[0]
            if readname in freadsfinal: last = True
            else: last = False
        if last: freads.write(line)
        linecount += 1
freads.close()

# if readsfasta.split('.')[-1] == 'fastq' or readsfasta.split('.')[-1] == 'fq':
#     out5 = open('.'.join(readsfasta.split('.')[:-1]) + '-isoSupport.fastq', 'w')
#     last = False
#     for line in open(readsfasta):
#         if line[0] == '@':
#             if line.lstrip('@').split(' ')[0] in freadsfinal:
#                 last = True
#             else:
#                 last = False
#         if last: out5.write(line)
# elif readsfasta.split('.')[-1] == 'fasta' or readsfasta.split('.')[-1] == 'fa':
#     out5 = open('.'.join(readsfasta.split('.')[:-1]) + '-isoSupport.fasta', 'w')
#     last = False
#     for line in open(readsfasta):
#         if line[0] == '>':
#             if line.rstrip().lstrip('>').split(' ')[0] in freadsfinal:
#                 last = True
#             else:
#                 last = False
#         if last: out5.write(line)

synthchrtoinfo = {}
for line in open(breakpointfile):
    line = line.rstrip().split('\t')
    synthchrtoinfo[line[0]] = '--'.join(line[-1].split('--')[1:])

out = open(outname, 'w')
for line in open(isoformsbed):
    line = line.split('\t')
    iso, start, esizes, estarts = line[3], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
    fusionchr = line[0]
    if isoreadsup[iso] >= 1:
        synthinfo = [x.split('..') for x in synthchrtoinfo[fusionchr].split('--')]
        synthinfo = [[y[0], y[1], int(y[2]), int(y[3])] for y in synthinfo]

        locuslen = [abs(x[3] - x[2]) for x in synthinfo]
        locusbounds = []
        laststart = 0
        for i in range(len(locuslen)):
            locusbounds.append((laststart, locuslen[i] + laststart))
            laststart += locuslen[i]

        numloci = len(synthinfo)
        introns, exons = [[] for x in range(numloci)], [[] for x in range(numloci)]
        starts = [None for x in range(numloci)]
        exonindexes = [[] for x in range(numloci)]
        lastexonend = 0
        if start < locusbounds[0][1] and locusbounds[-1][0] < int(line[2]):
            for i in range(len(esizes)):
                thisstart, thisend = start + estarts[i], start + estarts[i] + esizes[i]

                for order in range(numloci):
                    if locusbounds[order][0] < thisstart and thisend <= locusbounds[order][1]:
                        if starts[order] == None:
                            starts[order] = estarts[i]#thisstart #- locusbounds[order][0]
                        exonindexes[order].append(i)
                        # else: introns[order].append(estarts[i]-lastexonend)
        #                 exons[order].append(esizes[i])
        #
        #                 lastexonend = estarts[i] + esizes[i]
        #                 break
            if None in starts:
                # print(iso, 'not able to be converted to genomic coordinates')
                continue
        #     # print('exons', exons)
        #     # print('introns', introns)
            for order in range(numloci):
                genename, genomicchr, leftbound, rightbound = synthinfo[order]
                locusesizes = [esizes[i] for i in exonindexes[order]]
                locusestarts = [estarts[i] for i in exonindexes[order]]
                locusestarts = [x - starts[order] for x in locusestarts]
                totlen = locusestarts[-1] + locusesizes[-1]
                locusdir = '+' if leftbound < rightbound else '-'
                if leftbound > rightbound: #reverse direction
                    temp = []
                    for i in range(len(locusestarts)-1, -1, -1):
                        temp.append(totlen - (locusestarts[i] + locusesizes[i]))
                    locusestarts = temp
                    locusesizes = locusesizes[::-1]
                    outline = [genomicchr, str(leftbound - (((start+ starts[order])-locusbounds[order][0]) + totlen)), str(leftbound - ((start+ starts[order])-locusbounds[order][0])), 'gene' + str(order+1) + '_' + iso,
                               '1000', locusdir, str(leftbound - (start+totlen)), str(leftbound - start), '0',
                               str(len(exons)),
                               ','.join([str(x) for x in locusesizes]),
                               ','.join([str(x) for x in locusestarts])]
                else:
                    outline = [genomicchr, str(leftbound + ((start+ starts[order])-locusbounds[order][0])), str(leftbound + ((start+ starts[order])-locusbounds[order][0]) + totlen), 'gene' + str(order+1) + '_' + iso,
                           '1000', locusdir, str(leftbound + ((start+ starts[order])-locusbounds[order][0])), str(leftbound + ((start+ starts[order])-locusbounds[order][0]) + totlen), '0', str(len(exons)),
                                                  ','.join([str(x) for x in locusesizes]),
                                                  ','.join([str(x) for x in locusestarts])]
                out.write('\t'.join(outline) + '\n')
                # totlen, geneestarts = 0, [0, ]
                # print(len(introns[order]))
        #         if leftbound < rightbound:
        #             for i in range(len(introns[order])):
        #                 geneestarts.append(totlen + introns[order][i])
        #                 totlen += introns[order][i] + exons[order][i]
        #             outline = [genomicchr, str(leftbound+starts[order]), str(leftbound + starts[order] + totlen), iso, '1000', '+',
        #                        str(leftbound+starts[order]), str(leftbound + starts[order] + totlen), '0', str(len(exons)),
        #                        ','.join([str(x) for x in exons[order]]),
        #                        ','.join([str(x) for x in geneestarts])]
        #         else:
        #             for i in range(len(introns[order])-1, -1, -1):
        #                 geneestarts.append(totlen + introns[order][i])
        #                 totlen += introns[order][i] + exons[order][i]
        #             outline = [genomicchr, str(leftbound - (starts[order]+totlen)), str(leftbound - starts[order]), iso, '1000', '-',
        #                        str(leftbound - (starts[order]+totlen)), str(leftbound - starts[order]), '0', str(len(exons)),
        #                        ','.join([str(x) for x in exons[order][::-1]]),
        #                        ','.join([str(x) for x in geneestarts])]
        #         out.write('\t'.join(outline) + '\n')
out.close()
#
# out = open('test.txt', 'w')
# out.write('hi')
