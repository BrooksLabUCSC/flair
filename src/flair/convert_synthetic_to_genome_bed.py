
import sys, os


isoformsbed = sys.argv[1]
readmapfile = sys.argv[2]
file = sys.argv[3]
breakpointfile = sys.argv[4]
outname = sys.argv[5]


isoreadsup = {}
# freadsfinal = set()
for line in open(readmapfile):
    line = line.split('\t')
    reads = line[1].split(',')
    isoreadsup[line[0]] = set(reads)#len(reads)
    # if len(reads) > 1:
    #     freadsfinal.update(reads)
# print('done reading read map')

freadsfinal = set()
finalisos = set()
synthchrtoinfo = {}
for line in open(breakpointfile):
    line = line.rstrip().split('\t')
    synthchrtoinfo[line[0]] = '--'.join(line[-1].split('--')[1:])

grouptogenes = {}
genetoparalogs = {}
# genetoname = {}
for line in open(os.path.realpath(__file__).split('convert_')[0] + 'dgd_Hsa_all_v71.tsv'):
    line = line.split('\t')
    group = line[1]
    gname = line[6]
    hname = line[7]
    # genetoname[gname] = hname
    if group not in grouptogenes: grouptogenes[group] = set()
    grouptogenes[group].add(gname)
    #5	244	2	69321074	69338940	1	ENSG00000205572	SERF1B	small EDRK-rich factor 1B (centromeric) [Source:HGNC Symbol;Acc:10756]
for group in grouptogenes:
    gset = grouptogenes[group]
    for g in gset:
        genetoparalogs[g] = group

genetoname = {}
for line in open('/private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf'):
    if line[0] != '#':
        line = line.split('\t', 3)
        if line[2] == 'gene':
            geneid = line[-1].split('gene_id "')[1].split('"')[0]
            genename = line[-1].split('gene_name "')[1].split('"')[0]
            genetoname[geneid.split('.')[0]] = genename
            # print(geneid.split('.')[0], genename)



locustopartners = {}
for line in open(isoformsbed):
    fnames = set(line.split('\t',1)[0].split('--'))
    fnames = {x.split('.')[0] for x in fnames}
    for i in fnames:
        other = fnames - {i,}
        newother = frozenset([genetoparalogs[g] if g in genetoparalogs else g for g in other])
        # if any(g in genetoparalogs for g in other):
        #     print(i, newother)
        if i not in locustopartners: locustopartners[i] = set()
        locustopartners[i].add(newother)

maxpromiscuity = 4

out = open(outname, 'w')
for line in open(isoformsbed):
    line = line.split('\t')
    iso, start, esizes, estarts = line[3], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
    fusionchr = line[0]
    fgenes, ispromiscuous = set([x.split('.')[0] for x in fusionchr.split('--')]), False
    areparalogs, areig, allnames = [], [], set()
    for g in fgenes:
        if len(locustopartners[g]) > maxpromiscuity:
            ispromiscuous = True
            # print(g, locustopartners[g])
        if g in genetoname:
            allnames.add(genetoname[g][:4])
            if genetoname[g][:3] in {'IGK', 'IGH', 'IGL', 'IGF'}: areig.append(True)
            else: areig.append(False)
        else:
            allnames.add(g)
            if g[:3] == 'chr':areig.append(True)
            else: areig.append(False)
        if g in genetoparalogs:
            other = fgenes - {g,}
            other = {genetoparalogs[g2] if g2 in genetoparalogs else g2 for g2 in other}
            # if 'ENSG00000211640' in fgenes:  print(g, genetoparalogs[g], other)
            # if g == 'ENSG00000211640': print(g, genetoparalogs[g], other)
            if genetoparalogs[g] in other or len(other) == 0: areparalogs.append(True)
            else: areparalogs.append(False)
        else:
            areparalogs.append(False)
    # if 'ENSG00000211640' in fgenes: print(areparalogs) #, locustopartners['ENSG00000211640'])
        # if g in genetoparalogs:
        #     paraset = grouptogenes[genetoparalogs[g]]
        #     for g2 in fgenes - {g,}:
        #         if g2 in paraset:
        #             areparalogs.append(True)
        #             break
        # else: areparalogs.append(False)
    synthinfo = [x.split('..') for x in synthchrtoinfo[fusionchr].split('--')]
    synthinfo = [[y[0], y[1], int(y[2]), int(y[3])] for y in synthinfo]
    if len(isoreadsup[iso]) >= 1 and not ispromiscuous and any(x==False for x in areparalogs) and any(x==False for x in areig) and len(allnames) > 1 and 'chrM' not in [x[1] for x in synthinfo]:
        # synthinfo = [x.split('..') for x in synthchrtoinfo[fusionchr].split('--')]
        # synthinfo = [[y[0], y[1], int(y[2]), int(y[3])] for y in synthinfo]

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

            # finalisos.add(iso)
            ###VERY IMPORTANT


            genomicbounds = []
            outlines = []
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
                    genomicbounds.append((genomicchr, leftbound - ((start+ starts[order])-locusbounds[order][0]), leftbound - (((start+ starts[order])-locusbounds[order][0]) + totlen), locusdir))
                else:
                    outline = [genomicchr, str(leftbound + ((start+ starts[order])-locusbounds[order][0])), str(leftbound + ((start+ starts[order])-locusbounds[order][0]) + totlen), 'gene' + str(order+1) + '_' + iso,
                           '1000', locusdir, str(leftbound + ((start+ starts[order])-locusbounds[order][0])), str(leftbound + ((start+ starts[order])-locusbounds[order][0]) + totlen), '0', str(len(exons)),
                                                  ','.join([str(x) for x in locusesizes]),
                                                  ','.join([str(x) for x in locusestarts])]
                    genomicbounds.append((genomicchr, leftbound + ((start+ starts[order])-locusbounds[order][0]), leftbound + ((start+ starts[order])-locusbounds[order][0]) + totlen, locusdir))
                # out.write('\t'.join(outline) + '\n')
                outlines.append('\t'.join(outline) + '\n')
            istooclose = False
            for i in range(numloci - 1):
                if genomicbounds[i][0] == genomicbounds[i+1][0] and genomicbounds[i][3] == genomicbounds[i+1][3] and abs(genomicbounds[i][2] - genomicbounds[i+1][1]) < 350000:
                    # print(genomicbounds)
                    istooclose = True
            if istooclose: continue
            for l in outlines:
                out.write(l)
            print(iso)
            freadsfinal.update(isoreadsup[iso])
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
