import sys
from collections import Counter
import pysam

alignedbedfile = sys.argv[1]
referencegtffile = sys.argv[2]
outfilename = sys.argv[3]
refbpfile = sys.argv[4]
sjwiggle = int(sys.argv[5]) #15
readcov = int(sys.argv[6]) #2
refgenomefile = sys.argv[7]

def grouper(iterable):
    prev = None
    group = []
    for item in iterable:
        if prev is None or item - prev <= sjwiggle:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group

fusiontoannotsj = {}
for line in open(referencegtffile):
    line = line.rstrip().split('\t', 5)
    if line[2] == 'exon':
        if line[0] not in fusiontoannotsj: fusiontoannotsj[line[0]] = set()
        fusiontoannotsj[line[0]].add(int(line[3]))
        fusiontoannotsj[line[0]].add(int(line[4]))

fusiontobp = {}
for line in open(refbpfile):
    line = line.rstrip().split('\t')
    fusiontobp[line[0]] = int(line[1])


genome = pysam.FastaFile(refgenomefile)


splicejunctosupport = {}
chrtonovelss = {}
reconsideredss = []
c = 0
for line in open(alignedbedfile):
    line = line.rstrip().split('\t')
    thischr, iso, dir, start, esizes, estarts = line[0], line[3], line[5], int(line[1]), \
                            [int(x) for x in line[10].rstrip(',').split(',')], [ int(x) for x in line[11].rstrip(',').split(',')]
    for i in range(len(esizes) - 1):
        thisintron = [thischr, start + estarts[i] + esizes[i] , start + estarts[i + 1] + 1]
        ###only process around bp
        # print(iso, thisintron[1], genome.fetch(thischr, thisintron[1], thisintron[1]+2))
        # print(iso, thisintron[2], genome.fetch(thischr, thisintron[2]-2, thisintron[2]))
        if thisintron[1] <= fusiontobp[thischr] <= thisintron[2] or \
                genome.fetch(thischr, thisintron[1], thisintron[1]+2) == 'GT' and \
                genome.fetch(thischr, thisintron[2]-3, thisintron[2]-1) == 'AG':
            doreconsider = False
            if thischr in fusiontoannotsj:
                for i in range(1,3):#, 3):
                    closestdist, closestpos = 1000, None
                    for sj in fusiontoannotsj[thischr]:
                        thisdist = abs(thisintron[i] - sj)
                        if thisdist <= sjwiggle:#== 0:
                            # print(iso, thisintron[i], genome.fetch(thischr, thisintron[i], thisintron[i]+2))
                            # if i == 1: print(iso, thisintron[i], genome.fetch(thischr, thisintron[i], thisintron[i]+2))
                            # else: print(iso, thisintron[i], genome.fetch(thischr, thisintron[i]-3, thisintron[i]-1))
                            if thisdist < closestdist: closestdist, closestpos = thisdist, sj
                    if not closestpos: doreconsider = True
            if doreconsider:
                if thischr not in chrtonovelss: chrtonovelss[thischr] = []
                chrtonovelss[thischr].append(thisintron[1])
                chrtonovelss[thischr].append(thisintron[2])
                reconsideredss.append(thisintron)
                    # if thisdist <= sjwiggle:
#         if thisintron[1] <= fusiontobp[thischr] <= thisintron[2] or
#         # if True: #thisintron[1] <= fusiontobp[thischr] <= thisintron[2]:
#             doreconsider = False
#             newintron = [thischr, None, None]
#             for i in range(1,3): ##checking each splice site
#                 ###FIXME currently processing each splice site individually
#                 closestdist, closestpos = 1000, None
#                 if thischr in fusiontoannotsj:
#                     for sj in fusiontoannotsj[thischr]:
#                         thisdist = abs(thisintron[i]-sj)
#                         if thisdist <= sjwiggle:
#                             if thisdist < closestdist: closestdist, closestpos = thisdist, sj
#                 if closestpos: newintron[i] = closestpos
#                 else:
#                     # if thischr not in chrtonovelss: chrtonovelss[thischr] = []
#                     # chrtonovelss[thischr].append(thisintron[i])
#                     doreconsider = True
#             if doreconsider:
#                 if thischr not in chrtonovelss: chrtonovelss[thischr] = []
#                 chrtonovelss[thischr].append(thisintron[1])
#                 chrtonovelss[thischr].append(thisintron[2])
#                 reconsideredss.append(thisintron)
#         # if thischr == 'ENSG00000146872.19--ENSG00000172354.10': print(doreconsider, chrtonovelss[thischr])
#         c += 1
#         # if c < 5: print(doreconsider, thisintron, chrtonovelss[thischr])
#         # thisintron[2] -= 1
#         # thisintron = tuple(thisintron)
#         # if thisintron not in splicejunctosupport: splicejunctosupport[thisintron] = 0
#         # splicejunctosupport[thisintron] += 1
#
# # print('novel', chrtonovelss['ENSG00000146872.19--ENSG00000172354.10'])
#
chrtogoodss = {}
for chr in chrtonovelss:
    chrtogoodss[chr] = set()
    chrtonovelss[chr] = sorted(chrtonovelss[chr])
    for group in grouper(chrtonovelss[chr]):
        # if chr == 'ENSG00000211640.4--ENSG00000211677.2': print('group', group)
        if len(group) >= readcov:
            groupsize = len(group)
            finalpos = []
            for pos, count in Counter(group).most_common():
                # if chr == 'ENSG00000211640.4--ENSG00000211677.2': print(pos, count)
                if count >= readcov and count >= groupsize/10: ###needs to also be 1/10 of locus
                    isdifferent = True
                    for fp in finalpos:
                        if fp-sjwiggle <= pos <= fp+sjwiggle: isdifferent = False
                    if isdifferent: finalpos.append(pos)
            # if chr == 'ENSG00000211640.4--ENSG00000211677.2': print('final', finalpos)
            for fp in finalpos:
                chrtogoodss[chr].add(fp)
# print('good', chrtogoodss['ENSG00000146872.19--ENSG00000172354.10'])

for thisintron in reconsideredss:
    thischr = thisintron[0]
    for i in range(1, 3):  ##checking each splice site
        ###FIXME currently processing each splice site individually
        closestdist, closestpos = 1000, None
        if thischr in chrtogoodss:
            for sj in chrtogoodss[thischr]:
                thisdist = abs(thisintron[i] - sj)
                if thisdist <= sjwiggle:
                    if thisdist < closestdist: closestdist, closestpos = thisdist, sj
        if closestpos: thisintron[i] = closestpos
    thisintron[2] -= 1
    thisintron = tuple(thisintron)
    if thisintron not in splicejunctosupport: splicejunctosupport[thisintron] = 0
    splicejunctosupport[thisintron] += 1


goodsj, lowsupsj = {}, {}

out = open(outfilename, 'w')
for j in splicejunctosupport:
    if splicejunctosupport[j] >= 2:
        # goodsj[j] = splicejunctosupport[j]
        out.write('\t'.join([str(x) for x in j]) + '\t+\n')
    # else: lowsupsj[j] = splicejunctosupport[j]
out.close()

# for chr, start, end in lowsupsj:
#     for i in range(-2, 3):
#         for j in range(-2, 3):
#             if (chr, start+i, end+j) in goodsj:
#                 goodsj[(chr, start+i, end+j)] += lowsupsj[(chr, start, end)]


# out = open(outfilename, 'w')
# for j in goodsj:
#     # if j[0] == 'TP53*ENSG00000141510.18.chr17.7671131.7688547--chr6-30221000.chr6.30223661.30220450':
#     #     print(goodsj[j], j)
#     if goodsj[j] >= 2:
#         out.write('\t'.join([str(x) for x in j]) + '\t+\n')
# out.close()











# fusiontoannotsj = {}
# for line in open(referencegtffile):
#     line = line.rstrip().split('\t', 5)
#     if line[2] == 'exon':
#         if line[0] not in fusiontoannotsj: fusiontoannotsj[line[0]] = set()
#         fusiontoannotsj[line[0]].add(int(line[3]))
#         fusiontoannotsj[line[0]].add(int(line[4]))
#
# splicejunctosupport = {}
# for line in open(alignedbedfile):
#     line = line.rstrip().split('\t')
#     thischr, iso, dir, start, esizes, estarts = line[0], line[3], line[5], int(line[1]), \
#                             [int(x) for x in line[10].rstrip(',').split(',')], [ int(x) for x in line[11].rstrip(',').split(',')]
#     for i in range(len(esizes) - 1):
#         thisintron = [thischr, start + estarts[i] + esizes[i] , start + estarts[i + 1] + 1]
#         if thischr in fusiontoannotsj:
#             for i in range(1,3):
#                 for sj in fusiontoannotsj[thischr]:
#                     if abs(thisintron[i]-sj) <= 3:
#                         thisintron[i] = sj
#                         break
#         thisintron[2] -= 1
#         thisintron = tuple(thisintron)
#         if thisintron not in splicejunctosupport: splicejunctosupport[thisintron] = 0
#         splicejunctosupport[thisintron] += 1
#
# goodsj, lowsupsj = {}, {}
#
# for j in splicejunctosupport:
#     if splicejunctosupport[j] >= 2:
#         goodsj[j] = splicejunctosupport[j]
#     else: lowsupsj[j] = splicejunctosupport[j]
#
#
# for chr, start, end in lowsupsj:
#     for i in range(-2, 3):
#         for j in range(-2, 3):
#             if (chr, start+i, end+j) in goodsj:
#                 goodsj[(chr, start+i, end+j)] += lowsupsj[(chr, start, end)]
#
#
#
#
# out = open(outfilename, 'w')
# for j in goodsj:
#     # if j[0] == 'TP53*ENSG00000141510.18.chr17.7671131.7688547--chr6-30221000.chr6.30223661.30220450':
#     #     print(goodsj[j], j)
#     if goodsj[j] >= 2:
#         out.write('\t'.join([str(x) for x in j]) + '\t+\n')
# out.close()