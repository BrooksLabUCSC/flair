#!/usr/bin/env python3

import pysam, sys, argparse
from collections import defaultdict
from statistics import median

def def_value():
    return set()

def getCorrectGene(annot, geneannot, chr, readblocks, thisdir=None):
    # print(chr, readblocks)
    if chr not in annot:
        return chr + '-' + str(round(readblocks[0][0], -3))
    intToCheck, geneinttocheck = set(), set()
    for start, stop in readblocks:
        # for i in range(round(start, -1), round(stop, -1), 10):
        for i in range(round(start, -2), round(stop, -2) + 1, 100):
            intToCheck.add(i)
        # for i in range(round(start, -2), round(stop, -2), 100):
        #     geneinttocheck.add(i)
    # thisdir = '-' if is_reverse else '+'
    geneOptions = {}
    # secondaryGeneOptions = {}
    for i in intToCheck:
        if i in annot[chr]:
            # if is_reverse:
                # for g, dir in annot[chr][i]:
                #     if dir == thisdir:
                #         if g not in geneOptions: geneOptions[g] = 0
                #         geneOptions[g] += 1
            # else:
            for g, dir in annot[chr][i]:
                if g not in geneOptions:
                    if not thisdir: geneOptions[g] = 0
                    elif dir == thisdir: geneOptions[g] = 5 ##boost of 50bp for going the same direction
                    else: geneOptions[g] = 0
                geneOptions[g] += 1
    # print(len(intToCheck))
    # print(geneOptions)
    if len(geneOptions.keys()) > 0:
        bestgene, besttot = None, 0
        for g in geneOptions:
            if geneOptions[g] > besttot: bestgene, besttot = g, geneOptions[g]
        if besttot/len(intToCheck) > 0.5:
            return bestgene
    return chr + '-' + str(round(readblocks[0][0], -3))
    # for i in geneinttocheck:
    #     if i in geneannot[chr]:
    #         # if is_reverse:
    #         #     for g, dir in geneannot[chr][i]:
    #         #         if dir == thisdir:
    #         #             if g not in secondaryGeneOptions: secondaryGeneOptions[g] = 0
    #         #             secondaryGeneOptions[g] += 1
    #         # else:
    #         for g, dir in geneannot[chr][i]:
    #             if g not in secondaryGeneOptions: #secondaryGeneOptions[g] = 0
    #                 if not thisdir: secondaryGeneOptions[g] = 0
    #                 elif dir == thisdir: secondaryGeneOptions[g] = 1  ##boost of 100bp for going the same direction
    #                 else: secondaryGeneOptions[g] = 0
    #             secondaryGeneOptions[g] += 1
    # print(len(geneinttocheck))
    # print(secondaryGeneOptions)
    # if len(secondaryGeneOptions.keys()) > 0:
    #     bestgene, besttot = None, 0
    #     for g in secondaryGeneOptions:
    #         if secondaryGeneOptions[g] > besttot: bestgene, besttot = g, secondaryGeneOptions[g]
    #     if besttot / len(geneinttocheck) > 0.5:
    #         return bestgene
    #     else:
    #         return chr + '-' + str(round(readblocks[0][0], -3))
    # else:
    #     return chr + '-' + str(round(readblocks[0][0], -3))

# def gchimparsegtf(gtffilename):
#     annot = {}
#     geneannot = {}
#     genetoinfo = {}
#     for line in open(gtffilename):
#         if line[0] != '#':
#             line = line.split('\t')
#             if line[2] == 'exon':
#                 chr, start, stop, dir = line[0], int(line[3]), int(line[4]), line[6]
#                 genename = line[8].split('gene_id "')[1].split('"')[0]
#                 # genename += '*' + line[8].split('gene_id "')[1].split('"')[0]
#                 if chr not in annot: annot[chr] = defaultdict(def_value)
#                 for i in range(round(start, -1), round(stop, -1), 10):
#                     # if i not in annot[chr]: annot[chr][i] = set()
#                     annot[chr][i].add((genename, dir))
#             elif line[2] == 'gene':
#                 chr, start, stop, dir = line[0], int(line[3]), int(line[4]), line[6]
#                 genename = line[8].split('gene_id "')[1].split('"')[0]
#                 # genename += '*' + line[8].split('gene_id "')[1].split('"')[0]
#                 genetoinfo[genename] = (chr, start, stop, dir)
#                 if chr not in geneannot: geneannot[chr] = defaultdict(def_value)
#                 for i in range(round(start, -2), round(stop, -2), 100):
#                     geneannot[chr][i].add((genename, dir))
#     print('done loading annot')
#     return annot, geneannot, genetoinfo


def idGenomicChimeras(bam, annot, geneannot, genetoinfo, minsup, maxloci=10, reqdisttostart=None):
    isrevtosign = {True: '-', False: '+'}
    withsup = pysam.AlignmentFile(bam, "rb")
    c = 0
    readToAligns = {}
    for read in withsup:
        if read.is_mapped and not read.is_secondary and read.has_tag('SA'):
            rname = read.query_name
            if rname not in readToAligns:
                readToAligns[rname] = []
            # print(rname, read.get_blocks())
            refchr, refstart, refend, dir = read.reference_name, read.reference_start, read.reference_end, isrevtosign[read.is_reverse]
            genename = getCorrectGene(annot, geneannot, refchr, read.get_blocks())#, isrevtosign[read.is_reverse])  # + '|' + isrevtosign[read.is_reverse]
            qstart, qend = read.query_alignment_start, read.query_alignment_end
            readlen = read.infer_read_length()
            cigar = read.cigartuples
            if cigar[0][0] == 5:  ##just hard clipping
                qstart += cigar[0][1]
                qend += cigar[0][1]
            if dir == '+':
                readToAligns[rname].append([(qstart, refstart), (qend, refend), genename, dir, refchr])
            else:
                readToAligns[rname].append([(readlen - qend, refend), (readlen - qstart, refstart), genename, dir, refchr])

    interestingloci = {}
    for read in readToAligns:
        alignedloci = sorted(readToAligns[read])
        ###ADD: need to handle non-genic alignments
        g1annot = genetoinfo[alignedloci[0][2]] if alignedloci[0][2] in genetoinfo else None
        g2annot = genetoinfo[alignedloci[-1][2]] if alignedloci[-1][2] in genetoinfo else None
        end1strandalignswithtranscript = True if (g1annot and alignedloci[0][0][1] < alignedloci[0][1][1] and g1annot[3] == '+') or (g1annot and alignedloci[0][0][1] > alignedloci[0][1][1] and g1annot[3] == '-') else False
        end2strandalignswithtranscript = True if (g2annot and alignedloci[-1][1][1] < alignedloci[-1][0][1] and g2annot[3] == '+') or (g2annot and alignedloci[-1][1][1] > alignedloci[-1][0][1] and g2annot[3] == '-') else False
        if not end1strandalignswithtranscript and not end2strandalignswithtranscript: continue
        g15end = sorted(g1annot[-1], key=lambda x: abs(x - alignedloci[0][0][1]))[0] if g1annot else None  ###compare to all annot transcript ends
        g25end = sorted(g2annot[-1], key=lambda x: abs(x - alignedloci[-1][1][1]))[0] if g2annot else None
        if (not end1strandalignswithtranscript and end2strandalignswithtranscript) or (
                end2strandalignswithtranscript and end1strandalignswithtranscript and abs(
                alignedloci[-1][1][1] - g25end) < abs(
                alignedloci[0][0][1] - g15end)):  # alignedloci[-1][1][1] < alignedloci[0][0][1]):
            alignedloci = [[x[1], x[0]] + x[2:] for x in alignedloci][
                          ::-1]  ###reverse so smaller transcript pos is start
        readToAligns[read] = alignedloci
        # print(alignedloci)
        ##if end2tpos < end1tpos: alignedloci = [x[::-1] for x in alignedloci][::-1]  ###reverse so smaller transcript pos is start
        readToAligns[read] = alignedloci
        readgenes = [x[2] for x in alignedloci]
        if len(set(readgenes)) > 1:
            # print(readgenes)
            info = tuple(readgenes)
            if info not in interestingloci: interestingloci[info] = []
            interestingloci[info].append(read)

    # fusionsout.write('\t'.join(
    #     ['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']) + '\n')
    fusiontoinfo = {}
    for l in interestingloci:
        if 2 <= len(l) <= maxloci and len(interestingloci[l]) >= minsup:
            numloci = len(l)
            qdist, readsup = [[] for x in range(numloci-1)], 0
            alignblocks = [[[], []] for x in range(numloci)]
            aligngenes = [[] for x in range(numloci)]
            goodreads = []
            for r in interestingloci[l]:
                if len(readToAligns[r]) == numloci:
                    alignedloci = readToAligns[r]
                    goodreads.append(r)
                    readsup += 1
                    for i in range(numloci - 1):
                        qdist[i].append(alignedloci[i][1][0] - alignedloci[i + 1][0][0])

                    for i in range(numloci):
                        alignblocks[i][0].append(alignedloci[i][0][1])
                        alignblocks[i][1].append(alignedloci[i][1][1])
                        aligngenes[i].append((alignedloci[i][2], alignedloci[i][4])) ###gene, chr
            # print(l, readsup, [set(x) for x in aligngenes], alignblocks)
            shortgenes = [list(set([z[0] for z in x])) for x in aligngenes]
            # if aligngenes[0][0][0].split('.')[0] == 'ENSG00000141510': print(shortgenes, 'readsup', readsup >= minsup)
            if readsup >= minsup:
                consistentGenes = True
                # print(aligngenes)
                for i in range(numloci):
                    if len(set(aligngenes[i])) > 1: consistentGenes = False
                # print('consgenes', consistentGenes)
                # if aligngenes[0][0][0].split('.')[0] == 'ENSG00000141510': print(shortgenes, 'consgenes', consistentGenes)

                if consistentGenes:  # check that gene orders for all reads are consistent
                    ###to start, no clustering, take simple min/max
                    for i in range(numloci):
                        aligngenes[i] = list(set(aligngenes[i]))[0]

                        strand = '+' if median(alignblocks[i][1]) > median(alignblocks[i][0]) else '-'
                        for j in range(2):
                            poslist = sorted(alignblocks[i][j])
                            simplemed = median(poslist)
                            groups, g = [], [-500]
                            for p in poslist:
                                if p - g[-1] > 300:
                                    if g[0] != -500: groups.append(g)
                                    g = [p]
                                else:
                                    g.append(p)
                            groups.append(g)
                            goodpos = []
                            for g in groups:
                                if len(g) > 1: goodpos.extend(g)
                            if (strand == '+' and j == 0) or (strand == '-' and j == 1):
                                outpos = int(min(simplemed, min(goodpos))) - 1000 if len(goodpos) > 0 else int(simplemed) - 1000
                            else:
                                outpos = int(max(simplemed, max(goodpos))) + 1000 if len(goodpos) > 0 else int(simplemed) + 1000
                            outpos = max(0, outpos)
                            # print(outpos, simplemed, median(goodpos) if len(goodpos) > 0 else goodpos, goodpos, poslist)
                            alignblocks[i][j] = outpos

                        # if median(alignblocks[i][1]) > median(alignblocks[i][0]):
                        #     alignblocks[i][0] = int(median(alignblocks[i][0])) - 1000
                        #     alignblocks[i][1] = int(median(alignblocks[i][1])) + 1000
                        # else:
                        #     alignblocks[i][0] = int(median(alignblocks[i][0])) + 1000
                        #     alignblocks[i][1] = int(median(alignblocks[i][1])) - 1000

                    ###check that 5' gene is in the forward direction, implies plausible promoter
                    if aligngenes[0][0] in genetoinfo:
                        firstgenedir = genetoinfo[aligngenes[0][0]][3]
                        firstgenetstarts = genetoinfo[aligngenes[0][0]][-1]
                        end5 = alignblocks[0][0] #if firstgenedir == '+' else alignblocks[0][1]
                        if alignblocks[0][0] > alignblocks[0][1]: end5 -= 1000
                        else: end5 += 1000
                        mindisttostart = min([abs(end5-x) for x in firstgenetstarts])
                        # print(mindisttostart, end5, firstgenetstarts)
                        # if aligngenes[0][0].split('.')[0] == 'ENSG00000141510':
                        #     print(firstgenedir)
                        #     print(firstgenetstarts)
                        #     print(alignblocks[0])
                        #     print(end5)
                        #     print([abs(end5-x) for x in firstgenetstarts])
                        #     print(mindisttostart)
                        #     print(qdist)
                        if ((firstgenedir == '+' and alignblocks[0][0] < alignblocks[0][1]) or (
                                firstgenedir == '-' and alignblocks[0][0] > alignblocks[0][1]))\
                                and (reqdisttostart == None or mindisttostart <= reqdisttostart):
                            simscores = []
                            for qdistlist in qdist:
                                simscore = []
                                qdistlist = sorted(qdistlist)
                                for i in range(1, len(qdistlist)):
                                    simscore.append(qdistlist[i] - qdistlist[i-1])
                                simscores.append(median(simscore))
                            # qdist = [median(x) for x in qdist]
                            # if aligngenes[0][0].split('.')[0] == 'ENSG00000141510': print('simscores', simscores, max([abs(min(x)) for x in qdist]), max([abs(min(x)) for x in qdist]) <= 10, max(simscores) <= 3)
                            if max([abs(median(x)) for x in qdist]) <= 10 or (max([abs(min(x)) for x in qdist]) <= 10 and max(simscores) <= 3): ###alignments have to either have few gaps or be very consistent
                                ###['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']
                                # if aligngenes[0][0].split('.')[0] == 'ENSG00000141510': print(shortgenes, 'passes firstgenedir', alignblocks)
                                fname = '__'.join([x[0] for x in aligngenes])
                                fusiontoinfo[fname] = {'reads': set(goodreads), 'disttostart':[mindisttostart], 'qdist':qdist}
                                for i in range(numloci):
                                    # fusiontoinfo[fname][i] = [aligngenes[i][0], aligngenes[i][2], alignblocks[i][0],
                                    #                    alignblocks[i][1]]

                                    fusiontoinfo[fname][aligngenes[i][0]] = [aligngenes[i][1], alignblocks[i][0], alignblocks[i][1]]
                                # for i in range(numloci):
                                #     outline = ['__'.join([x[0] for x in aligngenes]), aligngenes[i][0],
                                #                "gene" + str(i), aligngenes[i][2], alignblocks[i][0], alignblocks[i][1],
                                #                readsup]
                                #     if i == 0: outline.append(','.join(goodreads))
                                #     fusionsout.write('\t'.join([str(x) for x in outline]) + '\n')
    return fusiontoinfo
    # fusionsout.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''for identifying chimeras from genomic alignment''')
    parser.add_argument('-f', '--gtf',
                        help='specify isoforms.gtf')
    parser.add_argument('-b', '--bam',
                        help='filtered bam file from alignment to genome')
    parser.add_argument('-o', '--output',
                        help='output file name')
    args = parser.parse_args()

    annot, geneannot, genetoinfo = gchimparsegtf(args.gtf)







# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='''for identifying chimeras from genomic alignment''')
#     parser.add_argument('-f', '--gtf',
#                         help='specify isoforms.gtf')
#     parser.add_argument('-b', '--bam',
#                         help='filtered bam file from alignment to genome')
#     parser.add_argument('-o', '--output',
#                         help='output file name')
#     args = parser.parse_args()
#
#     annot = {}
#     geneannot = {}
#     genetoinfo = {}
#     for line in open(args.gtf):
#         if line[0] != '#':
#             line = line.split('\t')
#             if line[2] == 'exon':
#                 chr, start, stop, dir = line[0], int(line[3]), int(line[4]), line[6]
#                 genename = line[8].split('gene_id "')[1].split('"')[0]
#                 # genename += '*' + line[8].split('gene_id "')[1].split('"')[0]
#                 if chr not in annot: annot[chr] = defaultdict(def_value)
#                 for i in range(round(start, -1), round(stop, -1), 10):
#                     # if i not in annot[chr]: annot[chr][i] = set()
#                     annot[chr][i].add((genename, dir))
#             elif line[2] == 'gene':
#                 chr, start, stop, dir = line[0], int(line[3]), int(line[4]), line[6]
#                 genename = line[8].split('gene_id "')[1].split('"')[0]
#                 # genename += '*' + line[8].split('gene_id "')[1].split('"')[0]
#                 genetoinfo[genename] = (chr, start, stop, dir)
#                 if chr not in geneannot: geneannot[chr] = defaultdict(def_value)
#                 for i in range(round(start, -2), round(stop, -2), 100):
#                     geneannot[chr][i].add((genename, dir))
#     print('done loading annot')
#
#     isrevtosign = {True: '-', False: '+'}
#     fusionsout = open(args.output, 'w')
#     withsup = pysam.AlignmentFile(args.bam, "rb")
#     c = 0
#     readToAligns = {}
#     for read in withsup:
#         rname = read.query_name
#         if rname not in readToAligns:
#             readToAligns[rname] = []
#         # print(rname, read.get_blocks())
#         genename = getCorrectGene(annot, read.reference_name, read.get_blocks(),
#                                   isrevtosign[read.is_reverse])  # + '|' + isrevtosign[read.is_reverse]
#         refchr, refstart, refend, dir = read.reference_name, read.reference_start, read.reference_end, isrevtosign[
#             read.is_reverse]
#         qstart, qend = read.query_alignment_start, read.query_alignment_end
#         readlen = read.infer_read_length()
#         cigar = read.cigartuples
#         if cigar[0][0] == 5:  ##just hard clipping
#             qstart += cigar[0][1]
#             qend += cigar[0][1]
#         if dir == '+':
#             readToAligns[rname].append(((qstart, refstart), (qend, refend), genename, dir, refchr))
#         else:
#             readToAligns[rname].append(((readlen - qend, refend), (readlen - qstart, refstart), genename, dir, refchr))
#
#     interestingloci = {}
#     for read in readToAligns:
#         readgenes = [x[2] for x in sorted(readToAligns[read])]
#         if len(set(readgenes)) > 1:
#             # print(readgenes)
#             info = tuple(readgenes)
#             if info not in interestingloci: interestingloci[info] = []
#             interestingloci[info].append(read)
#
#     fusionsout.write('\t'.join(
#         ['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']) + '\n')
#
#     for l in interestingloci:
#         if len(l) >= 2 and len(interestingloci[l]) >= 3:
#             qdist, readsup = [], 0
#             numloci = len(l)
#             alignblocks = [[[], []] for x in range(numloci)]
#             aligngenes = [[] for x in range(numloci)]
#             goodreads = []
#             for r in interestingloci[l]:
#                 if len(readToAligns[r]) == numloci:
#                     alignedloci = sorted(readToAligns[r])
#                     goodreads.append(r)
#                     readsup += 1
#                     for i in range(numloci - 1):
#                         qdist.append(alignedloci[i][0][0] - alignedloci[i + 1][1][0])
#
#                     for i in range(numloci):
#                         alignblocks[i][0].append(alignedloci[i][0][1])
#                         alignblocks[i][1].append(alignedloci[i][1][1])
#                         aligngenes[i].append(alignedloci[i][2:])
#             # print(l, readsup, [set(x) for x in aligngenes], alignblocks)
#             if readsup >= 3:
#                 consistentGenes = True
#                 for i in range(numloci):
#                     if len(set(aligngenes[i])) > 1: consistentGenes = False
#                 if consistentGenes:  # check that gene orders for all reads are consistent
#                     ###to start, no clustering, take simple min/max
#                     for i in range(numloci):
#                         aligngenes[i] = list(set(aligngenes[i]))[0]
#
#                     for i in range(numloci):
#                         if aligngenes[i][1] == '+':
#                             alignblocks[i][0] = min(alignblocks[i][0]) - 1000
#                             alignblocks[i][1] = max(alignblocks[i][1]) + 1000
#                         else:
#                             alignblocks[i][0] = max(alignblocks[i][0]) + 1000
#                             alignblocks[i][1] = min(alignblocks[i][1]) - 1000
#
#                     ###check that 5' gene is in the forward direction, implies plausible promoter
#                     if aligngenes[0][0] in genetoinfo:
#                         firstgenedir = genetoinfo[aligngenes[0][0]][-1]
#                         if (firstgenedir == '+' and alignblocks[0][0] < alignblocks[0][1]) or (
#                                 firstgenedir == '-' and alignblocks[0][0] > alignblocks[0][1]):
#                             ###['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']
#                             for i in range(numloci):
#                                 outline = ['__'.join([x[0] for x in aligngenes]), aligngenes[i][0],
#                                            "gene" + str(i), aligngenes[i][2], alignblocks[i][0], alignblocks[i][1],
#                                            readsup]
#                                 if i == 0: outline.append(','.join(goodreads))
#                                 fusionsout.write('\t'.join([str(x) for x in outline]) + '\n')
#
#     fusionsout.close()
