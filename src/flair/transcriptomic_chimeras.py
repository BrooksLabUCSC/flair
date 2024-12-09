#!/usr/bin/env python3

import pysam, sys, argparse
from statistics import median


def binarySearch(arr, t):
    if t <= arr[0]: return arr[0]
    if t >= arr[-1]: return arr[-1]
    i, j, mid = 0, len(arr) - 1, 0
    while i < j:
        mid = int((i + j) / 2)
        if arr[mid] == t:
            return arr[mid]
        elif t < arr[mid]:
            if mid > 0 and t > arr[mid - 1]:
                if abs(arr[mid] - t) < abs(arr[mid - 1] - t):
                    return arr[mid]
                else:
                    return arr[mid - 1]
            j = mid
        else:
            if mid < len(arr) - 1 and t < arr[mid + 1]:
                if abs(arr[mid] - t) < abs(arr[mid + 1] - t):
                    return arr[mid]
                else:
                    return arr[mid + 1]
            i = mid + 1


def getGenomicPreciseLoc(tname, bpCoord, genedir, intronLocs, intronToGenome):
    closestSS = binarySearch(intronLocs[tname], bpCoord)
    bpIntronEnds = intronToGenome[tname][closestSS]
    genomeSS = None
    diffFromSS = bpCoord-closestSS
    if closestSS == 0: ##start of gene
        if genedir == '+': genomeSS = bpIntronEnds[1]
        else: genomeSS = bpIntronEnds[0]
    elif closestSS == max(intronLocs[tname]): ###end of gene
        if genedir == '+': genomeSS = bpIntronEnds[0]
        else: genomeSS = bpIntronEnds[1]
    elif genedir == '+':
        if bpCoord >= closestSS: genomeSS = bpIntronEnds[1]
        else: genomeSS = bpIntronEnds[0]
    else:
        if bpCoord >= closestSS: genomeSS = bpIntronEnds[0]
        else: genomeSS = bpIntronEnds[1]
    genomepos = None
    if genedir == '+': genomepos = genomeSS + diffFromSS
    else: genomepos = genomeSS - diffFromSS
    return genomepos

def parsegtftranscriptomic(gtffilename):
    geneannot = {}
    genetoinfo = {}
    for line in open(gtffilename):
        if line.startswith('#'):
            continue
        line = line.rstrip().split('\t')
        chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
        if ty == 'gene':
            genename = line[8].split('gene_id "')[1].split('"')[0]
            geneannot[genename] = (chrom, start, end, strand)
        if ty != 'exon': continue
        transcript_id = line[8].split('transcript_id "')[1].split('"')[0]
        gene_id = line[8].split('gene_id "')[1].split('"')[0]
        if gene_id not in genetoinfo: genetoinfo[gene_id] = {}
        if transcript_id not in genetoinfo[gene_id]: genetoinfo[gene_id][transcript_id] = []
        genetoinfo[gene_id][transcript_id].append((start, end))

    intronLocs, intronToGenome = {}, {}
    for g in genetoinfo:
        chrom, start, end, strand = geneannot[g]
        for t in genetoinfo[g]:
            myexons = sorted(genetoinfo[g][t])
            first, last = myexons[0], myexons[-1]
            mylocs = [[0, myexons[0][0] - 500, myexons[0][0]]]  ##add start of transcript
            runningtot = 0
            for i in range(len(myexons) - 1):
                runningtot += myexons[i][1] - myexons[i][0]  ##add size of last exon
                mylocs.append([runningtot, myexons[i][1], myexons[i + 1][0]])  ##add intron
            runningtot += myexons[-1][1] - myexons[-1][0]
            mylocs.append([runningtot, myexons[-1][1], myexons[-1][1] + 500])
            if strand == '-':
                mylocs = [[runningtot - mylocs[x][0], mylocs[x][1], mylocs[x][2]] for x in range(len(mylocs))]
            intronLocs[t] = sorted([x[0] for x in mylocs])
            intronToGenome[t] = {x[0]: (x[1], x[2]) for x in mylocs}
    print('loaded annot')
    return geneannot, intronLocs, intronToGenome



def idTranscriptomicChimeras(bam, geneannot, intronLocs, intronToGenome, minsup):
    isrevtosign = {True: '-', False: '+'}
    withsup = pysam.AlignmentFile(bam, "rb")
    c = 0
    readToAligns = {}
    for read in withsup.fetch():
        rname = read.query_name
        if rname not in readToAligns:
            readToAligns[rname] = []
        # genename = read.reference_name.split('|')[1]
        # tname = read.reference_name.split('|')[4]
        ###if aligning to annotated_transcripts.fa ##need to add better flexibility for formatting here
        # tname, genename = read.reference_name.split('_')
        tname = read.reference_name.split('ENSG', 1)[0].rstrip('_')
        genename = 'ENSG' + read.reference_name.split('ENSG', 1)[1]
        genedir = geneannot[genename][3]
        refstart, refend, dir = read.reference_start, read.reference_end, isrevtosign[read.is_reverse]
        refstart = getGenomicPreciseLoc(tname, refstart, genedir, intronLocs, intronToGenome)
        refend = getGenomicPreciseLoc(tname, refend, genedir, intronLocs, intronToGenome)
        qstart, qend = read.query_alignment_start, read.query_alignment_end
        refchr = geneannot[genename][0]
        readlen = read.infer_read_length()
        cigar = read.cigartuples
        if cigar[0][0] == 5:  ##just hard clipping
            qstart += cigar[0][1]
            qend += cigar[0][1]
        if dir == '+':
            readToAligns[rname].append(((qstart, refstart), (qend, refend), genename, genedir, refchr))
        else:
            readToAligns[rname].append(
                ((readlen - qend, refend), (readlen - qstart, refstart), genename, genedir, refchr))
    print('processed bam file')
    interestingloci = {}
    for read in readToAligns:
        readgenes = [x[2] for x in sorted(readToAligns[read])]
        if len(set(readgenes)) > 1:
            info = tuple(readgenes)
            if info not in interestingloci: interestingloci[info] = []
            interestingloci[info].append(read)

    fusiontoinfo = {}
    # fusionsout.write('\t'.join(
    #     ['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']) + '\n')
    for l in interestingloci:
        if len(interestingloci[l]) >= minsup and len(l) >= 2:
            qdist, readsup = [], 0
            numloci = len(l)
            alignblocks = [[[], []] for x in range(numloci)]
            aligngenes = [[] for x in range(numloci)]
            mygeneannot = [geneannot[x] for x in l]
            genesep = True

            for i in range(numloci):
                for j in range(numloci):
                    if i != j:
                        if not (mygeneannot[i][0] != mygeneannot[j][0] or
                                max(0, min(mygeneannot[i][2], mygeneannot[j][2]) - max(mygeneannot[i][1],
                                                                                       mygeneannot[j][1])) == 0):
                            genesep = False
            goodreads = []
            if genesep:
                for r in interestingloci[l]:
                    if len(readToAligns[r]) == numloci:
                        alignedloci = sorted(readToAligns[r])
                        goodreads.append(r)
                        readsup += 1
                        for i in range(numloci - 1):
                            qdist.append(alignedloci[i][0][0] - alignedloci[i + 1][1][0])
                        for i in range(numloci):
                            alignblocks[i][0].append(alignedloci[i][0][1])
                            alignblocks[i][1].append(alignedloci[i][1][1])
                            aligngenes[i].append(alignedloci[i][2:])
                if readsup >= minsup:
                    consistentGenes = True
                    for i in range(numloci):
                        if len(set(aligngenes[i])) > 1: consistentGenes = False
                    if consistentGenes:  # check that the 5' gene and 3' gene are consistent
                        ###to start, no clustering, take simple min/max
                        for i in range(numloci):
                            aligngenes[i] = list(set(aligngenes[i]))[0]

                        for i in range(numloci):
                            if median(alignblocks[i][1]) > median(alignblocks[i][0]):
                                alignblocks[i][0] = min(alignblocks[i][0]) - 1000
                                alignblocks[i][1] = max(alignblocks[i][1]) + 1000
                            else:
                                alignblocks[i][0] = max(alignblocks[i][0]) + 1000
                                alignblocks[i][1] = min(alignblocks[i][1]) - 1000

                        ###check that 5' gene is in the forward direction, implies plausible promoter
                        firstgenedir = geneannot[aligngenes[0][0]][-1]
                        if (firstgenedir == '+' and alignblocks[0][0] < alignblocks[0][1]) or (
                                firstgenedir == '-' and alignblocks[0][0] > alignblocks[0][1]):
                            ###['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']
                            fname = '__'.join([x[0] for x in aligngenes])
                            fusiontoinfo[fname] = {'reads': set(goodreads)}
                            for i in range(numloci):
                                # fusiontoinfo[fname][i] = [aligngenes[i][0], aligngenes[i][2], alignblocks[i][0], alignblocks[i][1]]
                                fusiontoinfo[fname][aligngenes[i][0]] = [aligngenes[i][2], alignblocks[i][0], alignblocks[i][1]]
                                # outline = ['__'.join([x[0] for x in aligngenes]), aligngenes[i][0],
                                #            "gene" + str(i), aligngenes[i][2], alignblocks[i][0], alignblocks[i][1],
                                #            readsup]
                                # if i == 0: outline.append(','.join(goodreads))
                                # fusionsout.write('\t'.join([str(x) for x in outline]) + '\n')
    return fusiontoinfo



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''for identifying chimeras from transcriptomic alignment''')
    parser.add_argument('-f', '--gtf',
                        help='specify isoforms.gtf')
    parser.add_argument('-b', '--bam',
                        help='filtered bam file from alignment to genome')
    parser.add_argument('-o', '--output',
                        help='output file name')
    args = parser.parse_args()

    geneannot, intronLocs, intronToGenome = parsegtftranscriptomic(args.gtf)
    chim = idTranscriptomicChimeras(args.bam, geneannot, intronLocs, intronToGenome)







# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='''for identifying chimeras from transcriptomic alignment''')
#     parser.add_argument('-f', '--gtf',
#                         help='specify isoforms.gtf')
#     parser.add_argument('-b', '--bam',
#                         help='filtered bam file from alignment to genome')
#     parser.add_argument('-o', '--output',
#                         help='output file name')
#     args = parser.parse_args()
#
#     geneannot = {}
#     genetoinfo = {}
#     for line in open(args.gtf):
#         if line.startswith('#'):
#             continue
#         line = line.rstrip().split('\t')
#         chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
#         if ty == 'gene':
#             genename = line[8].split('gene_id "')[1].split('"')[0]
#             geneannot[genename] = (chrom, start, end, strand)
#         if ty != 'exon': continue
#         transcript_id = line[8].split('transcript_id "')[1].split('"')[0]
#         gene_id = line[8].split('gene_id "')[1].split('"')[0]
#         if gene_id not in genetoinfo: genetoinfo[gene_id] = {}
#         if transcript_id not in genetoinfo[gene_id]: genetoinfo[gene_id][transcript_id] = []
#         genetoinfo[gene_id][transcript_id].append((start, end))
#
#     intronLocs, intronToGenome = {}, {}
#     for g in genetoinfo:
#         chrom, start, end, strand = geneannot[g]
#         for t in genetoinfo[g]:
#             myexons = sorted(genetoinfo[g][t])
#             first, last = myexons[0], myexons[-1]
#             mylocs = [[0, myexons[0][0]-500, myexons[0][0]]] ##add start of transcript
#             runningtot = 0
#             for i in range(len(myexons)-1):
#                 runningtot += myexons[i][1] - myexons[i][0] ##add size of last exon
#                 mylocs.append([runningtot, myexons[i][1], myexons[i+1][0]]) ##add intron
#             runningtot += myexons[-1][1] - myexons[-1][0]
#             mylocs.append([runningtot, myexons[-1][1], myexons[-1][1]+500])
#             if strand == '-':
#                 mylocs = [[runningtot - mylocs[x][0], mylocs[x][1], mylocs[x][2]] for x in range(len(mylocs))]
#             intronLocs[t] = sorted([x[0] for x in mylocs])
#             intronToGenome[t] = {x[0]:(x[1], x[2]) for x in mylocs}
#     print('loaded annot')
#
#     isrevtosign = {True: '-', False: '+'}
#     fusionsout = open(args.output, 'w')
#     withsup = pysam.AlignmentFile(args.bam, "rb")
#     c = 0
#     readToAligns = {}
#     for read in withsup.fetch():
#         rname = read.query_name
#         if rname not in readToAligns:
#             readToAligns[rname] = []
#         # genename = read.reference_name.split('|')[1]
#         # tname = read.reference_name.split('|')[4]
#         ###if aligning to annotated_transcripts.fa ##need to add better flexibility for formatting here
#         tname, genename = read.reference_name.split('_')
#         genedir = geneannot[genename][3]
#         refstart, refend, dir = read.reference_start, read.reference_end, isrevtosign[read.is_reverse]
#         refstart = getGenomicPreciseLoc(tname, refstart, genedir)
#         refend = getGenomicPreciseLoc(tname, refend, genedir)
#         qstart, qend = read.query_alignment_start, read.query_alignment_end
#         refchr = geneannot[genename][0]
#         readlen = read.infer_read_length()
#         cigar = read.cigartuples
#         if cigar[0][0] == 5:  ##just hard clipping
#             qstart += cigar[0][1]
#             qend += cigar[0][1]
#         if dir == '+':
#             readToAligns[rname].append(((qstart, refstart), (qend, refend), genename, genedir, refchr))
#         else:
#             readToAligns[rname].append(
#                 ((readlen - qend, refend), (readlen - qstart, refstart), genename, genedir, refchr))
#     print('processed bam file')
#     interestingloci = {}
#     for read in readToAligns:
#         readgenes = [x[2] for x in sorted(readToAligns[read])]
#         if len(set(readgenes)) > 1:
#             info = tuple(readgenes)
#             if info not in interestingloci: interestingloci[info] = []
#             interestingloci[info].append(read)
#
#     fusionsout.write('\t'.join(
#         ['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']) + '\n')
#     for l in interestingloci:
#         if len(interestingloci[l]) >= 3 and len(l) >= 2:
#             qdist, readsup = [], 0
#             numloci = len(l)
#             alignblocks = [[[], []] for x in range(numloci)]
#             aligngenes = [[] for x in range(numloci)]
#             mygeneannot = [geneannot[x] for x in l]
#             genesep = True
#
#             for i in range(numloci):
#                 for j in range(numloci):
#                     if i != j:
#                         if not (mygeneannot[i][0] != mygeneannot[j][0] or
#                                 max(0, min(mygeneannot[i][2],mygeneannot[j][2]) - max(mygeneannot[i][1], mygeneannot[j][1])) == 0):
#                             genesep = False
#             goodreads = []
#             if genesep:
#                 for r in interestingloci[l]:
#                     if len(readToAligns[r]) == numloci:
#                         alignedloci = sorted(readToAligns[r])
#                         goodreads.append(r)
#                         readsup += 1
#                         for i in range(numloci - 1):
#                             qdist.append(alignedloci[i][0][0] - alignedloci[i + 1][1][0])
#                         for i in range(numloci):
#                             alignblocks[i][0].append(alignedloci[i][0][1])
#                             alignblocks[i][1].append(alignedloci[i][1][1])
#                             aligngenes[i].append(alignedloci[i][2:])
#                 if readsup >= 3:
#                     consistentGenes = True
#                     for i in range(numloci):
#                         if len(set(aligngenes[i])) > 1: consistentGenes = False
#                     if consistentGenes:  # check that the 5' gene and 3' gene are consistent
#                         ###to start, no clustering, take simple min/max
#                         for i in range(numloci):
#                             aligngenes[i] = list(set(aligngenes[i]))[0]
#
#                         for i in range(numloci):
#                             if median(alignblocks[i][1]) > median(alignblocks[i][0]):
#                                 alignblocks[i][0] = min(alignblocks[i][0]) - 1000
#                                 alignblocks[i][1] = max(alignblocks[i][1]) + 1000
#                             else:
#                                 alignblocks[i][0] = max(alignblocks[i][0]) + 1000
#                                 alignblocks[i][1] = min(alignblocks[i][1]) - 1000
#
#                         ###check that 5' gene is in the forward direction, implies plausible promoter
#                         firstgenedir = geneannot[aligngenes[0][0]][-1]
#                         if (firstgenedir == '+' and alignblocks[0][0] < alignblocks[0][1]) or (
#                                 firstgenedir == '-' and alignblocks[0][0] > alignblocks[0][1]):
#                             ###['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'leftCoord', 'rightCoord', 'readSupport']
#                             for i in range(numloci):
#                                 outline = ['__'.join([x[0] for x in aligngenes]), aligngenes[i][0],
#                                            "gene" + str(i), aligngenes[i][2], alignblocks[i][0], alignblocks[i][1],
#                                            readsup]
#                                 if i == 0: outline.append(','.join(goodreads))
#                                 fusionsout.write('\t'.join([str(x) for x in outline]) + '\n')
#     fusionsout.close()
