#!/usr/bin/env python3

import pysam
import argparse
from collections import Counter
from statistics import median
from flair.convert_synthetic_to_genome_bed import identify_fusion_problems


def binarySearch(arr, t):
    if t <= arr[0]:
        return arr[0]
    if t >= arr[-1]:
        return arr[-1]
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
    diffFromSS = bpCoord - closestSS
    if closestSS == 0:  # start of gene
        if genedir == '+':
            genomeSS = bpIntronEnds[1]
        else:
            genomeSS = bpIntronEnds[0]
    elif closestSS == max(intronLocs[tname]):  # end of gene
        if genedir == '+':
            genomeSS = bpIntronEnds[0]
        else:
            genomeSS = bpIntronEnds[1]
    elif genedir == '+':
        if bpCoord >= closestSS:
            genomeSS = bpIntronEnds[1]
        else:
            genomeSS = bpIntronEnds[0]
    else:
        if bpCoord >= closestSS:
            genomeSS = bpIntronEnds[0]
        else:
            genomeSS = bpIntronEnds[1]
    genomepos = None
    if genedir == '+':
        genomepos = genomeSS + diffFromSS
    else:
        genomepos = genomeSS - diffFromSS
    return genomepos


def def_value():
    return set()

def get_exon_intron_blocks(read):
    align_start = read.reference_start
    align_end = read.reference_end
    ref_pos = align_start
    intron_blocks = []
    has_match = False
    for block in read.cigartuples:
        if block[0] == 3:  # intron
            if has_match:
                intron_blocks.append([ref_pos, ref_pos + block[1]])
            # this fixes weird bug if there's an intron, then an insertion, then another intron???
            elif len(intron_blocks) > 0:
                intron_blocks[-1][1] += block[1]
            has_match = False
            ref_pos += block[1]
        elif block[0] in {0, 7, 8, 2}:  # consumes reference
            ref_pos += block[1]
            if block[0] in {0, 7, 8}:
                has_match = True
    intron_blocks = [tuple(x) for x in intron_blocks]
    if len(intron_blocks) == 0:
        exon_blocks = [(align_start, align_end),]
    else:
        exon_blocks = [(align_start, intron_blocks[0][0])] + [(intron_blocks[x][1], intron_blocks[x + 1][0]) for x in range(len(intron_blocks) - 1)] + [(intron_blocks[-1][1], align_end)]
    return intron_blocks, exon_blocks


def getCorrectGene(chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, chrom, readblocks, thisdir):  # noqa: C901 - FIXME: reduce complexity
    intron_blocks, exon_blocks = readblocks
    if chrom not in chrom_to_gene_pos:
        my_gene = chrom + ':' + str(round(exon_blocks[0][0], -4))
    else:
        found_genes = []
        for j in intron_blocks:
            if j in juncs_to_gene[chrom]:
                found_genes.extend(list(juncs_to_gene[chrom][j]))
        if len(found_genes) > 0:
            my_gene = Counter(found_genes).most_common()[0][0]
        else:
            s, e = exon_blocks[0][0], exon_blocks[-1][1]
            gene_overlaps = []
            for start, end, strand, gene in chrom_to_gene_pos[chrom]:
                if start > e:  # assumes sorted
                    break
                if strand == thisdir and min(e, end) > max(start, s):  # require strand match
                    totoverlap = 0
                    for es, ee in gene_to_all_exons[gene]:
                        overlap = min(e, ee) - max(es, s)
                        if overlap > 0:
                            totoverlap += overlap
                    gene_overlaps.append((totoverlap, gene))
            if len(gene_overlaps) > 0:
                gene_overlaps.sort(reverse=True)
                my_gene = gene_overlaps[0][1]
            else:
                my_gene = chrom + ':' + str(round(exon_blocks[0][0], -4))
    return my_gene


def id_chimeras(mode, bam, genetoinfo, chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, gene_to_paralogs,  # noqa: C901 - FIXME: reduce complexity
                genetoname, minsup, maxloci=10, reqdisttostart=None, maxpromiscuity=4, intronLocs=None, intronToGenome=None):
    isrevtosign = {True: '-', False: '+'}
    withsup = pysam.AlignmentFile(bam, "rb")
    readToAligns = {}
    for read in withsup:
        if read.is_mapped and not read.is_secondary and read.has_tag('SA') and (mode == 'genomic' or not read.is_reverse):  # if aligned to transcriptome, must match strand of transcript
            rname = read.query_name
            if rname not in readToAligns:
                readToAligns[rname] = []

            refstart, refend, readdir = read.reference_start, read.reference_end, isrevtosign[read.is_reverse]
            qstart, qend = read.query_alignment_start, read.query_alignment_end
            readlen = read.infer_read_length()
            cigar = read.cigartuples
            if cigar[0][0] == 5:  # just hard clipping
                qstart += cigar[0][1]
                qend += cigar[0][1]

            if mode == 'genomic':
                refchr = read.reference_name
                genename = getCorrectGene(chrom_to_gene_pos, gene_to_all_exons, juncs_to_gene, refchr, get_exon_intron_blocks(read), readdir)
                if readdir == '+':
                    readToAligns[rname].append([(qstart, refstart), (qend, refend), genename, readdir, refchr])
                else:
                    readToAligns[rname].append([(readlen - qend, refend), (readlen - qstart, refstart), genename, readdir, refchr])
            else:
                genename = read.reference_name.split('_')[-1].split('.')[0]
                tname = '_'.join(read.reference_name.split('_')[:-1])
                refchr, genedir = genetoinfo[genename][0], genetoinfo[genename][3]  # can do this because already required the read to be forward strand
                refstart = getGenomicPreciseLoc(tname, refstart, genedir, intronLocs, intronToGenome)
                refend = getGenomicPreciseLoc(tname, refend, genedir, intronLocs, intronToGenome)
                refstart, refend = min(refstart, refend), max(refstart, refend)
                if genedir == '+':
                    readToAligns[rname].append([(qstart, refstart), (qend, refend), genename, readdir, refchr])
                else:
                    readToAligns[rname].append([(qstart, refend), (qend, refstart), genename, genedir, refchr])
    withsup.close()

    interestingloci = {}
    for read in readToAligns:
        alignedloci = sorted(readToAligns[read])

        readgenes = [x[2] for x in alignedloci]
        if 2 <= len(set(readgenes)) <= maxloci:
            if readgenes[0] in gene_to_all_exons:  # 5' gene is annotated
                info = tuple(readgenes)
                if info not in interestingloci:
                    interestingloci[info] = []
                interestingloci[info].append(read)

    goodcov = []
    for l in interestingloci:
        if len(interestingloci[l]) >= minsup:
            goodcov.append(l)

    locustopartners = {}
    for fgenes in goodcov:
        for i in fgenes:
            other = set(fgenes) - {i, }
            newother = frozenset([gene_to_paralogs[g] if g in gene_to_paralogs else g for g in other])
            if i not in locustopartners:
                locustopartners[i] = set()
            locustopartners[i].add(newother)

    fusiontoinfo = {}
    for fgenes in goodcov:
        genomic_chroms = [genetoinfo[g][0] if g in genetoinfo else g.split(':')[0] for g in fgenes]
        is_good_fusion = identify_fusion_problems(fgenes, locustopartners, maxpromiscuity, genetoname, gene_to_paralogs, genomic_chroms, len(interestingloci[fgenes]))
        if is_good_fusion:
            numloci = len(fgenes)
            qdist = [[] for x in range(numloci - 1)]
            alignblocks = [[[], []] for x in range(numloci)]
            goodreads = []
            for r in interestingloci[fgenes]:
                if len(readToAligns[r]) == numloci:
                    alignedloci = sorted(readToAligns[r])
                    goodreads.append(r)
                    for i in range(numloci - 1):
                        qdist[i].append(alignedloci[i][1][0] - alignedloci[i + 1][0][0])
                    for i in range(numloci):
                        alignblocks[i][0].append(alignedloci[i][0][1])
                        alignblocks[i][1].append(alignedloci[i][1][1])

            # check that gene orders for all reads are consistent
            # to start, no clustering, take simple min/max
            for i in range(numloci):
                strand = '+' if median(alignblocks[i][1]) > median(alignblocks[i][0]) else '-'
                for j in range(2):
                    poslist = sorted(alignblocks[i][j])
                    simplemed = median(poslist)
                    groups, g = [], [-500]
                    for p in poslist:
                        if p - g[-1] > 300:
                            if g[0] != -500:
                                groups.append(g)
                            g = [p]
                        else:
                            g.append(p)
                    groups.append(g)
                    goodpos = []
                    for g in groups:
                        if len(g) > 1:
                            goodpos.extend(g)
                    if (strand == '+' and j == 0) or (strand == '-' and j == 1):
                        outpos = int(min(simplemed, min(goodpos))) - 1000 if len(goodpos) > 0 else int(simplemed) - 1000
                    else:
                        outpos = int(max(simplemed, max(goodpos))) + 1000 if len(goodpos) > 0 else int(simplemed) + 1000
                    outpos = max(0, outpos)
                    alignblocks[i][j] = outpos

            # check that 5' gene is in the forward direction, implies plausible promoter
            if fgenes[0] in genetoinfo:
                firstgenedir = genetoinfo[fgenes[0]][3]
                firstgenetstarts = genetoinfo[fgenes[0]][-1]
                end5 = alignblocks[0][0]  # if firstgenedir == '+' else alignblocks[0][1]
                if alignblocks[0][0] > alignblocks[0][1]:
                    end5 -= 1000
                else:
                    end5 += 1000
                mindisttostart = min([abs(end5 - x) for x in firstgenetstarts])

                if ((firstgenedir == '+' and alignblocks[0][0] < alignblocks[0][1])
                        or (firstgenedir == '-' and alignblocks[0][0] > alignblocks[0][1]))\
                        and (reqdisttostart is None or mindisttostart <= reqdisttostart):
                    simscores = []
                    for qdistlist in qdist:
                        simscore = []
                        qdistlist = sorted(qdistlist)
                        for i in range(1, len(qdistlist)):
                            simscore.append(qdistlist[i] - qdistlist[i - 1])
                        simscores.append(median(simscore))
                    if max([abs(median(x)) for x in qdist]) <= 10 \
                            or (max([abs(min(x)) for x in qdist]) <= 10 and max(simscores) <= 3):  # alignments have to either have few gaps or be very consistent
                        fname = '__'.join(fgenes)
                        fusiontoinfo[fname] = {'reads': set(goodreads), 'disttostart': [mindisttostart], 'qdist': qdist}
                        for i in range(numloci):
                            fusiontoinfo[fname][fgenes[i]] = [genomic_chroms[i], alignblocks[i][0], alignblocks[i][1]]

    return fusiontoinfo


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''for identifying chimeras from genomic alignment''')
    parser.add_argument('-f', '--gtf',
                        help='specify isoforms.gtf')
    parser.add_argument('-b', '--bam',
                        help='filtered bam file from alignment to genome')
    parser.add_argument('-o', '--output',
                        help='output file name')
    args = parser.parse_args()

    # annot, geneannot, genetoinfo = gchimparsegtf(args.gtf)
