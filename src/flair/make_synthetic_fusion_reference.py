import sys, os, argparse
from statistics import median,stdev
from datetime import date
import itertools
import pysam

parser = argparse.ArgumentParser(description='make synthetic fusion reference')
parser.add_argument('-c', '--chimbp', action='store', help='bed file of fusion breakpoints')
parser.add_argument('-g', '--genome', action='store', dest='g', help='path to genome')
parser.add_argument('-a', '--anno', action='store', dest='a', help='path to anno.gtf')
parser.add_argument('-o', '--output', action='store', help='output file prefix')
args = parser.parse_args()

prefix = args.output#'.'.join(args.chimbp.split('.')[:-2])

def revComp(seq):
    newseq = ''
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'}
    for char in seq[::-1].upper():
        newseq += comp[char]
    return newseq

#####DONELoad in transcriptome and genome breakpoints and process them into one list of breakpoint locations
####DONEGo through transcript reference and load in gene start/end locations
####DONE    figure out if any predicted breakpoints are in the same intron node and collapse them
####Cut genes at all predicted breakpoints, label with gene names and cut locations, make synthetic fasta
####make synthetic annotation for these sequences

###check if any reads map to wrong orientation of fusion loci

allBP = {}
fgenes = {}
transcripts = {}
###['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'breakpointCoord', 'outerEdgeCoord', 'readSupport']
for line in open(args.chimbp):#'31-01-2023DRR059313-transcriptomeChimericBreakpoints-correctDir.tsv'):
    chrom, start, end, name, score, strand = line.rstrip().split('\t')[:6]
    start, end = max(0, int(start)), max(0, int(end))
    gene, fusion = name.split('__')
    fusion = tuple(fusion.split('--'))
    if fusion not in allBP: allBP[fusion] = []
    if strand == '+': allBP[fusion].append((gene, chrom, start, end))
    else: allBP[fusion].append((gene, chrom, end, start))
    if gene not in transcripts: transcripts[gene] = {}
    if gene not in fgenes: fgenes[gene] = [chrom, start, end, strand]
    else:
        fgenes[gene][1] = min(start, fgenes[gene][1])
        fgenes[gene][2] = max(end, fgenes[gene][2])

genome=pysam.FastaFile(args.g)
print('loaded genome')

####To make synthetic transcriptome:
####DONE Get transcript/exon annotation for fusion genes
####    all annotation is recorded in plain left-right direction
####Filter this annotation to 5'/3' ends based on each breakpoint
####    Convert annotation values to be 0-based depending on start of gene (5' end) or breakpoint location (3' end)
####        Make sure to flip - strand values accordingly
####When making synthetic references, simulatneously make gtf annotation file - make sure to convert 3' side values based on

for line in open(args.a):#'/private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf'):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'gene' or line[2] == 'exon' or line[2] == 'start_codon':
            genename = line[8].split('gene_id "')[1].split('"')[0]
            genename = genename.replace('_','-')
            # genename += '*' + line[8].split('gene_id "')[1].split('"')[0]
            if genename in fgenes:
                if line[2] == 'gene':
                    ###learned that can't assume that transcript appears in anno only once - two diff ENSG can have same hugo name
                    fgenes[genename] = (line[0], min([int(line[3]) - 501, fgenes[genename][1]]), max([int(line[4])+500, fgenes[genename][2]]), line[6])
                elif line[2] == 'exon' or line[2] == 'start_codon':
                    tname = line[8].split('transcript_id "')[1].split('"')[0]
                    if tname not in transcripts[genename]: transcripts[genename][tname] = []
                    if line[6] == '+': transcripts[genename][tname].append((int(line[3])-1, int(line[4]), line[2]))
                    else: transcripts[genename][tname].insert(0,(int(line[3])-1, int(line[4]), line[2]))
# print(transcripts)
print('loaded transcripts')

out = open(prefix + '-syntheticFusionGenome.fa', 'w')#'syntheticFusionGenomeAttempt4.fa', 'w')
annoOut = open(prefix + '-syntheticReferenceAnno.gtf', 'w')#'syntheticReferenceAnnoAttempt1.gtf', 'w')
# sjOut = open(prefix + '-syntheticReferenceSJ.bed', 'w')#'syntheticReferenceAnnoAttempt1.gtf', 'w')
bpOut = open(prefix + '-syntheticBreakpointLoc.bed', 'w')#'syntheticFusionBreakpointLoc.bed', 'w')
c = 0
for fusion in allBP:
    # print(fusion)
    allisochunks = []
    labels, sequence = [], []
    allstartloc = []
    startLoc = 0
    for order in range(len(fusion)):
        isochunks = {}
        gene, thisChr, leftbound, rightbound = allBP[fusion][order]
        if order == 0: #5' gene, correct end to 5' end of gene
            if leftbound < rightbound: leftbound = fgenes[gene][1]
            else: leftbound = fgenes[gene][2]
        elif order == len(fusion) - 1:
            if leftbound < rightbound: rightbound = fgenes[gene][2]
            else: rightbound = fgenes[gene][1]

        if leftbound < rightbound: sequence.append(genome.fetch(thisChr, leftbound, rightbound))#genome[thisChr][leftbound:rightbound])
        else: sequence.append(revComp(genome.fetch(thisChr, rightbound, leftbound)))#genome[thisChr][rightbound:leftbound]))
        labels.append('..'.join([str(x) for x in [gene, thisChr, leftbound, rightbound]]))

        ###Add a justends isoform here
        chunksize = abs(leftbound-rightbound)
        percentchange = int(chunksize * 0.1)
        isochunks['justends-' + str(c)] = [(percentchange+startLoc, (chunksize-percentchange)+startLoc, 'exon'),]
        c += 1

        for tname in transcripts[gene]:
            isochunks[tname] = []
            if leftbound < rightbound: ###positive strand transcript
                for exon in transcripts[gene][tname]:
                    synthexon = None
                    ###check if exon is within bounds of fusion
                    if exon[0] >= leftbound and exon[1] <= rightbound and (order == 0 or exon[2] == 'exon'):
                        synthexon = [(exon[0]-leftbound)+startLoc, (exon[1]-leftbound) + startLoc, exon[2]]
                    if synthexon: isochunks[tname].append(tuple(synthexon))
            else: #negative strand transcript
                for exon in reversed(transcripts[gene][tname]):
                    synthexon = None
                    ###check if exon is within bounds of fusion
                    if exon[0] >= rightbound and exon[1] <= leftbound and (order == 0 or exon[2] == 'exon'):
                        synthexon = [(leftbound-exon[1]) + startLoc, (leftbound-exon[0]) + startLoc, exon[2]]
                    if synthexon: isochunks[tname].append(tuple(synthexon))
        allstartloc.append(startLoc)
        startLoc = startLoc + abs(rightbound-leftbound)
        seqlen = len(''.join(sequence))
        # print(order, isochunks)
        allisochunks.append(isochunks)
    finalisochunks = [[] for x in range(len(fusion))]
    for order in range(len(fusion)):
        seen = []
        for iso in list(allisochunks[order].keys()):
            if not(allisochunks[order][iso] == [] or allisochunks[order][iso] in seen):
                # allisochunks[order].pop(iso)
                finalisochunks[order].append([iso, allisochunks[order][iso]])
                seen.append(allisochunks[order][iso])
        # print(order, len(finalisochunks[order]), finalisochunks[order])

    fusionname = '--'.join(labels)
    geneid = '--'.join(fusion)
    fusionchrname = '--'.join(['..'.join(x.split('*')) for x in fusion])
    # print(fusion)
    # print(fusionchrname)
    out.write('>' + fusionchrname + '\n')
    out.write(''.join(sequence) + '\n')
    for s in range(1, len(fusion)):
        bpOut.write('\t'.join([fusionchrname, str(allstartloc[s]), str(allstartloc[s]), 'breakpoint-' + str(s) + '--' + fusionname]) + '\n')
    annoOut.write('\t'.join([fusionchrname, 'SYNTHFUSION', 'gene', '1', str(len(''.join(sequence))), '.', '+', '.','gene_id "' + geneid + '"']) + '\n')

    finalisocomb = list(itertools.product(*finalisochunks))
    # print(len(finalisocomb))
    for isocomb in finalisocomb:
        isonames = [x[0] for x in isocomb]
        isoexons = [x[1] for x in isocomb]
        transcriptid = '_'.join(isonames)
        annoOut.write('\t'.join(
            [fusionchrname, 'SYNTHFUSION', 'transcript', str(isoexons[0][0][0] + 1), str(isoexons[-1][-1][1]), '.', '+', '.', '; '.join(['gene_id "' + geneid + '"','transcript_id "' + transcriptid + '"'])]) + '\n')
        for exonset in isoexons:
            for exon in exonset:
                annoOut.write('\t'.join([fusionchrname, 'SYNTHFUSION', exon[2], str(exon[0] + 1), str(exon[1]), '.', '+', '.', '; '.join(['gene_id "' + geneid + '"','transcript_id "' + transcriptid + '"'])]) + '\n')


out.close()
annoOut.close()
bpOut.close()
