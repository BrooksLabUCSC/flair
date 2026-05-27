import argparse
import pysam
from flair.isoform_data import get_reverse_complement
from flair.gtf_io import gtf_record_parser, gtf_write_row, GtfAttrsSet
from flair.pycbio.hgdata.bed import BedReader

parser = argparse.ArgumentParser(description='make synthetic fusion reference')
parser.add_argument('-c', '--chimbp', action='store', help='bed file of fusion breakpoints')
parser.add_argument('-g', '--genome', action='store', dest='g', help='path to genome')
parser.add_argument('-a', '--anno', action='store', dest='a', help='path to anno.gtf')
parser.add_argument('-o', '--output', action='store', help='output file prefix')
args = parser.parse_args()

prefix = args.output  # '.'.join(args.chimbp.split('.')[:-2])


# DONE Load in transcriptome and genome breakpoints and process them into one list of breakpoint locations
# DONE Go through transcript reference and load in gene start/end locations
# DONE     figure out if any predicted breakpoints are in the same intron node and collapse them
# DONE Cut genes at all predicted breakpoints, label with gene names and cut locations, make synthetic fasta
# DONE make synthetic annotation for these sequences

# check if any reads map to wrong orientation of fusion loci

allBP = {}
fgenes = {}
transcripts = {}
# ['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'breakpointCoord', 'outerEdgeCoord', 'readSupport']
for bed in BedReader(args.chimbp, numStdCols=6, fixScores=True):
    chrom, name, score, strand = bed.chrom, bed.name, bed.score, bed.strand
    start, end = max(0, bed.chromStart), max(0, bed.chromEnd)
    gene, fusion = name.split('__')
    fusion = tuple(fusion.split('--'))
    if fusion not in allBP:
        allBP[fusion] = []
    if strand == '+':
        allBP[fusion].append((gene, chrom, start, end))
    else:
        allBP[fusion].append((gene, chrom, end, start))
    if gene not in transcripts:
        transcripts[gene] = {}
    if gene not in fgenes:
        fgenes[gene] = [chrom, start, end, strand]
    else:
        fgenes[gene][1] = min(start, fgenes[gene][1])
        fgenes[gene][2] = max(end, fgenes[gene][2])

genome = pysam.FastaFile(args.g)
# print('loaded genome')

# To make synthetic transcriptome:
# DONE Get transcript/exon annotation for fusion genes
#     all annotation is recorded in plain left-right direction
# Filter this annotation to 5'/3' ends based on each breakpoint
#     Convert annotation values to be 0-based depending on start of gene (5' end) or breakpoint location (3' end)
#         Make sure to flip - strand values accordingly
# When making synthetic references, simulatneously make gtf annotation file - make sure to convert 3' side values based on

for rec in gtf_record_parser(args.a, include_features={'gene', 'exon', 'start_codon'}, attrs=GtfAttrsSet.ALL):
    genename = rec.gene_id.replace('_', '-').split('.')[0]
    if genename in fgenes:
        if rec.feature == 'gene':
            # learned that can't assume that transcript appears in anno only once - two diff ENSG can have same hugo name
            fgenes[genename] = (rec.chrom, min([rec.start - 500, fgenes[genename][1]]), max([rec.end + 500, fgenes[genename][2]]), rec.strand)
        elif rec.feature in ('exon', 'start_codon'):
            tname = rec.transcript_id
            if tname not in transcripts[genename]:
                transcripts[genename][tname] = []
            if rec.strand == '+':
                transcripts[genename][tname].append((rec.start, rec.end, rec.feature))
            else:
                transcripts[genename][tname].insert(0, (rec.start, rec.end, rec.feature))
# print(transcripts)
# print('loaded transcripts')

out = open(prefix + '-syntheticFusionGenome.fa', 'w')  # 'syntheticFusionGenomeAttempt4.fa', 'w')
annoOut = open(prefix + '-syntheticReferenceAnno.gtf', 'w')  # 'syntheticReferenceAnnoAttempt1.gtf', 'w')
# sjOut = open(prefix + '-syntheticReferenceSJ.bed', 'w')  # 'syntheticReferenceAnnoAttempt1.gtf', 'w')
bpOut = open(prefix + '-syntheticBreakpointLoc.bed', 'w')  # 'syntheticFusionBreakpointLoc.bed', 'w')
c = 0
isocount = 1
for fusion in allBP:  # noqa: C901 - FIXME: reduce complexity
    # print(fusion)
    allisochunks = []
    labels, sequence = [], []
    allstartloc = []
    startLoc = 0
    for order in range(len(fusion)):
        isochunks = {}
        gene, thisChr, leftbound, rightbound = allBP[fusion][order]
        if order == 0:  # 5' gene, correct end to 5' end of gene
            if leftbound < rightbound:
                leftbound = fgenes[gene][1]
            else:
                leftbound = fgenes[gene][2]
        elif order == len(fusion) - 1:
            if leftbound < rightbound:
                rightbound = fgenes[gene][2]
            else:
                rightbound = fgenes[gene][1]

        if leftbound < rightbound:
            sequence.append(genome.fetch(thisChr, leftbound, rightbound))  # genome[thisChr][leftbound:rightbound])
        else:
            sequence.append(get_reverse_complement(genome.fetch(thisChr, rightbound, leftbound)))  # genome[thisChr][rightbound:leftbound]))
        labels.append('..'.join([str(x) for x in [gene, thisChr, leftbound, rightbound]]))

        # Add a justends isoform here
        chunksize = abs(leftbound - rightbound)
        percentchange = int(chunksize * 0.1)
        isochunks['justends-' + str(c)] = [(percentchange + startLoc, (chunksize - percentchange) + startLoc, 'exon'),]
        c += 1

        for tname in transcripts[gene]:
            isochunks[tname] = []
            if leftbound < rightbound:  # positive strand transcript
                for exon in transcripts[gene][tname]:
                    synthexon = None
                    # check if exon is within bounds of fusion
                    if exon[0] >= leftbound and exon[1] <= rightbound and (order == 0 or exon[2] == 'exon'):
                        synthexon = [(exon[0] - leftbound) + startLoc, (exon[1] - leftbound) + startLoc, exon[2]]
                    if synthexon:
                        isochunks[tname].append(tuple(synthexon))
            else:  # negative strand transcript
                for exon in reversed(transcripts[gene][tname]):
                    synthexon = None
                    # check if exon is within bounds of fusion
                    if exon[0] >= rightbound and exon[1] <= leftbound and (order == 0 or exon[2] == 'exon'):
                        synthexon = [(leftbound - exon[1]) + startLoc, (leftbound - exon[0]) + startLoc, exon[2]]
                    if synthexon:
                        isochunks[tname].append(tuple(synthexon))
        allstartloc.append(startLoc)
        startLoc = startLoc + abs(rightbound - leftbound)
        seqlen = len(''.join(sequence))
        allisochunks.append(isochunks)
    # finalisochunks = [[] for x in range(len(fusion))]
    finalisochunks = []
    for order in range(len(fusion)):
        seen = []
        for iso in list(allisochunks[order].keys()):
            if allisochunks[order][iso] != [] and len([x for x in allisochunks[order][iso] if x[2] == 'exon']) > 1 and allisochunks[order][iso] not in seen:
                # finalisochunks[order].append([iso, allisochunks[order][iso]])
                finalisochunks.append([iso, allisochunks[order][iso]])
                seen.append(allisochunks[order][iso])

    fusionname = '--'.join(labels)
    geneid = '--'.join(fusion)
    fusionchrname = '--'.join(['..'.join(x.split('*')) for x in fusion])
    out.write('>' + fusionchrname + '\n')
    out.write(''.join(sequence) + '\n')
    for s in range(1, len(fusion)):
        bpOut.write('\t'.join([fusionchrname, str(allstartloc[s]), str(allstartloc[s]), 'breakpoint-' + str(s) + '--' + fusionname]) + '\n')
    gtf_write_row(annoOut, fusionchrname, 'SYNTHFUSION', 'gene', 0, len(''.join(sequence)), None, '+', None,
                  gene_id=geneid)

    seen = set()
    for isocomb in finalisochunks:
        isoexons = isocomb[1]
        transcriptid = str(isocount) + '-' + isocomb[0]
        isocount += 1
        seen.add(transcriptid)
        gtf_write_row(annoOut, fusionchrname, 'SYNTHFUSION', 'transcript', isoexons[0][0], isoexons[-1][1], None, '+', None,
                      gene_id=geneid, transcript_id=transcriptid)
        for exon in isoexons:
            gtf_write_row(annoOut, fusionchrname, 'SYNTHFUSION', exon[2], exon[0], exon[1], None, '+', None,
                          gene_id=geneid, transcript_id=transcriptid)


out.close()
annoOut.close()
bpOut.close()
