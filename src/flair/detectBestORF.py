

import sys
import argparse
import os, glob, math
import pipettor
import pysam
import pybedtools
import numpy as np


class Isoform(object):
    '''
    Object to handle isoform related data.

    attributes:

    methods:

    '''

    def __init__(self, name=None):
        self.name = name
        # self.pro = "UNK"
        self.chrom = ""
        self.strand = None
        self.seqvariants = {} ##dictionary of sequence variants (variant aware isoforms) with this isoform bed structure
        self.exons     = set()
        self.starts    = set() ###annotated start coords on transcript
        # self.orfs      = list()
        # self.exonSizes = list()
        self.ptcpoint = 0 ##location in transcript where if stop codon is before this point, transcript has a premature termination codon
        self.startcodons = [] ###position of all ATG codons on transcript

class SeqVar(object):
    def __init__(self, name=None, seq=None):
        self.name = name
        self.pro = "UNK" ##productivity prediction
        self.sequence = seq
        self.orfs = list()
        self.bestorf = None ##only set if pro='PRO'
        self.aaseq = None ###amino acid sequence from best orf - only set if pro='PRO'

def parseGTFline(l):
    cols = l.rstrip().split("\t")
    chrom, c1, c2, strand = cols[0], int(cols[3]) - 1, int(cols[4]), cols[6]
    if cols[2] == "start_codon":
        gene = cols[8][cols[8].find('gene_id') + len('gene_id') + 2:]
        gene = gene[:gene.find('"')]
        return chrom,c1,c2,gene,".",strand
    else: return None

def getStarts(gtf):
    starts = []
    for l in open(gtf):
        if l[0] != "#":
            startinfo = parseGTFline(l)
            if startinfo: starts.append(startinfo)
    if (len(starts)) == 0:
        sys.stderr.write('ERROR, no start codons were found in', gtf)
        sys.exit(1)
    return starts


def getStartRelPos(genomicStartPos,exon, exons, isoObj):
    '''
    is handed a genomic position, the exon it occurs in, all exons,
    and returns the position relative to all exons
    '''
    exonNum = exons.index(exon)
    isoObj.exonSizes = [x[1]-x[0] for x in exons]
    relativeStart = None
    # First get start position relative to transcript sequence.
    if isoObj.strand == "+":
        relativeStart = genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])
    elif isoObj.strand == "-":
        relativeStart = len(isoObj.sequence) - (genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])) - 3

    return relativeStart

def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_prefix', default='flair',
                        help="output prefix. default: 'flair'")
    parser.add_argument('-i', '--isoforms',
                        help="path to transcriptome fasta file", required=True)
    parser.add_argument('-b', '--bedisoforms',
                        help="path to transcriptome bed file", required=True)
    parser.add_argument('-f', '--gtf',
                        type=str, required=True,
                        help='GTF annotation file')
    parser.add_argument('-l', '--logitmodel',
                        type=str, required=True, help='logimodel.RData file for cpat')
    parser.add_argument('-x', '--hexamer',
                        type=str, required=True,
                        help='hexamer.tsv file for cpat')
    args = parser.parse_args()
    return args



def addATGpos(isoDict, fastafile):
    ###get positions of ATGs on transcript
    isoname = None
    for line in open(fastafile):
        if line[0] == '>':
            isoname = line[1:].rstrip().split(' ')[0].upper()
        else:
            seq = line.rstrip()
            bedisoname = isoname  # '-'.join(isoname.split('-')[1:])
            if bedisoname not in isoDict:
                isoDict[bedisoname] = Isoform(bedisoname)
            isoDict[bedisoname].seqvariants[isoname] = SeqVar(isoname, seq)
            for i in range(len(seq) - 3):
                codon = seq[i:i + 3]
                if codon == 'ATG': isoDict[bedisoname].startcodons.append(i)
    return isoDict

def getBedInfo(line):
    esizes = [int(x) for x in line[10].rstrip(',').split(',')]
    isostart = int(line[1])
    estarts = [int(x) for x in line[11].rstrip(',').split(',')]
    if line[9] == '1': ptcpointont = 0
    else:
        if esizes[-2] < 55: ptcpointont = sum(esizes[:-2])
        else: ptcpointont = sum(esizes[:-1]) - 55
    exons = [(isostart + estarts[x], isostart + estarts[x] + esizes[x]) for x in range(len(estarts))]
    return exons, ptcpointont

def processBedFile(isoDict, bedobj):
    for line in bedobj:
        isoname = line[3].upper()
        if isoname in isoDict:
            exons, ptcpointont = getBedInfo(line.fields)
            isoDict[isoname].chrom = line.chrom
            isoDict[isoname].strand = line.strand
            isoDict[isoname].ptcpoint = ptcpointont
            isoDict[isoname].exons = exons
    return isoDict

def getRelativePos(exons, thisexon, genomestartpos, strand):
    exonNum = exons.index(thisexon)
    exonSizes = [x[1] - x[0] for x in exons]
    # First get start position relative to transcript sequence. ##NOTE This might be incorrect for transcripts that have insertions or deletions, which could impact sequence length beyond what is represented in the bed file
    relativeStart = None
    if strand == "+":
        relativeStart = (genomestartpos - exons[exonNum][0]) + sum([x for x in exonSizes[:exonNum]])
    elif strand == "-":
        relativeStart = (sum(exonSizes) - (genomestartpos - exons[exonNum][0]) + sum(
            [x for x in exonSizes[:exonNum]])) - 3
    return relativeStart

def addAnnotStarts(isoDict, gtffile, isobed):
    annotstarts = getStarts(gtffile)
    isoexons = isobed.bed6()
    startsbed = pybedtools.BedTool(annotstarts)
    bt_st = isoexons.intersect(startsbed, s=True, split=True, wao=True)
    for intersection in bt_st:
        read = intersection[3].upper()
        if read in isoDict:
            # iso,gene = read.split("_")
            overlap = intersection[-1]
            if overlap == "3": # check for full overlap of start codon with exon
                genomestartpos = int(intersection[-6])
                thisexon = (int(intersection[1]), int(intersection[2]))
                relativeStart = getRelativePos(isoDict[read].exons, thisexon, genomestartpos, isoDict[read].strand)
                isoDict[read].starts.add(relativeStart)
    return isoDict

def addNoOrf(isoDict, noorffile):
    for line in open(noorffile):
        isoname = line.rstrip().upper()
        bedisoname = isoname  # '-'.join(isoname.split('-')[1:])
        isoDict[bedisoname].seqvariants[isoname].pro = 'NGO'

def modifyOrfScore(codingprob, numpriorstarts, atg_growth, atg_shift):
    atgscore = 1 - 1 / (1 + np.exp(-atg_growth * (numpriorstarts - atg_shift)))
    return codingprob * atgscore


def calcOrfProb(line, atg_growth, atg_shift, isoobj, seqlen):
    isoid, isolen, orfstart, orfend, orflen, codingprob = line[0], int(line[1]), int(line[4]), int(line[5]), int \
        (line[6]), float(line[9])
    orfstart -= 1
    if codingprob > 0.35:
        hasstopcodon = orfend < seqlen - 1
        isnotPTC = not (orfend < isoobj.ptcpoint)
        hasannotstart = orfstart in isoobj.starts
        numpriorstarts = sum([1 if x < orfstart else 0 for x in isoobj.startcodons])
        orfscore = modifyOrfScore(codingprob, numpriorstarts, atg_growth, atg_shift)
        ###FIXME I think I want to modify this so having an annotated start bumps up your score instead of overwriting it
        return hasstopcodon, isnotPTC, hasannotstart, orfscore, orfstart, orfend, orflen
    else: return None


def addOrfProb(isoDict, orfprobfile, atg_growth, atg_shift):
    for line in open(orfprobfile):
        line = line.rstrip().split('\t')
        if line[0] != 'ID':
            isoid = '_'.join(line[0].split('_')[:-2])
            bedisoname = isoid  # '-'.join(isoid.split('-')[1:])
            orfinfo = calcOrfProb(line, atg_growth, atg_shift, isoDict[bedisoname],
                                  len(isoDict[bedisoname].seqvariants[isoid].sequence))
            if orfinfo: isoDict[bedisoname].seqvariants[isoid].orfs.append(orfinfo)
    return isoDict

def calcBestOrf(isoDict):
    for bedisoname in isoDict:
        for isoid in isoDict[bedisoname].seqvariants:
            orfs = sorted(isoDict[bedisoname].seqvariants[isoid].orfs, reverse=True)
            if len(orfs) == 0: isoDict[bedisoname].seqvariants[isoid].pro = 'NGO'
            else:
                bestorf = orfs[0]
                if not bestorf[0]: isoDict[bedisoname].seqvariants[isoid].pro = 'NST'
                elif not bestorf[1]: isoDict[bedisoname].seqvariants[isoid].pro = 'PTC'
                else:
                    isoDict[bedisoname].seqvariants[isoid].pro = 'PRO'
                    isoDict[bedisoname].seqvariants[isoid].bestorf = bestorf
                    orfstart, orfend = bestorf[4], bestorf[5]
                    orfseq = isoDict[bedisoname].seqvariants[isoid].sequence[orfstart:orfend]
                    isoDict[bedisoname].seqvariants[isoid].aaseq = translate(orfseq)
    return isoDict

def printoutfile(isoDict, outfilename):
    out = open(outfilename, 'w')
    for bedisoname in isoDict:
        for isoid in isoDict[bedisoname].seqvariants:
            outline = [isoid, isoDict[bedisoname].seqvariants[isoid].pro]
            if isoDict[bedisoname].seqvariants[isoid].pro == 'PRO':
                bestorf = isoDict[bedisoname].seqvariants[isoid].bestorf
                orfstart, orfend = bestorf[4], bestorf[5]
                aaseq = isoDict[bedisoname].seqvariants[isoid].aaseq
                outline.extend([str(orfstart), str(orfend), aaseq])
            out.write('\t'.join(outline) + '\n')
    out.close()


def getorfs(args):
    logitmodelpath = args.logitmodel
    hexamerpath = args.hexamer

    pipettor.run([('cpat', '-g', args.isoforms, '-o', args.output_prefix + '.isoforms.cpat', '-d', logitmodelpath, '-x', hexamerpath)])
    print('done with cpat')

    isoDict = {}
    isobed = pybedtools.BedTool(args.bedisoforms)
    isoDict = addATGpos(isoDict, args.isoforms) ##generate iso dict and add atg positions
    isoDict = processBedFile(isoDict, isobed) ##add exon pos and ptc pos
    isoDict = addAnnotStarts(isoDict, args.gtf, isobed)
    print('done loading annotation')

    ##these alter the way upstream start codons modify the orf score
    atg_growth = 0.5
    atg_shift = 10

    ##adding no orfs isn't necessary, this is checked at the end
    # isoDict = addNoOrf(isoDict, args.output_prefix + '.isoforms.cpat.no_ORF.txt')
    isoDict = addOrfProb(isoDict, args.output_prefix + '.isoforms.cpat.ORF_prob.tsv', atg_growth, atg_shift)
    isoDict = calcBestOrf(isoDict)
    print('done calculating best orfs, writing output')

    printoutfile(isoDict, args.output_prefix + '.isoforms.orfs.tsv')

if __name__ == "__main__":
    args = getargs()
    getorfs(args)
