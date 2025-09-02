#!/usr/bin/env python3

########################################################################
# File: predictProductivity.py
#  executable: predictProductivity.py
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 10/09/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import sys
import pipettor
import os
from flair import FlairInputDataError
########################################################################
# CommandLine
########################################################################


class CommandLine(object):
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    and a standard usage and help,
    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']

    methods:

    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description=' predictProductivity - a tool.')
        # Add args
        self.parser.add_argument('-i', "--input_isoforms", action='store', required=True, help='Input collapsed isoforms in bed12 format.')
        self.parser.add_argument('-g', "--gtf", action='store', required=True, help='Gencode annotation file.')
        self.parser.add_argument('-f', "--genome_fasta", action='store', required=True, help='Fasta file containing transcript sequences.')
        self.parser.add_argument('-o', "--output", action='store', required=True, help='prefix of output files')
        self.parser.add_argument("--quiet", action='store_false', required=False, default=True, help='Do not display progress')
        self.parser.add_argument("--append_column", action='store_true', required=False, default=False, help='Append prediction as an additional column in file')

        self.group = self.parser.add_mutually_exclusive_group(required=True)
        self.group.add_argument('--firstTIS', action='store_true', default=False, help='Defined ORFs by the first annotated TIS.')
        self.group.add_argument('--longestORF',action='store_true', default=False, help='Defined ORFs by the longest open reading frame.')

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# Isoform
########################################################################

class Isoform(object):
    '''
    Object to handle isoform related data.

    attributes:

    methods:

    '''

    def __init__(self, name=None, seq=None):
        self.name = name
        self.pro = "UNK"
        self.chrom = ""

        self.sequence  = seq
        self.exons     = set()
        self.starts    = set()
        self.orfs      = list()
        self.exonSizes = list()
        self.ptcpoint = ''
        self.allEsizes = list()
        self.allExons = list()


########################################################################
# MAIN
########################################################################


def getStarts(gtf):
    # starts = list()
    starts = 'tmp_predictProd_starts.bed'
    scount = 0
    out = open(starts, 'w')
    tnamenmdexcep = set()
    with open(gtf) as lines:
        for l in lines:
            if l[0] == "#": continue
            cols = l.rstrip().split("\t")
            chrom, c1, c2, strand = cols[0], int(cols[3])-1, int(cols[4]), cols[6]
            if cols[2] == "start_codon":
                gene = cols[8][cols[8].find('gene_id')+len('gene_id')+2:]
                gene = gene[:gene.find('"')]
                scount += 1
                out.write('\t'.join([str(x) for x in [chrom,c1,c2,gene,".",strand]]) + '\n')
                # starts.append((chrom,c1,c2,gene,".",strand))
            if cols[2] == 'transcript':
                if 'NMD_exception' in cols[8]:
                    transcript = cols[8].split('transcript_id "')[1].split('"')[0]
                    tnamenmdexcep.add(transcript)
    out.close()
    if scount == 0:
        raise FlairInputDataError(f'ERROR, no start codons were found in {gtf}')
    return starts, tnamenmdexcep


def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
        splitchar = '_chr'
    elif '_XM' in iso_gene:
        splitchar = '_XM'
    elif '_XR' in iso_gene:
        splitchar = '_XR'
    elif '_NM' in iso_gene:
        splitchar = '_NM'
    elif '_NR' in iso_gene:
        splitchar = '_NR'
    elif '_R2_' in iso_gene:
        splitchar = '_R2_'
    elif '_NC_' in iso_gene:
        splitchar = '_NC_'
    else:
        splitchar = '_'
    iso = iso_gene[:iso_gene.rfind(splitchar)]
    gene = iso_gene[iso_gene.rfind(splitchar)+1:]
    return iso, gene


def getSeqs(bed, genome):

    isoDict = dict()
    dr = pipettor.DataReader()
    bedtools_cmd = ('bedtools', 'getfasta', '-fi', genome, '-bed', bed, '-tab', '-split', '-s', '-name')
    pipettor.run([bedtools_cmd], stdout=dr)

    for entry in dr.data.split('\n'):
        if len(entry) > 0:
            read, seq = entry.split()
            # accommodate different bedtools versions - they use different separators

            iso = read.split('::')[0]
            iso = iso.split("(")[0]
            if iso[:10] == 'fusiongene': iso = '_'.join(iso.split('_')[1:])
            if iso not in isoDict:
                isoDict[iso] = Isoform(iso,seq)
            else:
                isoDict[iso].sequence = isoDict[iso].sequence + seq
    return isoDict


def getStartRelPos(genomicStartPos,exon, exons, isoObj):
    '''
    is handed a genomic position, the exon it occurs in, all exons,
    and returns the position relative to all exons
    '''
    exonNum = exons.index(exon)
    isoObj.exonSizes = [x[1]-x[0] for x in exons]

    # First get start position relative to transcript sequence.
    if isoObj.strand == "+":
        relativeStart = genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])
    elif isoObj.strand == "-":
        # print('calc', sum(isoObj.exonSizes), genomicStartPos - exons[exonNum][0], sum(isoObj.exonSizes[:exonNum]), sum(isoObj.exonSizes) - (genomicStartPos - exons[exonNum][0] + sum(isoObj.exonSizes[:exonNum])))
        relativeStart = sum(isoObj.exonSizes) - (genomicStartPos - exons[exonNum][0] + sum(isoObj.exonSizes[:exonNum])) #- 3

    return relativeStart


def checkPTC(orfEndPos, exonSizes, allExons, nmdexcep, isoname):
    '''
    takes a transcript sequence position, and list of exon sizes to detemine
    if that position occurs more than 55nucleotides away from a splice junction.
    ptc = True if yes, ptc = False if not.
    the genomic position is also reported.
    '''
    stopDistFromExon = None
    exonWithStop = None
    ptc  = None
    genomicPos = int()
    distance   = 0
    maxdistfromexonedge = 55

    for num,e in enumerate(exonSizes,0):

        distance += e

        # if the stop codon is in the last exon, then not ptc.
        if num == len(exonSizes)-1:
            ptc = False
            if exonWithStop is None:
                exonWithStop = num
                stopDistFromExon = distance - orfEndPos

        # if the distance is greater than the stop position, then check if
        # the difference in distance is more than 55nt
        # if yet then ptc = True
        # also, track which exon the stop codon is in to get genomic position
        elif orfEndPos < distance:
            distToJunc = distance - orfEndPos
            if exonWithStop is None:
                exonWithStop = num
                stopDistFromExon = int(distToJunc)

            # allow stop codon in second to last exon if fusion breakpoint between second to last and last exon
            if allExons[-2][3] != allExons[-1][3]:
                if num == len(exonSizes)-2:
                    ptc = False
                    break
                elif distToJunc > maxdistfromexonedge or num < len(exonSizes)-3:
                    ptc = True
                    break

            if distToJunc > maxdistfromexonedge or num < len(exonSizes)-2:
                ptc = True
                break

    if len(exonSizes) == 1: ptcpointont = 0
    elif allExons[-2][3] != allExons[-1][3]:
        if len(exonSizes) == 2: ptcpointont = 0
        elif exonSizes[-3] < maxdistfromexonedge: ptcpointont = sum(exonSizes[:-3])
        else: ptcpointont = sum(exonSizes[:-2]) - maxdistfromexonedge
    elif exonSizes[-2] < maxdistfromexonedge: ptcpointont = sum(exonSizes[:-2])
    else: ptcpointont = sum(exonSizes[:-1]) - maxdistfromexonedge

    isoname = '_'.join(isoname.split('_')[:-1])
    if isoname[-2] == '-': isoname = isoname[:-2]
    if isoname in nmdexcep:
        ptc, ptcpointont = False, 0

    exonsWithStop = allExons[exonWithStop]
    left,right,strand,fusionindex    = exonsWithStop

    genomicPos = right - stopDistFromExon if strand == "+" else left + stopDistFromExon

    return genomicPos, ptc, ptcpointont


def get_exons(bedline):
    gstart = int(bedline[1])
    esizes, estarts = [int(x) for x in bedline[10].rstrip(',').split(',')], \
                      [int(x) for x in bedline[11].rstrip(',').split(',')],
    return [(gstart+estarts[i], gstart+estarts[i]+esizes[i]) for i in range(len(esizes))]


def predict(bed, starts, isoDict, nmdexcep):
    for line in open(bed):
        line = line.rstrip().split('\t')
        read   = line[3]
        strand = line[5]
        for exonCoord in get_exons(line):
            elen = exonCoord[1]-exonCoord[0]
            if read[:10] == 'fusiongene':
                fusionindex = read.split('_')[0]
                read = '_'.join(read.split('_')[1:])
            else: fusionindex = "NA"
            if strand == '+':
                isoDict[read].allEsizes.append(elen)
                isoDict[read].allExons.append((exonCoord[0], exonCoord[1], strand, fusionindex))
            elif strand == '-':
                isoDict[read].allEsizes.insert(0, elen)
                isoDict[read].allExons.insert(0, (exonCoord[0], exonCoord[1], strand, fusionindex))

    dr = pipettor.DataReader()
    bedtools_cmd = ('bedtools', 'intersect', '-a', bed, '-b', starts, '-split', '-s', '-wao')
    pipettor.run([bedtools_cmd], stdout=dr)
    os.remove(starts)

    for intersection in dr.data.split('\n'):
        if len(intersection) > 0:
            intersection = intersection.split('\t')
            read   = intersection[3]
            if read[:10] == 'fusiongene' and read[10] != '1': continue # only getting starts for 5' genes
            if read[:10] == 'fusiongene': read = '_'.join(read.split('_')[1:])
            overlap  = intersection[-1]
            goStart  = int(intersection[-6]) if intersection[5] == '+' else int(intersection[-5])
            if intersection[5] != intersection[-2]: overlap = '0' ##if start is not on same strand, doesn't count
            isoDict[read].strand = intersection[5]
            isoDict[read].chrom = intersection[0]

            for exonCoord in get_exons(intersection):
                isoDict[read].exons.add(exonCoord)
                if overlap == "3" and exonCoord[0] <= goStart <= exonCoord[1]:
                    isoDict[read].starts.add((exonCoord,goStart))

    stops = set(['TAA','TGA','TAG'])
    for iso,o in isoDict.items():
        exons = list(o.exons)
        exons.sort()
        if len(o.starts) < 1:
            o.orfs.append(["NGO", exons[0][0], exons[0][0], 0, 0])

        else:
            for start in o.starts:
                exon,startPos = start
                relativeStart = getStartRelPos(startPos,exon,exons,o)
                fiveUTR,rest  = o.sequence[:relativeStart], o.sequence[relativeStart:].upper()
                # Next find first stop codon
                stopReached = False
                for i in range(0, len(rest), 3):
                    if rest[i:i+3] in stops:
                        stopReached = True
                        break

                # i is the last position after going through all codons and breaking at a stop
                # is a stop was never reached then i should represent the last NT in the entire seq
                # therefore, i+3 should be longer than the entire potential orf is a stop was never reached.
                # lets call these nonstop, or nst for now.
                if not stopReached:
                    orfEndPos = len(fiveUTR)+i
                    o.orfs.append(["NST", startPos, exons[-1][-1] if o.strand == "+" else exons[0][0], orfEndPos-relativeStart, relativeStart])

                # else if a stop was reached...
                else:
                    orfEndPos = len(fiveUTR)+i+3
                    genomicStopPos, ptc, ptcdecidingpoint = checkPTC(orfEndPos, o.allEsizes, o.allExons, nmdexcep, iso)
                    ptc = "PTC" if ptc else "PRO"
                    o.orfs.append([ptc, startPos, genomicStopPos, orfEndPos - relativeStart, relativeStart])
                    o.ptcpoint = ptcdecidingpoint

    return isoDict


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
    seq = seq.upper()
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


def main():
    '''
    maine
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()
    bed    = myCommandLine.args['input_isoforms']
    genome = myCommandLine.args['genome_fasta']
    gtf    = myCommandLine.args['gtf']
    output = myCommandLine.args['output']
    extra_col = myCommandLine.args['append_column']

    if myCommandLine.args['firstTIS']:
        defineORF = 'first'
    elif myCommandLine.args['longestORF']:
        defineORF = 'longest'
    else:
        raise FlairInputDataError('** ERR. Select method for ORF definition with --firstTIS or --longestORF')

    starts, nmdexcep      = getStarts(gtf)
    isoformObjs = getSeqs(bed, genome)
    isoformObjs = predict(bed, starts, isoformObjs, nmdexcep)

    beaut = {"PRO":"103,169,207", "PTC":"239,138,98", "NST":"0,0,0","NGO":"0,0,0"}

    bedout = open(output + '.bed', 'w')
    infoout = open(output + '.info.tsv', 'w')
    infoout.write('\t'.join(['#isoname', 'tstartont', 'tendont', 'ptcpointont', 'AAseq']) + '\n')
    with open(bed) as lines:
        for line in lines:
            bedCols = line.rstrip().split()
            if bedCols[3][:10] == 'fusiongene': isoname = '_'.join(bedCols[3].split('_')[1:])
            else: isoname = bedCols[3]
            isoObj = isoformObjs[isoname]

            if defineORF == 'longest':
                isoObj.orfs.sort(key=lambda x: x[3], reverse=True)
            elif defineORF == 'first':
                isoObj.orfs.sort(key=lambda x: x[4])
            pro,start,end,orfLen, tisPos = isoObj.orfs[0]

            if extra_col:
                bedCols += [pro]
            else:
                iso, gene = split_iso_gene(bedCols[3])
                bedCols[3] = "%s_%s_%s" % (iso, pro, gene)

            bedCols[8] = beaut[pro]
            if 'fusiongene' in bedCols[3]:
                if not int(bedCols[1]) < start < int(bedCols[2]) and not int(bedCols[1]) < end < int(bedCols[2]):
                    bedCols[6],bedCols[7] = bedCols[1], bedCols[1]
                else:
                    if isoObj.strand == "+":
                        if int(bedCols[1]) < start < int(bedCols[2]): bedCols[6] = str(start)
                        if int(bedCols[1]) < end < int(bedCols[2]): bedCols[7] = str(end)
                    else:
                        if int(bedCols[1]) < start < int(bedCols[2]): bedCols[7] = str(start)
                        if int(bedCols[1]) < end < int(bedCols[2]): bedCols[6] = str(end)
            else:
                if isoObj.strand == "+":
                    bedCols[6],bedCols[7] = str(start),str(end)
                else:
                    bedCols[7],bedCols[6] = str(start),str(end)
            bedout.write("\t".join(bedCols) + '\n')
            infoout.write('\t'.join([bedCols[3], str(tisPos), str(tisPos+orfLen), str(isoObj.ptcpoint), translate(isoObj.sequence[tisPos:tisPos+orfLen])]) + '\n')



if __name__ == "__main__":
    main()
