#! /usr/bin/env python3

import sys
import argparse


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
    if len(seq) % 3 != 0: seq = seq[:(len(seq)//3)*3]
    # print(len(seq), len(seq) % 3, (len(seq)%3)*3)
    # print(seq)
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if 'N' in codon: protein += '?'
            else: protein += table[codon]

    stoppos = protein.find('_')
    # print(protein)
    if stoppos >= 0: protein = protein[:stoppos+1]
    return protein


class TranscriptInfo:
    def __init__(self, aaseq, s, e, p, ptc):
        self.aaseq = aaseq
        self.origstart = s
        self.origstop = e
        self.origpro = p
        self.ptcpoint = ptc

referencetranscriptfile = sys.argv[1]
modtranscriptsfile = sys.argv[2]
outfile = sys.argv[3]


transcriptToInfo = {}
for line in open(referencetranscriptfile):
    if line[0] != '#':
        line = line.rstrip().split('\t')
        if len(line) == 3: line += ['','']
        elif len(line) == 4: line.append('')
        tinfo, orfstart, orfstop, ptcpoint, seq = line
        tinfo = tinfo.split('_')
        tname, propred, gname = '_'.join(tinfo[:-2]), tinfo[-2], tinfo[-1]
        if 'fusiongene' in tname: tname = '_'.join(tname.split('_')[1:])
        transcriptToInfo[(tname, gname)] = TranscriptInfo(seq, int(orfstart), int(orfstop), propred, ptcpoint)

cl = modtranscriptsfile.split('-')[0]
out = open(outfile, 'w')



last = None
for line in open(modtranscriptsfile):
    if line[0] == '>': last = line[1:].rstrip().split(' ')
    else:
        predProd = None
        modseq = line.rstrip()
        # tinfo, genomevars, seqvars, support = last
        if len(last) == 2: tinfo, seqvars = last
        else: tinfo, seqvars = last[0], ''
        if seqvars == 'nomuts': seqvars = ''
        temp = '-'.join(tinfo.split('-')[1:])
        tname = '_'.join(temp.split('_')[:-1])
        gname = temp.split('_')[-1]
        if last[0] == '1-flairiso23973-1_ENSG00000105173.14--chr19-27666000': print(tname, gname, (tname, gname) in transcriptToInfo)
        if (tname, gname) not in transcriptToInfo: continue
        thist = transcriptToInfo[(tname, gname)]
        newpredseq = thist.aaseq
        has5utrvars, has3utrvars = False, False
        # print(tname)
        if len(seqvars) == 0: ###no change from original isoform
            predProd = thist.origpro
        else:
            # chr2..-..65312507..2862..S..A
            #chr2..-..65312507..2862..S..A,chr2..-..65313771..1598..S..C
            seqvars = [x.split('..') for x in seqvars.split(',')]
            for i in range(len(seqvars)):
                #seqvars[i][0] = int(seqvars[i][0])
                seqvars[i] = [int(seqvars[i][3]), seqvars[i][4], seqvars[i][5]]
            seqvars.sort()

            if thist.origstart == thist.origstop:
                thist.origstart = len(modseq)
                thist.origstop = len(modseq)
            if thist.ptcpoint == '': thist.ptcpoint = 0
            else: thist.ptcpoint = int(thist.ptcpoint)

            posToVar = {}
            for pos, ref, alt in seqvars:
                if ref == 'I' or ref == 'D':
                    posToVar[pos] = (ref, alt)
            refseqpostomodseqpos = {}
            modseqpos, refseqpos = 0,0
            while modseqpos <= len(modseq):
                if modseqpos in posToVar:
                    ref, alt = posToVar[modseqpos]
                    if ref == 'I': modseqpos += len(alt)
                    elif ref == 'D':
                        for i in range(int(alt)):
                            refseqpostomodseqpos[refseqpos+i] = modseqpos
                        refseqpos += int(alt)
                refseqpostomodseqpos[refseqpos] = modseqpos
                modseqpos += 1
                refseqpos += 1
            newptcpoint = refseqpostomodseqpos[thist.ptcpoint]
            ###Is there a base change before the predicted start?
            prestartvars = [x for x in seqvars if x[0] < thist.origstart]
            hasnovelstart = False
            startposchange = 0
            newstart, newend = thist.origstart, thist.origstop

            if len(prestartvars) > 0:
                for pos, ref, alt in prestartvars:
                    checkstart, checkend = pos-2, pos+3
                    if ref == 'I':
                        checkend += len(alt)
                        startposchange += len(alt)
                    elif ref == 'D': startposchange -= int(alt)
                    checkseq = modseq[checkstart:checkend]
                    hasstart = checkseq.find('ATG')
                    if hasstart >= 0:
                        newstart=checkstart + hasstart
                        newpredseq = translate(modseq[newstart:])
                        hasnovelstart = True
                        if newpredseq[-1] != '_':
                            predProd = 'NST'
                        else:
                            newend = newstart + ((len(newpredseq))*3)
                            if newend < newptcpoint:
                                predProd = 'PTC'
                            else:  ##same length or longer
                                predProd = 'PRO'
                        continue
                        ####Add code to get new end and AA sequence and predict productivity
                        ###decide if want to look at all potential novel starts or just the first one
            if len([x for x in seqvars if x[0] < newstart]) > 0 and newstart != len(modseq): has5utrvars = True ####only if there is actually a start
            if not hasnovelstart and thist.origstart == len(modseq):
                predProd = 'NGO'
            elif not hasnovelstart and thist.origstart != len(modseq): #do not execute if there was no original AA sequence
                poststartvars = [x for x in seqvars if x[0] >= thist.origstart]
                if len(poststartvars) > 0:
                    newpredseq = translate(modseq[refseqpostomodseqpos[thist.origstart]:refseqpostomodseqpos[thist.origstop]])
                    if len(newpredseq) == 0:
                        predProd = 'NGO'
                    elif newpredseq[-1] != '_':
                        predProd = 'NST'
                    else:
                        newend = newstart + ((len(newpredseq)) * 3)
                        if newend < newptcpoint:
                            predProd = 'PTC'  ###Do we care about how early the stop is? Maybe as long as protein is still 50% of original length?
                        else:
                            predProd = 'PRO'
                else:
                    predProd = thist.origpro
            if len([x for x in seqvars if x[0] >= newend]) > 0 and newend != len(modseq): has3utrvars = True
        utrvars = []
        if has5utrvars: utrvars.append('5utr')
        if has3utrvars: utrvars.append('3utr')
        out.write('\t'.join([tinfo, predProd, ','.join(utrvars),  newpredseq]) + '\n')
