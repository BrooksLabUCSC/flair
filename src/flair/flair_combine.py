#! /usr/bin/env python3

import sys
import argparse
import os
import pipettor
import pysam
import math
from bed_to_gtf import bed_to_gtf
os.environ['OPENBLAS_NUM_THREADS'] = '1'


def bedReadToIntronChain(line): #line is a list of strings from a tab separated line
    dir, start, esizes, estarts = line[5], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
    introns = []
    for i in range(len(esizes) - 1):
        introns.append((start + estarts[i] + esizes[i], start + estarts[i + 1]))
    # if dir == '-': introns = [x[::-1] for x in introns[::-1]]
    return tuple(introns)

def intronChainToestarts(ichain, start, end):
    esizes, estarts = [], [0,]
    for i in ichain:
        esizes.append(i[0] - (start + estarts[-1]))
        estarts.append(i[1] - start)
    esizes.append(end - (start + estarts[-1]))
    return esizes, estarts

def getbestends2(isodata):
    # beststart, bestend = None, None
    # for info in isodata:
    #     start, end = info[0], info[1]
    #     if not beststart: beststart, bestend = start, end
    #     elif end-start > bestend-beststart: beststart, bestend = start, end
    # beststart, bestend, bestisousage, bestname = None, None, 0, None
    bestiso = (None, None, None, None, 0)
    for info in isodata:
        # start, end = info[0], info[1]
        # if info[-1] > bestisousage: beststart, bestend, bestisousage, bestname = start, end, info[-1]
        if info[-1] > bestiso[-1]: bestiso = info
        elif info[-1] == bestiso[-1]:
            if info[1]-info[0] > bestiso[1]-bestiso[0]: bestiso = info
    # return (beststart, bestend)
    return bestiso

def combineIsos2(isolist, endwindow):
    isolist.sort()
    isoendgroups = {}
    laststart, lastend = 0, 0
    currgroup = []
    for isoinfo in isolist:
        start, end = isoinfo[0], isoinfo[1]
        if start-laststart <= endwindow and end - lastend <= endwindow:
            currgroup.append(isoinfo)
        else:
            if len(currgroup) > 0:
                isoendgroups[getbestends2(currgroup)] = currgroup
            currgroup = [isoinfo]
        laststart, lastend = start, end
    if len(currgroup) > 0:
        isoendgroups[getbestends2(currgroup)] = currgroup
    return isoendgroups


revcomp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y',
	'Y':'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}


def revcomp(seq):
    rev_seq = []
    for i in reversed(range(len(seq))):
        rev_seq.append(revcomp_dict[seq[i]])
    return ''.join(rev_seq)


def combine():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--manifest', required=True, type=str,
                        help="path to manifest files that points to transcriptomes to combine. Each line of file should be tab separated with sample name, condition, batch, path/to/isoforms.bed, path/to/isoforms.fa, path/to/combined.isoform.read.map.txt. fa and read.map.txt files are not required, although if .fa files are not provided for each sample a .fa output will not be generated")
    parser.add_argument('-o', '--output_prefix', default='flair.combined.isoforms',
                        help="path to collapsed_output.bed file. default: 'collapsed_flairomes'")
    parser.add_argument('-w', '--endwindow', type=int, default=200,
                        help="window for comparing ends of isoforms with the same intron chain. Default:200bp")
    parser.add_argument('-p', '--minpercentusage', type=int, default=10,
                        help="minimum percent usage required in one sample to keep isoform in combined transcriptome. Default:10")
    parser.add_argument('-c', '--convert_gtf', action='store_true',
                        help="[optional] whether to convert the combined transcriptome bed file to gtf")
    parser.add_argument('-s', '--include_se', action='store_true',
                        help='whether to include single exon isoforms. Default: dont include')
    parser.add_argument('-f', '--filter', default='usageandlongest',
                        help='type of filtering. Options: usageandlongest(default), usageonly, none')

    args = parser.parse_args()
    manifest = args.manifest
    outprefix = args.output_prefix
    endwindow = int(args.endwindow)
    minpercentusage = int(args.minpercentusage) / 100.

    bedfiles, mapfiles, samples, fafiles = [], [], [], []
    for line in open(manifest):
        line = line.rstrip().split('\t')
        samples.append(line[0])
        bedfiles.append(line[1])
        if len(line) > 2: fafiles.append(line[2])
        else: fafiles.append('')
        if len(line) > 3: mapfiles.append(line[3])
        else: mapfiles.append('')

    generatefa = all([len(x) > 0 for x in fafiles])

    intronchaintoisos = {}
    sampletoseq = {}
    for i in range(len(samples)):
        sample = samples[i]
        genetoreads, isotoreads = {}, {}
        if mapfiles[i] != '':
            for line in open(mapfiles[i]):
                iso, reads = line.split('\t', 1)
                gene = iso.split('_')[-1]
                numreads = len(reads.split(','))
                if gene not in genetoreads:
                    genetoreads[gene] = numreads
                else:
                    genetoreads[gene] += numreads
                isotoreads[iso] = numreads
        for line in open(bedfiles[i]):
            line = line.rstrip().split('\t')
            chr, start, end, strand, isoname = line[0], int(line[1]), int(line[2]), line[5], line[3]
            gene = isoname.split('_')[-1]
            if int(line[9]) > 1:  ###removing single exon isoforms, may want to add this as a user input option later - although how am I handling single exon isoforms? Are they all getting stored as the same empty intron chain? that seems bad
                ichain = bedReadToIntronChain(line)
                ichainid = (chr, strand, gene, ichain)
                if mapfiles[i] != '':
                    isousage = isotoreads[isoname] / genetoreads[gene]
                else:
                    isousage = 1
                if ichainid not in intronchaintoisos: intronchaintoisos[ichainid] = []
                # intronchaintoisos[ichainid].append([start, end, sample, isoname, isotoreads[isoname], isousage])
                intronchaintoisos[ichainid].append((start, end, sample, isoname, isousage))
            elif args.include_se:
                ichain = chr + '-' + str(int(round(start, -4)))
                ichainid = (chr, strand, gene, ichain)
                if mapfiles[i] != '':
                    isousage = isotoreads[isoname] / genetoreads[gene]
                else:
                    isousage = 1
                if ichainid not in intronchaintoisos: intronchaintoisos[ichainid] = []
                intronchaintoisos[ichainid].append((start, end, sample, isoname, isousage))
        if generatefa:
            last = None
            sampletoseq[sample] = {}
            for line in open(fafiles[i]):
                if line[0] == '>': last = line[1:].rstrip()
                else: sampletoseq[sample][last] = line.rstrip()

    finalisostosupport = {}

    chromtobedinfo = {}
    isocount = 1
    outbed, outcounts = open(outprefix + '.bed', 'w'), open(outprefix + '.counts.tsv', 'w')
    if generatefa: outfa = open(outprefix + '.fa', 'w')
    for ichainid in intronchaintoisos:
        chr, strand, gene, ichain = ichainid
        collapsedIsos = combineIsos2(intronchaintoisos[ichainid], endwindow)
        longestEnds = (None, None)
        biggestdiff = 0
        maxintronchainusage = 0
        ichainendscount = 1
        for start, end, sample, isoname, isousage in collapsedIsos:
            if abs(end - start) > biggestdiff: longestEnds = (start, end)
            maxisousage = max([x[-1] for x in collapsedIsos[(start, end, sample, isoname, isousage)]])
            if maxisousage > maxintronchainusage: maxintronchainusage = maxisousage
        if args.filter == 'none' or maxintronchainusage > minpercentusage:
            for start, end, sample, isoname, isousage in collapsedIsos:
                theseisos = collapsedIsos[(start, end, sample, isoname, isousage)]
                theseisos.sort(key=lambda x: x[1] - x[0], reverse=True)  ##longest first
                maxisousage = max([x[-1] for x in theseisos])
                if args.filter == 'none' or maxisousage > minpercentusage or ((start, end) == longestEnds and type(
                        ichain) != str and args.filter == 'usageandlongest'):  # True:#
                    outname = None
                    for i in theseisos:
                        if i[3][:4] == 'ENST':
                            outname = str(isocount) + '-' + str(ichainendscount) + '_' + i[3]
                            break
                    if not outname: outname = 'flairiso' + str(isocount) + '-' + str(ichainendscount) + '_' + gene

                    ###output bed line
                    if type(ichain) == str:
                        esizes, estarts = [end - start], [0]
                    else:
                        esizes, estarts = intronChainToestarts(ichain, start, end)
                    outbed.write(
                        '\t'.join([chr, str(start), str(end), outname, '1000', strand, str(start), str(end), '0',
                                   str(len(esizes)), ','.join([str(x) for x in esizes]) + ',',
                                   ','.join([str(x) for x in estarts]) + ',']) + '\n')
                    if chr not in chromtobedinfo: chromtobedinfo[chr] = []
                    chromtobedinfo[chr].append([start, end, outname, strand, esizes, estarts])

                    ##output sequence
                    if generatefa:
                        isoseq = sampletoseq[sample][isoname]
                        outfa.write('>' + outname + '\n' + isoseq + '\n')
                else:
                    outname = 'lowexpiso_' + gene

                ##get counts
                if outname not in finalisostosupport: finalisostosupport[outname] = {x: 0 for x in samples}
                for isoinfo in theseisos:
                    finalisostosupport[outname][isoinfo[2]] += isoinfo[4]
                ichainendscount += 1
            isocount += 1
    outbed.close()
    if generatefa: outfa.close()

    outcounts.write('\t'.join(['ids'] + samples) + '\n')
    for name in finalisostosupport:
        outline = [name]
        for s in samples:
            outline.append(str(finalisostosupport[name][s]))
        outcounts.write('\t'.join(outline) + '\n')
    outcounts.close()

    if args.convert_gtf:
        bed_to_gtf(query=outprefix + '.bed', outputfile=outprefix + '.gtf')



if __name__ == "__main__":
    combine()
