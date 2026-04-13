#!/usr/bin/env python3
import sys
import csv
import os
from flair import FlairInputDataError
from flair.pycbio.hgdata.bed import BedReader

try:
    bedfh = open(sys.argv[1])
    outfilename = sys.argv[2]
    txtout = sys.argv[3]
except Exception:
    raise FlairInputDataError('usage: mark_intron_retention in.bed out_isoforms.bed out_introns.txt')


def overlap(coords0, coords1):
    return coords1[0] >= coords0[0] and coords1[0] <= coords0[1] or \
        coords1[1] >= coords0[0] and coords1[1] <= coords0[1]


isoforms = {}
for bed in BedReader(bedfh, fixScores=True):
    chrom, name, start, end, strand = bed.chrom, bed.name, bed.chromStart, bed.chromEnd, bed.strand
    blockstarts = [blk.start for blk in bed.blocks]
    blocksizes = [len(blk) for blk in bed.blocks]
    if chrom not in isoforms:
        isoforms[chrom] = {}
    isoforms[chrom][name] = {}
    isoforms[chrom][name]['entry'] = bed.toRow()
    isoforms[chrom][name]['sizes'] = blocksizes
    isoforms[chrom][name]['strand'] = strand
    isoforms[chrom][name]['starts'] = blockstarts
    isoforms[chrom][name]['range'] = start, end
    isoforms[chrom][name]['ir'] = False  # detection of intron retention event

introncoords = set()
allcoords = set()

for chrom in isoforms:
    for iname0 in isoforms[chrom]:
        for iname1 in isoforms[chrom]:  # compare with all other isoforms to find IR
            if iname0 == iname1:
                continue
            if not overlap(isoforms[chrom][iname0]['range'], isoforms[chrom][iname1]['range']):
                continue
            starts0, sizes0 = isoforms[chrom][iname0]['starts'], isoforms[chrom][iname0]['sizes']
            starts1, sizes1 = isoforms[chrom][iname1]['starts'], isoforms[chrom][iname1]['sizes']
            prev5 = starts1[0] + sizes1[0]  # previous 5' end of isoform1's intron
            for start1, size1 in zip(starts1[1:], sizes1[1:]):
                if start1 - prev5 < 100:  # do not count exons spanning introns smaller than 100 bp
                    prev5 = start1 + size1
                    continue
                for start0, size0 in zip(starts0, sizes0):
                    allcoords.add((chrom, str(prev5), (start1), isoforms[chrom][iname1]['strand']))
                    if start0 < prev5 + 10 and start0 + size0 > start1 - 10:  # if isoform 0 has exon where isoform 1 has intron
                        isoforms[chrom][iname0]['ir'] = True
                        introncoords.add((chrom, str(prev5), (start1), isoforms[chrom][iname0]['strand']))
                prev5 = start1 + size1

with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    for chrom in isoforms:
        for name in isoforms[chrom]:
            if isoforms[chrom][name]['ir']:
                writer.writerow(isoforms[chrom][name]['entry'] + [1])
            else:
                writer.writerow(isoforms[chrom][name]['entry'] + [0])

with open(txtout, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
    for intron in introncoords:
        writer.writerow(intron)
