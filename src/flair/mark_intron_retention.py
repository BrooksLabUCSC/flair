#!/usr/bin/env python3
import argparse
import csv
import os
from flair.pycbio.sys import cli
from flair.pycbio.hgdata.bed import BedReader

def overlap(coords0, coords1):
    """Do two closed ranges share any position.  The old test asked only whether an
    end of coords1 fell inside coords0, so a coords1 containing coords0 came back
    False, which is the very shape an intron retention pair has."""
    return (coords1[0] <= coords0[1]) and (coords0[0] <= coords1[1])


def build_parser():
    desc = "Mark isoforms that retain an intron of another isoform of the same locus"
    parser = argparse.ArgumentParser(prog='mark_intron_retention', description=desc)
    parser.add_argument('isoform_bed', help='isoforms in bed format')
    parser.add_argument('marked_bed', help='output bed, with intron retention marked in an extra column')
    parser.add_argument('intron_txt', help='output text file of the retained introns')
    return parser

def mark_intron_retention(isoform_bed, marked_bed, intron_txt):  # noqa C901
    bedfh = open(isoform_bed)
    outfilename = marked_bed
    txtout = intron_txt
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


def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        mark_intron_retention(args.isoform_bed, args.marked_bed, args.intron_txt)


if __name__ == "__main__":
    main()
