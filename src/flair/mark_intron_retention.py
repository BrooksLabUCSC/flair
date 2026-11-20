#!/usr/bin/env python3
import argparse
import logging
from flair.cli_help import FlairHelpFormatter
from flair.marked_isoforms_bed import MarkedIsoformsWriter, bed_column_count
from flair.interval_index import IntervalIndex
from flair.io_utils import bed_row_sort_key
from flair.pycbio.sys import cli, loggingOps
from flair.retained_introns_tsv import RetainedIntronsWriter
from flair.pycbio.hgdata.bed import BedReader
from flair.run_checks import check_input_files
from flair.run_logging import configure_logging

def build_parser():
    desc = "Mark isoforms that retain an intron of another isoform of the same locus"
    parser = argparse.ArgumentParser(prog='mark_intron_retention', formatter_class=FlairHelpFormatter, description=desc)
    parser.add_argument('isoform_bed', help='isoforms in bed format')
    parser.add_argument('marked_bed', help='output bed, with intron retention marked in an extra column')
    parser.add_argument('intron_txt', help='output text file of the retained introns')
    loggingOps.addCmdOptions(parser, defaultLevel=logging.INFO)
    return parser


# an intron shorter than this is not counted: an exon spanning it is an
# alignment difference, not intron retention
MIN_RETAINED_INTRON = 100
# how far an exon may fall short of an intron's ends and still retain it
RETAINED_INTRON_SLOP = 10


def _introns(isoform):
    "the introns of an isoform long enough to count, as (start, end) pairs"
    starts, sizes = isoform['starts'], isoform['sizes']
    ends = [start + size for start, size in zip(starts, sizes)]
    return [(ends[i], starts[i + 1]) for i in range(len(starts) - 1)
            if starts[i + 1] - ends[i] >= MIN_RETAINED_INTRON]


def _retains_intron(isoform, intron_start, intron_end):
    "does an exon of the isoform span the intron"
    return any((start < intron_start + RETAINED_INTRON_SLOP) and (start + size > intron_end - RETAINED_INTRON_SLOP)
               for start, size in zip(isoform['starts'], isoform['sizes']))


def _mark_retained_introns(chrom, retaining, other, introncoords):
    "mark retaining if one of its exons spans an intron of other"
    for intron_start, intron_end in _introns(other):
        if _retains_intron(retaining, intron_start, intron_end):
            retaining['ir'] = True
            introncoords.add((chrom, intron_start, intron_end, retaining['strand']))


def _chrom_isoform_index(chrom_isoforms):
    """The isoforms of one chromosome indexed by their span.  A scan of every pair
    was quadratic in the isoforms of a chromosome, where only the overlapping ones
    can retain each other's introns."""
    index = IntervalIndex()
    for name, iso in chrom_isoforms.items():
        index.add(iso['range'][0], iso['range'][1], (name, iso))
    index.build()
    return index


def _mark_chrom_intron_retention(chrom, chrom_isoforms, introncoords):
    "the isoforms of one chromosome against those overlapping them on the same strand"
    index = _chrom_isoform_index(chrom_isoforms)
    for name0, iso0 in chrom_isoforms.items():
        for name1, iso1 in index.overlap(iso0['range'][0], iso0['range'][1]):
            if (name0 != name1) and (iso0['strand'] == iso1['strand']):
                _mark_retained_introns(chrom, iso0, iso1, introncoords)


def mark_intron_retention(isoform_bed, marked_bed, intron_txt):
    outfilename = marked_bed
    txtout = intron_txt
    isoforms = {}
    for bed in BedReader(isoform_bed, fixScores=True):
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
    for chrom, chrom_isoforms in isoforms.items():
        _mark_chrom_intron_retention(chrom, chrom_isoforms, introncoords)

    marked = sorted([(isoforms[chrom][name]['entry'], isoforms[chrom][name]['ir'])
                     for chrom in isoforms for name in isoforms[chrom]],
                    key=lambda entry_ir: bed_row_sort_key(entry_ir[0]))
    with MarkedIsoformsWriter(outfilename, bed_column_count([e for e, _ in marked], isoform_bed)) as writer:
        for entry, retains_intron in marked:
            writer.writeIsoform(entry, retains_intron)

    # sorted, since introncoords is a set and an unordered file cannot be diffed
    with RetainedIntronsWriter(txtout) as writer:
        for intron in sorted(introncoords):
            writer.writeRow(intron)


def main():
    args = cli.parseArgsWithLogging(build_parser())
    configure_logging()
    with cli.ErrorHandler():
        check_input_files({'isoform BED': args.isoform_bed})
        mark_intron_retention(args.isoform_bed, args.marked_bed, args.intron_txt)


if __name__ == "__main__":
    main()
