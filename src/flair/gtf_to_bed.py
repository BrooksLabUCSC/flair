#!/usr/bin/env python3
import argparse
from flair.gtf_io import gtf_record_parser, GtfAttrsSet
from flair.pycbio.sys import cli
from flair.pycbio.hgdata.bed import Bed, BedBlock
from flair.io_utils import open_sorted_bed

# with --include_gene, the column after the twelve standard ones holds the gene
GENE_ID_EXTRA_COL = 0


def build_parser():
    desc = ('converts a gtf to a bed, depending on the output filename extension; '
            'gtf exons need to be grouped by transcript and sorted by coordinate within a transcript')
    parser = argparse.ArgumentParser(prog='gtf_to_bed', description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('gtf', type=str, help='annotated gtf')
    required.add_argument('bed', type=str, help='bed file')
    parser.add_argument('--include_gene', action='store_true',
                        help='add a column naming the gene of each isoform')
    return parser

def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        gtf_to_bed(args.bed, args.gtf, args.include_gene)

def write_bed_row(include_gene, iso_to_cds, prev_transcript, blockstarts, blocksizes, prev_gene, prev_chrom, prev_strand, fh):
    blockcount = len(blockstarts)
    if blockcount > 1 and blockstarts[0] > blockstarts[1]:  # need to reverse exons
        blocksizes = blocksizes[::-1]
        blockstarts = blockstarts[::-1]

    tstart, tend = blockstarts[0], blockstarts[-1] + blocksizes[-1]  # target (e.g. chrom)
    if prev_transcript in iso_to_cds:
        cds_start, cds_end = iso_to_cds[prev_transcript]
    else:
        cds_start, cds_end = tstart, tend
    blocks = [BedBlock(blockstarts[i], blockstarts[i] + blocksizes[i]) for i in range(blockcount)]
    # the gene goes in a column of its own rather than into the name, which could only
    # be split again with a heuristic that is wrong for a gene id containing '_'.  One
    # extra column rather than a full FLAIR BED: bedtools rejects the empty fields of
    # the columns a FLAIR BED would add
    extra_cols = (prev_gene,) if include_gene else ()
    Bed(prev_chrom, tstart, tend, name=prev_transcript, score=1000, strand=prev_strand,
        thickStart=cds_start, thickEnd=cds_end, itemRgb='0', blocks=blocks,
        extraCols=extra_cols).write(fh)


def get_iso_info(gtf):
    iso_to_cds = {}
    iso_to_exons = {}
    iso_to_info = {}
    for rec in gtf_record_parser(gtf, include_features={'exon', 'CDS'}, attrs=GtfAttrsSet.ALL):
        if rec.feature == 'CDS':
            # both ends: GTF need not be coordinate sorted and the parser does not
            # sort CDS records, so taking the start from the first record seen put
            # thickStart inside the CDS whenever the lines were out of order
            if rec.transcript_id not in iso_to_cds:
                iso_to_cds[rec.transcript_id] = [rec.start, rec.end]
            else:
                cds = iso_to_cds[rec.transcript_id]
                cds[0] = min(cds[0], rec.start)
                cds[1] = max(cds[1], rec.end)
        elif rec.feature == 'exon':
            if rec.transcript_id not in iso_to_exons:
                iso_to_exons[rec.transcript_id] = []
                gene_id = rec.gene_id
                iso_to_info[rec.transcript_id] = (rec.chrom, rec.strand, gene_id)
            iso_to_exons[rec.transcript_id].append((rec.start, rec.end))
    return iso_to_info, iso_to_exons, iso_to_cds

def gtf_to_bed(outputfile, gtf, include_gene=False):

    with open_sorted_bed(outputfile) as outfile:
        iso_to_info, iso_to_exons, iso_to_cds = get_iso_info(gtf)

        for this_transcript in iso_to_exons:
            chrom, strand, gene = iso_to_info[this_transcript]
            exons = sorted(iso_to_exons[this_transcript])
            blockstarts = [x[0] for x in exons]
            blocksizes = [x[1] - x[0] for x in exons]

            write_bed_row(include_gene, iso_to_cds, this_transcript, blockstarts, blocksizes, gene, chrom, strand, outfile)


if __name__ == "__main__":
    main()
