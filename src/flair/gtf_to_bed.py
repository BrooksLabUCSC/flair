#!/usr/bin/env python3
import argparse
from flair.gtf_io import gtf_record_parser, GtfAttrsSet
from flair.pycbio.hgdata.bed import Bed, BedBlock


def main():
    parser = argparse.ArgumentParser(description='''converts a gtf to a bed, depending on the output filename extension;
            gtf exons need to be grouped by transcript and sorted by coordinate w/in a transcript''')
    required = parser.add_argument_group('required named arguments')
    required.add_argument('gtf', type=str, help='annotated gtf')
    required.add_argument('bed', type=str, help='bed file')
    parser.add_argument('--include_gene', action='store_true', dest='include_gene', required=False,
                        help='''Include gene name in the isoform name''')
    args = parser.parse_args()

    gtf_to_bed(args.bed, args.gtf, args.include_gene)

def write_bed_row(include_gene, iso_to_cds, prev_transcript, blockstarts, blocksizes, prev_gene, prev_chrom, prev_strand, fh):
    blockcount = len(blockstarts)
    if blockcount > 1 and blockstarts[0] > blockstarts[1]:  # need to reverse exons
        blocksizes = blocksizes[::-1]
        blockstarts = blockstarts[::-1]

    tstart, tend = blockstarts[0], blockstarts[-1] + blocksizes[-1]  # target (e.g. chrom)
    if include_gene:
        qname = prev_transcript + '_' + prev_gene
    else:
        qname = prev_transcript

    if qname in iso_to_cds:
        cds_start, cds_end = iso_to_cds[qname]
    else:
        cds_start, cds_end = tstart, tend
    blocks = [BedBlock(blockstarts[i], blockstarts[i] + blocksizes[i]) for i in range(blockcount)]
    Bed(prev_chrom, tstart, tend, name=qname, score=1000, strand=prev_strand,
        thickStart=cds_start, thickEnd=cds_end, itemRgb='0', blocks=blocks).write(fh)


def get_iso_info(gtf):
    iso_to_cds = {}
    iso_to_exons = {}
    iso_to_info = {}
    for rec in gtf_record_parser(gtf, include_features={'exon', 'CDS'}, attrs=GtfAttrsSet.ALL):
        if rec.feature == 'CDS':
            if rec.transcript_id not in iso_to_cds:
                iso_to_cds[rec.transcript_id] = [rec.start, rec.end]
            elif rec.end > iso_to_cds[rec.transcript_id][1]:
                iso_to_cds[rec.transcript_id][1] = rec.end
        elif rec.feature == 'exon':
            if rec.transcript_id not in iso_to_exons:
                iso_to_exons[rec.transcript_id] = []
                gene_id = rec.gene_id.replace('_', '-')
                iso_to_info[rec.transcript_id] = (rec.chrom, rec.strand, gene_id)
            iso_to_exons[rec.transcript_id].append((rec.start, rec.end))
    return iso_to_info, iso_to_exons, iso_to_cds

def gtf_to_bed(outputfile, gtf, include_gene=False):

    with open(outputfile, 'wt') as outfile:
        iso_to_info, iso_to_exons, iso_to_cds = get_iso_info(gtf)

        for this_transcript in iso_to_exons:
            chrom, strand, gene = iso_to_info[this_transcript]
            exons = sorted(iso_to_exons[this_transcript])
            blockstarts = [x[0] for x in exons]
            blocksizes = [x[1] - x[0] for x in exons]

            write_bed_row(include_gene, iso_to_cds, this_transcript, blockstarts, blocksizes, gene, chrom, strand, outfile)


if __name__ == "__main__":
    main()
