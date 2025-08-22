#!/usr/bin/env python3
import sys
import csv
import os
import argparse

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

def write_bed_row(include_gene, iso_to_cds, prev_transcript, blockstarts, blocksizes, prev_gene, prev_chrom, prev_strand, writer):
    blockcount = len(blockstarts)
    if blockcount > 1 and blockstarts[0] > blockstarts[1]:  # need to reverse exons
        blocksizes = blocksizes[::-1]
        blockstarts = blockstarts[::-1]

    tstart, tend = blockstarts[0], blockstarts[-1] + blocksizes[-1]  # target (e.g. chrom)
    qsize = sum(blocksizes)  # query (e.g. transcript)
    if include_gene:
        qname = prev_transcript+'_'+prev_gene
    else:
        qname = prev_transcript

    blocksizes = ','.join([str(b) for b in blocksizes]) + ','

    relblockstarts = [block - tstart for block in blockstarts]
    relblockstarts = ','.join([str(b) for b in relblockstarts]) + ','
    if qname in iso_to_cds:
        cds_start, cds_end = iso_to_cds[qname]
    else:
        cds_start, cds_end = tstart, tend
    writer.writerow([prev_chrom, tstart, tend, qname, 1000, prev_strand, cds_start,
                     cds_end, 0, blockcount, blocksizes, relblockstarts])


def get_iso_info(gtf, adjustpos = True):
    iso_to_cds = {}
    iso_to_exons = {}
    iso_to_info = {}
    # prev_transcript, blockstarts, blocksizes, prev_gene, prev_chrom, prev_strand = [None, None, None, None, None, None]
    for line in open(gtf):  # extract all exons from the gtf, keep exons grouped by transcript
        if line.startswith('#') or (len(line.rstrip()) == 0):
            continue
        line = line.rstrip().split('\t')
        chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
        if adjustpos: start -= 1
        if ty == 'CDS' or ty == 'exon':
            this_transcript = line[8][line[8].find('transcript_id') + 15:]
            this_transcript = this_transcript[:this_transcript.find('"')]
            if ty == 'CDS':
                if this_transcript not in iso_to_cds:
                    iso_to_cds[this_transcript] = [start, end]
                elif end > iso_to_cds[this_transcript][1]:
                    iso_to_cds[this_transcript][1] = end
            if ty == 'exon':
                if this_transcript not in iso_to_exons:
                    iso_to_exons[this_transcript] = []
                    prev_gene = line[8][line[8].find('gene_id') + 9:]
                    prev_gene = prev_gene[:prev_gene.find('"')]
                    prev_gene = prev_gene.replace('_', '-')
                    iso_to_info[this_transcript] = (chrom, strand, prev_gene)
                iso_to_exons[this_transcript].append((start, end))
    return iso_to_info, iso_to_exons, iso_to_cds

def gtf_to_bed(outputfile, gtf, include_gene=False):


    with open(outputfile, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        iso_to_info, iso_to_exons, iso_to_cds = get_iso_info(gtf)


        for this_transcript in iso_to_exons:
            chrom, strand, gene = iso_to_info[this_transcript]
            exons = sorted(iso_to_exons[this_transcript])
            blockstarts = [x[0] for x in exons]
            blocksizes = [x[1]-x[0] for x in exons]

            write_bed_row(include_gene, iso_to_cds, this_transcript, blockstarts, blocksizes, gene, chrom, strand, writer)





if __name__ == "__main__":
    main()
