#!/usr/bin/env python3
import argparse
from flair import FlairInputDataError
from flair.gtf_io import gtf_write_row
from flair.pycbio.hgdata.bed import BedReader


def main():
    parser = argparse.ArgumentParser(description='options')
    parser.add_argument('inputfile', type=str,
                        action='store', help='isoforms in bed format')
    parser.add_argument('--force', action='store_true', dest='force',
                        help='specify to not split isoform name by underscore into isoform and gene ids')
    parser.add_argument('--add_reference_transcript_id', action='store_true', dest='reference_transcript_id',
                        help='specify to add reference_transcript_id attribute')
    parser.add_argument('--noCDS', action='store_true',
                        help='do not carry forward CDS from bed file (thickstart and thickend) to gtf file')
    parser.add_argument('--gene_col', type=int, help='if bed file is bed12+, specifiy index of column with gene ID in context of extra columns - if gene_id is in column 13, enter index 0, since it is the first extra column')
    parser.add_argument('--as_file', help='if input is bed12 and want to add extra columns from bed to gtf (besides gene ID), specify path to as file with column definitions')
    parser.add_argument('--extra_cols_to_use', help='comma separated list of additional bed columns to carry over to gtf in context of extra columns (see gene_col)')
    args = parser.parse_args()
    if args.as_file:
        extra_cols = [int(x) for x in args.extra_cols_to_use.split(',')]
        c = 0
        my_fields = []
        for line in open(args.as_file):
            c += 1
            if c > 3 and line[0] != ')':
                line = line.rstrip().split('\t')
                this_data = (line[0], line[1].rstrip(';'), line[2].strip('"'))
                my_fields.append(this_data)
        extracolindexnames = []
        for i in range(12, len(my_fields)):
            if i - 12 in extra_cols:
                extracolindexnames.append((i - 12, my_fields[i][1]))
        args.extra_cols_to_use = extracolindexnames
    bed_to_gtf(query=args.inputfile, force=args.force, outputfile='/dev/stdout',
               reference_transcript_id=args.reference_transcript_id, useCDS=not args.noCDS, genecol=args.gene_col, extracolindexnames=args.extra_cols_to_use)

def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_chr')]
        gene = iso_gene[iso_gene.rfind('_chr') + 1:]
    elif '_XM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XM')]
        gene = iso_gene[iso_gene.rfind('_XM') + 1:]
    elif '_XR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XR')]
        gene = iso_gene[iso_gene.rfind('_XR') + 1:]
    elif '_NM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NM')]
        gene = iso_gene[iso_gene.rfind('_NM') + 1:]
    elif '_NR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NR')]
        gene = iso_gene[iso_gene.rfind('_NR') + 1:]
    elif '_R2_' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_R2_')]
        gene = iso_gene[iso_gene.rfind('_R2_') + 1:]
    else:
        iso = iso_gene[:iso_gene.rfind('_')]
        gene = iso_gene[iso_gene.rfind('_') + 1:]
    return iso, gene

def _make_attrs(gene_id, transcript_id, reference_transcript_id, extraCols, extracolindexnames, exon_number=None):
    """Build attrs dict for GTF output."""
    attrs = {'gene_id': gene_id, 'transcript_id': transcript_id}
    if exon_number is not None:
        attrs['exon_number'] = str(exon_number)
    if reference_transcript_id is not None:
        attrs['reference_transcript_id'] = reference_transcript_id
    if extracolindexnames is not None:
        for index, name in extracolindexnames:
            if extraCols[index] != 'NA':
                attrs[name] = extraCols[index]
    return attrs

def bed_to_gtf(query, outputfile, force=False, reference_transcript_id=False, useCDS=True, genecol=None, extracolindexnames=None):  # noqa: C901 - FIXME: reduce complexity
    outfile = open(outputfile, 'w')
    gene_to_records = {}
    gene_to_chrom_strand = {}
    for bed in BedReader(query):
        name = bed.name
        if '_' not in name and not force and genecol is None:
            raise FlairInputDataError('Entry name should contain underscore-delimited transcriptid and geneid like so: \n'
                                      'ENST00000318842.11_ENSG00000156313.12 or a4bab8a3-1d28_chr8:232000\n'
                                      'So no GTF conversion was done. Please run identify_gene_isoform first\n'
                                      'for best results, or run with --force')

        if ';' in name:
            name = name.replace(';', ':')

        if genecol is not None:
            transcript_id = name
            gene_id = bed.extraCols[genecol]
        elif force is True:
            transcript_id, gene_id = name, name
        else:
            transcript_id, gene_id = split_iso_gene(name)

        ref_tid = None
        if reference_transcript_id and '-referencetranscript' in transcript_id:
            ref_tid = transcript_id[:transcript_id.find('-referencetranscript')]
            transcript_id = ref_tid

        if gene_id not in gene_to_records:
            gene_to_records[gene_id] = []
            gene_to_chrom_strand[gene_id] = (bed.chrom, bed.strand)

        attrs = _make_attrs(gene_id, transcript_id, ref_tid, bed.extraCols, extracolindexnames)
        gene_to_records[gene_id].append(('transcript', bed.chrom, bed.chromStart, bed.chromEnd, bed.strand, attrs))

        for b, blk in enumerate(bed.blocks):
            exon_attrs = _make_attrs(gene_id, transcript_id, ref_tid, bed.extraCols, extracolindexnames, exon_number=b)
            gene_to_records[gene_id].append(('exon', bed.chrom, blk.start, blk.end, bed.strand, exon_attrs))
            if bed.thickStart != bed.thickEnd and (bed.thickStart != bed.chromStart or bed.thickEnd != bed.chromEnd) and useCDS:
                if bed.thickStart < blk.start and blk.end < bed.thickEnd:  # in CDS
                    gene_to_records[gene_id].append(('CDS', bed.chrom, blk.start, blk.end, bed.strand, exon_attrs))
                elif blk.end < bed.thickStart:  # fully left of CDS
                    if bed.strand == '+':
                        gene_to_records[gene_id].append(('5UTR', bed.chrom, blk.start, blk.end, bed.strand, exon_attrs))
                    elif bed.strand == '-':
                        gene_to_records[gene_id].append(('3UTR', bed.chrom, blk.start, blk.end, bed.strand, exon_attrs))
                elif blk.start > bed.thickEnd:  # fully right of CDS
                    if bed.strand == '+':
                        gene_to_records[gene_id].append(('3UTR', bed.chrom, blk.start, blk.end, bed.strand, exon_attrs))
                    elif bed.strand == '-':
                        gene_to_records[gene_id].append(('5UTR', bed.chrom, blk.start, blk.end, bed.strand, exon_attrs))
                elif blk.start <= bed.thickStart <= blk.end:  # left end of CDS in exon
                    if bed.strand == '+':
                        gene_to_records[gene_id].append(('5UTR', bed.chrom, blk.start, bed.thickStart, bed.strand, exon_attrs))
                        gene_to_records[gene_id].append(('start_codon', bed.chrom, bed.thickStart, bed.thickStart + 3, bed.strand, attrs))
                        gene_to_records[gene_id].append(('CDS', bed.chrom, bed.thickStart, blk.end, bed.strand, exon_attrs))
                    elif bed.strand == '-':
                        gene_to_records[gene_id].append(('3UTR', bed.chrom, blk.start, bed.thickStart, bed.strand, exon_attrs))
                        gene_to_records[gene_id].append(('CDS', bed.chrom, bed.thickStart, blk.end, bed.strand, exon_attrs))
                elif blk.start <= bed.thickEnd <= blk.end:  # right end of CDS in exon
                    if bed.strand == '+':
                        gene_to_records[gene_id].append(('CDS', bed.chrom, blk.start, bed.thickEnd, bed.strand, exon_attrs))
                        gene_to_records[gene_id].append(('3UTR', bed.chrom, bed.thickEnd, blk.end, bed.strand, exon_attrs))
                    elif bed.strand == '-':
                        gene_to_records[gene_id].append(('CDS', bed.chrom, blk.start, bed.thickEnd, bed.strand, exon_attrs))
                        gene_to_records[gene_id].append(('start_codon', bed.chrom, bed.thickEnd - 3, bed.thickEnd, bed.strand, attrs))
                        gene_to_records[gene_id].append(('5UTR', bed.chrom, bed.thickEnd, blk.end, bed.strand, exon_attrs))

    for gene_id, records in gene_to_records.items():
        chrom, strand = gene_to_chrom_strand[gene_id]
        gene_start = min(r[2] for r in records)
        gene_end = max(r[3] for r in records)
        gtf_write_row(outfile, chrom, 'FLAIR', 'gene', gene_start, gene_end, None, strand, None,
                      gene_id=gene_id)
        for feature, chrom, start, end, strand, attrs in records:
            gtf_write_row(outfile, chrom, 'FLAIR', feature, start, end, None, strand, None,
                          attrs=attrs)
    outfile.close()


if __name__ == "__main__":
    main()
