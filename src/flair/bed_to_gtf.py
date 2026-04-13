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
    args = parser.parse_args()
    bed_to_gtf(query=args.inputfile, force=args.force, outputfile='/dev/stdout',
               reference_transcript_id=args.reference_transcript_id, useCDS=not args.noCDS)


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


def _make_attrs(gene_id, transcript_id, reference_transcript_id, exon_number=None):
    """Build attrs dict for GTF output."""
    attrs = {'gene_id': gene_id, 'transcript_id': transcript_id}
    if exon_number is not None:
        attrs['exon_number'] = str(exon_number)
    if reference_transcript_id is not None:
        attrs['reference_transcript_id'] = reference_transcript_id
    return attrs


def bed_to_gtf(query, outputfile, force=False, reference_transcript_id=False, useCDS=True):  # noqa: C901 - FIXME: reduce complexity
    outfile = open(outputfile, 'w')
    gene_to_records = {}
    gene_to_chrom_strand = {}
    for bed in BedReader(query):
        name = bed.name
        if '_' not in name and not force:
            raise FlairInputDataError('Entry name should contain underscore-delimited transcriptid and geneid like so: \n'
                                      'ENST00000318842.11_ENSG00000156313.12 or a4bab8a3-1d28_chr8:232000\n'
                                      'So no GTF conversion was done. Please run identify_gene_isoform first\n'
                                      'for best results, or run with --force')

        if ';' in name:
            name = name.replace(';', ':')

        if force is True:
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

        attrs = _make_attrs(gene_id, transcript_id, ref_tid)
        gene_to_records[gene_id].append(('transcript', bed.chrom, bed.chromStart, bed.chromEnd, bed.strand, attrs))

        if bed.thickStart != bed.thickEnd and (bed.thickStart != bed.chromStart or bed.thickEnd != bed.chromEnd) and useCDS:
            gene_to_records[gene_id].append(('CDS', bed.chrom, bed.thickStart, bed.thickEnd, bed.strand, attrs))
            if bed.strand == '+':
                gene_to_records[gene_id].append(('start_codon', bed.chrom, bed.thickStart, bed.thickStart + 3, bed.strand, attrs))
                gene_to_records[gene_id].append(('5UTR', bed.chrom, bed.chromStart, bed.thickStart, bed.strand, attrs))
                gene_to_records[gene_id].append(('3UTR', bed.chrom, bed.thickEnd, bed.chromEnd, bed.strand, attrs))
            elif bed.strand == '-':
                gene_to_records[gene_id].append(('start_codon', bed.chrom, bed.thickEnd - 3, bed.thickEnd, bed.strand, attrs))
                gene_to_records[gene_id].append(('3UTR', bed.chrom, bed.chromStart, bed.thickStart, bed.strand, attrs))
                gene_to_records[gene_id].append(('5UTR', bed.chrom, bed.thickEnd, bed.chromEnd, bed.strand, attrs))

        for b, blk in enumerate(bed.blocks):
            exon_attrs = _make_attrs(gene_id, transcript_id, ref_tid, exon_number=b)
            gene_to_records[gene_id].append(('exon', bed.chrom, blk.start, blk.end, bed.strand, exon_attrs))

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
