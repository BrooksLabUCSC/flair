#!/usr/bin/env python3
import sys
import argparse
from flair import FlairInputDataError

def main():
    parser = argparse.ArgumentParser(description='options')
    parser.add_argument('inputfile', type=str,
            action='store', help='isoforms in bed format')
    parser.add_argument('--force', action='store_true', dest='force',
            help='specify to not split isoform name by underscore into isoform and gene ids')
    parser.add_argument('--add_reference_transcript_id', action='store_true', dest='reference_transcript_id',
            help='specify to add reference_transcript_id attribute')
    parser.add_argument('--dontuseCDS', action='store_true')
    args = parser.parse_args()
    bed_to_gtf(query=args.inputfile, force=args.force, outputfile='/dev/stdout',
               reference_transcript_id=args.reference_transcript_id, useCDS= not args.dontuseCDS)


def split_iso_gene(iso_gene):
    if '_chr' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_chr')]
        gene = iso_gene[iso_gene.rfind('_chr')+1:]
    elif '_XM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XM')]
        gene = iso_gene[iso_gene.rfind('_XM')+1:]
    elif '_XR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_XR')]
        gene = iso_gene[iso_gene.rfind('_XR')+1:]
    elif '_NM' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NM')]
        gene = iso_gene[iso_gene.rfind('_NM')+1:]
    elif '_NR' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_NR')]
        gene = iso_gene[iso_gene.rfind('_NR')+1:]
    elif '_R2_' in iso_gene:
        iso = iso_gene[:iso_gene.rfind('_R2_')]
        gene = iso_gene[iso_gene.rfind('_R2_')+1:]
    else:
        iso = iso_gene[:iso_gene.rfind('_')]
        gene = iso_gene[iso_gene.rfind('_')+1:]
    return iso, gene

def bed_to_gtf(query, outputfile, force=False, reference_transcript_id=False, useCDS=True):
    outfile = open(outputfile, 'w')
    for line in open(query):
        line = line.rstrip().split('\t')
        start = int(line[1])
        chrom, strand, score, name, start = line[0], line[5], line[4], line[3], int(line[1])
        tstarts = [int(n) + start for n in line[11].rstrip(',').split(',')]
        bsizes = [int(n) for n in line[10].rstrip(',').split(',')]
        end, thick_start, thick_end = int(line[2]), int(line[6]), int(line[7])

        if '_' not in name and not force:
            raise FlairInputDataError('Entry name should contain underscore-delimited transcriptid and geneid like so: \n'
                             'ENST00000318842.11_ENSG00000156313.12 or a4bab8a3-1d28_chr8:232000\n'
                             'So no GTF conversion was done. Please run identify_gene_isoform first\n'
                             'for best results, or run with --force')

        if ';' in name:
            name = name.replace(';', ':')

        if force == True:
            transcript_id, gene_id = name, name
        else:
            transcript_id, gene_id = split_iso_gene(name)

        attributes = 'gene_id \"{}\"; transcript_id \"{}\";'\
                                .format(gene_id, transcript_id)
        if reference_transcript_id and '-referencetranscript' in transcript_id:
            trimmed_transcript_id = transcript_id[:transcript_id.find('-referencetranscript')]
            attributes = 'gene_id \"{}\"; transcript_id \"{}\"; reference_transcript_id \"{}\";'\
                    .format(gene_id, trimmed_transcript_id, trimmed_transcript_id)
        print('\t'.join([chrom, 'FLAIR', 'transcript', str(start+1), str(tstarts[-1]+bsizes[-1]), '.', strand, '.',
                attributes]), file=outfile)
        if thick_start != thick_end and (thick_start != start or thick_end != end) and useCDS:
            print('\t'.join([chrom, 'FLAIR', 'CDS', str(thick_start+1), str(thick_end), '.', strand, '.',
            attributes]), file=outfile)
            if strand == '+':
                print('\t'.join([chrom, 'FLAIR', 'start_codon', str(thick_start+1), str(thick_start+3), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '5UTR', str(start+1), str(thick_start+1), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '3UTR', str(thick_end), str(tstarts[-1]+bsizes[-1]), '.', strand, '.',
                attributes]), file=outfile)
            elif strand == '-':
                print('\t'.join([chrom, 'FLAIR', 'start_codon', str(thick_end-2), str(thick_end), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '3UTR', str(start+1), str(thick_start+1), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '5UTR', str(thick_end), str(tstarts[-1]+bsizes[-1]), '.', strand, '.',
                attributes]), file=outfile)
        # if strand == '-':  # to list exons in 5'->3'
        #       for b in range(len(tstarts)):  # exon number
        #               bi = len(tstarts) - 1 - b  # block index
        #               attributes = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";'\
        #                                               .format(gene_id, transcript_id, b)
        #               print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[bi]+1), \
        #                       str(tstarts[bi]+bsizes[bi]), '.', strand, '.', attributes]))
        # else:
        for b in range(len(tstarts)):
            attributes = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";'\
                    .format(gene_id, transcript_id, b)
            if reference_transcript_id and '-referencetranscript' in transcript_id:
                attributes = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\"; reference_transcript_id \"{}\";'\
                        .format(gene_id, trimmed_transcript_id, b, trimmed_transcript_id)
            print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[b]+1),
                    str(tstarts[b]+bsizes[b]), '.', strand, '.', attributes]), file=outfile)

if __name__ == "__main__":
    main()
