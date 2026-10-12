import argparse
from flair.flair_variantquant import (get_bedisoform_info, combine_vcf_files,
                                      group_annotated_ref_vars)
from flair.region_vars_tsv import RegionVarsWriter


def parse_args():
    parser = argparse.ArgumentParser(description='''This script is for annotating VCF files before passing them to FLAIR variantquant. Having pre-annotated files speeds FLAIR variantquant up considerably''')
    parser.add_argument('-i', '--bed_isoforms', required=True, help='bed file of isoforms, can be built from the reference GTF')
    parser.add_argument('-v', '--vcf', required=True, help='Vcf file with variant positions - can use dbgap, rediportal, or called variants from sample. Only requires the first few vcf columns, no header rows')
    parser.add_argument('-o', '--output', required=True, help='output name - Txt or TSV file, not standard format')
    args = parser.parse_args()
    return args

def annotate_vars_in_region(vcf_vars_for_region, chrom, region, out):
    for pos in vcf_vars_for_region:
        ref, alts, name = vcf_vars_for_region[pos]
        out.writeRow((chrom, region, pos, ref, alts, name))

def main():
    args = parse_args()

    print('retrieving gene info')
    iso_to_blocks, gene_to_iso, chr_region_to_genes, genes_to_boundaries = get_bedisoform_info(args.bed_isoforms)
    print('parsing vcf')
    var_to_alt = combine_vcf_files([args.vcf, ])
    print('annotating and grouping variants')
    vcfvars = group_annotated_ref_vars(var_to_alt, chr_region_to_genes, genes_to_boundaries, gene_to_iso, iso_to_blocks)
    print('outputting annotated variants')
    with RegionVarsWriter(args.output) as out:
        for chrom, region in vcfvars:
            annotate_vars_in_region(vcfvars[(chrom, region)], chrom, region, out)


if __name__ == "__main__":
    main()
