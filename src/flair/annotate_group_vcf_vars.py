import flair.flair_variantmodels as fv
from flair.flair_variantquant import _add_vcf_var, _get_correct_vcf_vars, group_annotated_ref_vars
import argparse


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
        outline = [chrom, region, pos, ref, ','.join(alts), name]
        out.write('\t'.join([str(x) for x in outline]) + '\n')

def main():
    args = parse_args

    fv.add_vcf_var = _add_vcf_var
    fv.get_correct_vcf_vars = _get_correct_vcf_vars
    print('retrieving gene info')
    isotoblocks, genetoiso, chrregiontogenes, genestoboundaries = fv.get_bedisoform_info(args.bed_isoforms)
    print('parsing vcf')
    vartoalt = fv.combine_vcf_files([args.vcf, ])
    print('annotating and grouping variants')
    vcfvars = group_annotated_ref_vars(vartoalt, chrregiontogenes, genestoboundaries, genetoiso, isotoblocks)
    print('outputting annotated variants')
    out = open(args.output, 'w')
    for chrom, region in vcfvars:
        annotate_vars_in_region(vcfvars[(chrom, region)], chrom, region, out)
        

if __name__ == "__main__":
    main()
