import sys
import flair.flair_variantmodels as fv
from flair.flair_variantquant import _add_vcf_var, _get_correct_vcf_vars, group_annotated_ref_vars

bedfile = sys.argv[1] ###Isoforms, can be built from the reference GTF
vcffile = sys.argv[2] ##Only requires the first few vcf columns, no header rows
outfile = sys.argv[3] ##Txt or TSV file, not standard format

###This script is for annotating VCF files before passing them to FLAIR variantquant
###Having pre-annotated files speeds FLAIR variantquant up considerably


fv.add_vcf_var = _add_vcf_var
fv.get_correct_vcf_vars = _get_correct_vcf_vars
print('retrieving gene info')
isotoblocks, genetoiso, chrregiontogenes, genestoboundaries = fv.get_bedisoform_info(bedfile)
print('parsing vcf')
vartoalt = fv.combine_vcf_files([vcffile, ])
print('annotating and grouping variants')
vcfvars = group_annotated_ref_vars(vartoalt, chrregiontogenes, genestoboundaries, genetoiso, isotoblocks)
print('outputting annotated variants')
out = open(outfile, 'w')
for chrom, region in vcfvars:
    for pos in vcfvars[(chrom, region)]:
        ref, alts, name = vcfvars[(chrom, region)][pos]
        outline = [chrom, region, pos, ref, ','.join(alts), name]
        out.write('\t'.join([str(x) for x in outline]) + '\n')
