"""Write a --pos_ref file of flair variantquant: the SNV sites of a vcf, each with the
genes of a flair isoform bed or a GTF it is in.  variantquant loads it much faster than the vcf
and bed, for a database of millions of sites, such as the RNA editing sites of
REDIportal, counted again and again"""
import argparse
from flair import FlairInputDataError
from flair.flair_variantquant import KEEP_ALL_FLAGS, name_sites_by_genes, read_vcf_sites, site_exons
from flair.pycbio.sys import cli


def build_parser():
    parser = argparse.ArgumentParser(prog='annotate_group_vcf_vars', description=__doc__.split('\n\n')[0])
    parser.add_argument('-i', '--isoform_bed', '--bed_isoforms',
                        help="bed file of isoforms whose exons name the sites by their genes, each its 13th column, as "
                             "a flair isoform bed's, or, without it, the --gtf gene it overlaps most; as flair "
                             'variantquant --isoform_bed')
    parser.add_argument('-f', '--gtf',
                        help='GTF annotation file, plain or gzipped, whose exons name the sites by their gene_ids, or, '
                             'with --isoform_bed, the genes of its isoforms without them')
    parser.add_argument('-v', '--vcf', required=True,
                        help='vcf file, plain or gzipped, of the sites: of a database such as REDIportal or dbSNP, or '
                             'called in the samples; only its first five columns are needed.  Only SNVs of records '
                             'whose FILTER is PASS or . are written, or one of --keep_flags')
    parser.add_argument('--keep_flags', type=lambda flags: set(flags.split(',')), default=set(), metavar='FLAGS',
                        help=f'comma-separated vcf FILTER flags whose records are written too, or {KEEP_ALL_FLAGS}, '
                             'every record')
    parser.add_argument('-o', '--output', required=True,
                        help="the --pos_ref file to write, of flair variantquant's format")
    return parser


def write_pos_ref(vcfvars, path):
    """write the sites of vcfvars, (chrom, bin) -> pos -> (ref, alts, genes), a line each:
    chromosome, bin, position, ref, comma-separated alts, comma-separated genes"""
    with open(path, 'w') as out:
        for (chrom, region), sites in vcfvars.items():
            for pos, (ref, alts, genes) in sorted(sites.items()):
                out.write('\t'.join([chrom, str(region), str(pos), ref, ','.join(alts), genes]) + '\n')


def annotate_group_vcf_vars(vcf, output, *, isoform_bed=None, gtf=None, keep_flags=frozenset()):
    if isoform_bed is None and gtf is None:
        raise FlairInputDataError('give --isoform_bed, --gtf or both, the genes to name the sites with')
    write_pos_ref(name_sites_by_genes(read_vcf_sites([vcf], keep_flags), site_exons(isoform_bed, gtf)), output)


def main():
    args = build_parser().parse_args()
    with cli.ErrorHandler():
        annotate_group_vcf_vars(args.vcf, args.output, isoform_bed=args.isoform_bed, gtf=args.gtf,
                                keep_flags=args.keep_flags)


if __name__ == '__main__':
    main()
