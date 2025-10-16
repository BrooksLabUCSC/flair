#! /usr/bin/env python3

import argparse
import os
import shutil
import math
import flair.flair_variantmodels as fv
os.environ['OPENBLAS_NUM_THREADS'] = '1'

compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def parse_var_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--manifest', type=str, #required=True, 
                        help="[USED INSTEAD OF input_bam AND vcf] path to manifest files that points to sample names, bam files aligned to transcriptome, "
                             "and vcf vars for that sample called on the genome. "
                             "Each line of file should be tab separated. "
                             "If you are using just one reference vcf file, "
                             "just include it in the third column for the first sample "
                             "and leave that column blank for the rest.")
    parser.add_argument('-i', '--input_bam', help='[USE INSTEAD OF MANIFEST] Path to bam file for individual sample')
    parser.add_argument('-r', '--pos_ref', help='[EITHER THIS, VCF, OR MANIFEST] Path to reference file of sites to check for variants')
    parser.add_argument('-v', '--vcf', help='[EITHER THIS, POS_REF, OR MANIFEST] Path to reference vcf file of sites to check for variants')
    parser.add_argument('-o', '--output_prefix', default='flair',
                        help="path to collapsed_output.bed file. default: 'flair'")
    parser.add_argument('-b', '--bedisoforms',
                        help="path to transcriptome bed file")
    parser.add_argument('-t', '--threshold', type=int, default=5, 
                        help='specify minimum total read coverage threshold to output a site')
    parser.add_argument('-k', '--output_all', action='store_true',
                        help="specify this option if you want to output read counts for all putative RNA editing sites that pass the coverage threshold, regardless of whether any reads are edited")
    parser.add_argument('--mode', default='genomic', help='type of alignment used - default "genomic", '
                                                          'but can specify "transcriptomic" if bam files in manifest '
                                                          'are aligned to the transcriptome (recommended to use FLAIR '
                                                          'quantify with --output_bam to do this)')
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.''')
    args = parser.parse_args()
    return args

# this is called by extract_vcf_vars
def _add_vcf_var(vcfvars, chrom, ref, alts, tpos2, name):
    roundedpos = math.floor(tpos2 / (10 ** 6)) * 10 ** 6
    poskey = (chrom, roundedpos)
    if poskey not in vcfvars:
        vcfvars[poskey] = {}
    vcfvars[poskey][tpos2] = (ref, alts, name)
    return vcfvars

def _get_correct_vcf_vars(vcfvars, refname, startpos, endpos):
    myvcfvars = {}
    for roundedpos in range(math.floor(startpos / (10 ** 6)) * 10 ** 6, math.floor(endpos / (10 ** 6)) * 10 ** 6, 10**6):
        if (refname, roundedpos) in vcfvars:
            myvcfvars.update(vcfvars[(refname, roundedpos)])
    return myvcfvars


def _parse_cigar(cigar, alignstart, transcriptvars):
    ref, quer = 0, 0
    coveredvars = {}
    for block in cigar:
        if block[0] in {0, 7, 8}:  # match, consumes both
            for pos in range(alignstart + ref + 1, alignstart + ref + block[1] + 1):
                if pos in transcriptvars:
                    coveredvars[pos] = 0
            ref += block[1]
            quer += block[1]
        elif block[0] in {1, 4}:  # consumes query ###1 is insertion
            quer += block[1]
        elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
            ref += block[1]
    return coveredvars


def _parse_single_bam_read(s, tempdir, vcfvars, sampleindex, tempfilename): ##add mode
    """for each read, figure out what variants it overlaps with. Then figure out whether it's modified or not at that variant"""
    # check for which vars are covered
    coveredvars = _parse_cigar(s.cigartuples, s.reference_start, vcfvars)

    # check for which covered vars are actually mod
    readseq = s.query_sequence
    alignedbases = s.get_aligned_pairs(with_seq=True, matches_only=True)
    for i in alignedbases:
        if i[2] and i[2].islower():
            if i[1] + 1 in coveredvars:
                ref, alts = vcfvars[i[1] + 1][:2]
                for alt in alts:  # could be multiple alt alleles
                    # restrict to same mutation from vcf, but not strand specific
                    if (i[2].upper() == ref and readseq[i[0]] == alt) \
                            or (compbase[i[2].upper()] == ref and compbase[readseq[i[0]]] == alt):
                        coveredvars[i[1] + 1] = 1

    if coveredvars:  # for now, only outputting reads that cover var pos
        coveredvarstrings = [','.join([str(v), vcfvars[v][2], str(coveredvars[v])]) for v in coveredvars]
        with open(tempdir + tempfilename + '.txt', 'a') as tempvarout:
            tempvarout.write(
                '\t'.join([s.reference_name, str(sampleindex) + '__' + s.query_name, ';'.join(coveredvarstrings)]) + '\n')


def read_vars_to_genome_pos_counts(tempfilenames, tempdir, outprefix, mode, sampledata, threshold, output_all):
    samplenames = [x[0] for x in sampledata]

    with open(f'{outprefix}.{mode}.var.counts.tsv', 'w') as out:
        out.write('\t'.join(['varpos', 'gene', 'transcript'] + samplenames) + '\n')
        vartocounts = {}
        for tf in tempfilenames:
            with open(tempdir + tf + '.txt', 'r') as mutstringfile:
                for line in mutstringfile:
                    refname, readname, mutstring = line.rstrip('\n').split('\t')
                    allmuts = [x.split(',') for x in mutstring.split(';')]
                    # each mut: position, varname (may be chrom:pos or gene), varstatus
                    sampleindex = int(readname.split('__')[0])
                    for m in allmuts:
                        if mode == 'genomic':
                            varpos = refname + ':' + m[0]
                            gene = m[1]
                            transcript = ''
                        else:
                            varpos = m[1]
                            gene = refname.split('_')[-1]
                            transcript = '_'.join(refname.split('_')[:-1])
                        var = (varpos, gene, transcript)
                        if var not in vartocounts:
                            vartocounts[var] = [[0, 0] for x in range(len(samplenames))]  # [unmod counts, mod counts]
                        vartocounts[var][sampleindex][int(m[-1])] += 1
        for var in vartocounts:
            if any([x[0] + x[1] >= threshold for x in vartocounts[var]]) and (any([x[1] > 0 for x in vartocounts[var]]) or output_all):  # any modified reads in any sample
                varcounts = [f'{x[0]};{x[1]}' for x in vartocounts[var]]
                outline = list(var) + varcounts
                out.write('\t'.join(outline) + '\n')


fv.add_vcf_var = _add_vcf_var
fv.parse_single_bam_read = _parse_single_bam_read
fv.get_correct_vcf_vars = _get_correct_vcf_vars

def quantvarpos():
    args = parse_var_args()
    # Load reference data
    if args.manifest:
        sampledata = fv.extract_sample_data(args.manifest)
    elif args.input_bam and (args.pos_ref or args.vcf):
        if args.vcf:
            sampledata = [['sample', args.input_bam, args.vcf]]
        else:
            sampledata = [['sample', args.input_bam]]
    else:
        raise ValueError("please provide either manifest or bam and vcf")
    isotoblocks, genetoiso, chrregiontogenes, genestoboundaries = fv.get_bedisoform_info(args.bedisoforms)
    print('done loading annot')

    vcfvars = {}
    if args.manifest or args.vcf:
        vartoalt = fv.combine_vcf_files([x[2] for x in sampledata if len(x) > 2])
        print('done combining vcfs')
        if args.mode == 'transcriptomic':
            vcfvars = fv.convert_vars_to_tpos(vartoalt, isotoblocks, genetoiso, chrregiontogenes, genestoboundaries)
        else:
            for refinfo in vartoalt:
                chrom, gpos, ref, alts = fv.extract_varinfo(refinfo, vartoalt[refinfo])

                potgenes = fv.get_potential_genes(chrom, gpos, chrregiontogenes)

                overlapgenes = set()
                ##THIS MAY BE THE BOTTLENECK
                for gene, _, _ in fv.retrieve_good_iso_pos(potgenes, genestoboundaries, gpos, genetoiso, isotoblocks):
                    overlapgenes.add(gene)
                vcfvars = _add_vcf_var(vcfvars, chrom, ref, alts, gpos, ','.join(overlapgenes))
    else:
        for line in open(args.pos_ref):
            line = line.rstrip('\n').split('\t')
            chrom, region, pos, ref, alt, name = line
            region, pos = int(region), int(pos)
            alts = alt.split(',')
            poskey = (chrom, region)
            if poskey not in vcfvars:
                vcfvars[poskey] = {}
            vcfvars[poskey][pos] = (ref, alts, name)

    print('done combining vcf variants')
    tempdir = fv.make_temp_dir(args.output_prefix)
    fv.parse_all_bam_files(sampledata, tempdir, vcfvars, args.mode)  # parses to intermediate files with read name to all vars
    print('parsed all reads')
    genenames = fv.get_genes_from_tempdir(tempdir)
    read_vars_to_genome_pos_counts(genenames, tempdir, args.output_prefix, args.mode, sampledata, args.threshold, args.output_all)

    if not args.keep_intermediate:
        shutil.rmtree(tempdir)


if __name__ == "__main__":
    quantvarpos()
