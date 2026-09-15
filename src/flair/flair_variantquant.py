#! /usr/bin/env python3

import argparse
import os
import shutil
import math
from contextlib import ExitStack
import pysam
from flair.pycbio.hgdata.bed import BedReader
from flair.flair_bed import FlairBed
from flair.io_utils import make_temp_dir
os.environ['OPENBLAS_NUM_THREADS'] = '1'

compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W',
            'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}

# genes are binned by the start coordinate of the gene for the variant lookup
GENE_BIN_SIZE = 10 ** 7

# a temp record is pos, gene list, status; the gene list is itself a list, so the two
# need separators that cannot be confused, or the field count varies with gene count
VAR_FIELD_SEP = '|'
GENE_LIST_SEP = ','


def parse_var_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--manifest', type=str,  # required=True,
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
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.''')
    args = parser.parse_args()
    return args

def extract_sample_data(manifestfile):
    sampledata = []
    for line in open(manifestfile):
        line = line.rstrip().split()
        sampledata.append(line)
    return sampledata

def extract_varinfo(refinfo, alt):
    chrom, gpos, id, ref = refinfo
    alt = list(alt)
    gpos = int(gpos)
    return chrom, gpos, ref, alt

def get_potential_genes(chrom, gpos, chrregiontogenes):
    """Genes that might contain gpos.  chrregiontogenes is keyed by the bin of each
    gene's start, and a gene containing gpos starts at or before it, so the bin that
    matters beside gpos's own is the one before, not the one after.  A gene longer
    than the bin is still missed; the longest human gene is a fifth of one."""
    this_bin = math.floor(gpos / GENE_BIN_SIZE) * GENE_BIN_SIZE
    potgenes = set()
    for bin_start in (this_bin - GENE_BIN_SIZE, this_bin):
        potgenes.update(chrregiontogenes.get((chrom, bin_start), ()))
    return potgenes

def process_bedline(bed):
    thischr, iso, dir, start, end = bed.chrom, bed.name, bed.strand, bed.chromStart, bed.chromEnd
    gene = bed.gene_id
    esizes = [len(blk) for blk in bed.blocks]
    estarts = [blk.start - start for blk in bed.blocks]
    exonblocks = []  # block is gstart, tstart, len
    if dir == '-':
        esizes = esizes[::-1]
        estarts = estarts[::-1]
    currtstart = 0
    for i in range(len(esizes)):
        exonblocks.append((currtstart, estarts[i] + start, esizes[i], thischr, dir))
        currtstart += esizes[i]
    exonblocks.sort()

    return thischr, iso, dir, start, esizes, estarts, end, exonblocks, currtstart, gene


def add_gene_to_boundaries(genestoboundaries, gene, thischr, dir, start, end):
    if gene not in genestoboundaries:
        genestoboundaries[gene] = [thischr, dir, start, end]
    else:
        if start < genestoboundaries[gene][2]:
            genestoboundaries[gene][2] = start
        if end > genestoboundaries[gene][3]:
            genestoboundaries[gene][3] = end
    return genestoboundaries

def add_iso_to_blocks(isotoblocks, iso, exonblocks):
    if iso not in isotoblocks:
        isotoblocks[iso] = exonblocks
    else:
        tstartadj = isotoblocks[iso][-1][0] + isotoblocks[iso][-1][2]
        exonblocks = [(x[0] + tstartadj,) + x[1:] for x in exonblocks]
        isotoblocks[iso].extend(exonblocks)
    return isotoblocks

def get_bedisoform_info(bedisofile):
    isotoblocks = {}
    genetoiso = {}
    chrregiontogenes, genestoboundaries = {}, {}
    for bed in BedReader(bedisofile, bedClass=FlairBed):
        thischr, iso, dir, start, esizes, estarts, end, exonblocks, currtstart, gene = process_bedline(bed)

        if iso[:10] == 'fusiongene':
            iso = '_'.join(iso.split('_')[1:])

        isotoblocks = add_iso_to_blocks(isotoblocks, iso, exonblocks)

        if gene not in genetoiso:
            genetoiso[gene] = set()
        genetoiso[gene].add(iso)

        chromregion = (thischr, math.floor(start / GENE_BIN_SIZE) * GENE_BIN_SIZE)
        if chromregion not in chrregiontogenes:
            chrregiontogenes[chromregion] = set()
        chrregiontogenes[chromregion].add(gene)

        genestoboundaries = add_gene_to_boundaries(genestoboundaries, gene, thischr, dir, start, end)
    return isotoblocks, genetoiso, chrregiontogenes, genestoboundaries

# this is called by extract_vcf_vars
def add_vcf_var(vcfvars, chrom, ref, alts, tpos2, name):
    roundedpos = math.floor(tpos2 / (10 ** 6)) * 10 ** 6
    poskey = (chrom, roundedpos)
    if poskey not in vcfvars:
        vcfvars[poskey] = {}
    vcfvars[poskey][tpos2] = (ref, alts, name)
    return vcfvars

def get_correct_vcf_vars(vcfvars, refname, startpos, endpos):
    myvcfvars = {}
    for roundedpos in range(math.floor(startpos / (10 ** 6)) * 10 ** 6, math.ceil(endpos / (10 ** 6)) * 10 ** 6, 10**6):
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


def parse_single_bam_read(s, temp_files, vcfvars, sampleindex):
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
        coveredvarstrings = [VAR_FIELD_SEP.join([str(v), vcfvars[v][2], str(coveredvars[v])]) for v in coveredvars]
        tempvarout = temp_files.open(s.reference_name)
        tempvarout.write(
            '\t'.join([s.reference_name, str(sampleindex) + '__' + s.query_name, ';'.join(coveredvarstrings)]) + '\n')


def read_vars_to_genome_pos_counts(tempfilenames, tempdir, outprefix, sampledata, threshold, output_all):  # noqa: C901 - FIXME: reduce complexity
    samplenames = [x[0] for x in sampledata]

    with open(f'{outprefix}.var.counts.tsv', 'w') as out, open(f'{outprefix}.vargroup.counts.tsv', 'w') as out2:
        out.write('\t'.join(['varpos', 'gene', 'transcript'] + samplenames) + '\n')
        vartocounts = {}
        vargroup_to_data = {}
        for tf in tempfilenames:
            with open(tempdir + tf + '.txt', 'r') as mutstringfile:
                for line in mutstringfile:
                    refname, readname, mutstring = line.rstrip('\n').split('\t')
                    # three fields per mut now, whatever the gene count: the gene list
                    # used to be comma joined into a comma joined record
                    allmuts = [x.split(VAR_FIELD_SEP) for x in mutstring.split(';')]

                    allgenes = []
                    for x in allmuts:
                        allgenes.extend([g for g in x[1].split(GENE_LIST_SEP) if g != ''])

                    if len(allgenes) > 0:  # only use reads that overlap annotated genes
                        mygene = max(set(allgenes), key=allgenes.count)

                        varpos = [x[0] for x in allmuts]
                        varkey = ','.join(varpos)
                        key = (refname, mygene, varkey)
                        # print(allmuts)
                        modinfo = tuple([int(x[-1]) for x in allmuts])

                        if key not in vargroup_to_data:
                            vargroup_to_data[key] = []
                        vargroup_to_data[key].append(modinfo)

                    # each mut: position, varname (may be chrom:pos or gene), varstatus
                    sampleindex = int(readname.split('__')[0])
                    for m in allmuts:
                        varpos = refname + ':' + m[0]
                        # first of the sorted gene list, as before, when a variant falls
                        # in more than one gene
                        gene = m[1].split(GENE_LIST_SEP)[0]
                        transcript = ''
                        var = (varpos, gene, transcript)
                        if var not in vartocounts:
                            vartocounts[var] = [[0, 0] for x in range(len(samplenames))]  # [unmod counts, mod counts]
                        vartocounts[var][sampleindex][int(m[-1])] += 1

        for chrom, gene, varpos in vargroup_to_data:
            readinfo = vargroup_to_data[(chrom, gene, varpos)]
            outmods = []
            totpos = len(varpos.split(','))
            for readmods in readinfo:
                totmods = len([x for x in readmods if x == 1])
                outmods.append(str(totmods))
            outline = [chrom, gene, str(len(outmods)), str(totpos), ','.join(outmods), varpos]
            out2.write('\t'.join(outline) + '\n')

        for var in vartocounts:
            if any([x[0] + x[1] >= threshold for x in vartocounts[var]]) and (any([x[1] > 0 for x in vartocounts[var]]) or output_all):  # any modified reads in any sample
                varcounts = [f'{x[0]};{x[1]}' for x in vartocounts[var]]
                outline = list(var) + varcounts
                out.write('\t'.join(outline) + '\n')

def retrieve_good_iso_pos(potgenes, genestoboundaries, gpos, genetoiso, isotoblocks):
    for gene in potgenes:
        genechr, genedir, genestart, geneend = genestoboundaries[gene]
        if genestart <= gpos < geneend:
            for iso2 in genetoiso[gene]:
                blocks = isotoblocks[iso2]
                tpos2 = None
                for tstart, gstart, bsize, thischr, dir in blocks:
                    if gstart <= gpos < gstart + bsize:
                        if dir == '+':
                            tpos2 = tstart + (gpos - gstart)
                        else:
                            # mirror within the block, covering the same
                            # tstart .. tstart + bsize - 1 the plus branch does
                            tpos2 = tstart + (bsize - 1 - (gpos - gstart))
                        break
                if tpos2 is not None:
                    yield gene, iso2, tpos2

def group_annotated_ref_vars(vartoalt, chrregiontogenes, genestoboundaries, genetoiso, isotoblocks):
    vcfvars = {}
    for refinfo in vartoalt:
        chrom, gpos, ref, alts = extract_varinfo(refinfo, vartoalt[refinfo])

        potgenes = get_potential_genes(chrom, gpos, chrregiontogenes)

        overlapgenes = set()
        # THIS MAY BE THE BOTTLENECK
        for gene, _, _ in retrieve_good_iso_pos(potgenes, genestoboundaries, gpos, genetoiso, isotoblocks):
            overlapgenes.add(gene)
        # sorted: overlapgenes is a set, so the order, and therefore which gene the
        # counts are attributed to below, varied between runs
        vcfvars = add_vcf_var(vcfvars, chrom, ref, alts, gpos, GENE_LIST_SEP.join(sorted(overlapgenes)))
    return vcfvars

def combine_vcf_files(vcffilelist):
    vartoalt = {}
    for samplevcf in vcffilelist:
        for line in open(samplevcf):
            if line[0] != '#':
                line = line.rstrip().split('\t')
                refinfo, alt = tuple(line[:4]), line[4]
                if refinfo not in vartoalt:
                    vartoalt[refinfo] = set()
                vartoalt[refinfo].add(alt)
    return vartoalt

class TempVarFiles:
    """Append handles for the per-reference temp files, opened once each.  The write
    path used to open, append and close the file for every read covering a variant."""
    def __init__(self, tempdir, stack):
        self.tempdir = tempdir
        self.stack = stack
        self.by_refname = {}

    def open(self, refname):
        out = self.by_refname.get(refname)
        if out is None:
            out = self.stack.enter_context(open(self.tempdir + refname + '.txt', 'a'))
            self.by_refname[refname] = out
        return out


def _parse_one_bam_file(bamfile, temp_files, vcfvars, sindex):
    with pysam.AlignmentFile(bamfile, 'rb') as samfile:
        c = 0
        for s in samfile:
            if s.is_mapped:  # and not s.is_supplementary: ##not s.is_secondary and
                c += 1
                if c % 100000 == 0:
                    print(c, 'reads checked')
                myvcfvars = get_correct_vcf_vars(vcfvars, s.reference_name, s.reference_start, s.reference_end)
                parse_single_bam_read(s, temp_files, myvcfvars, sindex)


def parse_all_bam_files(sampledata, tempdir, vcfvars):
    with ExitStack() as stack:
        temp_files = TempVarFiles(tempdir, stack)
        for sindex in range(len(sampledata)):
            sample, bamfile = sampledata[sindex][0], sampledata[sindex][1]
            _parse_one_bam_file(bamfile, temp_files, vcfvars, sindex)
            print('done parsing reads for', sample)

def get_genes_from_tempdir(tempdir):
    genenames = set()
    for f in os.listdir(tempdir):
        if f[0] != '.' and 'processed' not in f:
            genenames.add(f.split('.txt')[0])
    return genenames

def quantvarpos():
    args = parse_var_args()
    # Load reference data
    if args.manifest:
        sampledata = extract_sample_data(args.manifest)
    elif args.input_bam and (args.pos_ref or args.vcf):
        if args.vcf:
            sampledata = [['sample', args.input_bam, args.vcf]]
        else:
            sampledata = [['sample', args.input_bam]]
    else:
        raise ValueError("please provide either manifest or bam and vcf")

    print('done loading annot')

    vcfvars = {}
    if args.manifest or args.vcf:
        isotoblocks, genetoiso, chrregiontogenes, genestoboundaries = get_bedisoform_info(args.bedisoforms)
        vartoalt = combine_vcf_files([x[2] for x in sampledata if len(x) > 2])
        print('done combining vcfs')
        vcfvars = group_annotated_ref_vars(vartoalt, chrregiontogenes, genestoboundaries, genetoiso, isotoblocks)
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
    tempdir = make_temp_dir(args.output_prefix)
    parse_all_bam_files(sampledata, tempdir, vcfvars)  # parses to intermediate files with read name to all vars
    print('parsed all reads')
    genenames = get_genes_from_tempdir(tempdir)
    read_vars_to_genome_pos_counts(genenames, tempdir, args.output_prefix, sampledata, args.threshold, args.output_all)

    if not args.keep_intermediate:
        shutil.rmtree(tempdir)


if __name__ == "__main__":
    quantvarpos()
