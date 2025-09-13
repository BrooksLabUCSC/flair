#! /usr/bin/env python3

import argparse
import os
import math
import pipettor
import pysam
import shutil
from flair import FlairInputDataError
os.environ['OPENBLAS_NUM_THREADS'] = '1'


compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def getStarts(gtf):
    starts = list()
    with open(gtf) as lines:
        for l in lines:
            if l[0] == "#":
                continue
            cols = l.rstrip().split("\t")
            chrom, c1, c2, strand = cols[0], int(cols[3]) - 1, int(cols[4]), cols[6]
            if cols[2] == "start_codon":
                gene = cols[8][cols[8].find('gene_id') + len('gene_id') + 2:]
                gene = gene[:gene.find('"')]

                starts.append((chrom, c1, c2, gene, ".", strand))
    if (len(starts)) == 0:
        raise FlairInputDataError(f'ERROR, no start codons were found in {gtf}')
    return starts

class Isoform(object):
    '''
    Object to handle isoform related data.
    '''

    def __init__(self, name=None):
        self.name = name
        self.chrom = ""
        self.seqvariants = {}
        self.exons = set()
        self.starts = set()
        self.exonSizes = list()
        self.ptcpoint = 0
        self.startcodons = []

class SeqVar(object):
    def __init__(self, name=None, seq=None):
        self.name = name
        self.pro = "UNK"
        self.sequence = seq
        self.orfs = list()

def getStartRelPos(genomicStartPos, exon, exons, isoObj):
    '''
    is handed a genomic position, the exon it occurs in, all exons,
    and returns the position relative to all exons
    '''
    exonNum = exons.index(exon)
    isoObj.exonSizes = [x[1] - x[0] for x in exons]
    relativeStart = None
    # First get start position relative to transcript sequence.
    if isoObj.strand == "+":
        relativeStart = genomicStartPos - exons[exonNum][0] + sum([x for x in isoObj.exonSizes[:exonNum]])
    elif isoObj.strand == "-":
        relativeStart = len(isoObj.sequence) - (genomicStartPos - exons[exonNum][0] +
                                                sum([x for x in isoObj.exonSizes[:exonNum]])) - 3

    return relativeStart

def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--manifest', required=True, type=str,
                        help="path to manifest files that points to sample names, bam files aligned to transcriptome, "
                             "and vcf vars for that sample called on the genome. "
                             "Each line of file should be tab separated. "
                             "If you are using just one reference vcf file, "
                             "just include it in the third column for the first sample "
                             "and leave that column blank for the rest.")
    parser.add_argument('-o', '--output_prefix', default='flair',
                        help="path to collapsed_output.bed file. default: 'flair'")
    parser.add_argument('-i', '--isoforms',
                        help="path to transcriptome fasta file")
    parser.add_argument('-b', '--bedisoforms',
                        help="path to transcriptome bed file")
    parser.add_argument('-g', '--genome',
                        type=str, required=True, help='FastA of reference genome')
    parser.add_argument('-f', '--gtf',
                        type=str, required=True,
                        help='GTF annotation file')
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

def process_bedline(line):
    line = line.rstrip().split('\t')
    thischr, iso, dir, start, esizes, estarts, end = line[0], line[3], line[5], int(line[1]), \
        [int(x) for x in line[10].split(',')[:-1]], \
        [int(x) for x in line[11].split(',')[:-1]], int(line[2])
    exonblocks = []  # block is gstart, tstart, len
    if dir == '-':
        esizes = esizes[::-1]
        estarts = estarts[::-1]
    currtstart = 0
    for i in range(len(esizes)):
        exonblocks.append((currtstart, estarts[i] + start, esizes[i], thischr, dir))
        currtstart += esizes[i]
    exonblocks.sort()

    return thischr, iso, dir, start, esizes, estarts, end, exonblocks, currtstart


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
    for line in open(bedisofile):
        thischr, iso, dir, start, esizes, estarts, end, exonblocks, currtstart = process_bedline(line)

        if iso[:10] == 'fusiongene':
            iso = '_'.join(iso.split('_')[1:])

        isotoblocks = add_iso_to_blocks(isotoblocks, iso, exonblocks)

        gene = iso.split('_')[-1]
        if gene not in genetoiso:
            genetoiso[gene] = set()
        genetoiso[gene].add(iso)

        chromregion = (thischr, math.floor(start / (10 ** 7)) * 10 ** 7)
        if chromregion not in chrregiontogenes:
            chrregiontogenes[chromregion] = set()
        chrregiontogenes[chromregion].add(gene)

        genestoboundaries = add_gene_to_boundaries(genestoboundaries, gene, thischr, dir, start, end)
    return isotoblocks, genetoiso, chrregiontogenes, genestoboundaries

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


def add_vcf_var(vcfvars, iso2, ref, alt, tpos2, name=None):  # chrom and gpos required for other script
    if iso2 not in vcfvars:
        vcfvars[iso2] = set()
    for a in alt:
        if len(a) == len(ref):
            thistype = 'S'
        elif len(a) < len(ref):
            thistype = 'D'
        else:
            thistype = 'I'
        vcfvars[iso2].add((tpos2, thistype))
    return vcfvars


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
                            tpos2 = (tstart + bsize + 1) - (gpos - gstart)
                        break
                if tpos2:
                    yield gene, iso2, tpos2

def extract_varinfo(refinfo, alt):
    chrom, gpos, id, ref = refinfo
    alt = list(alt)
    gpos = int(gpos)
    return chrom, gpos, ref, alt

def get_potential_genes(chrom, gpos, chrregiontogenes):
    chromregion1, chromregion2 = (chrom, math.floor(gpos / (10 ** 7)) * 10 ** 7), \
                                 (chrom, (math.floor(gpos / (10 ** 7)) + 1) * 10 ** 7)
    potgenes = set()
    if chromregion1 in chrregiontogenes:
        potgenes.update(chrregiontogenes[chromregion1])
    if chromregion2 in chrregiontogenes:
        potgenes.update(chrregiontogenes[chromregion2])
    return potgenes

def convert_vars_to_tpos(vartoalt, isotoblocks, genetoiso, chrregiontogenes, genestoboundaries):
    # load vcf variants - could skip this and instead load variants dict directly when combining vcfs above
    # for longshot variant calling, can force it to output ‘genotype’ snvs, but it refuses to output alt allele
    # current workaround is to make the longshot vcf report positions only, but that doesn’t seem good long-term
    vcfvars = {}
    for refinfo in vartoalt:
        chrom, gpos, ref, alt = extract_varinfo(refinfo, vartoalt[refinfo])
        # Handles - chrom region - what if start and end of gene are in different rounding blocks??
        potgenes = get_potential_genes(chrom, gpos, chrregiontogenes)

        # First need to get what genes overlap with this genomic pos
        # Then go through transcripts in that gene and save the transcript pos of the variant
        # Don't forget to correct for strand for alt!!
        for gene, iso2, tpos2 in retrieve_good_iso_pos(potgenes, genestoboundaries, gpos, genetoiso, isotoblocks):
            vcfvars = add_vcf_var(vcfvars, iso2, ref, alt, tpos2, chrom + ':' + str(gpos))

    return vcfvars

def parse_isoform_fa(isoformfile):
    isotoseq = {}
    isoname = None
    for line in open(isoformfile):
        if line[0] == '>':
            isoname = line[1:].rstrip()
        else:
            isotoseq[isoname] = line.rstrip()
    return isotoseq

def make_temp_dir(outprefix):
    dir = outprefix + '_tempvarfiles/'
    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.makedirs(dir)
    return dir

def parse_cigar(cigar, alignstart, readseq):
    thesemutations = []
    ref, quer = 0, 0
    for block in cigar:
        if block[0] in {0, 7, 8}:  # match, consumes both
            ref += block[1]
            quer += block[1]
        elif block[0] in {1, 4}:  # consumes query, 1 is insertion
            if block[0] == 1:
                # insertion has the inserted base sequence
                thesemutations.append((ref + alignstart + 1, 'I', readseq[quer:quer + block[1]]))
            quer += block[1]
        elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
            if block[0] == 2:
                thesemutations.append((ref + alignstart + 1, 'D', block[1]))  # deletion has a number
            ref += block[1]
    return thesemutations


def parse_single_bam_read(s, tempdir, vcfvars, sampleindex, tempfilename):
    readseq = s.query_sequence
    thesemutations = parse_cigar(s.cigartuples, s.reference_start, readseq)
    alignedbases = s.get_aligned_pairs(with_seq=True, matches_only=True)
    for i in alignedbases:
        if i[2] and i[2].islower():
            thesemutations.append((i[1] + 1, 'S', readseq[i[0]]))
    filteredmuts = []
    for m in thesemutations:
        if m[:2] in vcfvars or (m[1] == 'D' and m[2] >= 10) or (m[1] == 'I' and len(m[2]) >= 10):
            filteredmuts.append(m)
    filteredmuts = ['.'.join([str(z) for z in x]) for x in filteredmuts]
    with open(tempdir + tempfilename + '.txt', 'a') as tempvarout:
        tempvarout.write('\t'.join([s.reference_name, str(sampleindex) + '__' + s.query_name, ','.join(filteredmuts)]) + '\n')

def gettempfilename(refname, mode):
    if mode and mode == 'genomic':
        return refname
    else:
        return refname.split('_')[-1].split('.')[0][:-2]  # partial gene name

def get_correct_vcf_vars(vcfvars, refname, startpos, endpos):
    if refname in vcfvars:
        return vcfvars[refname]
    else:
        return None


def parse_all_bam_files(sampledata, tempdir, vcfvars, mode=None):  # mode is only necessary for varquant mode
    for sindex in range(len(sampledata)):
        sample, bamfile = sampledata[sindex][0], sampledata[sindex][1]
        samfile = pysam.AlignmentFile(bamfile, 'rb')
        c = 0
        for s in samfile:
            # trusting FLAIR quantify to only output one alignment per read
            if s.is_mapped:  # and not s.is_supplementary: ##not s.is_secondary and
                c += 1
                if c % 100000 == 0:
                    print(c, 'reads checked')
                tempfilename = gettempfilename(s.reference_name, mode)
                myvcfvars = get_correct_vcf_vars(vcfvars, s.reference_name, s.reference_start, s.reference_end)
                if myvcfvars:
                    parse_single_bam_read(s, tempdir, myvcfvars, sindex, tempfilename)
        samfile.close()
        print('done parsing reads for', sample)

def get_genes_from_tempdir(tempdir):
    genenames = set()
    for f in os.listdir(tempdir):
        if f[0] != '.' and 'processed' not in f:
            genenames.add(f.split('.txt')[0])
    return genenames

def get_mutstrings_for_iso_group(isofile):
    isotomutstringtoreads = {}
    with open(isofile, 'r') as mutstringfile:
        for line in mutstringfile:
            isoname, readname, mutstring = line.rstrip('\n').split('\t')
            if mutstring == '':
                thesemutations = tuple()
            else:
                thesemuts = [x.split('.') for x in mutstring.split(',')]
                thesemutations = []
                for m in thesemuts:
                    m[0] = int(m[0])
                    thesemutations.append(tuple(m))
                thesemutations = tuple(thesemutations)
            if isoname not in isotomutstringtoreads:
                isotomutstringtoreads[isoname] = {}
            if thesemutations not in isotomutstringtoreads[isoname]:
                isotomutstringtoreads[isoname][thesemutations] = set()
            isotomutstringtoreads[isoname][thesemutations].add(readname)
    return isotomutstringtoreads


def write_mutsets_for_iso(iso, blocks, mutsetstoreads, outfile, genomepostosupport):
    for mutset in mutsetstoreads:
        rcount = len(mutsetstoreads[mutset])
        newmutpos = []
        for m in mutset:  # pos, type, str
            tpos = m[0]
            genomeset = None
            for tstart, gstart, bsize, thischr, dir in blocks:
                if (dir == '+' and tstart <= tpos < tstart + bsize) or (
                        dir == '-' and tstart <= tpos - 1 < tstart + bsize):
                    if dir == '+':
                        genomePos = gstart + (tpos - tstart)
                    else:
                        genomePos = gstart + ((tstart + bsize) - (tpos - 1))
                    genomeset = (thischr, dir, str(genomePos), str(tpos), m[1], m[2])
                    break
            if genomeset:
                if genomeset not in genomepostosupport:
                    genomepostosupport[genomeset] = rcount
                else:
                    genomepostosupport[genomeset] += rcount
                newmutpos.append('..'.join(genomeset))
        outfile.write(
            '\t'.join([iso, ','.join(newmutpos), ','.join(mutsetstoreads[mutset])]) + '\n')
    return genomepostosupport


def read_vars_to_var_groups(genenames, tempdir, isotoblocks):
    genomepostosupport = {}
    for gene in genenames:
        isotomutstringtoreads = get_mutstrings_for_iso_group(tempdir + gene + '.txt')
        # go through mutstrings for each transcript
        # convert all vars to genomic positions
        # get support for each var at the genomic level
        with open(tempdir + gene + '_processed.txt', 'w') as mutstringfile2:
            for iso in isotomutstringtoreads:
                write_mutsets_for_iso(iso, isotoblocks[iso], isotomutstringtoreads[iso], mutstringfile2, genomepostosupport)
    return genomepostosupport


def extract_readmuts_from_temp(tempdir, gene):
    isotomutstringtoreads = {}
    with open(tempdir + gene + '_processed.txt', 'r') as mutstringfile3:
        for line in mutstringfile3:
            iso, genomemutpos, readnames = line.rstrip('\n').split('\t')
            if iso not in isotomutstringtoreads:
                isotomutstringtoreads[iso] = {}
            if genomemutpos == '':
                genomemutpos = frozenset()
            else:
                genomemutpos = frozenset((tuple(x.split('..')) for x in genomemutpos.split(',')))
            isotomutstringtoreads[iso][genomemutpos] = set(readnames.split(','))
    return isotomutstringtoreads

def get_goodsupport_muts(mutsetstoreads, genomepostosupport):
    allmuts = set()
    for mutset in mutsetstoreads:
        allmuts.update(mutset)
    goodmuts = set()
    for m in allmuts:
        if genomepostosupport[m] >= 3:
            goodmuts.add(m)

    goodsupport = {}
    for mutset in mutsetstoreads:
        filtset = mutset & goodmuts
        if filtset not in goodsupport:
            goodsupport[filtset] = mutsetstoreads[mutset]
        else:
            goodsupport[filtset].update(mutsetstoreads[mutset])
    return goodsupport

def convert_muts_to_genome_pos(mutpos, thisseq):
    mutinfo = sorted([x[:3] + (int(x[3]) - 1,) + x[4:] for x in mutpos], reverse=True, key=lambda x: x[3])
    realmutpos = []
    for thischr, strand, genomepos, pos, muttype, alt in mutinfo:
        if muttype == 'S':
            thisseq[pos] = alt
        elif muttype == 'D':
            count = int(alt)
            thisseq = thisseq[:pos] + thisseq[pos + count:]
            realmutpos = [x[:3] + (x[3] - count,) + x[4:] for x in realmutpos]
        elif muttype == 'I':
            count = len(alt)
            thisseq = thisseq[:pos] + list(alt) + thisseq[pos:]
            realmutpos = [x[:3] + (x[3] + count,) + x[4:] for x in realmutpos]
        realmutpos.append((thischr, strand, genomepos, pos, muttype, alt))
    realmutpos.sort()
    return realmutpos


def get_good_support_varisos(genenames, tempdir, isotoseq, numsamples, genomepostosupport):
    isovarcounts, isovarfa = {}, {}
    for gene in genenames:
        isotomutstringtoreads = extract_readmuts_from_temp(tempdir, gene)
        for iso in isotomutstringtoreads:
            # for each isoform, get the count of total reads supporting each individual variant
            # also get total reads for isoform
            # for each individual mutation, figure out which ones pass threshold of support
            goodsupport = get_goodsupport_muts(isotomutstringtoreads[iso], genomepostosupport)

            for mutpos in goodsupport:
                thisseq = list(isotoseq[iso])

                realmutpos = convert_muts_to_genome_pos(mutpos, thisseq)

                if len(realmutpos) == 0:
                    muttotext = 'nomuts'
                else:
                    muttotext = ','.join(['..'.join([str(y) for y in x]) for x in realmutpos])
                outname = iso + '__' + muttotext

                if iso not in isovarcounts:
                    isovarcounts[iso] = {}
                if outname not in isovarcounts[iso]:
                    isovarcounts[iso][outname] = [0 for s in range(numsamples)]
                    isovarfa[outname] = ''.join(thisseq)

                for r in goodsupport[mutpos]:
                    s = int(r.split('__')[0])
                    isovarcounts[iso][outname][s] += 1
    return isovarcounts, isovarfa


def write_genomevars(outprefix, genomepostosupport):
    genomevarsout = open(outprefix + '.isovars.genomicpos.bed', 'w')
    for gset in genomepostosupport:
        if genomepostosupport[gset] >= 3:
            genomevarsout.write('\t'.join([gset[0], str(int(gset[2]) - 1), gset[2],
                                           str(genomepostosupport[gset]) + '-' + '..'.join(gset[3:]),
                                           '0', gset[1]]) + '\n')
    genomevarsout.close()

def write_varisos(outprefix, samplenames, isovarcounts, isovarfa):
    isomutsout = open(outprefix + '.isoswithvars.counts.tsv', 'w')
    isomutsfa = open(outprefix + '.isoswithvars.fa', 'w')
    isomutsout.write('\t'.join(['isoname', 'varsontranscript', 'varsongenome'] + samplenames) + '\n')

    for iso in sorted(list(isovarcounts.keys())):
        varperisocount = 1
        for isovar in sorted(list(isovarcounts[iso].keys())):
            i, m = isovar.split('__')
            counts = isovarcounts[iso][isovar]
            outseq = isovarfa[isovar]
            outname = str(varperisocount) + '-' + iso
            isomutsout.write('\t'.join([outname, m] + [str(x) for x in counts]) + '\n')
            isomutsfa.write('>' + outname + ' ' + m + '\n')
            isomutsfa.write(outseq + '\n')
            varperisocount += 1

    isomutsout.close()
    isomutsfa.close()

def align_variso_models(outprefix, genome):
    # align to transcriptome with --secondary=no
    mm2_cmd = ('minimap2', '-ax', 'splice', '-s', '40', '-t', '12', '--secondary=no', genome, outprefix + '.isoswithvars.fa')
    samtools_sort_cmd = ('samtools', 'sort', '-o', outprefix + '.isoswithvars.bam', '-')
    samtools_index_cmd = ('samtools', 'index', outprefix + '.isoswithvars.bam')
    pipettor.run([mm2_cmd, samtools_sort_cmd])
    pipettor.run([samtools_index_cmd])


def get_iso_prod_info(filename):
    isotoprodaaseq, aaseqtoinfo = {}, {}
    aaseqcount = 1
    for line in open(filename):
        line = line.rstrip().split('\t')
        isoname, prod = line[0], line[1]
        if len(line) == 2:
            aaseq = ''
        else:
            aaseq = line[-1]
        gname = isoname.split('_')[-1]
        if prod != 'PRO':
            aaseq = 'NOTPRO'
        else:
            if (gname, aaseq) not in aaseqtoinfo:
                aaseqtoinfo[(gname, aaseq)] = [aaseqcount, []]
                aaseqcount += 1
            aaseqtoinfo[(gname, aaseq)][1].append(isoname)
            aaseq = 'aaseq' + str(aaseqtoinfo[(gname, aaseq)][0])
        isotoprodaaseq[isoname] = (gname, prod, aaseq)
    return isotoprodaaseq, aaseqtoinfo

def get_iso_var_counts(filename, isotoprodaaseq):
    genetoaaseqtocounts, genetoprodtocounts = {}, {}
    countssamples = None
    for line in open(filename):
        line = line.rstrip().split('\t')
        if line[0] != 'isoname':
            isoname = line[0]  # '__'.join(line[0].split('__')[:-1])
            counts = [int(x) for x in line[2:]]
            if isoname in isotoprodaaseq:
                gname, prod, aaseq = isotoprodaaseq[isoname]
                if gname not in genetoprodtocounts:
                    genetoprodtocounts[gname] = {}
                    genetoaaseqtocounts[gname] = {}

                if prod not in genetoprodtocounts[gname]:
                    genetoprodtocounts[gname][prod] = counts
                else:
                    genetoprodtocounts[gname][prod] = [genetoprodtocounts[gname][prod][x]
                                                       + counts[x] for x in range(len(counts))]
                if aaseq not in genetoaaseqtocounts[gname]:
                    genetoaaseqtocounts[gname][aaseq] = counts
                else:
                    genetoaaseqtocounts[gname][aaseq] = [genetoaaseqtocounts[gname][aaseq][x]
                                                         + counts[x] for x in range(len(counts))]
        else:
            countssamples = line[3:]
    return countssamples, genetoaaseqtocounts, genetoprodtocounts


def write_aaseq_counts_tables(outprefix):
    # need to get counts table for productivity and aaseq
    isotoprodaaseq, aaseqtoinfo = get_iso_prod_info(outprefix + '.isoswithvars.productivity.info.tsv')
    countssamples, genetoaaseqtocounts, genetoprodtocounts = get_iso_var_counts(outprefix + '.isoswithvars.counts.tsv',
                                                                                isotoprodaaseq)

    out = open(outprefix + '.aaseq.key.tsv', 'w')
    for gname, aaseq in aaseqtoinfo:
        aaseqid, tlist = aaseqtoinfo[(gname, aaseq)]
        out.write('\t'.join([gname, 'aaseq' + str(aaseqid), ','.join(tlist), aaseq]) + '\n')
    out.close()

    out = open(outprefix + '.aaseq.counts.tsv', 'w')
    out.write('\t'.join(['aaseqid_gene'] + countssamples) + '\n')
    for gname in genetoaaseqtocounts:
        for id in genetoaaseqtocounts[gname]:
            counts = genetoaaseqtocounts[gname][id]
            out.write('\t'.join([id + '_' + gname] + [str(x) for x in counts]) + '\n')
    out.close()

    prodtypes = ['PRO', 'PTC', 'NST', 'NGO']
    out = open(outprefix + '.productivity.counts.tsv', 'w')
    out.write('\t'.join(['productivity_gene'] + countssamples) + '\n')
    for gname in genetoprodtocounts:
        for p in prodtypes:
            if p in genetoprodtocounts[gname]:
                counts = genetoprodtocounts[gname][p]
                out.write('\t'.join([p + '_' + gname] + [str(x) for x in counts]) + '\n')
    out.close()


def getvariants():
    args = parse_args()

    # Load reference data
    sampledata = extract_sample_data(args.manifest)
    isotoblocks, genetoiso, chrregiontogenes, genestoboundaries = get_bedisoform_info(args.bedisoforms)
    print('done parsing iso blocks')
    vartoalt = combine_vcf_files([x[2] for x in sampledata if len(x) > 2])  # combines vcf files, does not output combined file
    vcfvars = convert_vars_to_tpos(vartoalt, isotoblocks, genetoiso, chrregiontogenes, genestoboundaries)
    print('done combining vcf variants')
    isotoseq = parse_isoform_fa(args.isoforms)
    tempdir = make_temp_dir(args.output_prefix)

    parse_all_bam_files(sampledata, tempdir, vcfvars)  # parses to intermediate files with read name to all vars
    genenames = get_genes_from_tempdir(tempdir)
    genomepostosupport = read_vars_to_var_groups(genenames, tempdir, isotoblocks)
    print('done preprocessing iso muts')
    isovarcounts, isovarfa = get_good_support_varisos(genenames, tempdir, isotoseq, len(sampledata), genomepostosupport)
    write_genomevars(args.output_prefix, genomepostosupport)
    write_varisos(args.output_prefix, [x[0] for x in sampledata], isovarcounts, isovarfa)

    # align_variso_models(args.output_prefix, args.genome) ###TAKES TOO LONG, ADD OPTION TO TOGGLE

    # get productivity for transcripts without variants
    prodcmd = ('predictProductivity.py',
               '-i', args.bedisoforms,
               '-o', args.output_prefix + '.isoforms.productivity',
               '--gtf', args.gtf,
               '--genome_fasta', args.genome,
               '--longestORF')
    pipettor.run([prodcmd])

    # adjust productivity prediction to account for variants
    prodcmd = ('predict_aaseq_withvar.py',
               args.output_prefix + '.isoforms.productivity.info.tsv',
               args.output_prefix + '.isoswithvars.fa',
               args.output_prefix + '.isoswithvars.productivity.info.tsv')
    pipettor.run([prodcmd])

    write_aaseq_counts_tables(args.output_prefix)

    if not args.keep_intermediate:
        shutil.rmtree(tempdir)


if __name__ == "__main__":
    getvariants()
