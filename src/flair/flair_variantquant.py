#! /usr/bin/env python3

import os
import shutil
import logging
import math
import multiprocessing as mp
from collections import Counter
import pysam
from flair.pycbio.hgdata.bed import BedReader
from flair.flair_bed import FlairBed
from flair.io_utils import make_run_temp_dir
from flair import FlairInputDataError
os.environ['OPENBLAS_NUM_THREADS'] = '1'

compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W',
            'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}

# the bases of the ref and alt of an SNV, the only variants counted
SNV_BASES = frozenset('ACGT')

# genes are binned by the start coordinate of the gene for the variant lookup
GENE_BIN_SIZE = 10 ** 7
# variants are binned by position, for the lookup of the variants a read covers
VAR_BIN_SIZE = 10 ** 6

# a read's allele at a site where it shows none of the site's alts
NO_ALT = '.'

# a temp record is pos, gene list, status; the gene list is itself a list, so the two
# need separators that cannot be confused, or the field count varies with gene count
VAR_FIELD_SEP = '|'
GENE_LIST_SEP = ','


def add_subparser(subparsers):
    desc = "Quantify variants at genome positions from reads aligned to the transcriptome"
    parser = subparsers.add_parser('variantquant', help="Quantify variants at genome positions",
                                   description=desc)
    parser.add_argument('--manifest', type=str,
                        help="used instead of --transcriptome_bam and --vcf: path to manifest file that points to sample names, "
                             "bam files aligned to transcriptome, "
                             "and vcf vars for that sample called on the genome. "
                             "Each line of file should be tab separated. "
                             "If you are using just one reference vcf file, "
                             "just include it in the third column for the first sample "
                             "and leave that column blank for the rest.")
    parser.add_argument('--transcriptome_bam',
                        help='used instead of --manifest: path to bam file for an individual sample')
    parser.add_argument('--pos_ref',
                        help='either this, --vcf or --manifest: path to reference file of sites to check for variants')
    parser.add_argument('--vcf',
                        help='either this, --pos_ref or --manifest: path to reference vcf file of sites to check for variants')
    parser.add_argument('-o', '--output', default='flair',
                        help="prefix for output files (default: %(default)s)")
    parser.add_argument('--isoform_bed',
                        help="path to transcriptome bed file")
    parser.add_argument('--min_coverage', type=int, default=5,
                        help='minimum total read coverage required to output a site (default: %(default)s)')
    parser.add_argument('--output_all', action='store_true',
                        help="output read counts for all putative RNA editing sites that pass the coverage threshold, "
                             "regardless of whether any reads are edited")
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='number of samples to read at once, each in its own process (default: %(default)s)')
    parser.add_argument('--temp_dir',
                        help='directory for temporary files; each run makes its own directory, named for the '
                             'output, in it, with a file per sample of its reads covering the sites '
                             '(default: $TMPDIR or the system temporary directory)')
    parser.add_argument('--keep_intermediate', action='store_true',
                        help="keep intermediate and temporary files for debugging, in the run's directory in "
                             '--temp_dir, which must be given')
    parser.set_defaults(entry=quantvarpos_cmd)

def quantvarpos_cmd(args):
    if args.keep_intermediate and args.temp_dir is None:
        raise FlairInputDataError('--keep_intermediate requires --temp_dir, the directory to keep them in')
    quantvarpos(manifest=args.manifest, transcriptome_bam=args.transcriptome_bam,
                pos_ref=args.pos_ref, vcf=args.vcf, output=args.output,
                isoform_bed=args.isoform_bed, min_coverage=args.min_coverage,
                output_all=args.output_all, keep_intermediate=args.keep_intermediate,
                temp_dir=args.temp_dir, threads=args.threads)

def extract_sample_data(manifestfile):
    sampledata = []
    for line in open(manifestfile):
        line = line.rstrip().split()
        sampledata.append(line)
    return sampledata

def extract_varinfo(refinfo, alt):
    chrom, gpos, ref = refinfo
    return chrom, int(gpos), ref, sorted(alt)

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

def var_bin(pos):
    "the start of the bin of variants pos is in"
    return math.floor(pos / VAR_BIN_SIZE) * VAR_BIN_SIZE


def add_vcf_var(vcfvars, chrom, ref, alts, tpos2, name):
    vcfvars.setdefault((chrom, var_bin(tpos2)), {})[tpos2] = (ref, alts, name)
    return vcfvars


def site_alleles(vcfvars, chrom, pos):
    "the ref and alts of the site at chrom:pos"
    ref, alts, _ = vcfvars[(chrom, var_bin(pos))][pos]
    return ref, alts


def get_correct_vcf_vars(vcfvars, refname, startpos, endpos):
    myvcfvars = {}
    for roundedpos in range(var_bin(startpos), math.ceil(endpos / VAR_BIN_SIZE) * VAR_BIN_SIZE, VAR_BIN_SIZE):
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
                    coveredvars[pos] = NO_ALT
            ref += block[1]
            quer += block[1]
        elif block[0] in {1, 4}:  # consumes query ###1 is insertion
            quer += block[1]
        elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
            ref += block[1]
    return coveredvars


def parse_single_bam_read(s, out, vcfvars, sampleindex):
    """Write the variant sites a read covers to out, its sample's temp file, each with
    its genes and the alt the read shows there, or NO_ALT.  vcfvars are the sites
    near the read (get_correct_vcf_vars); most reads cover none, and finding the
    read's bases at its sites is by far the slowest part, so it is done only for
    the reads covering one"""
    if not vcfvars:
        return
    # check for which vars are covered
    coveredvars = _parse_cigar(s.cigartuples, s.reference_start, vcfvars)
    if not coveredvars:
        return

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
                        coveredvars[i[1] + 1] = alt

    coveredvarstrings = [VAR_FIELD_SEP.join([str(v), vcfvars[v][2], coveredvars[v]]) for v in coveredvars]
    out.write('\t'.join([s.reference_name, str(sampleindex) + '__' + s.query_name, ';'.join(coveredvarstrings)]) + '\n')


def read_vars_to_genome_pos_counts(temp_paths, outprefix, sampledata, threshold, output_all, vcfvars):
    """Write the counts of each allele of each site, a row of OUTPUT.var.counts.tsv:
    per sample, the reads covering the site not showing the alt, and those showing
    it, so they add up to the site's coverage whatever its other alts.  A row is
    written if a sample has threshold reads covering the site, and if a sample has a
    read showing the alt, or output_all.  temp_paths are the samples' temp files,
    from parse_all_bam_files"""
    samplenames = [x[0] for x in sampledata]

    with open(f'{outprefix}.var.counts.tsv', 'w') as out:
        out.write('\t'.join(['varpos', 'ref', 'alt', 'gene', 'transcript'] + samplenames) + '\n')
        vartocounts = {}
        for temp_path in temp_paths:
            with open(temp_path) as mutstringfile:
                for line in mutstringfile:
                    refname, readname, mutstring = line.rstrip('\n').split('\t')
                    # three fields per mut now, whatever the gene count: the gene list
                    # used to be comma joined into a comma joined record
                    allmuts = [x.split(VAR_FIELD_SEP) for x in mutstring.split(';')]

                    # each mut: position, varname (may be chrom:pos or gene), varstatus
                    sampleindex = int(readname.split('__')[0])
                    for m in allmuts:
                        # first of the sorted gene list, as before, when a variant falls
                        # in more than one gene
                        var = (refname, int(m[0]), m[1].split(GENE_LIST_SEP)[0])
                        if var not in vartocounts:
                            # per sample, the reads covering the site, by the alt they show
                            vartocounts[var] = [Counter() for x in range(len(samplenames))]
                        vartocounts[var][sampleindex][m[-1]] += 1

        transcript = ''
        for (refname, pos, gene), sample_alleles in vartocounts.items():
            coverage = [sum(alleles.values()) for alleles in sample_alleles]
            if not any(cov >= threshold for cov in coverage):
                continue
            ref, alts = site_alleles(vcfvars, refname, pos)
            for alt in alts:
                modified = [alleles[alt] for alleles in sample_alleles]
                if any(mod > 0 for mod in modified) or output_all:
                    varcounts = [f'{cov - mod};{mod}' for cov, mod in zip(coverage, modified)]
                    out.write('\t'.join([f'{refname}:{pos}', ref, alt, gene, transcript] + varcounts) + '\n')

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


def is_snv(ref, alt):
    "whether a vcf ref and one of its alts are an SNV, which one read base can show"
    return ref in SNV_BASES and alt in SNV_BASES and alt != ref


def combine_vcf_files(vcffilelist):
    """(chrom, pos, ref) -> the alt bases of the SNVs of the vcf files, of each site's
    alts in all of them, whatever their IDs.  Each alt of a multi-allelic site is its
    own SNV, a read showing any of them modified.  Other variants, indels and
    substitutions of more than one base, are left out with a warning.  An SNV site
    given different refs is an error: the vcfs aren't of one genome"""
    vartoalt = {}
    site_refs = {}  # (chrom, pos) -> (ref, the vcf file first giving it)
    nskipped = 0
    for samplevcf in vcffilelist:
        for line in open(samplevcf):
            if line[0] != '#':
                chrom, pos, _, ref, alts = line.rstrip('\n').split('\t')[:5]
                ref = ref.upper()
                alts = {alt.upper() for alt in alts.split(',')}
                snv_alts = {alt for alt in alts if is_snv(ref, alt)}
                nskipped += len(alts) - len(snv_alts)
                if snv_alts:
                    # of SNVs only: an indel at the site has a longer ref, as TA of a deletion of A
                    first_ref, first_vcf = site_refs.setdefault((chrom, pos), (ref, samplevcf))
                    if ref != first_ref:
                        raise FlairInputDataError(f'vcf files disagree on the reference base of {chrom}:{pos}: {first_ref} '
                                                  f'in {first_vcf}, {ref} in {samplevcf}; are they called on the same genome?')
                    vartoalt.setdefault((chrom, pos, ref), set()).update(snv_alts)
    if nskipped > 0:
        logging.warning(f'{nskipped} vcf alleles that are not SNVs (indels or substitutions of more than '
                        'one base) left out: variantquant counts SNVs only')
    return vartoalt

def _parse_one_bam_file(bamfile, out, vcfvars, sindex):
    with pysam.AlignmentFile(bamfile, 'rb') as samfile:
        c = 0
        for s in samfile:
            if s.is_mapped:  # and not s.is_supplementary: ##not s.is_secondary and
                c += 1
                if c % 100000 == 0:
                    logging.info(f'{c} reads checked')
                myvcfvars = get_correct_vcf_vars(vcfvars, s.reference_name, s.reference_start, s.reference_end)
                parse_single_bam_read(s, out, myvcfvars, sindex)


# the sites of parse_all_bam_files, for its forked processes, which share them
# rather than each task being sent a copy of what may be millions of sites
_pool_vcfvars = None


def _parse_sample(sindex, sample, bamfile, temp_path, vcfvars=None):
    "parse a sample's reads to its temp file, of the sites, or the pool's"
    with open(temp_path, 'w') as out:
        _parse_one_bam_file(bamfile, out, _pool_vcfvars if vcfvars is None else vcfvars, sindex)
    logging.info(f'done parsing reads for {sample}')


def parse_all_bam_files(sampledata, tempdir, vcfvars, threads):
    """Parse each sample's reads covering the sites to a temp file of its own,
    tempdir/<sample index>.txt, threads samples at once.  The files' paths, by sample"""
    global _pool_vcfvars
    tasks = [(sindex, row[0], row[1], f'{tempdir}{sindex}.txt') for sindex, row in enumerate(sampledata)]
    if threads <= 1 or len(tasks) <= 1:
        for task in tasks:
            _parse_sample(*task, vcfvars=vcfvars)
    else:
        _pool_vcfvars = vcfvars
        try:
            mp.set_start_method('fork', force=True)
            with mp.Pool(min(threads, len(tasks))) as pool:
                pool.starmap(_parse_sample, tasks)
        finally:
            _pool_vcfvars = None
    return [task[-1] for task in tasks]


def quantvarpos(*, manifest, transcriptome_bam, pos_ref, vcf, output, isoform_bed,
                min_coverage, output_all, keep_intermediate, temp_dir=None, threads=1):
    # Load reference data
    if manifest:
        sampledata = extract_sample_data(manifest)
    elif transcriptome_bam and (pos_ref or vcf):
        if vcf:
            sampledata = [['sample', transcriptome_bam, vcf]]
        else:
            sampledata = [['sample', transcriptome_bam]]
    else:
        raise FlairInputDataError("specify --manifest, or --transcriptome_bam with one of --vcf or --pos_ref")

    logging.info('done loading annot')

    vcfvars = {}
    if manifest or vcf:
        isotoblocks, genetoiso, chrregiontogenes, genestoboundaries = get_bedisoform_info(isoform_bed)
        vartoalt = combine_vcf_files([x[2] for x in sampledata if len(x) > 2])
        logging.info('done combining vcfs')
        vcfvars = group_annotated_ref_vars(vartoalt, chrregiontogenes, genestoboundaries, genetoiso, isotoblocks)
    else:
        nskipped = 0
        for line in open(pos_ref):
            line = line.rstrip('\n').split('\t')
            # the file's bin of the position, its second column, is that of add_vcf_var
            chrom, _, pos, ref, alt, name = line
            alts = [a for a in alt.split(',') if is_snv(ref, a)]
            nskipped += len(alt.split(',')) - len(alts)
            if alts:
                add_vcf_var(vcfvars, chrom, ref, alts, int(pos), name)
        if nskipped > 0:
            logging.warning(f'{nskipped} --pos_ref alleles that are not SNVs left out: variantquant counts SNVs only')

    logging.info('done combining vcf variants')
    # a directory of this run's own, so no earlier run's files are read
    tempdir = make_run_temp_dir(output, temp_dir)
    logging.info(f'temporary files in {tempdir}')
    try:
        temp_paths = parse_all_bam_files(sampledata, tempdir, vcfvars, threads)
        logging.info('parsed all reads')
        read_vars_to_genome_pos_counts(temp_paths, output, sampledata, min_coverage, output_all, vcfvars)
    finally:
        # also on failure, so runs don't leave their files in $TMPDIR
        if keep_intermediate:
            logging.info(f'intermediate files kept in {tempdir}')
        else:
            shutil.rmtree(tempdir, ignore_errors=True)
