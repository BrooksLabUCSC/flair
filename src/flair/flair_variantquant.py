#! /usr/bin/env python3

import os
import shutil
import logging
import math
import multiprocessing as mp
from bisect import bisect_right
from collections import Counter
import pysam
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.sys import fileOps
from flair.flair_bed import FlairBed
from flair.io_utils import make_run_temp_dir
from flair import FlairInputDataError
os.environ['OPENBLAS_NUM_THREADS'] = '1'

# cigar operations aligning read bases to the genome (M, =, X), only of the read
# (I, S), and only of the genome (D, N)
ALIGNED_OPS = (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF)
QUERY_OPS = (pysam.CINS, pysam.CSOFT_CLIP)
REF_OPS = (pysam.CDEL, pysam.CREF_SKIP)

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
    desc = ('count the reads showing each SNV allele at sites given by a vcf or a --pos_ref file, as the '
            'sites of RNA editing in a database such as REDIportal, or variants called in the samples, per '
            'sample, from reads aligned to the genome')
    parser = subparsers.add_parser('variantquant', help="Count the reads showing variants at genome positions",
                                   description=desc)
    reads = parser.add_mutually_exclusive_group(required=True)
    reads.add_argument('--bam',
                       help='bam file of one sample, aligned to the genome; or give --manifest')
    reads.add_argument('--manifest',
                       help='tab-separated manifest of the samples, a line each: sample name and bam file, '
                            'aligned to the genome, with --vcf or --pos_ref for the sites, or sample name, bam '
                            "file and vcf file of the sample's own variants, on every line, whose sites are "
                            'merged, every sample counted at all of them')
    sites = parser.add_mutually_exclusive_group()
    sites.add_argument('--vcf',
                       help='vcf file, plain or gzipped, of the sites to count, for all samples.  Only SNVs of '
                            'records whose FILTER is PASS or . are counted, or one of --keep_flags')
    sites.add_argument('--pos_ref',
                       help='file of the sites to count, for all samples, as annotate_group_vcf_vars writes, '
                            'faster to load than a vcf with --isoform_bed: tab-separated chromosome, the start of '
                            'the position\'s 1 Mb bin, 1-based position, ref base, comma-separated alt bases, '
                            'comma-separated genes')
    parser.add_argument('--keep_flags', type=lambda flags: set(flags.split(',')), default=set(), metavar='FLAGS',
                        help='comma-separated vcf FILTER flags whose records are counted too, as RefCall, or all, '
                             'every record; a record of more than one flag is counted if all of them are given')
    parser.add_argument('--sample_name',
                        help="the sample's name, with --bam (default: the bam file's name, less .bam)")
    parser.add_argument('-o', '--output', default='flair',
                        help="prefix for output files (default: %(default)s)")
    parser.add_argument('--isoform_bed',
                        help='flair isoform bed file, as flair transcriptome writes, of the genes to name the '
                             'sites of a vcf with: a site is of the genes with an isoform exon over it.  Without '
                             'it, or with --pos_ref, whose sites have their genes, the gene column is empty')
    parser.add_argument('--min_coverage', type=int, default=5,
                        help='minimum reads covering a site, in at least one sample, for its alleles to be '
                             'written (default: %(default)s)')
    parser.add_argument('--output_all', action='store_true',
                        help='write the alleles of every site with --min_coverage, whether or not any read '
                             'shows them')
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
    if args.sample_name is not None and args.bam is None:
        raise FlairInputDataError("--sample_name names --bam's sample; a manifest names its samples")
    quantvarpos(manifest=args.manifest, bam=args.bam, sample_name=args.sample_name,
                pos_ref=args.pos_ref, vcf=args.vcf, keep_flags=args.keep_flags, output=args.output,
                isoform_bed=args.isoform_bed, min_coverage=args.min_coverage,
                output_all=args.output_all, keep_intermediate=args.keep_intermediate,
                temp_dir=args.temp_dir, threads=args.threads)


def _check_exists(path, what, where):
    if not os.path.exists(path):
        raise FlairInputDataError(f'{where}: {what} {path} does not exist')


def read_manifest(path):
    """A manifest's samples, [name, bam] or [name, bam, vcf] each, a tab-separated line
    each but comments and blank lines.  Either every sample has its vcf or none does:
    a vcf on some lines, as one for all samples on the first, is for --vcf"""
    samples = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip('\r\n')
            if line.strip() == '' or line.startswith('#'):
                continue
            fields = line.split('\t')
            where = f'{path}:{lineno}'
            if len(fields) not in (2, 3) or '' in fields:
                raise FlairInputDataError(f'{where}: a manifest line is a sample name, its bam file and optionally '
                                          f'its own vcf file, tab-separated; for one vcf of all samples give --vcf: '
                                          f'{line!r}')
            _check_exists(fields[1], 'bam file', where)
            if len(fields) == 3:
                _check_exists(fields[2], 'vcf file', where)
            samples.append(fields)
    if len(samples) == 0:
        raise FlairInputDataError(f'{path}: manifest has no samples')
    if len({len(sample) for sample in samples}) > 1:
        raise FlairInputDataError(f'{path}: either every sample has its own vcf in the manifest or none does; '
                                  'for one vcf of all samples, leave them out and give --vcf')
    names = [sample[0] for sample in samples]
    duplicated = sorted({name for name in names if names.count(name) > 1})
    if duplicated:
        raise FlairInputDataError(f'{path}: samples named more than once: {", ".join(duplicated)}')
    return samples


def bam_sample_name(bam):
    "a sample's default name, its bam file's name less .bam"
    name = os.path.basename(bam)
    return name[:-len('.bam')] if name.endswith('.bam') else name


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


def site_positions(vcfvars):
    "the sorted positions of each bin's sites, by bin, to find the sites of each read (sites_in_span)"
    return {key: sorted(sites) for key, sites in vcfvars.items()}


def sites_in_span(positions, chrom, start, end):
    """the sorted positions, 1-based, of the sites in chrom:start-end, 0-based, of
    positions (site_positions), found by bisection, so a read's lookup costs little
    with millions of sites, as of a database of RNA editing sites"""
    sites = []
    for bin_start in range(var_bin(start + 1), var_bin(end) + 1, VAR_BIN_SIZE):
        bin_sites = positions.get((chrom, bin_start))
        if bin_sites:
            sites.extend(bin_sites[bisect_right(bin_sites, start):bisect_right(bin_sites, end)])
    return sites


def read_bases_at(read, sites):
    """(site, base) of each of sites, sorted 1-based positions in the read's span, that
    the read aligns a base to, with that base: one walk of its cigar, rather than
    decoding all its aligned bases.  Sites in its deletions and introns have none"""
    seq = read.query_sequence
    ref, query, i = read.reference_start, 0, 0
    for op, length in read.cigartuples:
        if i == len(sites):
            return
        if op in ALIGNED_OPS:
            end = ref + length
            # the sites up to the block's end; those before it were passed already
            while i < len(sites) and sites[i] <= end:
                yield sites[i], seq[query + sites[i] - 1 - ref]
                i += 1
            ref, query = end, query + length
        elif op in QUERY_OPS:
            query += length
        elif op in REF_OPS:
            ref += length
            while i < len(sites) and sites[i] <= ref:
                i += 1


def parse_single_bam_read(s, out, vcfvars, positions, sampleindex):
    """Write the sites a read aligns a base to, of vcfvars, to out, its sample's temp
    file, each with its genes and the alt the read shows there, or NO_ALT.  positions
    are vcfvars' site_positions.  Most reads cover no site, and return at once"""
    sites = sites_in_span(positions, s.reference_name, s.reference_start, s.reference_end)
    if not sites or s.query_sequence is None:
        return
    records = []
    for pos, base in read_bases_at(s, sites):
        _, alts, name = vcfvars[(s.reference_name, var_bin(pos))][pos]
        records.append(VAR_FIELD_SEP.join([str(pos), name, base if base in alts else NO_ALT]))
    if records:
        out.write('\t'.join([s.reference_name, str(sampleindex) + '__' + s.query_name, ';'.join(records)]) + '\n')


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


# the FILTER of the vcf records whose sites are counted; a vcf of just its first
# five columns, as a database of sites may be, has none, and all are counted
PASSING_FILTERS = ('PASS', '.')
# the --keep_flags of every record
KEEP_ALL_FLAGS = 'all'


def passes_filters(vcf_filter, keep_flags):
    """whether a vcf record of FILTER vcf_filter is counted: it passes, or all its flags,
    separated by ;, are of keep_flags, or keep_flags has KEEP_ALL_FLAGS"""
    return (vcf_filter in PASSING_FILTERS or KEEP_ALL_FLAGS in keep_flags
            or all(flag in keep_flags for flag in vcf_filter.split(';')))


def combine_vcf_files(vcffilelist, keep_flags=frozenset()):
    """(chrom, pos, ref) -> the alt bases of the SNVs of the vcf files, plain or
    gzipped, of each site's alts in all of them, whatever their IDs, of the records
    passing their filters or flagged only by keep_flags (passes_filters).  Each alt of a multi-allelic site is its
    own SNV, a read showing any of them modified.  Other variants, indels and
    substitutions of more than one base, are left out with a warning.  An SNV site
    given different refs is an error: the vcfs aren't of one genome"""
    vartoalt = {}
    site_refs = {}  # (chrom, pos) -> (ref, the vcf file first giving it)
    nskipped = nfiltered = 0
    for samplevcf in vcffilelist:
        with fileOps.opengz(samplevcf) as fh:
            for line in fh:
                if line[0] != '#':
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) > 6 and not passes_filters(fields[6], keep_flags):
                        nfiltered += 1
                        continue
                    chrom, pos, _, ref, alts = fields[:5]
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
    if nfiltered > 0:
        logging.info(f'{nfiltered} vcf records not passing their filters left out; --keep_flags counts them')
    if nskipped > 0:
        logging.warning(f'{nskipped} vcf alleles that are not SNVs (indels or substitutions of more than '
                        'one base) left out: variantquant counts SNVs only')
    return vartoalt


def unannotated_vars(vartoalt):
    "the sites of vartoalt (combine_vcf_files) as group_annotated_ref_vars gives them, of no genes"
    vcfvars = {}
    for refinfo in vartoalt:
        chrom, gpos, ref, alts = extract_varinfo(refinfo, vartoalt[refinfo])
        add_vcf_var(vcfvars, chrom, ref, alts, gpos, '')
    return vcfvars


def load_pos_ref(pos_ref):
    "the sites of a --pos_ref file, as group_annotated_ref_vars gives them, of its SNVs"
    vcfvars = {}
    nskipped = 0
    for line in open(pos_ref):
        # the file's bin of the position, its second column, is that of add_vcf_var
        chrom, _, pos, ref, alt, name = line.rstrip('\n').split('\t')
        alts = [a for a in alt.split(',') if is_snv(ref, a)]
        nskipped += len(alt.split(',')) - len(alts)
        if alts:
            add_vcf_var(vcfvars, chrom, ref, alts, int(pos), name)
    if nskipped > 0:
        logging.warning(f'{nskipped} --pos_ref alleles that are not SNVs left out: variantquant counts SNVs only')
    return vcfvars

def _parse_one_bam_file(bamfile, out, vcfvars, positions, sindex):
    with pysam.AlignmentFile(bamfile, 'rb') as samfile:
        c = 0
        for s in samfile:
            if s.is_mapped:  # and not s.is_supplementary: ##not s.is_secondary and
                c += 1
                if c % 100000 == 0:
                    logging.info(f'{c} reads checked')
                parse_single_bam_read(s, out, vcfvars, positions, sindex)


# the sites of parse_all_bam_files and their site_positions, for its forked
# processes, which share them rather than each task being sent a copy of what may
# be millions of sites
_pool_sites = None


def _parse_sample(sindex, sample, bamfile, temp_path, sites=None):
    "parse a sample's reads to its temp file, of sites, (vcfvars, positions), or the pool's"
    with open(temp_path, 'w') as out:
        _parse_one_bam_file(bamfile, out, *(_pool_sites if sites is None else sites), sindex)
    logging.info(f'done parsing reads for {sample}')


def parse_all_bam_files(sampledata, tempdir, vcfvars, threads):
    """Parse each sample's reads covering the sites to a temp file of its own,
    tempdir/<sample index>.txt, threads samples at once.  The files' paths, by sample"""
    global _pool_sites
    sites = (vcfvars, site_positions(vcfvars))
    tasks = [(sindex, row[0], row[1], f'{tempdir}{sindex}.txt') for sindex, row in enumerate(sampledata)]
    if threads <= 1 or len(tasks) <= 1:
        for task in tasks:
            _parse_sample(*task, sites=sites)
    else:
        _pool_sites = sites
        try:
            mp.set_start_method('fork', force=True)
            with mp.Pool(min(threads, len(tasks))) as pool:
                pool.starmap(_parse_sample, tasks)
        finally:
            _pool_sites = None
    return [task[-1] for task in tasks]


def samples_and_site_vcfs(manifest, bam, sample_name, pos_ref, vcf):
    """The samples, [name, bam] each, and the vcf files of the sites, or none for
    --pos_ref: --vcf, or each sample's own of the manifest"""
    if (manifest is None) == (bam is None):
        raise FlairInputDataError('give one of --bam, a sample, or --manifest')
    if pos_ref is not None and vcf is not None:
        raise FlairInputDataError('give one of --vcf or --pos_ref, the sites to count')
    if manifest is not None:
        samples = read_manifest(manifest)
        if len(samples[0]) == 3:
            if vcf is not None or pos_ref is not None:
                raise FlairInputDataError(f'{manifest} gives each sample its vcf, so --vcf and --pos_ref can not be given')
            return [sample[:2] for sample in samples], [sample[2] for sample in samples]
    else:
        _check_exists(bam, 'bam file', '--bam')
        samples = [[sample_name or bam_sample_name(bam), bam]]
    if vcf is None and pos_ref is None:
        raise FlairInputDataError('give --vcf or --pos_ref, the sites to count, or a manifest with the vcf of each sample')
    return samples, ([vcf] if vcf is not None else None)


def quantvarpos(*, manifest=None, bam=None, sample_name=None, pos_ref=None, vcf=None, keep_flags=frozenset(), output,
                isoform_bed=None, min_coverage, output_all, keep_intermediate, temp_dir=None, threads=1):
    sampledata, site_vcfs = samples_and_site_vcfs(manifest, bam, sample_name, pos_ref, vcf)
    if site_vcfs is None and keep_flags:
        logging.warning('--keep_flags is not used with --pos_ref, which has no vcf filters')
    if site_vcfs is None:
        if isoform_bed is not None:
            logging.warning('--isoform_bed is not used with --pos_ref, whose sites have their genes')
        vcfvars = load_pos_ref(pos_ref)
    else:
        vartoalt = combine_vcf_files(site_vcfs, keep_flags)
        logging.info('done combining vcfs')
        if isoform_bed is None:
            vcfvars = unannotated_vars(vartoalt)
        else:
            isotoblocks, genetoiso, chrregiontogenes, genestoboundaries = get_bedisoform_info(isoform_bed)
            vcfvars = group_annotated_ref_vars(vartoalt, chrregiontogenes, genestoboundaries, genetoiso, isotoblocks)
    logging.info(f'{sum(len(sites) for sites in vcfvars.values())} sites to count')

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
