"""flair variantquant: the reads of each sample showing each SNV allele at sites of a
vcf or a --pos_ref file, as RNA editing sites of a database like REDIportal, from
reads aligned to the genome.  Each sample's reads covering a site go to a temp file
of its own, in parallel, a record per read: its reference, sample index and read
name, and per site covered its position, genes and the alt it shows, or NO_ALT;
then the counts of every sample are written from them"""

import os
import shutil
import heapq
import logging
import math
import multiprocessing as mp
from bisect import bisect_right
from collections import Counter
import pysam
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.sys import fileOps
from flair.interval_index import IntervalIndex
from flair.io_utils import make_run_temp_dir
from flair import FlairInputDataError

# cigar operations aligning read bases to the genome (M, =, X), only of the read
# (I, S), and only of the genome (D, N)
ALIGNED_OPS = (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF)
QUERY_OPS = (pysam.CINS, pysam.CSOFT_CLIP)
REF_OPS = (pysam.CDEL, pysam.CREF_SKIP)

# the bases of the ref and alt of an SNV, the only variants counted
SNV_BASES = frozenset('ACGT')

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
                        help='bed file of isoforms whose exons name the sites of a vcf by their genes: a site is '
                             'of the genes of the exons over it.  An isoform\'s gene is its 13th column, as of a '
                             'flair isoform bed, or, for a BED12 or an isoform with it empty, that of the --gtf '
                             'gene its exons overlap most, else it names none.  Without it or --gtf the gene '
                             'column is empty; with --pos_ref, whose sites have their genes, neither is used')
    parser.add_argument('-f', '--gtf',
                        help='GTF annotation file, plain or gzipped, whose exons name the sites of a vcf by their '
                             'gene_ids, or, with --isoform_bed, the genes of its isoforms without them')
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
                isoform_bed=args.isoform_bed, gtf=args.gtf, min_coverage=args.min_coverage,
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


def var_bin(pos):
    "the start of the bin of variants pos is in"
    return math.floor(pos / VAR_BIN_SIZE) * VAR_BIN_SIZE


def add_vcf_var(vcfvars, chrom, ref, alts, pos, name):
    """add the site at 1-based pos to vcfvars, (chrom, bin) -> pos -> (ref, alts, name),
    name its comma-separated genes"""
    vcfvars.setdefault((chrom, var_bin(pos)), {})[pos] = (ref, alts, name)
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


def read_vcf_sites(vcffilelist, keep_flags=frozenset()):
    """The SNV sites of the vcf files, plain or gzipped, as add_vcf_var keeps them, of
    no genes yet (name_sites_by_genes): each site's alts of all the vcfs, whatever its
    IDs, of the records passing their filters or flagged only by keep_flags
    (passes_filters).  Each alt of a multi-allelic site is its own SNV, a read showing
    any of them modified.  Other variants, indels and substitutions of more than one
    base, are left out with a warning.  An SNV site given different refs is an error:
    the vcfs aren't of one genome.  Built as the files are read, its refs and alts
    shared, as millions of sites, of a database such as REDIportal, may be read"""
    vcfvars = {}
    shared = {}  # each ref and alts once
    nskipped = nfiltered = 0
    for vcf_index, samplevcf in enumerate(vcffilelist):
        with fileOps.opengz(samplevcf) as fh:
            for line in fh:
                if line[0] == '#':
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) > 6 and not passes_filters(fields[6], keep_flags):
                    nfiltered += 1
                    continue
                chrom, pos, _, ref, alts = fields[:5]
                ref = ref.upper()
                alts = set(alts.upper().split(','))
                snv_alts = {alt for alt in alts if is_snv(ref, alt)}
                nskipped += len(alts) - len(snv_alts)
                if not snv_alts:
                    continue
                pos = int(pos)
                sites = vcfvars.setdefault((chrom, var_bin(pos)), {})
                site = sites.get(pos)
                if site is not None:
                    # of SNVs only: an indel at the site has a longer ref, as TA of a deletion of A
                    if site[0] != ref:
                        first_vcf = _vcf_giving_ref(vcffilelist[:vcf_index + 1], chrom, pos, site[0])
                        raise FlairInputDataError(f'vcf files disagree on the reference base of {chrom}:{pos}: {site[0]} '
                                                  f'in {first_vcf}, {ref} in {samplevcf}; are they called on the same genome?')
                    snv_alts |= set(site[1])
                alts = tuple(sorted(snv_alts))
                sites[pos] = (shared.setdefault(ref, ref), shared.setdefault(alts, alts), '')
    if nfiltered > 0:
        logging.info(f'{nfiltered} vcf records not passing their filters left out; --keep_flags counts them')
    if nskipped > 0:
        logging.warning(f'{nskipped} vcf alleles that are not SNVs (indels or substitutions of more than '
                        'one base) left out: variantquant counts SNVs only')
    return vcfvars


def _vcf_giving_ref(vcffilelist, chrom, pos, ref):
    "the first of the vcf files giving the SNV site chrom:pos ref, for the error of read_vcf_sites"
    for samplevcf in vcffilelist:
        with fileOps.opengz(samplevcf) as fh:
            for line in fh:
                fields = line.split('\t', 5)
                if line[0] != '#' and fields[0] == chrom and int(fields[1]) == pos and fields[3].upper() == ref:
                    return samplevcf
    return vcffilelist[-1]


def _sorted_exons(exons):
    "chrom -> (start, end, gene) of each exon, sorted, of sets of them by chrom"
    return {chrom: sorted(chrom_exons) for chrom, chrom_exons in exons.items()}


# the column, 0-based, of a bed line's gene, as flair isoform beds have it, after the
# 12 of BED12; a bed without it, or with it empty, has no genes
BED_GENE_COLUMN = 12


def isoforms_of_bed(isoform_bed):
    """The isoforms of a bed file, (chrom, strand, exons, gene) each, exons 0-based,
    gene that of its 13th column, BED_GENE_COLUMN, as a flair isoform bed gives it, or
    None for a bed of 12 columns or an isoform with it empty"""
    isoforms = []
    for bed in BedReader(isoform_bed):
        extra = bed.extraCols
        gene = extra[BED_GENE_COLUMN - bed.numStdCols] if bed.numStdCols == 12 and len(extra) > 0 else ''
        exons = [(block.start, block.end) for block in bed.blocks] if bed.blocks else [(bed.chromStart, bed.chromEnd)]
        isoforms.append((bed.chrom, bed.strand, exons, gene or None))
    return isoforms


def _read_gtf_exons(gtf):
    """chrom -> the set of exons, (start, end, gene_id) 0-based, of a GTF file, plain or
    gzipped, and gene_id -> strand: only its exon lines and their gene_ids are read,
    several times faster than parsing its transcripts"""
    exons, gene_strand = {}, {}
    with fileOps.opengz(gtf) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            fields = line.split('\t', 9)
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            attrs = fields[8]
            start = attrs.find('gene_id "')
            if start < 0:
                raise FlairInputDataError(f'{gtf}: exon without a gene_id: {line.rstrip()}')
            start += len('gene_id "')
            gene = attrs[start:attrs.index('"', start)]
            exons.setdefault(fields[0], set()).add((int(fields[3]) - 1, int(fields[4]), gene))
            gene_strand[gene] = fields[6]
    return exons, gene_strand


def exons_of_gtf(gtf):
    "chrom -> the exons, (start, end, gene_id) 0-based and sorted, of a GTF file, each exon of a gene once"
    return _sorted_exons(_read_gtf_exons(gtf)[0])


def gtf_gene_of_isoform(chrom, strand, exons, gtf_index, gene_strand):
    """the GTF gene of an isoform without one: of the most bases of the isoform's exons
    in its exons, of the isoform's strand, unless it has none, of the gtf_index
    (chrom -> IntervalIndex of the GTF exons), or None if none overlaps it"""
    chrom_index = gtf_index.get(chrom)
    if chrom_index is None:
        return None
    overlap = Counter()
    for start, end in exons:
        for exon_start, exon_end, gene in chrom_index.overlap(start, end):
            if strand not in ('+', '-') or gene_strand[gene] == strand:
                overlap[gene] += min(end, exon_end) - max(start, exon_start)
    if not overlap:
        return None
    # the lowest gene id of those tied, so a run's names don't depend on set order
    return min(overlap, key=lambda gene: (-overlap[gene], gene))


def site_exons(isoform_bed, gtf):
    """chrom -> the exons, (start, end, gene) 0-based and sorted, to name the sites by
    (name_sites_by_genes): of the isoform bed, the genes of its isoforms its own, as a
    flair bed's, or else those of the GTF they overlap most (gtf_gene_of_isoform);
    else of the GTF alone.  Isoforms without a gene, and no GTF, name no sites"""
    if isoform_bed is None:
        return exons_of_gtf(gtf)
    isoforms = isoforms_of_bed(isoform_bed)
    exons = {}
    without_gene = [isoform for isoform in isoforms if isoform[3] is None]
    if without_gene and gtf is None:
        logging.warning(f'{len(without_gene)} of the {len(isoforms)} isoforms of {isoform_bed} have no gene, '
                        'its 13th column, and name no sites; give --gtf to name them by the genes they overlap')
    elif without_gene:
        gtf_exons, gene_strand = _read_gtf_exons(gtf)
        gtf_index = {}
        for chrom, chrom_exons in gtf_exons.items():
            gtf_index[chrom] = IntervalIndex()
            for exon_start, exon_end, gene in chrom_exons:
                gtf_index[chrom].add(exon_start, exon_end, (exon_start, exon_end, gene))
            gtf_index[chrom].build()
        named = [(chrom, strand, iso_exons, gtf_gene_of_isoform(chrom, strand, iso_exons, gtf_index, gene_strand))
                 for chrom, strand, iso_exons, _ in without_gene]
        logging.info(f'{sum(gene is not None for *_, gene in named)} of {len(without_gene)} isoforms of '
                     f'{isoform_bed} without a gene given that of {gtf}')
        isoforms = [isoform for isoform in isoforms if isoform[3] is not None] + named
    elif gtf is not None:
        logging.info(f'the isoforms of {isoform_bed} have their genes, so {gtf} is not used')
    for chrom, _, iso_exons, gene in isoforms:
        if gene is not None:
            exons.setdefault(chrom, set()).update((start, end, gene) for start, end in iso_exons)
    return _sorted_exons(exons)


def name_sites_by_genes(vcfvars, exons):
    """Name each site of vcfvars (add_vcf_var) by the genes, comma-separated and sorted,
    of the exons over it, of exons (site_exons).  One sweep
    of each chromosome's sorted sites and exons, holding the exons begun and not yet
    ended, rather than a search for each of millions of sites"""
    bins_of_chrom = {}
    for chrom, bin_start in vcfvars:
        bins_of_chrom.setdefault(chrom, []).append(bin_start)
    for chrom, bin_starts in bins_of_chrom.items():
        chrom_exons = exons.get(chrom, [])
        next_exon, active = 0, []  # the next exon, and a heap of (end, gene) of those begun
        for bin_start in sorted(bin_starts):
            sites = vcfvars[(chrom, bin_start)]
            for pos in sorted(sites):
                # an exon, 0-based [start, end), is over 1-based pos when start < pos <= end
                while next_exon < len(chrom_exons) and chrom_exons[next_exon][0] < pos:
                    heapq.heappush(active, chrom_exons[next_exon][1:])
                    next_exon += 1
                while active and active[0][0] < pos:
                    heapq.heappop(active)
                if active:
                    ref, alts, _ = sites[pos]
                    sites[pos] = (ref, alts, GENE_LIST_SEP.join(sorted({gene for _, gene in active})))
    return vcfvars


def load_pos_ref(pos_ref):
    "the SNV sites of a --pos_ref file, as add_vcf_var keeps them, with their genes, refs and alts shared"
    vcfvars = {}
    shared = {}
    nskipped = 0
    for line in open(pos_ref):
        # the file's bin of the position, its second column, is that of add_vcf_var
        chrom, _, pos, ref, alt, name = line.rstrip('\n').split('\t')
        alts = tuple(a for a in alt.split(',') if is_snv(ref, a))
        nskipped += len(alt.split(',')) - len(alts)
        if alts:
            add_vcf_var(vcfvars, chrom, shared.setdefault(ref, ref), shared.setdefault(alts, alts), int(pos),
                        shared.setdefault(name, name))
    if nskipped > 0:
        logging.warning(f'{nskipped} --pos_ref alleles that are not SNVs left out: variantquant counts SNVs only')
    return vcfvars

def _parse_one_bam_file(bamfile, out, vcfvars, positions, sindex):
    """parse the reads of a bam file: its primary and supplementary alignments, the
    supplementary of other parts of a read, as of a chimera, but not its secondary,
    other alignments of the same bases, which would count them again"""
    with pysam.AlignmentFile(bamfile, 'rb') as samfile:
        c = 0
        for s in samfile:
            if s.is_mapped and not s.is_secondary:
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
                isoform_bed=None, gtf=None, min_coverage, output_all, keep_intermediate, temp_dir=None, threads=1):
    sampledata, site_vcfs = samples_and_site_vcfs(manifest, bam, sample_name, pos_ref, vcf)
    if site_vcfs is None and keep_flags:
        logging.warning('--keep_flags is not used with --pos_ref, which has no vcf filters')
    if site_vcfs is None:
        if isoform_bed is not None or gtf is not None:
            logging.warning('--isoform_bed and --gtf are not used with --pos_ref, whose sites have their genes')
        vcfvars = load_pos_ref(pos_ref)
    else:
        vcfvars = read_vcf_sites(site_vcfs, keep_flags)
        logging.info('done reading vcfs')
        if isoform_bed is not None or gtf is not None:
            name_sites_by_genes(vcfvars, site_exons(isoform_bed, gtf))
            logging.info('done finding the genes of the sites')
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
