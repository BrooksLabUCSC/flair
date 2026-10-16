#! /usr/bin/env python3

import logging
import pysam
from flair import FlairError
from flair.bed_to_gtf import bed_to_gtf
from flair.pycbio.hgdata.bed import BedReader
from flair.flair_bed import FlairBed
from statistics import median
from flair.isoform_data import make_big_bed, get_sequence_for_exons, Gene, Isoform, Junc
from flair.predictProductivity import predict_productivity
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import AnnotsByChrom
from flair.gene_identification import get_gene_name_firstpass, _group_novel_by_splicing

def add_subparser(subparsers):
    desc = "Combine FLAIR transcriptomes or annotation transcriptomes into one"
    parser = subparsers.add_parser('combine', help="Combine transcriptomes from multiple samples",
                                   description=desc)
    mutexc = parser.add_mutually_exclusive_group(required=True)
    mutexc.add_argument('--manifest', type=str,
                        help="path to manifest file that has a list of flair bed files paths, eg path/to/isoforms.bed")
    mutexc.add_argument('--prefixes', type=str,
                        help="comma separated list of file prefixes to combine, assuming command is being run "
                        "in folder that contains prefix.isoform.bed")
    mutexc.add_argument('--isoform_beds', type=str,
                        help="comma separated list of paths to flair bed files to combine")
    parser.add_argument('-g', '--genome', required=True,
                        help='genome fasta, required to generate an isoform fasta file and/or a bigbed file')
    parser.add_argument('-f', '--gtf', default=None,
                        help='an annotation GTF.  With it, each combined isoform is given its gene, reference transcript, '
                        'CDS, productivity and protein from this annotation, as flair transcriptome gives them, whatever '
                        'annotation the samples were run with, so samples run with different annotations or FLAIR versions '
                        'can be combined, and the proteins are written to combined.aaseq.tsv.  Fusion isoforms keep the '
                        'samples\' genes.  Without it, the samples\' genes, CDS and productivity are used, isoforms from '
                        'different samples with different genes are an error, and isoforms have no protein ids')
    parser.add_argument('-o', '--output', default='flair.combined.isoforms',
                        help="prefix for output file (default: %(default)s)")
    parser.add_argument('--end_window', type=int, default=200,
                        help="window, in bases, for comparing ends of isoforms with the same intron chain (default: %(default)s)")
    parser.add_argument('--min_percent_usage', type=int, default=5,
                        help="minimum percent usage required in one sample to keep isoform in combined transcriptome (default: %(default)s)")
    parser.add_argument('--remove_single_exon', action='store_true',
                        help='remove all single exon isoforms')
    parser.add_argument('--max_ends', type=int, default=1,
                        help='maximum number of TSS/TES picked per isoform; make higher for more precise end detection. '
                             'If max_ends is 1, all isoforms/read support for a splice junction chain will be condensed into one isoform. '
                             'max_ends > 1 introduces greater stringency for assigning og to new isoforms, more read support may be dropped '
                             '(default: %(default)s)')
    parser.add_argument('--min_reads', type=int, default=3,
                        help='min reads from all samples to call isoform (default: %(default)s)')
    parser.set_defaults(entry=combine_cmd)

def combine_cmd(args):
    # the option is a percentage, the code compares against fractional usage
    combine(manifest=args.manifest, prefixes=args.prefixes, isoform_beds=args.isoform_beds,
            genome=args.genome, annotation=args.gtf, output=args.output, end_window=args.end_window,
            min_frac_usage=args.min_percent_usage / 100.0,
            remove_single_exon=args.remove_single_exon, max_ends=args.max_ends,
            min_reads=args.min_reads)

# def bedReadToIntronChain(bed):
#     introns = []
#     for i in range(len(bed.blocks) - 1):
#         introns.append((bed.blocks[i].end, bed.blocks[i + 1].start))
#     # strand is not accounted for here, all intron chains will be left to right
#     return tuple(introns)

# def intronChainToestarts(ichain, start, end):
#     esizes, estarts = [], [0,]
#     for i in ichain:
#         esizes.append(i[0] - (start + estarts[-1]))
#         estarts.append(i[1] - start)
#     esizes.append(end - (start + estarts[-1]))
#     return esizes, estarts

def _longest_first(bed_list):
    """a sort key, reversed, for a sample's isoform: longest, then one without an
    allele's variants (flair isoalleles), whose CDS and productivity are the
    reference's, then most reads, then by sample"""
    bed = bed_list[0]
    return (bed.end - bed.start, bed.transcript_class != 'with_vars', bed.read_support, bed.samples, bed.start)


def combineIsos(beds, endwindow):
    beds.sort(key=lambda x: ([(y.start, y.end) for y in x], x[0].samples))
    isoendgroups = []
    last_vals = [(0, 0) for x in beds[0]]  # accounting for fusions, getting start + end of each fragment
    currgroup = []
    for bed_list in beds:
        my_edges = [(x.start, x.end) for x in bed_list]
        my_diffs = [(last_vals[x][0] - my_edges[x][0], last_vals[x][1] - my_edges[x][1]) for x in range(len(last_vals))]
        if not all((abs(x[0]) <= endwindow and abs(x[1]) <= endwindow for x in my_diffs)):
            if len(currgroup) > 0:
                isoendgroups.append(currgroup)
            currgroup = []
        currgroup.append(bed_list)
        last_vals = my_edges
    if len(currgroup) > 0:
        isoendgroups.append(currgroup)
    return isoendgroups


def parse_input(manifest, prefixes, input_beds):
    # assumes that only one input is not None
    sampledata = []
    if manifest is not None:
        for line in open(manifest):
            sampledata.append(line.rstrip())
    elif prefixes is not None:
        for prefix in prefixes.rstrip(',').split(','):
            sampledata.append(prefix + '.isoforms.bed')
    elif input_beds is not None:
        for bed_path in input_beds.rstrip(',').split(','):
            sampledata.append(bed_path)
    return sampledata


def _has_allele_info(bed):
    "is this an isoform of an allele, from flair isoalleles"
    return bed.transcript_class == 'with_vars' or bed.allele_group is not None or bed.source_isoform is not None


def load_isos_by_id(bedfile):
    isotoinfo = {}
    mysamples = set()
    num_alleles = 0
    for bed in BedReader(bedfile, bedClass=FlairBed):
        num_alleles += _has_allele_info(bed)
        # by strand too: a single-exon isoform and one antisense to it, with
        # similar ends, aren't the same isoform
        juncchaininfo = (bed.pos_in_fusion, bed.chrom, bed.strand, bed.getGaps())
        if bed.name not in isotoinfo:
            isotoinfo[bed.name] = []
        isotoinfo[bed.name].append((juncchaininfo, bed))
        mysamples.update(set(bed.samples))
    if num_alleles > 0:
        logging.warning(f"{bedfile}: {num_alleles} isoforms of alleles, from flair isoalleles.  flair combine drops "
                        f"the allele information, combining the alleles of an isoform into one isoform, with the CDS and "
                        f"productivity of an allele without variants where there is one (give -f to recompute them from "
                        f"the reference).  Their support is per allele, so each is a smaller fraction of its gene's reads "
                        f"than its isoform: set --min_percent_usage lower for them")
    return isotoinfo, mysamples

def organize_isos_by_junc_chain(isotoinfo, intronchaintoisos):
    for iso_id, bed_records in isotoinfo.items():
        bed_records.sort()  # organize fusion partners in order if necessary
        juncchaininfo = tuple([x[0] for x in bed_records])
        beds = tuple([x[1] for x in bed_records])
        if juncchaininfo not in intronchaintoisos:
            intronchaintoisos[juncchaininfo] = []
        intronchaintoisos[juncchaininfo].append(beds)

def load_isoforms_by_junc_chain(sampledata):
    intronchaintoisos = {}
    allsamples = set()
    for bed_path in sampledata:
        isotoinfo, mysamples = load_isos_by_id(bed_path)
        allsamples.update(mysamples)
        organize_isos_by_junc_chain(isotoinfo, intronchaintoisos)
    return intronchaintoisos, sorted(list(allsamples))

def filter_groups(ends_to_iso_groups, max_ends, minpercentusage, min_reads):
    if max_ends == 1:
        groups_after_end_filtering = [ends_to_iso_groups[0], ]
        for group in ends_to_iso_groups[1:]:
            groups_after_end_filtering[0].extend(group)  # adding to main group
    else:
        groups_after_end_filtering = [x for x in ends_to_iso_groups if x[0][0].frac_support > minpercentusage and x[0][0].read_support >= min_reads]
        if len(groups_after_end_filtering) == 0:
            groups_after_end_filtering = [ends_to_iso_groups[0], ]
    return groups_after_end_filtering

def _is_fusion(bed_list):
    return len(bed_list) > 1 or bed_list[0].pos_in_fusion is not None


def _combined_sort_key(bed_list_group):
    "a combined isoform's place in the output: by the coordinates of its isoform, the group's first"
    return tuple((bed.chrom, bed.start, bed.end, bed.strand, tuple((blk.start, blk.end) for blk in bed.blocks))
                 for bed in bed_list_group[0])


def _genes_str(ref_gene_mappings):
    "genes for a message, a fusion partner without a gene as 'none'"
    return ','.join(g or 'none' for g in ref_gene_mappings) or 'novel'


def _ref_genes(bed_list_group):
    """the annotated genes the samples gave a combined isoform, () for none, an
    error if they gave different ones"""
    ref_gene_ids = {tuple(x[0].ref_gene_mappings) for x in bed_list_group if len(x[0].ref_gene_mappings) > 0}
    if len(ref_gene_ids) > 1:
        primary = bed_list_group[0][0]
        samples_genes = '; '.join(sorted(f"{','.join(x[0].samples)}: {_genes_str(x[0].ref_gene_mappings)}"
                                         for x in bed_list_group))
        if _is_fusion(bed_list_group[0]):
            advice = "These are fusion isoforms, whose genes flair fusion chose, which flair combine can't reconcile."
        else:
            advice = ("Give flair combine an annotation (-f) to give each isoform the gene flair transcriptome would, "
                      "or check that all samples were run with the same annotation.")
        raise FlairError(f"isoforms of the same splice junctions have different genes, at {primary.chrom}:{primary.start}-"
                         f"{primary.end} ({primary.strand}): {samples_genes}.  {advice}")
    return ref_gene_ids.pop() if ref_gene_ids else ()


def _check_ref_transcripts(bed_list_group):
    "an error if the samples gave a combined isoform different reference transcripts"
    transcripts = {x[0].ref_transcript_id for x in bed_list_group if x[0].ref_transcript_id is not None}
    if len(transcripts) > 1:
        primary = bed_list_group[0][0]
        samples_transcripts = '; '.join(sorted(f"{','.join(x[0].samples)}: {x[0].ref_transcript_id or 'none'}"
                                               for x in bed_list_group))
        raise FlairError(f"isoforms of the same splice junctions have different reference transcripts, at "
                         f"{primary.chrom}:{primary.start}-{primary.end} ({primary.strand}): {samples_transcripts}.  "
                         f"Give flair combine an annotation (-f) to give each isoform the reference transcript flair "
                         f"transcriptome would, or check that all samples were run with the same annotation and FLAIR version.")


def _novel_gene_keys(combined, novel):
    """gene keys for the combined isoforms without an annotated gene, by their
    indexes in combined: grouped into novel genes as flair transcriptome groups
    them (_group_novel_by_splicing), on each chromosome and strand"""
    by_chrom_strand = {}
    for i in novel:
        bed = combined[i][0][0]
        by_chrom_strand.setdefault((bed.chrom, bed.strand), {})[i] = \
            Isoform(bed.chrom, bed.strand, tuple(Junc(gap.start, gap.end) for gap in bed.getGaps()), bed.start, bed.end)
    keys = {}
    for (chrom, strand), isoforms in by_chrom_strand.items():
        for group in _group_novel_by_splicing(sorted(isoforms), isoforms):
            key = ('novel', min(group))
            for i in group:
                keys[i] = key
    return keys


def assign_gene_keys(combined):
    """The gene of each combined isoform, and of each partner of a fusion, as keys
    standing for combined genes, which depend only on the combined isoforms, not on
    the order of the samples.  An isoform with annotated genes is in the gene of
    them, and the others are grouped into novel genes as flair transcriptome
    groups them.  A fusion partner without an annotated gene is in the combined
    gene of its sample genes' other isoforms, the first of them if more than one.
    Returns a list of (gene key, partner gene keys), and a map of gene key -> set
    of (samples, sample gene id)"""
    gene_keys = [None] * len(combined)
    novel = []
    for i, bed_list_group in enumerate(combined):
        ref_genes = _ref_genes(bed_list_group)
        _check_ref_transcripts(bed_list_group)
        if ref_genes != ():
            gene_keys[i] = ('ref', ref_genes)
        elif not _is_fusion(bed_list_group[0]):
            novel.append(i)
    for i, key in _novel_gene_keys(combined, novel).items():
        gene_keys[i] = key

    # sample genes -> the combined genes of their isoforms, the first by output order
    og_to_key = {}
    for i, bed_list_group in enumerate(combined):
        if not _is_fusion(bed_list_group[0]):
            for x in bed_list_group:
                og_to_key.setdefault((x[0].samples, x[0].gene_id), gene_keys[i])

    key_to_og = {}
    result = []
    for i, bed_list_group in enumerate(combined):
        og_genes = sorted((x[0].samples, x[0].gene_id) for x in bed_list_group)
        if gene_keys[i] is None:  # a fusion without annotated genes
            gene_keys[i] = _from_sample_genes(og_to_key, og_genes, ('fusion', i))
        key_to_og.setdefault(gene_keys[i], set()).update(og_genes)
        result.append((gene_keys[i], _partner_gene_keys(i, bed_list_group, og_to_key, key_to_og)))
    return result, key_to_og


def _from_sample_genes(og_to_key, og_genes, default):
    "the combined gene of the first of the sample genes with one"
    return next((og_to_key[g] for g in og_genes if g in og_to_key), default)


def _partner_gene_keys(i, bed_list_group, og_to_key, key_to_og):
    """gene keys of the partners of fusion isoform i: its annotated gene, or the
    combined gene of its sample genes"""
    ref_genes = _ref_genes(bed_list_group)
    partner_keys = []
    for p in range(len(bed_list_group[0][0].fused_genes)):
        og_partner = sorted((x[0].samples, x[0].fused_genes[p]) for x in bed_list_group)
        if p < len(ref_genes) and ref_genes[p] is not None:
            partner_key = ('ref', (ref_genes[p],))
        else:
            partner_key = _from_sample_genes(og_to_key, og_partner, ('fusion', i, p))
        key_to_og.setdefault(partner_key, set()).update(og_partner)
        partner_keys.append(partner_key)
    return partner_keys


def set_genes_from_annotation(bed_list_group, annots_by_chrom, genome):
    """Give the isoforms of a combined isoform the gene, reference transcript, CDS
    and productivity flair transcriptome gives an isoform of its splice junctions,
    strand and ends (those of the combined isoform, the group's first), whatever
    the samples gave them, with its protein sequence as aaseq_id, for
    write_combined to number.  Fusion isoforms keep theirs.  Returns whether a
    sample's genes differed"""
    primary = bed_list_group[0]
    if _is_fusion(primary):
        return False
    bed = primary[0]
    annots = annots_by_chrom.get(bed.chrom)
    isoform = Isoform(bed.chrom, bed.strand, tuple(Junc(gap.start, gap.end) for gap in bed.getGaps()), bed.start, bed.end)
    gene_ids, transcript_id = get_gene_name_firstpass(isoform, annots)
    genes = tuple(gene_ids) if gene_ids else ()
    Gene(None, genes, bed.chrom, bed.strand).add_isoform(isoform)
    isoform.ref_transcript_id = transcript_id
    thick_start, thick_end, productivity, aaseq = predict_productivity(
        isoform, annots.start_codon_count, annots.gene_to_cds_starts, annots.transcript_to_nmd_except, genome)
    differing = sorted({_genes_str(x[0].ref_gene_mappings) for x in bed_list_group if tuple(x[0].ref_gene_mappings) != genes})
    if differing:
        logging.debug(f"combined isoform at {bed.chrom}:{bed.start}-{bed.end} ({bed.strand}) given gene "
                      f"{','.join(genes) or 'novel'}, samples gave {'; '.join(differing)}")
    for bed_list in bed_list_group:
        sample_bed = bed_list[0]
        sample_bed.ref_gene_mappings = genes
        sample_bed.ref_transcript_id = transcript_id
        if sample_bed.transcript_class in ('basic', 'readthrough'):
            sample_bed.transcript_class = 'basic' if len(genes) <= 1 else 'readthrough'
    # as flair transcriptome writes them (convert_to_flair_bed), on the combined isoform
    bed.thickStart, bed.thickEnd, bed.productivity, bed.aaseq_id = thick_start, thick_end, productivity, aaseq
    return len(differing) > 0


def get_iso_counts_per_sample(bed_list_group, new_iso_id, iso_to_samples_to_counts):
    for bed_list in bed_list_group:
        for sample in bed_list[0].samples:
            if sample not in iso_to_samples_to_counts[new_iso_id]:
                iso_to_samples_to_counts[new_iso_id][sample] = 0
            iso_to_samples_to_counts[new_iso_id][sample] += bed_list[0].read_support

def correct_bed_fields_write_out(bed_list_group, new_iso_id, new_gene_id, new_fused_genes, new_iso_to_og,
                                 genome, bed_fh, fa_fh, iso_to_samples_to_counts):
    # need to adjust: new FLAIR isoform ID, new FLAIR gene ID, sum read support, median frac support, samples,
    # also need to adjust fused_gene_ids if relevant
    # keep dictionary of new isos to OG isos (with sample)

    tot_read_sup = sum([x[0].read_support for x in bed_list_group])
    med_usage = median([x[0].frac_support for x in bed_list_group])
    all_samples = [x[0].samples for x in bed_list_group]
    all_samples = tuple(sorted(list(set([x for xs in all_samples for x in xs]))))

    new_iso_to_og[new_iso_id] = set([(x[0].samples, x[0].name) for x in bed_list_group])
    iso_to_samples_to_counts[new_iso_id] = {}
    get_iso_counts_per_sample(bed_list_group, new_iso_id, iso_to_samples_to_counts)

    primary_bed_list = bed_list_group[0]
    my_sequence = ''
    for bed in primary_bed_list:
        bed.samples = all_samples
        bed.score = tot_read_sup
        bed.read_support = tot_read_sup
        bed.frac_support = med_usage
        bed.name = new_iso_id
        bed.gene_id = new_gene_id
        bed.fused_genes = tuple(new_fused_genes)
        # a sample's own isoform name and allele, from flair isoalleles, which don't
        # carry over to the combined isoform
        bed.source_isoform = bed.allele_group = None
        if bed.transcript_class == 'with_vars':
            bed.transcript_class = 'readthrough' if len([g for g in bed.ref_gene_mappings if g]) > 1 else 'basic'
        my_sequence += get_sequence_for_exons(genome, bed.chrom, bed.strand, bed.blocks)
        bed.write(bed_fh)
    fa_fh.write('>' + new_iso_id + '\n' + my_sequence + '\n')


def write_combined(combined, genome, output, number_aaseqs=False):
    """Write the combined isoforms, numbered in coordinate order, as are their
    genes, so the ids don't depend on the order of the samples.  With
    number_aaseqs, their aaseq_ids are protein sequences (set_genes_from_annotation),
    numbered in the same order and written to the aaseq file; otherwise they have
    none.  Returns the maps of new isoform and gene ids to the samples' ones, and
    the counts by sample"""
    combined.sort(key=_combined_sort_key)
    gene_keys, key_to_og = assign_gene_keys(combined)
    key_to_new = {}

    def new_gene(key):
        if key not in key_to_new:
            key_to_new[key] = f'FLG{len(key_to_new) + 1:08d}'
        return key_to_new[key]

    new_iso_to_og, iso_to_samples_to_counts, aaseq_to_id = {}, {}, {}
    with open(output + '.combined.isoforms.bed', 'w') as bed_fh, open(output + '.combined.isoforms.fa', 'w') as fa_fh:
        for i, (bed_list_group, (gene_key, partner_keys)) in enumerate(zip(combined, gene_keys)):
            new_gene_id = new_gene(gene_key)
            new_fused_genes = [new_gene(key) for key in partner_keys]
            for bed in bed_list_group[0]:
                if number_aaseqs and not _is_fusion(bed_list_group[0]) and bed.aaseq_id is not None:
                    bed.aaseq_id = aaseq_to_id.setdefault(bed.aaseq_id, f'FLP{len(aaseq_to_id) + 1:08d}')
                else:
                    # a sample's protein id: each sample numbers its own, so they
                    # don't name the same protein across samples
                    bed.aaseq_id = None
            correct_bed_fields_write_out(bed_list_group, f'FLT{i + 1:08d}', new_gene_id, new_fused_genes, new_iso_to_og,
                                         genome, bed_fh, fa_fh, iso_to_samples_to_counts)
    if number_aaseqs:
        with open(output + '.combined.aaseq.tsv', 'w') as fh:
            fh.write('aaseq_id\taaseq\n')
            for aaseq, aaseq_id in aaseq_to_id.items():
                fh.write(f'{aaseq_id}\t{aaseq}\n')
    new_gene_to_og = {new_id: key_to_og.get(key, set()) for key, new_id in key_to_new.items()}
    return new_iso_to_og, new_gene_to_og, iso_to_samples_to_counts


def write_map_files(output, new_iso_to_og, new_gene_to_og):
    with open(output + '.combined.isoform.map.txt', 'w') as fh:
        for new_id in new_iso_to_og:
            og_names = sorted(','.join(x[0]) + ':' + x[1] for x in new_iso_to_og[new_id])
            fh.write(new_id + '\t' + '; '.join(og_names) + '\n')
    with open(output + '.combined.gene.map.txt', 'w') as fh:
        for new_id in new_gene_to_og:
            og_names = sorted(','.join(x[0]) + ':' + x[1] for x in new_gene_to_og[new_id])
            fh.write(new_id + '\t' + '; '.join(og_names) + '\n')

def write_counts_file(output, allsamples, iso_to_samples_to_counts):
    with open(output + '.combined.isoform.counts.tsv', 'w') as fh:
        fh.write('\t'.join(['isoform_id'] + allsamples) + '\n')
        for new_id in iso_to_samples_to_counts:
            outline = [new_id]
            for sample in allsamples:
                if sample in iso_to_samples_to_counts[new_id]:
                    outline.append(str(iso_to_samples_to_counts[new_id][sample]))
                else:
                    outline.append('0')
            fh.write('\t'.join(outline) + '\n')

def combine(*, manifest, prefixes, isoform_beds, genome, output, end_window,
            min_frac_usage, remove_single_exon, max_ends, min_reads, annotation=None):
    logging.info('parsing manifest')
    sampledata = parse_input(manifest, prefixes, isoform_beds)
    annots_by_chrom = None
    if annotation is not None:
        logging.info('loading annotation GTF')
        annots_by_chrom = AnnotsByChrom(gtf_data_parser(annotation, attrs=GtfAttrsSet.FLAIR,
                                                        include_features=TRANSCRIPT_EXON_FEATURES))
    num_regened = 0

    logging.info('loading isoforms from individual samples')
    intronchaintoisos, allsamples = load_isoforms_by_junc_chain(sampledata)

    logging.info('combining isoforms')

    genome_fa = pysam.FastaFile(genome)
    combined = []
    for juncchaininfo in intronchaintoisos:
        all_juncs = [x[3] for x in juncchaininfo]
        og_isos = intronchaintoisos[juncchaininfo]
        if not remove_single_exon or (len(all_juncs) > 1 or len(all_juncs[0]) > 0):
            ends_to_iso_groups = combineIsos(og_isos, end_window)
            max_usage, tot_counts = 0, 0
            for bed_lists in ends_to_iso_groups:
                # For transcripts grouped by ends, sort by longest ends, ties by
                # support then sample, not by the order of the samples
                bed_lists.sort(key=_longest_first, reverse=True)
                max_usage = max((max_usage, max([x[0].frac_support for x in bed_lists])))
                tot_counts += sum([x[0].read_support for x in bed_lists])
            # sort all ends groups with longest ends first
            ends_to_iso_groups.sort(key=lambda x: _longest_first(x[0]), reverse=True)

            if all(len(juncs) == 0 for juncs in all_juncs):
                # single-exon isoforms all have the same, empty, junction chain: each
                # group of them with similar ends is an isoform of its own, which
                # max_ends mustn't merge, filtered on its own
                groups_after_end_filtering = [group for group in ends_to_iso_groups
                                              if max(x[0].frac_support for x in group) > min_frac_usage
                                              and sum(x[0].read_support for x in group) >= min_reads]
            elif max_usage > min_frac_usage and tot_counts >= min_reads:
                # sum all reads from all files, but do median usage from original usages
                # that way high counts files don't throw off the usage estimates
                groups_after_end_filtering = filter_groups(ends_to_iso_groups, max_ends, min_frac_usage, min_reads)
            else:
                groups_after_end_filtering = []
            for bed_list_group in groups_after_end_filtering:
                if annots_by_chrom is not None:
                    num_regened += set_genes_from_annotation(bed_list_group, annots_by_chrom, genome_fa)
                combined.append(bed_list_group)

    new_iso_to_og, new_gene_to_og, iso_to_samples_to_counts = write_combined(
        combined, genome_fa, output, number_aaseqs=annots_by_chrom is not None)

    if annots_by_chrom is not None:
        logging.info(f'{num_regened} combined isoforms given a gene from the annotation different from '
                     f'one or more samples\' (listed with flair --log-level DEBUG)')
    write_map_files(output, new_iso_to_og, new_gene_to_og)

    write_counts_file(output, allsamples, iso_to_samples_to_counts)

    bed_to_gtf(output + '.combined.isoforms.bed', output + '.combined.isoforms.gtf', is_flair_bed=True)

    make_big_bed(genome_fa, output + '.chrom.sizes', output + '.combined.isoforms')
    genome_fa.close()
