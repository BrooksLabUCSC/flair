#! /usr/bin/env python3

import logging
import pysam
from flair import FlairError
from flair.bed_to_gtf import bed_to_gtf
from flair.pycbio.hgdata.bed import BedReader
from flair.flair_bed import FlairBed
from statistics import median
from flair.isoform_data import make_big_bed, get_sequence_for_exons

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
            genome=args.genome, output=args.output, end_window=args.end_window,
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

def combineIsos(beds, endwindow):
    beds.sort(key=lambda x: [(y.start, y.end) for y in x])
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


def load_isos_by_id(bedfile):
    isotoinfo = {}
    mysamples = set()
    for bed in BedReader(bedfile, bedClass=FlairBed):
        juncchaininfo = (bed.pos_in_fusion, bed.chrom, bed.getGaps())
        if bed.name not in isotoinfo:
            isotoinfo[bed.name] = []
        isotoinfo[bed.name].append((juncchaininfo, bed))
        mysamples.update(set(bed.samples))
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

def get_gene_id(ref_gene_ids, curr_gene_ids, ref_gene_to_new, new_gene_to_og, gene_count, og_flair_gene_to_new):
    if len(ref_gene_ids) > 0:  # has annotated gene name
        if ref_gene_ids[0] in ref_gene_to_new:
            new_gene_id = ref_gene_to_new[ref_gene_ids[0]]
        else:
            new_gene_id = f'FLG{gene_count:08d}'
            gene_count += 1
            ref_gene_to_new[ref_gene_ids[0]] = new_gene_id
            new_gene_to_og[new_gene_id] = set()
    else:
        new_gene_id = None
        for g in curr_gene_ids:
            if g in og_flair_gene_to_new:
                new_gene_id = og_flair_gene_to_new[g]
                break
        if new_gene_id is None:
            new_gene_id = f'FLG{gene_count:08d}'
            gene_count += 1
            new_gene_to_og[new_gene_id] = set()
    new_gene_to_og[new_gene_id].update(set(curr_gene_ids))
    for g in curr_gene_ids:
        og_flair_gene_to_new[g] = new_gene_id
    return new_gene_id, gene_count

def get_new_gene_ids(bed_list_group, ref_gene_to_new, new_gene_to_og, og_flair_gene_to_new, gene_count):
    curr_gene_ids = [(x[0].samples, x[0].gene_id) for x in bed_list_group]
    ref_gene_ids = list(set([x[0].ref_gene_mappings for x in bed_list_group if x[0].ref_gene_mappings != ()]))
    if len(ref_gene_ids) > 1:
        raise FlairError("Matching isoforms from different samples have different genes - did you use consistent annotation files for all samples?")
    new_gene_id, gene_count = get_gene_id(ref_gene_ids, curr_gene_ids, ref_gene_to_new, new_gene_to_og, gene_count, og_flair_gene_to_new)

    curr_fused_genes = [(x[0].samples, x[0].fused_genes) for x in bed_list_group]
    new_fused_genes = [None for x in range(len(curr_fused_genes[0][1]))]
    for i in range(len(curr_fused_genes[0][1])):
        og_fused_genes = [(x[0], x[1][i]) for x in curr_fused_genes]
        fused_ref_gene = list(set([(x[i],) for x in ref_gene_ids if x[i] is not None]))
        new_fused_gene_id, gene_count = get_gene_id(fused_ref_gene, og_fused_genes, ref_gene_to_new, new_gene_to_og, gene_count, og_flair_gene_to_new)
        new_fused_genes[i] = new_fused_gene_id
    return new_gene_id, new_fused_genes, gene_count

def get_iso_counts_per_sample(bed_list_group, new_iso_id, iso_to_samples_to_counts):
    for bed_list in bed_list_group:
        for sample in bed_list[0].samples:
            if sample not in iso_to_samples_to_counts[new_iso_id]:
                iso_to_samples_to_counts[new_iso_id][sample] = 0
            iso_to_samples_to_counts[new_iso_id][sample] += bed_list[0].read_support

def correct_bed_fields_write_out(bed_list_group, new_iso_to_og, ref_gene_to_new, new_gene_to_og, og_flair_gene_to_new,
                                 genome, bed_fh, fa_fh, iso_count, gene_count, iso_to_samples_to_counts):
    # need to adjust: new FLAIR isoform ID, new FLAIR gene ID, sum read support, median frac support, samples,
    # also need to adjust fused_gene_ids if relevant
    # keep dictionary of new isos to OG isos (with sample), new genes to OG genes (with sample)

    tot_read_sup = sum([x[0].read_support for x in bed_list_group])
    med_usage = median([x[0].frac_support for x in bed_list_group])
    all_samples = [x[0].samples for x in bed_list_group]
    all_samples = tuple(sorted(list(set([x for xs in all_samples for x in xs]))))

    new_iso_id = f'FLT{iso_count:08d}'
    iso_count += 1
    new_iso_to_og[new_iso_id] = set([(x[0].samples, x[0].name) for x in bed_list_group])
    iso_to_samples_to_counts[new_iso_id] = {}

    new_gene_id, new_fused_genes, gene_count = get_new_gene_ids(bed_list_group, ref_gene_to_new, new_gene_to_og, og_flair_gene_to_new, gene_count)

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
        my_sequence += get_sequence_for_exons(genome, bed.chrom, bed.strand, bed.blocks)
        bed.write(bed_fh)
    fa_fh.write('>' + new_iso_id + '\n' + my_sequence + '\n')
    return iso_count, gene_count

def write_map_files(output, new_iso_to_og, new_gene_to_og):
    with open(output + '.combined.isoform.map.txt', 'w') as fh:
        for new_id in new_iso_to_og:
            og_names = [','.join(x[0]) + ':' + x[1] for x in new_iso_to_og[new_id]]
            fh.write(new_id + '\t' + '; '.join(og_names) + '\n')
    with open(output + '.combined.gene.map.txt', 'w') as fh:
        for new_id in new_gene_to_og:
            og_names = [','.join(x[0]) + ':' + x[1] for x in new_gene_to_og[new_id]]
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
            min_frac_usage, remove_single_exon, max_ends, min_reads):
    logging.info('parsing manifest')
    sampledata = parse_input(manifest, prefixes, isoform_beds)

    logging.info('loading isoforms from individual samples')
    intronchaintoisos, allsamples = load_isoforms_by_junc_chain(sampledata)

    logging.info('combining isoforms')

    new_iso_to_og, new_gene_to_og, ref_gene_to_new, og_flair_gene_to_new = {}, {}, {}, {}
    iso_to_samples_to_counts = {}
    iso_count, gene_count = 1, 1
    genome_fa = pysam.FastaFile(genome)
    with open(output + '.combined.isoforms.bed', 'w') as bed_fh, open(output + '.combined.isoforms.fa', 'w') as fa_fh:
        for juncchaininfo in intronchaintoisos:
            all_juncs = [x[2] for x in juncchaininfo]
            og_isos = intronchaintoisos[juncchaininfo]
            if not remove_single_exon or (len(all_juncs) > 1 or len(all_juncs[0]) > 0):
                ends_to_iso_groups = combineIsos(og_isos, end_window)
                max_usage, tot_counts = 0, 0
                for bed_lists in ends_to_iso_groups:
                    # For transcripts grouped by ends, sort by longest ends
                    bed_lists.sort(key=lambda x: (x[0].end - x[0].start, x[0].read_support), reverse=True)
                    max_usage = max((max_usage, max([x[0].frac_support for x in bed_lists])))
                    tot_counts += sum([x[0].read_support for x in bed_lists])
                # sort all ends groups with longest ends first
                ends_to_iso_groups.sort(key=lambda x: (x[0][0].end - x[0][0].start, x[0][0].read_support), reverse=True)

                if max_usage > min_frac_usage and tot_counts >= min_reads:
                    # sum all reads from all files, but do median usage from original usages
                    # that way high counts files don't throw off the usage estimates
                    groups_after_end_filtering = filter_groups(ends_to_iso_groups, max_ends, min_frac_usage, min_reads)
                    for bed_list_group in groups_after_end_filtering:
                        iso_count, gene_count = correct_bed_fields_write_out(bed_list_group, new_iso_to_og, ref_gene_to_new, new_gene_to_og,
                                                                             og_flair_gene_to_new, genome_fa, bed_fh, fa_fh, iso_count, gene_count, iso_to_samples_to_counts)

    write_map_files(output, new_iso_to_og, new_gene_to_og)

    write_counts_file(output, allsamples, iso_to_samples_to_counts)

    bed_to_gtf(output + '.combined.isoforms.bed', output + '.combined.isoforms.gtf', is_flair_bed=True)

    make_big_bed(genome_fa, output + '.chrom.sizes', output + '.combined.isoforms')
    genome_fa.close()
