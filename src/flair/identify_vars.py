#! /usr/bin/env python3

import argparse
import logging
import pysam
import vcfpy
import shutil
import multiprocessing as mp
from flair.io_utils import make_temp_dir
from flair.pycbio.hgdata.bed import BedReader

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, type=str,
                        help="bam file (this is your primary bam input, required)")
    parser.add_argument('-o', '--output', default='flair',
                        help="output prefix. default: 'flair'")
    parser.add_argument('-g', '--genome',
                        type=str, required=True, help='FastA of reference genome')
    parser.add_argument('-t', '--threads', type=int,
                        action='store', default=4, help='minimap2 number of threads (4)')

    parser.add_argument('--min_cov', type=int, default=6,
                        help='Minimum total read coverage of a site to call indels')
    parser.add_argument('--min_var_reads', type=int, default=3,
                        help='Minimum number of variant reads required to call an indel')
    parser.add_argument('--min_vaf_pct', type=int, default=10,
                        help='Minimum percentage of reads required to call indel (int 1-100, default=10)')
    parser.add_argument('--min_sj_dist', type=int, default=3,
                        help='Minimum number of basepairs an indel must be from a splice junction in the same read in order to be considered')
    parser.add_argument('--min_read_end_dist', type=int, default=6,
                        help='Minimum number of basepairs an indel must be from the end of the read in order to be considered')
    parser.add_argument('--density_dist', type=int, default=20,
                        help='distance from variant considered for density calculation')
    parser.add_argument('--density_count', type=int, default=4,
                        help='number of nearby variants required to be flagged as a dense region')
    parser.add_argument('--num_repeats_for_filtering', type=int, default=1,
                        help='Required number of repeats to flag indels in tandem repeat regions as sequencing errors. '
                        'Set this to 1 to remove all calls in tandem repeat regions, set it high if you are interested in getting calls for these regions.')
    # parser.add_argument('--max_depth', type=int, default=100000,
    #                     help='Maximum read depth allowed to search for indels at locus (lower values will miss some genes, but will run substantially faster)')
    parser.add_argument('--region_bed', type=str,
                        help='bed (3+) file of regions on which to call indels')
    parser.add_argument('--identify_indels', type=str, default='yes',
                        help="[yes/no], default:yes - specify yes if you want this caller to identify SNVs")
    parser.add_argument('--identify_snvs', type=str, default='no',
                        help="[yes/no], default:no - specify yes if you want this caller to identify SNVs")

    # parser.add_argument('--keep_intermediate', default=False, action='store_true',
    #                     help='''specify if intermediate and temporary files are to be kept for debugging.''')
    args = parser.parse_args()
    args.identify_indels = True if args.identify_indels == 'yes' else False
    args.identify_snvs = True if args.identify_snvs == 'yes' else False
    if not (args.identify_indels or args.identify_snvs):
        raise ValueError('Please tell this caller to identify at least 1 of: indels, snvs')
    return args


def add_block_cov(ref_pos, block_len, r_start, r_end, pos_to_cov, is_del):
    for pos in range(ref_pos, ref_pos + block_len):
        if r_start <= pos <= r_end:
            if pos not in pos_to_cov:
                pos_to_cov[pos] = [0, 0]
            pos_to_cov[pos][1] += 1
            if not is_del:
                pos_to_cov[pos][0] += 1

def add_covered_pos(cigar, alignstart, pos_to_cov, r_start, r_end):
    # pos_to_cov is keyed by the 1-based position of the covered base
    ref_pos, quer_pos = alignstart + 1, 0
    introns = []
    for block in cigar:
        if block[0] in {0, 7, 8}:  # match, consumes both
            add_block_cov(ref_pos, block[1], r_start, r_end, pos_to_cov, False)
            ref_pos += block[1]
            quer_pos += block[1]
        elif block[0] in {1, 4}:  # consumes query ###1 is insertion
            quer_pos += block[1]
        elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
            if block[0] == 3:  # intron
                introns.append((ref_pos, ref_pos + block[1]))
            elif block[0] == 2:  # actually count deletions towards coverage
                add_block_cov(ref_pos, block[1], r_start, r_end, pos_to_cov, True)
            ref_pos += block[1]
    return introns, quer_pos

def check_dist_from_ends(min_read_end_dist, ref_start, ref_end, align_start, align_end):
    # NOTE: after testing, the distance from the end of the alignment tends to be more important than the distance from the actual read ends
    # it's not about the quality of the terminal read bases, it's about the quality of the end portion of the alignment (frequently it's about the combo of both)
    if ref_start >= align_start + min_read_end_dist and ref_end <= align_end - min_read_end_dist:
        return True
    return False

def check_dist_from_sj(min_sj_dist, ref_start, ref_end, introns):
    for i_start, i_end in introns:
        if i_start - min_sj_dist <= ref_end <= i_start + min_sj_dist or i_end - min_sj_dist <= ref_start <= i_end + min_sj_dist:
            return False
    return True

def add_var_to_dict(pos_to_var, ref_pos, indel_key, flagged_bad):
    if ref_pos not in pos_to_var:
        pos_to_var[ref_pos] = {}
    if indel_key not in pos_to_var[ref_pos]:
        pos_to_var[ref_pos][indel_key] = [0, 0]
    pos_to_var[ref_pos][indel_key][0] += 1
    if flagged_bad:
        pos_to_var[ref_pos][indel_key][1] += 1

def add_indel_info(cigar, align_start, align_end, region, readseq, genome, introns, read_length, min_read_end_dist, min_sj_dist, pos_to_var):
    chrom, r_start, r_end = region
    # ref_pos counts 0-based, so at an indel it holds the 1-based position of the
    # anchor base before it, which is what VCF POS wants and what is written.  SNVs
    # key pos_to_var by the 1-based position of the changed base instead, so the two
    # kinds of entry sit one base apart; get_cov_for_var_pos reconciles them
    ref_pos, quer_pos = align_start, 0
    # here is where we can filter out indels that are too close to splice junctions or read ends
    for block in cigar:
        if block[0] in {0, 7, 8}:  # match, consumes both
            ref_pos += block[1]
            quer_pos += block[1]
        elif block[0] in {1, 4}:  # consumes query ###1 is insertion
            if block[0] == 1:  # insertion
                if r_start <= ref_pos <= r_end:
                    flagged_bad = not (check_dist_from_ends(min_read_end_dist, ref_pos, ref_pos, align_start, align_end) and check_dist_from_sj(min_sj_dist, ref_pos, ref_pos, introns))
                    insertion_seq = readseq[quer_pos - 1:quer_pos + block[1]].upper()
                    add_var_to_dict(pos_to_var, ref_pos, ('INS', insertion_seq[0], insertion_seq), flagged_bad)
            quer_pos += block[1]
        elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
            if block[0] == 2:
                if r_start <= ref_pos <= r_end:
                    flagged_bad = not (check_dist_from_ends(min_read_end_dist, ref_pos, ref_pos + block[1], align_start, align_end) and check_dist_from_sj(min_sj_dist, ref_pos, ref_pos + block[1], introns))
                    deletion_seq = genome.fetch(chrom, ref_pos - 1, ref_pos + block[1]).upper()
                    add_var_to_dict(pos_to_var, ref_pos, ('DEL', deletion_seq, deletion_seq[0]), flagged_bad)
            ref_pos += block[1]

def len_same_base(seq):
    # get number of chars that match first char
    if len(seq) == 1:
        return 1
    else:
        m = 1
        while m < len(seq) and seq[m] == seq[0]:
            m += 1
        return m

def check_homopoly(indel_type, indel_seq, indel_vaf, genome, chrom, pos):
    indel_seq = indel_seq[1:]
    if len(indel_seq) == len_same_base(indel_seq):
        # ONLY NEED TO CHECK RIGHT SIDE DUE TO MINIMAP ALIGNMENT CONVENTION
        if indel_type == 'DEL':
            rightseq = genome.fetch(chrom, pos + len(indel_seq), pos + len(indel_seq) + 5).upper()
        else:
            rightseq = genome.fetch(chrom, pos + 1, pos + 6).upper()
        rightlen = len_same_base(rightseq)
        if rightseq[0] == indel_seq[0]:
            if rightlen >= 3:
                return True
            elif rightlen >= 2 and indel_vaf < 0.2:
                return True
    return False

def check_repeat(indel_type, indel_seq, genome, chrom, pos, num_repeats_for_filtering):
    if len(indel_seq) > 2:  # actual indel is more than 1bp, otherwise keep the preceding base for this check
        indel_seq = indel_seq[1:]
    elif indel_type == 'DEL':
        pos -= 1
    indel_len = len(indel_seq)
    if indel_type == 'DEL':
        rightseq = genome.fetch(chrom, pos + indel_len, pos + (indel_len * (num_repeats_for_filtering + 1))).upper()
    else:
        rightseq = genome.fetch(chrom, pos, pos + (indel_len * num_repeats_for_filtering)).upper()
    matches_repeats = []
    for i in range(num_repeats_for_filtering):
        if indel_seq == rightseq[indel_len * i:indel_len * (i + 1)]:
            matches_repeats.append(True)
        else:
            matches_repeats.append(False)
    if all(matches_repeats):
        return True
    return False

def add_ss_to_ref(introns, donor_counts, acceptor_counts):
    for i_start, i_end in introns:
        if i_start not in donor_counts:
            donor_counts[i_start] = 0
        donor_counts[i_start] += 1
        if i_end not in acceptor_counts:
            acceptor_counts[i_end] = 0
        acceptor_counts[i_end] += 1


def check_ss_dist_against_ref(pos, tot_cov, donor_counts, acceptor_counts, min_sj_dist, min_var_reads, min_vaf):
    # allowing min_var_reads and min_vaf to sub in for sj counts to try to get a sense of the user's error tolerance
    for i in range(pos - min_sj_dist, pos + min_sj_dist + 1):
        if (i in donor_counts and donor_counts[i] >= min_var_reads and donor_counts[i] / tot_cov > min_vaf) \
                or (i in acceptor_counts and acceptor_counts[i] >= min_var_reads and acceptor_counts[i] / tot_cov > min_vaf):
            return True
    return False

def get_snvs_from_bam(aligned_bases, read_seq, chrom, align_start, align_end, r_start, r_end, genome, read_length, introns, pos_to_var, min_read_end_dist, min_sj_dist):
    for quer_pos, ref_pos, base in aligned_bases:
        if base.islower():
            if r_start <= ref_pos <= r_end:
                flagged_bad = not (check_dist_from_ends(min_read_end_dist, ref_pos, ref_pos, align_start, align_end) and check_dist_from_sj(min_sj_dist, ref_pos, ref_pos, introns))
                snv_alt = read_seq[quer_pos].upper()
                snv_ref = genome.fetch(chrom, ref_pos, ref_pos + 1).upper()
                add_var_to_dict(pos_to_var, ref_pos + 1, ('SNV', snv_ref, snv_alt), flagged_bad)

def get_indels_from_bam(region, bam_file, genome, min_read_end_dist, min_sj_dist, identify_indels, identify_snvs):
    pos_to_cov, pos_to_var, donor_counts, acceptor_counts = {}, {}, {}, {}
    c = 0
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        chrom, r_start, r_end = region
        for a in bam.fetch(chrom, r_start, r_end):
            if not a.is_secondary:
                c += 1
                introns, read_length = add_covered_pos(a.cigartuples, a.reference_start, pos_to_cov, region[1], region[2])
                add_ss_to_ref(introns, donor_counts, acceptor_counts)
                if identify_indels:
                    add_indel_info(a.cigartuples, a.reference_start, a.reference_end, region, a.query_sequence, genome, introns, read_length, min_read_end_dist, min_sj_dist, pos_to_var)
                if identify_snvs:
                    get_snvs_from_bam(a.get_aligned_pairs(with_seq=True, matches_only=True), a.query_sequence, chrom, a.reference_start, a.reference_end, r_start, r_end, genome, read_length, introns, pos_to_var, min_read_end_dist, min_sj_dist)
                # inside the not-secondary block, and logging rather than print: c only
                # advances here, so every secondary alignment after a multiple of 10000
                # reprinted the same line, on the stdout of a worker process
                if c % 10000 == 0:
                    logging.info(f"{region[0]}: {c} reads, at {a.reference_start}")
    return pos_to_cov, pos_to_var, donor_counts, acceptor_counts

def get_cov_for_var_pos(pos, pos_to_cov):
    """Coverage at a variant position.  Both pos and pos + 1 are probed because an
    indel's key is its anchor base while an SNV's is the changed base, so which of the
    two carries the coverage depends on the kind of variant.  Verified against the
    alignments: a 2 bp deletion of chr12 25205729-25205730, 1-based, is written as
    POS 25205728 REF CTT ALT C, which is what these keys produce.

    A position can be missing entirely: an insertion in the middle of an intron, which
    minimap does produce, has no covered base."""
    # NOTE: account for insertion/deletion at edge of exon (introns have no coverage)
    tot_cov = []
    if pos in pos_to_cov:
        tot_cov.append(pos_to_cov[pos])
    if pos + 1 in pos_to_cov:
        tot_cov.append(pos_to_cov[pos + 1])
    if len(tot_cov) == 0:
        # raise ValueError(f'No coverage at indel locus: {pos}')
        # NOTE: this can happen if there is an insertion in the middle of an intron (strange but minimap does it for some reason)
        return [0, 0]
    elif len(tot_cov) == 1:
        return tot_cov[0]
    else:
        if abs(tot_cov[0][1] - tot_cov[1][1]) > 5 or abs(tot_cov[0][1] - tot_cov[1][1]) / tot_cov[0][1] > 0.1:
            return max(tot_cov, key=lambda x: x[1])
        else:
            return tot_cov[0]

def get_indel_filter_labels(indel_type, indel_seq, tot_cov, indel_reads, indel_vaf, genome, chrom, pos, donor_counts, acceptor_counts, min_sj_dist, min_var_reads, min_vaf, num_repeats_for_filtering):
    filters = []
    if not (indel_reads >= min_var_reads and indel_vaf >= min_vaf):
        filters.append('lq')
    if indel_type != 'SNV':
        if check_homopoly(indel_type, indel_seq, indel_vaf, genome, chrom, pos):
            filters.append('hp')
        elif check_repeat(indel_type, indel_seq, genome, chrom, pos, num_repeats_for_filtering):
            filters.append('tr')
    if check_ss_dist_against_ref(pos, tot_cov, donor_counts, acceptor_counts, min_sj_dist, min_var_reads, min_vaf):
        filters.append('sj')
    return filters

def filter_output_indel(pos_to_var_filtered, pos, pos_var_info, indel_type, ref_seq, var_seq, tot_cov_nodel, tot_cov_withdel,
                        genome, chrom, donor_counts, acceptor_counts, min_var_reads, min_vaf, min_sj_dist, num_repeats_for_filtering):
    tot_cov = tot_cov_nodel if indel_type == 'SNV' else tot_cov_withdel
    indel_reads, flagged_bad_reads = pos_var_info[(indel_type, ref_seq, var_seq)]
    indel_good_reads = indel_reads - flagged_bad_reads
    indel_vaf = indel_reads / tot_cov
    indel_good_vaf = indel_good_reads / tot_cov
    if indel_reads >= min_var_reads and indel_vaf >= min_vaf:
        filter_seq = var_seq if indel_type != 'DEL' else ref_seq
        filters = get_indel_filter_labels(indel_type, filter_seq, tot_cov, indel_good_reads, indel_good_vaf,
                                          genome, chrom, pos, donor_counts, acceptor_counts, min_sj_dist, min_var_reads, min_vaf, num_repeats_for_filtering)
        # out.write('\t'.join([str(x) for x in (chrom, pos, indel_type, ref_seq, var_seq, indel_reads, tot_cov, ','.join(filters))]) + '\n')
        if pos not in pos_to_var_filtered:
            pos_to_var_filtered[pos] = {}
        pos_to_var_filtered[pos][(ref_seq, var_seq)] = [chrom, pos, indel_type, ref_seq, var_seq, indel_reads, tot_cov, filters]

def process_var_pos(pos, pos_to_var_filtered, pos_to_var, pos_to_cov, genome, chrom, donor_counts, acceptor_counts, min_cov, min_var_reads, min_vaf, min_sj_dist, num_repeats_for_filtering):
    tot_cov_nodel, tot_cov_withdel = get_cov_for_var_pos(pos, pos_to_cov)
    if tot_cov_withdel >= min_cov:
        for indel_type, ref_seq, var_seq in pos_to_var[pos]:
            filter_output_indel(pos_to_var_filtered, pos, pos_to_var[pos], indel_type, ref_seq, var_seq, tot_cov_nodel, tot_cov_withdel,
                                genome, chrom, donor_counts, acceptor_counts, min_var_reads, min_vaf, min_sj_dist, num_repeats_for_filtering)


def process_reads_for_region(args):
    temp_dir, bam_file, region, genome_file, min_cov, min_var_reads, min_vaf, min_read_end_dist, min_sj_dist, num_repeats_for_filtering, identify_indels, identify_snvs, density_dist, density_count = args
    with pysam.FastaFile(genome_file) as genome, open(temp_dir + '-'.join([str(x) for x in region]) + '.txt', 'w') as out:
        pos_to_cov, pos_to_var, donor_counts, acceptor_counts = get_indels_from_bam(region, bam_file, genome, min_read_end_dist, min_sj_dist, identify_indels, identify_snvs)
        pos_to_var_filtered = {}
        for pos in pos_to_var:
            process_var_pos(pos, pos_to_var_filtered, pos_to_var, pos_to_cov, genome, region[0], donor_counts, acceptor_counts, min_cov, min_var_reads, min_vaf, min_sj_dist, num_repeats_for_filtering)
        for pos in pos_to_var_filtered:
            other_vars = 0
            for p2 in range(pos - density_dist, pos + density_dist + 1):
                if p2 != pos and p2 in pos_to_var_filtered:
                    other_vars += len(pos_to_var_filtered[p2])
            if other_vars >= density_count:
                for k in pos_to_var_filtered[pos]:
                    pos_to_var_filtered[pos][k][-1].append('dn')
            for k in pos_to_var_filtered[pos]:
                chrom, pos, indel_type, ref_seq, var_seq, indel_reads, tot_cov, filters = pos_to_var_filtered[pos][k]
                out.write('\t'.join([str(x) for x in (chrom, pos, indel_type, ref_seq, var_seq, indel_reads, tot_cov, ','.join(filters))]) + '\n')

def get_regions_from_bed(region_bed):
    regions = []
    for bed in BedReader(region_bed):
        regions.append((bed.chrom, bed.chromStart, bed.chromEnd))
    return regions

def get_regions_by_chrom(genome_file):
    regions = []
    with pysam.FastaFile(genome_file) as genome:
        for chrom in sorted(list(genome.references)):
            regions.append((chrom, 0, genome.get_reference_length(chrom)))
    return regions

def generate_new_vcf_header(genome):
    new_header = vcfpy.Header()
    new_header.add_line(vcfpy.HeaderLine('fileformat', 'VCFv4.2'))
    new_header.add_line(vcfpy.HeaderLine('source', 'FLAIR indels'))
    for chrom in genome.references:
        new_header.add_line(vcfpy.ContigHeaderLine.from_mapping({'ID': chrom, 'length': genome.get_reference_length(chrom)}))
    new_header.add_line(vcfpy.FilterHeaderLine.from_mapping({'ID': 'sj', 'Description': 'too close to splice junction'}))
    new_header.add_line(vcfpy.FilterHeaderLine.from_mapping({'ID': 'hp', 'Description': 'homopolymer region'}))
    new_header.add_line(vcfpy.FilterHeaderLine.from_mapping({'ID': 'tr', 'Description': 'tandem repeat of adjacent sequence'}))
    new_header.add_line(vcfpy.FilterHeaderLine.from_mapping({'ID': 'dn', 'Description': 'dense repeat region, likely misalignment and/or low complexity region'}))
    new_header.add_line(vcfpy.FilterHeaderLine.from_mapping({'ID': 'lq', 'Description': 'variant is in a low quality region of supporting reads (near read ends or splice junctions)'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'GT', 'Number': 1, 'Type': 'String', 'Description': 'Genotype'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'DP', 'Number': 1, 'Type': 'Integer', 'Description': 'Read depth'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'AD', 'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths for the ref and alt alleles in the order listed'}))
    new_header.samples = vcfpy.SamplesInfos(['SAMPLE'])
    return new_header

def write_output_vcf(genome_file, output, regions, temp_dir):
    with pysam.FastaFile(genome_file) as genome:
        vcf_header = generate_new_vcf_header(genome)
        with vcfpy.Writer.from_path(output + '.vcf', vcf_header) as writer:
            format_strings = ['GT', 'DP', 'AD']
            for region in regions:
                for line in open(temp_dir + '-'.join([str(x) for x in region]) + '.txt'):
                    chrom, start_pos, indel_type, ref_seq, alt_seq, indel_reads, tot_cov, filters = line.rstrip('\n').split('\t')
                    indel_reads, tot_cov = int(indel_reads), int(tot_cov)
                    filters = [x for x in filters.split(',') if x != '']
                    if len(filters) == 0:
                        filters.append('PASS')
                    alt_desc = vcfpy.Substitution(indel_type, alt_seq)
                    indel_vaf = round(indel_reads / tot_cov, 3)
                    gt = '0/1' if indel_vaf < 0.95 else '1/1'
                    these_calls = [vcfpy.Call('SAMPLE', {'GT': gt, 'DP': tot_cov, 'AD': [tot_cov - indel_reads, indel_reads]})]
                    new_record = vcfpy.Record(chrom, int(start_pos), [], ref_seq, [alt_desc], 100, filters, {}, format_strings, these_calls)
                    writer.write_record(new_record)

def get_indels():
    args = parse_args()
    temp_dir = make_temp_dir(args.output)
    if args.region_bed is not None:
        regions = get_regions_from_bed(args.region_bed)
    else:
        regions = get_regions_by_chrom(args.genome)

    packed = []
    for region in regions:
        packed.append((temp_dir, args.bam, region, args.genome, args.min_cov, args.min_var_reads, args.min_vaf_pct / 100,
                       args.min_read_end_dist, args.min_sj_dist, args.num_repeats_for_filtering,
                       args.identify_indels, args.identify_snvs, args.density_dist, args.density_count))

    if args.threads == 1:
        for p in packed:
            process_reads_for_region(p)
    else:
        mp.set_start_method('fork', force=True)
        with mp.Pool(args.threads) as pool:
            pool.map(process_reads_for_region, packed)
    write_output_vcf(args.genome, args.output, regions, temp_dir)

    shutil.rmtree(temp_dir)


if __name__ == "__main__":
    get_indels()
