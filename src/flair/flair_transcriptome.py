#! /usr/bin/env python3

import argparse
import os
import pipettor
import uuid
import shutil
import pysam
import logging
from flair.flair_align import inferMM2JuncStrand, intron_chain_to_exon_starts
from flair.ssUtils import addOtherJuncs, gtfToSSBed
from flair.ssPrep import buildIntervalTree, ssCorrect
import multiprocessing as mp
from collections import Counter
from flair import FlairInputDataError


def get_args():
    parser = argparse.ArgumentParser(description='generates confident transcript models directly from a bam file '
                                                 'of aligned long rna-seq reads')
    parser.add_argument('-b', '--genome_aligned_bam', required=True,
                        help='Sorted and indexed bam file aligned to the genome')
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help='FastA of reference genome, can be minimap2 indexed')
    parser.add_argument('-o', '--output', default='flair',
                        help='output file name base for FLAIR isoforms (default: flair)')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='number of threads to run with - related to parallel_mode')
    parser.add_argument('-f', '--gtf', default='',
                        help='GTF annotation file, used for identifying annotated isoforms')

    mutexc = parser.add_mutually_exclusive_group(required=False)
    mutexc.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    mutexc.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from long-read alignment with intronProspector)')
    parser.add_argument('--junction_support', type=int, default=1,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support.')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')
    parser.add_argument('-w', '--end_window', type=int, default=100,
                        help='window size for comparing TSS/TES (100)')

    parser.add_argument('--sjc_support', type=int, default=1,
                        help='''minimum number of supporting reads for a spliced isoform (1)''')
    parser.add_argument('--se_support', type=int, default=3,
                        help='''minimum number of supporting reads for a single exon isoform (3)''')
    parser.add_argument('--frac_support', type=float, default=0.05,
                        help='''minimum fraction of gene locus support for isoform to be called
                         default: 0.05, only isoforms that make up more than 5 percent of the gene locus are reported. Set to 0 for max recall''')


    parser.add_argument('--no_stringent', default=False, action='store_true',
                        help='''specify if all supporting reads don't need to be full-length
                (aligned to first and last exons of transcript). Use this for fragmented libraries, but understand that it will impact precision.''')
    parser.add_argument('--no_check_splice', default=False, action='store_true',
                        help='''don't enforce accurate alignment around splice site. Specify this for libraries with high error rates, but it will reduce precision''')


    parser.add_argument('--no_align_to_annot', default=False, action='store_true',
                        help='''related to old annotation_reliant, now specify if you don't want
                        an initial alignment to the annotated sequences and only want transcript
                        detection from the genomic alignment.
                         Will be slightly faster but less accurate if the annotation is good''')
    parser.add_argument('-n', '--no_redundant', default='none',
                        help='''For each unique splice junction chain, report options include:
            none-- multiple supported TSSs/TESs chosen for each set of splice junctions (modulated by max_ends);
            longest--single TSS/TES chosen to maximize length;
            best_only--single most supported TSS/TES used in conjunction chosen (none)''')
    parser.add_argument('--max_ends', type=int, default=1,
                        help='maximum number of TSS/TES picked per isoform (1) make higher for more precise end detection')
    parser.add_argument('--filter', default='nosubset',
                        help='''Report options include:
            nosubset--any isoforms that are a proper set of another isoform are removed;
            bysupport--subset isoforms are removed based on support;
            comprehensive--default set + all subset isoforms;
            ginormous--comprehensive set + single exon subset isoforms''')
    parser.add_argument('--quality', default=1, type=int,
                        help='minimum mapping quality threshold to consider genomic alignments for defining transcripts')
    parser.add_argument('--parallel_mode', default='auto:1GB',
                        help='''parallelization mode. Default: "auto:1GB" This indicates an automatic threshold where
                            if the file is less than 1GB, parallelization is done by chromosome, but if it's larger,
                            parallelization is done by region of non-overlapping reads. Other modes: bychrom, byregion,
                            auto:xGB - for setting the auto threshold, it must be in units of GB.''')
    parser.add_argument('--predict_cds', default=False, action='store_true',
                        help='specify if you want to predict the CDS of the final isoforms. '
                             'Will be output in the final bed file but not the gtf file. '
                             'Productivity annotation is also added in the name field, '
                             'which is detailed further in the predictProductivity documentation')
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.
                Intermediate files include: promoter-supported reads file,
                read assignments to firstpass isoforms''')
    parser.add_argument('--keep_sup', default=False, action='store_true',
                        help='''specify if you want to keep supplementary alignments to define isoforms''')
    parser.add_argument('--end_norm_dist',
                        help='specify the number of basepairs to extend transcript ends if you want to '
                             'normalize them across transcripts in a gene and extend them')
    parser.add_argument('--output_endpos', default=False, action='store_true',
                        help='specify if you want to ouptut a setparate file with corrected read end positions. '
                             'For development purposes')
    parser.add_argument('--output_bam', default=False, action='store_true',
                        help='output intermediate bams aligned to the transcriptome. '
                             'Only works with --keep_intermediate, for debugging')


    args = parser.parse_args()

    check_file_paths(args)
    args = add_preset_args(args)
    # if args.mm2_args:
    #     args.mm2_args = args.mm2_args.split(',')
    return args


def check_file_paths(args):
    if not args.genome_aligned_bam:
        raise FlairInputDataError('Please include the --genome_aligned_bam option')
    if not args.genome:
        raise FlairInputDataError('Please include the --genome option\n')
    if not os.path.exists(args.genome_aligned_bam):
        raise FlairInputDataError(f'Aligned reads file path does not exist: {args.genome_aligned_bam}')
    if not os.path.exists(args.genome):
        raise FlairInputDataError(f'Genome file path does not exist: {args.genome}')
    if not (args.parallel_mode in {'bychrom', 'byregion'}
            or (args.parallel_mode[:5] == 'auto:'
                and ((args.parallel_mode[-2:] == 'GB' and args.parallel_mode[5:-2].replace(".", "").isnumeric())
                     or args.parallel_mode[5:].replace(".", "").isnumeric()))):
        raise FlairInputDataError(
            f'parallel_mode {args.parallel_mode} is not in an allowed format. See docs for allowed formats')


def add_preset_args(args):
    args.mm2_args = []
    args.trust_ends = False
    args.remove_internal_priming = False
    args.isoformtss = True
    return args


def make_correct_temp_dir():
    temp_dir_name = str(uuid.uuid4())
    try:
        current_directory = os.getcwd()
        temp_dir = os.path.join(current_directory, temp_dir_name)
        os.mkdir(temp_dir)
    except OSError:
        raise OSError(f"Creation of the directory {temp_dir_name} failed")
    return temp_dir + '/'


def binary_search(query, data):
    """ Query is a coordinate interval. Binary search for the query in sorted data,
        which is a list of coordinates. Finishes when an overlapping value of query and
        data exists and returns the index in data. """
    i = int(round(len(data) / 2))  # binary search prep
    lower, upper = 0, len(data)
    while True:
        if upper - lower < 2:  # stop condition but not necessarily found
            break
        if data[i][1] < query[0]:
            lower = i
            i = int(round((i + upper) / 2))
        elif data[i][0] > query[1]:
            upper = i
            i = int(round((lower + i) / 2))
        else:  # found
            break
    return i


def generate_known_SS_database(args, temp_dir):
    # Convert gtf to bed and split by chromosome.
    # logging.info(f'Building ss database: {temp_dir}')
    juncs, chromosomes, knownSS = dict(), set(), dict()  # initialize juncs for adding to db

    if args.gtf:
        juncs, chromosomes, knownSS = gtfToSSBed(args.gtf, knownSS, False, False, False)

    # Do the same for the other juncs file.
    if args.junction_tab or args.junction_bed:
        if args.junction_tab:
            shortread, type = args.junction_tab, 'tab'
        else:
            shortread, type = args.junction_bed, 'bed'
        juncs, chromosomes, addFlag, hasNovelJuncs = addOtherJuncs(juncs, type, shortread, args.junction_support, chromosomes,
                                                    False, knownSS, False, False)
        if not addFlag:
            logging.info(f'WARNING: No junctions found in {shortread} that passed filters')
        if not hasNovelJuncs:
            logging.info(f'WARNING: {shortread} did not have any additional junctions that passed filters and were not in {args.gtf}')

    # added to allow annotations not to be used.
    if len(list(juncs.keys())) < 1:
        raise FlairInputDataError("No junctions from GTF or junctionsBed to correct with")

    annotation_files = dict()
    for chrom in chromosomes:
        annotation_files[chrom] = os.path.join(temp_dir, "%s_known_juncs.bed" % chrom)
        with open(os.path.join(temp_dir, "%s_known_juncs.bed" % chrom), "w") as out:
            if chrom in juncs:
                data = juncs[chrom]
                sortedData = sorted(list(data.keys()), key=lambda item: item[0])
                for k in sortedData:
                    annotation = data[k]
                    c1, c2, strand = k
                    print(chrom, c1, c2, annotation, ".", strand, sep="\t", file=out)
    return chromosomes, annotation_files



def correct_single_read(bed_read, intervalTree, junctionBoundaryDict):
    juncs = bed_read.juncs
    strand = bed_read.strand
    c1Type, c2Type = ("donor", "acceptor") if strand == "+" else ("acceptor", "donor")

    newJuncs = list()
    ssStrands = set()

    for x in juncs:
        c1, c2 = x[0], x[1]
        if c1 not in junctionBoundaryDict:
            junctionBoundaryDict = ssCorrect(c1, strand, c1Type, intervalTree, junctionBoundaryDict, False)
        if c2 not in junctionBoundaryDict:
            junctionBoundaryDict = ssCorrect(c2, strand, c2Type, intervalTree, junctionBoundaryDict, False)

        c1Corr = junctionBoundaryDict[c1].ssCorr.coord
        c2Corr = junctionBoundaryDict[c2].ssCorr.coord
        # don't allow junctions outside or near the ends of the reads
        ends_slop = 8
        if not ((bed_read.start + ends_slop) <= c1Corr < (bed_read.end - ends_slop)):
            return None
        if not ((bed_read.start + ends_slop) <= c2Corr < (bed_read.end - ends_slop)):
            return None

        ssTypes = [junctionBoundaryDict[c1].ssCorr.ssType, junctionBoundaryDict[c2].ssCorr.ssType]

        ssStrands.add(junctionBoundaryDict[c1].ssCorr.strand)
        ssStrands.add(junctionBoundaryDict[c2].ssCorr.strand)

        if None in ssTypes:  # or ssTypes[0] == ssTypes[1]: # Either two donors or two acceptors or both none.
            return None
        newJuncs.append((c1Corr, c2Corr))

    starts, sizes = get_exons_from_juncs(newJuncs, bed_read.start, bed_read.end)
    # 0 length exons, remove them.
    if min(sizes) == 0:
        return None

    else:
        bed_read.juncs = newJuncs
        bed_read.exon_sizes = sizes
        bed_read.exon_starts = starts
        bed_read.set_exons()
        return bed_read


def get_rgb(name, strand, junclen):
    if name[:4] == 'ENST':
        return '3,28,252'
    elif junclen == 0:
        return "99,99,99"
    elif strand == '+':
        return "27,158,119"
    else:
        return "217,95,2"


def get_exons_from_juncs(juncs, start, end):
    if len(juncs) == 0:
        exon_starts = [0]
        exon_sizes = [end - start]
    else:
        exon_starts = [0] + [x[1] - start for x in juncs]
        exon_sizes = [juncs[0][0] - start] + [juncs[i + 1][0] - juncs[i][1] for i in range(len(juncs) - 1)] + [
            end - juncs[-1][1]]
    return exon_starts, exon_sizes


class BedRead(object):
    def __init__(self):
        self.name = None

    def generate_from_cigar(self, align_start, is_reverse, cigar_tuples, read_name, reference_chrom, map_qual_score,
                            junc_direction):
        ref_pos = align_start
        intron_blocks = []
        has_match = False
        for block in cigar_tuples:
            if block[0] == 3 and has_match:  # intron, pay attention
                intron_blocks.append([ref_pos, ref_pos + block[1]])
                ref_pos += block[1]
            elif block[0] in {0, 7, 8, 2}:  # consumes reference
                ref_pos += block[1]
                if block[0] in {0, 7, 8}:
                    has_match = True  # match
        # dirtowrite = '-' if is_reverse else '+'
        # chr1  476363  497259  ENST00000455464.7_ENSG00000237094.12    1000    -
        # 476363  497259  0       3       582,169,151,    0,8676,20745,
        exon_sizes, exon_starts = intron_chain_to_exon_starts(intron_blocks, align_start, ref_pos)
        if junc_direction not in {'+', '-'}:
            junc_direction = "-" if is_reverse else "+"

        junctions = []
        for i in range(len(exon_starts) - 1):
            junctions.append((align_start + exon_starts[i] + exon_sizes[i], align_start + exon_starts[i + 1]))

        self.ref_chrom = reference_chrom
        self.start = align_start
        self.end = ref_pos
        self.name = read_name
        self.score = map_qual_score
        self.strand = junc_direction
        # self.blockcount = len(intron_blocks) + 1
        self.exon_sizes = exon_sizes
        self.exon_starts = exon_starts
        self.juncs = tuple(junctions)
        self.set_exons()

    def set_exons(self):
        self.exons = [(self.start + self.exon_starts[i], self.start + self.exon_starts[i] + self.exon_sizes[i]) for i in
                      range(len(self.exon_starts))]

    def reset_from_exons(self, exons):
        self.exons = exons
        self.start = exons[0][0]
        self.end = exons[-1][1]
        self.juncs = tuple([(exons[x][1], exons[x + 1][0]) for x in range(len(exons) - 1)])
        self.exon_sizes = [x[1] - x[0] for x in exons]
        self.exon_starts = [x[0] - self.start for x in exons]

    def get_sequence(self, genome):
        exons = self.exons
        if self.strand == '-':
            exons = exons[::-1]
        exon_seq = []
        for i in range(len(exons)):
            this_exon_seq = genome.fetch(self.ref_chrom, exons[i][0], exons[i][1])
            if self.strand == '-':
                this_exon_seq = get_reverse_complement(this_exon_seq)
            exon_seq.append(this_exon_seq)
        return ''.join(exon_seq)

    def generate_from_vals(self, chrom, start, end, name, score, strand, juncs):
        self.ref_chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.juncs = juncs
        self.exon_starts, self.exon_sizes = get_exons_from_juncs(juncs, start, end)
        self.set_exons()

    def get_bed_line(self):
        rgb_color = get_rgb(self.name, self.strand, len(self.juncs))
        bed_line = [self.ref_chrom, self.start, self.end, self.name, self.score, self.strand,
                   self.start, self.end, rgb_color, len(self.exon_starts), ','.join([str(x) for x in self.exon_sizes]),
                   ','.join([str(x) for x in self.exon_starts])]
        bed_line = [str(x) for x in bed_line]
        return bed_line


class AnnotData(object):
    def __init__(self):
        self.transcript_to_exons = {}
        self.transcript_to_info = {}
        self.all_transcripts = []
        self.juncchain_to_transcript = {}
        self.junc_to_gene = {}
        self.all_annot_SE = []
        self.all_spliced_exons = {'+': {}, '-': {}}
        self.gene_to_annot_juncs = {}
        self.gene_to_strand = {}

    def return_data(self):
        return self.juncchain_to_transcript, self.junc_to_gene, self.all_annot_SE, self.all_spliced_exons, self.gene_to_annot_juncs, \
            self.gene_to_strand, self.transcript_to_exons, self.all_transcripts

def generate_region_dict(all_regions):
    chrom_to_regions, regions_to_annot_data = {}, {}
    for chrom, region_start, region_end in all_regions:
        if chrom not in chrom_to_regions:
            chrom_to_regions[chrom] = []
        chrom_to_regions[chrom].append((region_start, region_end))
        regions_to_annot_data[(chrom, region_start, region_end)] = AnnotData()
    return chrom_to_regions, regions_to_annot_data

def get_t_name_to_exons(gtf):
    chrom_to_transcript_to_exons = {}
    for line in open(gtf):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
            if ty == 'transcript' or ty == 'exon':
                if chrom not in chrom_to_transcript_to_exons:
                    chrom_to_transcript_to_exons[chrom] = {}
                this_transcript = line[8][line[8].find('transcript_id') + 15:]
                this_transcript = this_transcript[:this_transcript.find('"')]
                this_gene = line[8].split('gene_id "')[1].split('"')[0].replace('_', '-')
                if (this_transcript, this_gene) not in chrom_to_transcript_to_exons[chrom]:
                    chrom_to_transcript_to_exons[chrom][(this_transcript, this_gene)] = [(None, None, strand), []]

                if ty == 'transcript':
                    chrom_to_transcript_to_exons[chrom][(this_transcript, this_gene)][0] = (start, end, strand)
                elif ty == 'exon':
                    chrom_to_transcript_to_exons[chrom][(this_transcript, this_gene)][1].append((start, end))
    return chrom_to_transcript_to_exons

def get_annot_t_ends(tinfo):
    t_start, t_end, strand = tinfo[0]
    if t_start is None:
        t_start = min([x[0] for x in tinfo[1]])
        t_end = max([x[1] for x in tinfo[1]])
    return t_start, t_end, strand

def save_transcript_annot_to_region(transcript, gene, region, regions_to_annot_data, t_start, t_end, strand, t_exons):
    sorted_exons = sorted(t_exons)
    regions_to_annot_data[region].transcript_to_exons[(transcript, gene)] = tuple(sorted_exons)
    juncs = []
    for i in range(len(sorted_exons) - 1):
        juncs.append((sorted_exons[i][1], sorted_exons[i + 1][0]))
    regions_to_annot_data[region].all_transcripts.append((transcript, gene, strand))
    if gene not in regions_to_annot_data[region].gene_to_strand:
        regions_to_annot_data[region].gene_to_strand[gene] = strand
    if len(juncs) == 0:
        regions_to_annot_data[region].all_annot_SE.append((t_start, t_end, strand, gene))
    else:
        if gene not in regions_to_annot_data[region].all_spliced_exons[strand]:
            regions_to_annot_data[region].all_spliced_exons[strand][gene] = set()
        regions_to_annot_data[region].all_spliced_exons[strand][gene].update(set(sorted_exons))
        regions_to_annot_data[region].juncchain_to_transcript[tuple(juncs)] = (transcript, gene)
        if gene not in regions_to_annot_data[region].gene_to_annot_juncs:
            regions_to_annot_data[region].gene_to_annot_juncs[gene] = set()
        for j in juncs:
            if j not in regions_to_annot_data[region].junc_to_gene:
                regions_to_annot_data[region].junc_to_gene[j] = set()
            regions_to_annot_data[region].junc_to_gene[j].add((transcript, gene))
            regions_to_annot_data[region].gene_to_annot_juncs[gene].add(j)
    regions_to_annot_data[region].all_annot_SE = sorted(regions_to_annot_data[region].all_annot_SE)
    return regions_to_annot_data

def get_annot_for_chrom(chromregions, region_chrom, regions_to_annot_data, chrom_transcript_to_exons):
    for transcript, gene in chrom_transcript_to_exons:
        tinfo = chrom_transcript_to_exons[(transcript, gene)]
        t_start, t_end, strand = get_annot_t_ends(tinfo)
        for region_start, region_end in chromregions:
            if region_start < t_start < region_end or region_start < t_end < region_end:
                region = (region_chrom, region_start, region_end)
                regions_to_annot_data = save_transcript_annot_to_region(transcript, gene, region, regions_to_annot_data,
                                                                     t_start, t_end, strand, tinfo[1])
    return regions_to_annot_data


def get_annot_info(gtf, all_regions):
    chrom_to_regions, regions_to_annot_data = generate_region_dict(all_regions)
    chrom_to_transcript_to_exons = get_t_name_to_exons(gtf)

    for region_chrom in chrom_to_transcript_to_exons:
        if region_chrom in chrom_to_regions:  # only get annot for regions that exist in reads
            regions_to_annot_data = get_annot_for_chrom(chrom_to_regions[region_chrom], region_chrom, regions_to_annot_data,
                                                     chrom_to_transcript_to_exons[region_chrom])
    return regions_to_annot_data


def get_filter_tome_align_cmd(args, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound):
    # count sam transcripts ; the dash at the end means STDIN
    count_cmd = ['filter_transcriptome_align.py', '--sam', '-',
                 '-o', output_name, '-t', 1,  # feeding 1 thread in because this is already multithreaded here
                 ]

    if clipping_file:
        count_cmd.extend(['--trimmedreads', clipping_file])
    if map_file:
        count_cmd.extend(['--generate_map', map_file])
    if args.end_norm_dist:
        count_cmd.extend(['--output_endpos', output_name.split('.counts.tsv')[0] + '.ends.tsv',
                          '--end_norm_dist', args.end_norm_dist])
    if not args.no_stringent or is_annot:
        count_cmd.extend(['--stringent', '--allow_UTR_indels'])
    if args.output_bam:
        count_cmd.extend(['--output_bam', output_name.split('.counts.tsv')[0] + '.bam'])
    if not args.no_check_splice:
        count_cmd.append('--check_splice')
    if not args.no_check_splice or not args.no_stringent or is_annot:
        count_cmd.extend(['-i', ref_bed])  # annotated isoform bed file
    if args.trust_ends:
        count_cmd.append('--trust_ends')
    if unique_bound and (not args.no_stringent or is_annot):
        count_cmd.extend(['--unique_bound', unique_bound])
    if args.remove_internal_priming:
        count_cmd.extend(['--remove_internal_priming',
                          '--intprimingthreshold', str(args.intprimingthreshold),
                          '--intprimingfracAs', str(args.intprimingfracAs),
                          '--transcriptomefasta', args.transcriptfasta])
    if args.remove_internal_priming and is_annot:
        count_cmd.append('--permissive_last_exons')
    return tuple(count_cmd)


def transcriptome_align_and_count(args, input_reads, align_ref_fasta, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound):
    # minimap (results are piped into count_sam_transcripts.py)
    # '--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    if isinstance(input_reads, str):
        input_reads = [input_reads]
    mm2_cmd = tuple(
        ['minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD'] + args.mm2_args + [align_ref_fasta] + input_reads)
    # FIXME add in step to filter out chimeric reads here
    # FIXME really need to go in and check on how count_sam_transcripts is working
    count_cmd = get_filter_tome_align_cmd(args, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound)
    pipettor.run([mm2_cmd, count_cmd])

# START METHODS TO EDIT - HARRISON

def get_best_ends(curr_group, end_window):
    best_ends = []
    if len(curr_group) > int(end_window):
        all_starts = Counter([x[0] for x in curr_group])
        all_ends = Counter([x[1] for x in curr_group])
        for start1, end1, strand1, name1 in curr_group:
            weighted_score = all_starts[start1] + all_ends[end1]
            best_ends.append((weighted_score, start1, end1, strand1, name1))
    else:
        # take most common non-ambiguous strand for group
        groupStrands = Counter([x[2] for x in curr_group]).most_common()
        myStrand = 'ambig'
        for i in range(len(groupStrands)):
            if groupStrands[i][0] != 'ambig':
                myStrand = groupStrands[i][0]
                break
        
        for start1, end1, strand1, name1 in curr_group:
            score, weighted_score = 0, 0
            for start2, end2, strand2, name2 in curr_group:
                if abs(start1 - start2) <= end_window and abs(end1 - end2) <= end_window:
                    score += 2
                    weighted_score += ((end_window - abs(start1 - start2)) / end_window) + \
                                     ((end_window - abs(end1 - end2)) / end_window)
            best_ends.append((weighted_score, start1, end1, myStrand, name1))
    best_ends.sort(reverse=True)
    # DO I WANT TO ADD CORRECTION TO NEARBY ANNOTATED TSS/TTS????
    return best_ends[0]


def combine_final_ends(curr_group):
    if len(curr_group) == 1:
        return curr_group[0]
    else:
        curr_group.sort(key=lambda x: x[2])  # sort by marker
        all_reads = [y for x in curr_group for y in x[-1]]
        if curr_group[0] != 'a':  # if no annotated iso, sort further
            curr_group.sort(key=lambda x: len(x[-1]), reverse=True)
        best_iso = curr_group[0]
        best_iso[-1] = all_reads
        return best_iso


def group_reads_by_ends(read_info_list, sort_index, end_window):
    sorted_ends = sorted(read_info_list, key=lambda x: x[sort_index])
    new_groups, group = [], []
    last_edge = 0
    for iso_info in sorted_ends:
        edge = iso_info[sort_index]
        if edge - last_edge <= end_window:
            group.append(iso_info)
        else:
            if len(group) > 0:
                new_groups.append(group)
            group = [iso_info]
        last_edge = edge
    if len(group) > 0:
        new_groups.append(group)
    return new_groups

# MAIN METHOD - CALLS OTHERS IN GROUP
# read_ends is a list containing elements with: (read.start, read.end, read.strand, read.name)
# If the reads are spliced, the group will contain only the info for reads with a shared splice junction
# if the reads are unspliced, the group will contain info for all unspliced reads in a given chromosome/region,
# depending on how you're parallelizing
# The output is a list containing elements with:
# [weighted.end.score, group.start, group.end, group.strand, representative.read.name,
# [list of all read names in group]]
# weighted.end.score (represents how many reads have ends similar to this exact position)
# You can rewrite this to redo the grouping method, but please maintain the inputs and outputs.
def collapse_end_groups(end_window, read_ends, do_get_best_ends=True):
    start_groups = group_reads_by_ends(read_ends, 0, end_window)
    all_end_groups, iso_end_groups = [], []
    for start_group in start_groups:
        all_end_groups.extend(group_reads_by_ends(start_group, 1, end_window))
    for end_group in all_end_groups:
        if do_get_best_ends:
            iso_end_groups.append(list(get_best_ends(end_group, end_window)) + [[x[3] for x in end_group]])
        else:
            iso_end_groups.append(combine_final_ends(end_group))
    return iso_end_groups

# END METHODS TO EDIT


def filter_spliced_iso(filter_type, support, juncs, exons, name, score, junc_to_gene,
                       annot_transcript_to_exons, firstpass_junc_to_name, firstpass_unfiltered,  
                       sup_annot_transcript_to_juncs, strand):#, check_term_exons):
    isos_with_similar_juncs = set()
    for j in juncs:
        if firstpass_junc_to_name:
            isos_with_similar_juncs.update(firstpass_junc_to_name[j])
        if j in junc_to_gene:
            isos_with_similar_juncs.update(junc_to_gene[j])  # annot isos
    terminal_exon_is_subset = [0, 0]  # first exon is a subset, last exon is a subset
    first_exon, last_exon = exons[0], exons[-1]
    exon_lengths = [x[1]-x[0] for x in exons]
    superset_support = []
    unique_seq_bound = []
    for otheriso_name in isos_with_similar_juncs:
        if otheriso_name != name:
            if isinstance(otheriso_name, tuple):  # annotated isoform
                if sup_annot_transcript_to_juncs:  # using only supported annotated
                    if otheriso_name in sup_annot_transcript_to_juncs:
                        otheriso_score, otheriso_juncs = sup_annot_transcript_to_juncs[otheriso_name]
                        otheriso_exons = annot_transcript_to_exons[otheriso_name]
                    else:
                        continue
                elif otheriso_name in annot_transcript_to_exons:  # using all annotated
                    otheriso_exons = annot_transcript_to_exons[otheriso_name]
                    otheriso_juncs = [(otheriso_exons[i][1], otheriso_exons[i+1][0]) for i in range(len(otheriso_exons)-1)]  # using all novel
                    otheriso_score = 0
                else:
                    continue
            else:  # firstpass isoform
                otheriso = firstpass_unfiltered[otheriso_name]
                otheriso_score, otheriso_juncs, otheriso_exons = otheriso.score, otheriso.juncs, otheriso.exons

            if len(juncs) < len(otheriso_juncs):
                iso_juncs_str, otheriso_juncs_str = str(juncs)[1:-1].rstrip(','), str(otheriso_juncs)[1:-1]
                if iso_juncs_str in otheriso_juncs_str:  # if junctions are subset
                    # if not check_term_exons:
                    #     return False
                    # else:
                    # check whether first + last exon overlap
                    for i in range(len(otheriso_exons)):
                        other_exon = otheriso_exons[i]
                        if i == 0 or i == len(otheriso_exons) - 1:  # is first or last exon of other transcript - just care that it shares the same terminal ss
                            if first_exon[1] == other_exon[1] or last_exon[0] == other_exon[0]:
                                if first_exon[1] == other_exon[1]:
                                    terminal_exon_is_subset[0] = 1
                                elif last_exon[0] == other_exon[0]:
                                    terminal_exon_is_subset[1] = 1
                                superset_support.append(otheriso_score)
                        else: ## is internal exon of other transcript - see if it contains additional sequence
                            if first_exon[1] == other_exon[1]:
                                unique_seq_bound.append((0, first_exon[1]-other_exon[0]))#other_exon[0]-first_exon[0]))
                                if first_exon[0] >= other_exon[0] - 20:
                                    terminal_exon_is_subset[0] = 1
                                    superset_support.append(otheriso_score)
                            if last_exon[0] == other_exon[0]:
                                unique_seq_bound.append((1, other_exon[1]-last_exon[0]))#sum(exon_lengths) - last_exon[1]-other_exon[1]))
                                if last_exon[1] <= other_exon[1] + 20:
                                    terminal_exon_is_subset[1] = 1
                                    superset_support.append(otheriso_score)
    ##unique_seq is pegged at distance from first/last splice junction
    unique_seq_bound = list(set(unique_seq_bound))
    if strand == '-':
        for i in range(len(unique_seq_bound)):
            # unique_seq_bound[i] = (abs(unique_seq_bound[i][0]-1), sum(exon_lengths) - unique_seq_bound[i][1])
            unique_seq_bound[i] = f'{abs(unique_seq_bound[i][0]-1)}_{unique_seq_bound[i][1]}' ##just invert the indexes
    else:
        for i in range(len(unique_seq_bound)):
            unique_seq_bound[i] = f'{unique_seq_bound[i][0]}_{unique_seq_bound[i][1]}'

    # if not check_term_exons:
    #     return True
    # else:
    if sum(terminal_exon_is_subset) < 2:  # both first and last exon have to overlap
        return True, unique_seq_bound
    elif filter_type != 'nosubset':
        if score >= support and score > max(superset_support) * 1.2:
            return True, unique_seq_bound
    return False, None




def generate_transcriptome_reference(temp_prefix, all_transcripts, annot_transcript_to_exons, chrom, genome, junc_to_gene,
                                     normalize_ends=False, add_length_at_ends=0):
    transcript_to_strand = {}
    transcript_to_new_exons = {}
    with open(temp_prefix + '.annotated_transcripts.bed', 'w') as annot_bed, \
         open(temp_prefix + '.annotated_transcripts.fa', 'w') as annot_fa, \
         open(temp_prefix + '.annotated_transcripts_uniquebound.txt', 'w') as annot_uniqueseq:
        gene_to_terminal_junction_specific_ends = {}
        if normalize_ends:
            for transcript, gene, strand in all_transcripts:
                exons = annot_transcript_to_exons[(transcript, gene)]
                if (gene, strand) not in gene_to_terminal_junction_specific_ends:
                    gene_to_terminal_junction_specific_ends[(gene, strand)] = {'left': {}, 'right': {}}
                if len(exons) > 1:  # don't normalize ends for single exon transcripts
                    if exons[0][1] not in gene_to_terminal_junction_specific_ends[(gene, strand)]['left']:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]] = exons[0][0]
                    else:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]] = \
                            min((gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]], exons[0][0]))
                    if exons[-1][0] not in gene_to_terminal_junction_specific_ends[(gene, strand)]['right']:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]] = exons[-1][1]
                    else:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]] = \
                            max((gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]], exons[-1][1]))

        for transcript, gene, strand in all_transcripts:
            transcript_to_strand[(transcript, gene)] = strand
            exons = list(annot_transcript_to_exons[(transcript, gene)])
            juncs = [(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)]
            # if not check_junction_subset(juncs, junc_to_gene, annot_transcript_to_exons):  #NEW removing subset isos from annot
            # if filter_spliced_iso('nosubset', 0, juncs, exons, (transcript, gene), 0, junc_to_gene, annot_transcript_to_exons, None, None, None, False):
            # if filter_spliced_iso('nosubset', 0, juncs, exons, (transcript, gene), 0, junc_to_gene,
            #                     annot_transcript_to_exons, None, None, None, True):
            is_not_subset, unique_seq = filter_spliced_iso('nosubset', 0, juncs, exons, (transcript, gene), 0, junc_to_gene,
                                annot_transcript_to_exons, None, None, None, strand)
            if is_not_subset:
                if normalize_ends and len(exons) > 1:  # don't normalize ends for single exon transcripts
                    exons[0] = (gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]] - add_length_at_ends,
                                exons[0][1])
                    exons[-1] = (exons[-1][0],
                                 gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]] + add_length_at_ends)
                    transcript_to_new_exons[(transcript, gene)] = tuple(exons)
                exons = tuple(exons)
                start, end = exons[0][0], exons[-1][1]

                exon_starts = [x[0] - start for x in exons]
                exon_sizes = [x[1] - x[0] for x in exons]
                bed_line = [chrom, start, end, transcript + '_' + gene, '.', strand, start, end, '0', len(exons),
                           ','.join([str(x) for x in exon_sizes]), ','.join([str(x) for x in exon_starts])]
                if strand == '-':
                    exons = exons[::-1]
                exon_seq = []
                for i in range(len(exons)):
                    this_exon_seq = genome.fetch(chrom, exons[i][0], exons[i][1])
                    if strand == '-':
                        this_exon_seq = get_reverse_complement(this_exon_seq)
                    exon_seq.append(this_exon_seq)
                annot_bed.write('\t'.join([str(x) for x in bed_line]) + '\n')
                annot_fa.write('>' + transcript + '_' + gene + '\n')
                annot_fa.write(''.join(exon_seq) + '\n')
                if len(unique_seq) > 0:
                    annot_uniqueseq.write(transcript + '_' + gene + '\t' + ','.join(unique_seq) + '\n')
    return transcript_to_strand, transcript_to_new_exons


def identify_good_match_to_annot(args, temp_prefix, chrom, annot_transcript_to_exons, all_transcripts, genome, junc_to_gene):
    good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs = [], set(), {}
    if not args.no_align_to_annot and len(all_transcripts) > 0:
        logging.info('generating transcriptome reference')
        if args.end_norm_dist:
            transcript_to_strand, transcript_to_new_exons = \
                generate_transcriptome_reference(temp_prefix,
                                                 all_transcripts,
                                                 annot_transcript_to_exons,
                                                 chrom,
                                                 genome,
                                                 junc_to_gene,
                                                 normalize_ends=True,
                                                 add_length_at_ends=int(args.end_norm_dist))
        else:
            transcript_to_strand, transcript_to_new_exons = generate_transcriptome_reference(temp_prefix, all_transcripts,
                                                                                        annot_transcript_to_exons,
                                                                                        chrom, genome, junc_to_gene)
        clipping_file = temp_prefix + '.reads.genomicclipping.txt'  # if args.trimmedreads else None
        logging.info('aligning to transcriptome reference')
        transcriptome_align_and_count(args, temp_prefix + '.reads.fasta',
                                   temp_prefix + '.annotated_transcripts.fa',
                                   temp_prefix + '.annotated_transcripts.bed',
                                   temp_prefix + '.matchannot.counts.tsv',
                                   temp_prefix + '.matchannot.read.map.txt', True, clipping_file, temp_prefix + '.annotated_transcripts_uniquebound.txt')
        logging.info('processing good matches')
        with open(temp_prefix + '.matchannot.bed', 'w') as annot_bed:
            for line in open(temp_prefix + '.matchannot.read.map.txt'):
                striso, reads = line.rstrip().split('\t', 1)
                reads = reads.split(',')
                if len(reads) >= args.sjc_support:
                    good_align_to_annot.extend(reads)
                    transcript = '_'.join(striso.split('_')[:-1])
                    gene = striso.split('_')[-1]
                    if (transcript, gene) in annot_transcript_to_exons:
                        if (transcript, gene) in transcript_to_new_exons:
                            exons = transcript_to_new_exons[(transcript, gene)]
                        else:
                            exons = annot_transcript_to_exons[(transcript, gene)]
                        start, end = exons[0][0], exons[-1][1]
                        exon_starts = [x[0] - start for x in exons]
                        exon_sizes = [x[1] - x[0] for x in exons]
                        strand = transcript_to_strand[(transcript, gene)]
                        bed_line = [chrom, start, end, transcript + '_' + gene, len(reads), strand, start, end, '0',
                                   len(exons), ','.join([str(x) for x in exon_sizes]), ','.join([str(x) for x in exon_starts])]
                        annot_bed.write('\t'.join([str(x) for x in bed_line]) + '\n')
                        firstpass_SE.update(set(exons))
                        annot_juncs = tuple([(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)])
                        sup_annot_transcript_to_juncs[(transcript, gene)] = (len(reads), annot_juncs)
    else:
        with open(temp_prefix + '.matchannot.counts.tsv', 'w') as _, \
                open(temp_prefix + '.matchannot.read.map.txt', 'w') as _, \
                open(temp_prefix + '.matchannot.bed', 'w') as _:
            pass
        if args.output_endpos:
            with open(temp_prefix + '.ends.tsv', 'w') as _:
                pass
    good_align_to_annot = set(good_align_to_annot)
    return good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs


def filter_correct_group_reads(args, temp_prefix, region_chrom, region_start, region_end, bam_file, good_align_to_annot, intervalTree,
                            junctionBoundaryDict, generate_fasta=True, sj_to_ends=None,
                            return_used_reads=False, allow_secondary=False):
    if not sj_to_ends:
        sj_to_ends = {}
    if generate_fasta:
        out_fasta = open(temp_prefix + 'reads.notannotmatch.fasta', 'w')
    c = 0
    used_reads = set()
    for read in bam_file.fetch(region_chrom, int(region_start), int(region_end)):
        if (not read.is_secondary or allow_secondary) and (not read.is_supplementary or args.keep_sup):
            if read.reference_name == region_chrom \
                    and int(region_start) <= read.reference_start \
                    and read.reference_end <= int(region_end):
                if read.query_name not in good_align_to_annot:
                    c += 1
                    if generate_fasta:
                        out_fasta.write('>' + read.query_name + '\n')
                        out_fasta.write(read.get_forward_sequence() + '\n')
                    if read.mapping_quality >= args.quality:  # TODO: test this more rigorously
                        used_reads.add(read.query_name)
                        bed_read = BedRead()
                        read_strand = '-' if read.is_reverse else '+'
                        bed_read.generate_from_cigar(read.reference_start, read.is_reverse, read.cigartuples,
                                                    read.query_name,
                                                    read.reference_name, read.mapping_quality, read_strand)
                        if len(bed_read.juncs) > 0:
                            new_strand = inferMM2JuncStrand(read)
                            if new_strand != 'ambig':
                                bed_read.strand = new_strand
                        corrected_read = correct_single_read(bed_read, intervalTree, junctionBoundaryDict)
                        if corrected_read:
                            junc_key = tuple(sorted(corrected_read.juncs))
                            if junc_key not in sj_to_ends:
                                sj_to_ends[junc_key] = []
                            sj_to_ends[junc_key].append((corrected_read.start, corrected_read.end,
                                                      corrected_read.strand, corrected_read.name))
    if generate_fasta:
        out_fasta.close()
    if return_used_reads:
        return sj_to_ends, used_reads
    else:
        return sj_to_ends


# def group_firstpass_single_exon(good_ends_with_sup_reads):
#     last_end, furthergroups, thisgroup = 0, [], []
#     good_ends_with_sup_reads.sort(key=lambda x: [x[1], x[2]])
#     for weighted_score, start, end, strand, name, endsupport in good_ends_with_sup_reads:
#         if (start >= last_end or (last_end - start) / (end - start) > 0.5) and len(thisgroup) > 0:
#             furthergroups.append(thisgroup)
#             thisgroup = []
#         thisgroup.append([weighted_score, start, end, strand, name, endsupport])
#         if end > last_end:
#             last_end = end
#     if len(thisgroup) > 0:
#         furthergroups.append(thisgroup)
#     logging.info(f'single exon groups: {len(furthergroups)}')
#     return furthergroups


def filter_ends_by_subset_and_support(args, good_ends_with_sup_reads):
    best_ends = []
    for i in range(len(good_ends_with_sup_reads)):
        good_ends_with_sup_reads[i][-1] = len(good_ends_with_sup_reads[i][-1])  # last val is read support
    good_ends_with_sup_reads.sort(key=lambda x: [x[0], x[2] - x[1]],
                              reverse=True)  # first by weighted score, then by length
    junc_support = sum([x[-1] for x in good_ends_with_sup_reads])
    if junc_support >= args.sjc_support:
        if args.no_redundant == 'none':
            if good_ends_with_sup_reads[0][-1] < args.sjc_support:
                good_ends_with_sup_reads = [good_ends_with_sup_reads[0]]
                good_ends_with_sup_reads[0][-1] = junc_support
            else:
                good_ends_with_sup_reads = [x for x in good_ends_with_sup_reads if x[-1] >= args.sjc_support]
                good_ends_with_sup_reads = good_ends_with_sup_reads[:args.max_ends]  # select only top most supported ends
            for these_ends in good_ends_with_sup_reads:
                best_ends.append(these_ends)
        else:
            # best_only uses the default sorting
            if args.no_redundant == 'longest':
                good_ends_with_sup_reads.sort(reverse=True, key=lambda x: x[2] - x[1])
            this_best = good_ends_with_sup_reads[0]
            this_best[-1] = junc_support  # all reads for junc are counted as support
            best_ends.append(this_best)
    return best_ends


def process_juncs_to_firstpass_isos(args, temp_prefix, chrom, sj_to_ends, firstpass_SE):
    firstpass_unfiltered, firstpass_junc_to_name = {}, {}
    with open(temp_prefix + '.firstpass.unfiltered.bed', 'w') as iso_out, \
            open(temp_prefix + '.firstpass.reallyunfiltered.bed', 'w') as iso_out_unfilt:
        for juncs in sj_to_ends:
            good_ends_with_sup_reads = collapse_end_groups(args.end_window, sj_to_ends[juncs])
            for best_score, best_start, best_end, best_strand, best_name, sup_reads in good_ends_with_sup_reads:
                this_score = len(sup_reads)
                this_iso = BedRead()
                this_iso.generate_from_vals(chrom, best_start, best_end, best_name, this_score, best_strand, juncs)
                iso_out_unfilt.write('\t'.join(this_iso.get_bed_line()) + '\n')
            if juncs == ():
                best_ends = [x[:-1] + [len(x[-1])] for x in good_ends_with_sup_reads if len(x[-1]) >= args.sjc_support]
            else:
                best_ends = filter_ends_by_subset_and_support(args, good_ends_with_sup_reads)
            for best_score, best_start, best_end, best_strand, best_name, this_score in best_ends:
                this_iso = BedRead()
                this_iso.generate_from_vals(chrom, best_start, best_end, best_name, this_score, best_strand, juncs)
                firstpass_unfiltered[best_name] = this_iso
                iso_out.write('\t'.join(this_iso.get_bed_line()) + '\n')
                if juncs == ():
                    firstpass_SE.add((this_iso.exons[0][0], this_iso.exons[0][1], this_iso.name))
                else:
                    for j in juncs:
                        if j not in firstpass_junc_to_name:
                            firstpass_junc_to_name[j] = set()
                        firstpass_junc_to_name[j].add(best_name)
                    firstpass_SE.update(this_iso.exons)

    firstpass_SE = sorted(list(firstpass_SE))
    return firstpass_unfiltered, firstpass_junc_to_name, firstpass_SE



def filter_single_exon_iso(args, grouped_iso, curr_group, firstpass_unfiltered):
    this_iso = firstpass_unfiltered[grouped_iso[2]]
    expression_comp_with_superset = []
    is_contained = False
    for comp_iso in curr_group:
        if comp_iso != grouped_iso:
            if comp_iso[0] - 10 <= grouped_iso[0] and grouped_iso[1] <= comp_iso[1] + 10:
                if len(comp_iso) == 2 or args.filter == 'nosubset':  # is exon from spliced transcript
                    is_contained = True
                    break  # filter out
                else:  # is other single exon - check relative expression
                    other_score = firstpass_unfiltered[comp_iso[2]].score
                    this_score = this_iso.score
                    if this_score >= args.sjc_support and other_score * 1.2 < this_score:
                        expression_comp_with_superset.append(True)
                    else:
                        expression_comp_with_superset.append(False)
    if not is_contained and all(expression_comp_with_superset):
        return True
    else:
        return False


def filter_single_exon_group(args, curr_group, firstpass_unfiltered, firstpass):
    for grouped_iso in curr_group:
        if len(grouped_iso) == 3:  # is single exon with name
            if filter_single_exon_iso(args, grouped_iso, curr_group, firstpass_unfiltered):
                firstpass[grouped_iso[2]] = firstpass_unfiltered[grouped_iso[2]]
    return firstpass


def filter_all_single_exon(args, firstpass_SE, firstpass_unfiltered, firstpass):
    # group_start = 0
    last_end = 0
    curr_group = []
    for iso_info in firstpass_SE:
        start, end = iso_info[0], iso_info[1]
        if start < last_end:
            curr_group.append(iso_info)
        else:
            if len(curr_group) > 0:
                firstpass = filter_single_exon_group(args, curr_group, firstpass_unfiltered, firstpass)
            curr_group = [iso_info]
            # group_start = start
        if end > last_end:
            last_end = end
    if len(curr_group) > 0:
        firstpass = filter_single_exon_group(args, curr_group, firstpass_unfiltered, firstpass)

    return firstpass


def filter_firstpass_isos(args, firstpass_unfiltered, firstpass_junc_to_name, firstpass_SE, junc_to_gene,
                        sup_annot_transcript_to_juncs, annot_transcript_to_exons):
    iso_to_unique_bound = {}
    if args.filter == 'ginormous':
        firstpass = firstpass_unfiltered
    else:
        firstpass = {}
        for iso_name in firstpass_unfiltered:
            this_iso = firstpass_unfiltered[iso_name]
            if this_iso.juncs != ():
                if args.filter == 'comprehensive':
                    firstpass[iso_name] = this_iso
                else:
                    # passesfiltering = filter_spliced_iso(args, this_iso, firstpass_junc_to_name, firstpass_unfiltered,
                    #                                    junc_to_gene, sup_annot_transcript_to_juncs, annot_transcript_to_exons)
                    # if passesfiltering:
                    is_not_subset, unique_seq = filter_spliced_iso(args.filter, args.sjc_support, this_iso.juncs, this_iso.exons, 
                                          this_iso.name, this_iso.score, junc_to_gene, annot_transcript_to_exons, 
                                          firstpass_junc_to_name, firstpass_unfiltered,  
                                          sup_annot_transcript_to_juncs, this_iso.strand)#, True):
                    if is_not_subset:
                        firstpass[iso_name] = this_iso
                        if len(unique_seq) > 0:
                            iso_to_unique_bound[iso_name] = ','.join(unique_seq)
        # HANDLE SINGLE EXONS SEPARATELY - group first - one traversal of list
        firstpass = filter_all_single_exon(args, firstpass_SE, firstpass_unfiltered, firstpass)

    return firstpass, iso_to_unique_bound


def combine_temp_files_by_suffix(args, temp_prefixes, suffixes):
    for filesuffix in suffixes:
        with open(args.output + filesuffix, 'wb') as combined_file:
            for temp_prefix in temp_prefixes:
                with open(temp_prefix + filesuffix, 'rb') as fd:
                    shutil.copyfileobj(fd, combined_file, 1024 * 1024 * 10)


def get_genes_with_shared_juncs(juncs, junc_to_gene, gene_to_annot_juncs):
    gene_hits = {}
    if juncs != ():
        for j in juncs:
            if j in junc_to_gene:
                for transcript, gene in junc_to_gene[j]:
                    if gene not in gene_hits:
                        gene_hits[gene] = [0, -1 * len(gene_to_annot_juncs[gene])]
                    gene_hits[gene][0] += 1
    return gene_hits


def get_single_exon_gene_overlaps(this_iso, all_annot_SE):
    gene_hits = {}
    this_exon = this_iso.exons[0]
    index = binary_search(this_exon, all_annot_SE)
    for annot_exon_info in all_annot_SE[index - 2:index + 2]:  # start, end, strand, gene
        overlap = min(this_exon[1], annot_exon_info[1]) - max(this_exon[0], annot_exon_info[0])
        if overlap > 0:
            # base coverage of long-read isoform by the annotated isoform
            frac_of_this_iso = float(overlap) / (this_exon[1] - this_exon[0])
            # base coverage of the annotated isoform by the long-read isoform
            frac_of_annot = float(overlap) / (annot_exon_info[1] - annot_exon_info[0])
            if frac_of_this_iso > 0.5 and frac_of_annot > 0.8:
                if annot_exon_info[3] not in gene_hits or frac_of_this_iso > gene_hits[annot_exon_info[3]][0]:
                    gene_hits[annot_exon_info[3]] = [frac_of_this_iso, frac_of_annot]
    return gene_hits


def get_spliced_exon_overlaps(mystrand, myexons, all_spliced_exons, gene_hits):
    for annot_gene in all_spliced_exons[mystrand]:
        annot_exons = sorted(list(all_spliced_exons[mystrand][annot_gene]))
        # check if there is overlap in the genes
        if min((annot_exons[-1][1], myexons[-1][1])) \
                > max((annot_exons[0][0], myexons[0][0])):
            covered_pos = set()
            for s, e in myexons:
                for ast, ae in annot_exons:
                    for p in range(max((ast, s)), min((ae, e))):
                        covered_pos.add(p)
            if len(covered_pos) > sum([x[1] - x[0] for x in myexons]) * 0.5:
                gene_hits.append([len(covered_pos), annot_gene, mystrand])
    return gene_hits


def get_gene_names_and_write_firstpass(temp_prefix, chrom, firstpass, juncchain_to_transcript, junc_to_gene, all_annot_SE,
                                  gene_to_annot_juncs, gene_to_strand, genome, all_spliced_exons,
                                  normalize_ends=False, add_length_at_ends=0, unique_bound=None):
    with open(temp_prefix + '.firstpass.bed', 'w') as iso_out, open(temp_prefix + '.firstpass.fa', 'w') as seq_out, \
         open(temp_prefix + '.firstpass.uniquebound.txt', 'w') as unique_out:
        # THIS IS WHERE WE CAN GET GENES AND ADJUST NAMES
        annot_name_to_used_counts = {}
        gene_to_terminal_junction_specific_ends = {}
        iso_to_info = {}
        novel_gene_isos_to_group = {'+': [], '-': []}
        for iso_name in firstpass:
            this_iso = firstpass[iso_name]
            # Adjust name based on annotation
            this_transcript, this_gene = this_iso.name, None
            if this_iso.juncs != () and this_iso.juncs in juncchain_to_transcript:
                this_transcript, this_gene = juncchain_to_transcript[this_iso.juncs]
                if this_transcript in annot_name_to_used_counts:
                    annot_name_to_used_counts[this_transcript] += 1
                    this_transcript = this_transcript + '-endvar' + str(annot_name_to_used_counts[this_transcript])
                else:
                    annot_name_to_used_counts[this_transcript] = 1
            else:
                if this_iso.juncs != ():
                    gene_hits = get_genes_with_shared_juncs(this_iso.juncs, junc_to_gene, gene_to_annot_juncs)
                else:
                    gene_hits = get_single_exon_gene_overlaps(this_iso, all_annot_SE)  # single exon
                if gene_hits:
                    sorted_genes = sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)
                    this_gene = sorted_genes[0][0]
                else:
                    # look for exon overlap
                    gene_hits = []
                    if this_iso.strand != 'ambig':
                        gene_hits = get_spliced_exon_overlaps(this_iso.strand, this_iso.exons, all_spliced_exons, gene_hits)
                    else:
                        gene_hits = get_spliced_exon_overlaps('+', this_iso.exons, all_spliced_exons, gene_hits)
                        gene_hits = get_spliced_exon_overlaps('-', this_iso.exons, all_spliced_exons, gene_hits)
                    if len(gene_hits) > 0:
                        gene_hits.sort(reverse=True)
                        this_gene = gene_hits[0][1]
                        if this_iso.strand == 'ambig':
                            this_iso.strand = gene_hits[0][2]
            if this_gene is not None:
                strand = gene_to_strand[this_gene]
            else:
                strand = this_iso.strand
                novel_gene_isos_to_group[strand].append((this_iso.start, this_iso.end, iso_name))
                # this_gene = chrom.replace('_', '-') + ':' + str(round(this_iso.start, -3))
            iso_to_info[iso_name] = [this_gene, this_transcript, strand, this_iso.exons]

        # generating non-gene iso groups
        for strand in novel_gene_isos_to_group:
            transcripts_to_group = sorted(novel_gene_isos_to_group[strand])
            last_end = 0
            group_start = 0
            curr_group = []
            for start, end, t_name in transcripts_to_group:
                if start < last_end:
                    curr_group.append((start, end, t_name))
                else:
                    if len(curr_group) > 0:
                        group_name = f'{chrom}:{group_start}-{last_end}:{strand}'
                        for s, e, t in curr_group:
                            iso_to_info[t][0] = group_name
                    curr_group = [(start, end, t_name)]
                    group_start = start
                if end > last_end:
                    last_end = end
            if len(curr_group) > 0:
                group_name = f'{chrom}:{group_start}-{last_end}:{strand}'
                for s, e, t in curr_group:
                    iso_to_info[t][0] = group_name

        # generating standardized set of ends for gene
        if normalize_ends:
            for iso_name in iso_to_info:
                gene, this_transcript, strand, exons = iso_to_info[iso_name]
                if (gene, strand) not in gene_to_terminal_junction_specific_ends:
                    gene_to_terminal_junction_specific_ends[(gene, strand)] = {'left': {}, 'right': {}}
                if len(exons) > 1:
                    if exons[0][1] not in gene_to_terminal_junction_specific_ends[(gene, strand)]['left']:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]] = exons[0][0]
                    else:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]] = min(
                            (gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]], exons[0][0]))
                    if exons[-1][0] not in gene_to_terminal_junction_specific_ends[(gene, strand)]['right']:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]] = exons[-1][1]
                    else:
                        gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]] = max(
                            (gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]], exons[-1][1]))

        for iso_name in iso_to_info:
            this_iso = firstpass[iso_name]
            gene, this_transcript, strand, exons = iso_to_info[iso_name]
            if normalize_ends and len(this_iso.juncs) > 0:
                exons[0] = (gene_to_terminal_junction_specific_ends[(gene, strand)]['left'][exons[0][1]] - add_length_at_ends,
                            exons[0][1])
                exons[-1] = (exons[-1][0],
                             gene_to_terminal_junction_specific_ends[(gene, strand)]['right'][exons[-1][0]] + add_length_at_ends)
                this_iso.reset_from_exons(exons)
            this_iso.strand = strand
            this_iso.name = this_transcript + '_' + gene
            
            if unique_bound and iso_name in unique_bound:
                unique_out.write(this_iso.name + '\t' + unique_bound[iso_name] + '\n')
            
            iso_out.write('\t'.join(this_iso.get_bed_line()) + '\n')
            seq_out.write('>' + this_iso.name + '\n')
            seq_out.write(this_iso.get_sequence(genome) + '\n')


def decode_name_to_iso_gene(name):
    iso = '_'.join(name.split('_')[:-1])
    gene = name.split('_')[-1]
    return iso, gene


def process_detected_isos(args, map_file, bed_file, marker, gene_to_juncs_to_ends, ends_file):
    og_iso_to_reads = {}
    for line in open(map_file):
        name, reads = line.rstrip().split('\t', 1)
        reads = reads.split(',')
        # support = len(reads)
        iso, gene = decode_name_to_iso_gene(name)
        iso_id = (marker, iso)
        og_iso_to_reads[iso_id] = reads

    iso_to_ends = {}
    if args.end_norm_dist and ends_file:
        for line in open(ends_file):
            read_name, transcript, start, end = line.rstrip().split('\t')
            start, end = int(start), int(end)
            if transcript not in iso_to_ends:
                iso_to_ends[transcript] = []
            iso_to_ends[transcript].append((start, end, None, None))
        for iso in iso_to_ends:
            # (weighted_score, start1, end1, strand1, name1)
            new_ends = get_best_ends(iso_to_ends[iso], args.end_window)[1:3]
            iso_to_ends[iso] = new_ends

    for line in open(bed_file):
        line = line.rstrip().split('\t')
        chrom, start, end, name, score, strand = line[:6]
        iso, gene = decode_name_to_iso_gene(name)
        iso_id = (marker, iso)
        start, end = int(start), int(end)
        exon_starts = [int(x) for x in line[11].rstrip(',').split(',')]
        exon_sizes = [int(x) for x in line[10].rstrip(',').split(',')]
        if iso_id in og_iso_to_reads and ((len(og_iso_to_reads[iso_id]) >= args.se_support and len(exon_sizes) == 1) or (len(og_iso_to_reads[iso_id]) >= args.sjc_support and len(exon_sizes) > 1)):
            juncs = []
            for i in range(len(exon_starts) - 1):
                juncs.append((start + exon_starts[i] + exon_sizes[i], start + exon_starts[i + 1]))
            junc_key = (chrom, strand, tuple(juncs))

            if args.end_norm_dist:
                if ends_file:
                    if name in iso_to_ends:
                        start, end = iso_to_ends[name]
                elif len(juncs) > 0:
                    start += int(args.end_norm_dist)
                    end -= int(args.end_norm_dist)

            if gene not in gene_to_juncs_to_ends:
                gene_to_juncs_to_ends[gene] = {}
            if junc_key not in gene_to_juncs_to_ends[gene]:
                gene_to_juncs_to_ends[gene][junc_key] = []
            gene_to_juncs_to_ends[gene][junc_key].append([start, end, iso_id, og_iso_to_reads[iso_id]])
    return gene_to_juncs_to_ends


def get_reverse_complement(seq):
    compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    seq = seq.upper()
    new_seq = []
    for base in seq:
        new_seq.append(compbase[base])
    return ''.join(new_seq[::-1])



def get_bed_gtf_out_from_info(end_info, chrom, strand, juncs, gene, genome):
    start, end, iso_id, read_names = end_info
    score = min(len(read_names), 1000)
    marker, iso = iso_id
    exon_starts, exon_sizes = get_exons_from_juncs(juncs, start, end)
    bed_line = [chrom, start, end, iso + '_' + gene, score, strand, start, end,
               get_rgb(iso, strand, juncs), len(exon_starts), ','.join([str(x) for x in exon_sizes]),
               ','.join([str(x) for x in exon_starts])]
    gtf_lines = []
    gtf_lines.append([chrom, 'FLAIR', 'transcript', start + 1, end, score, strand, '.',
                     'gene_id "' + gene + '"; transcript_id "' + iso + '";'])
    exons = [(start + exon_starts[i], start + exon_starts[i] + exon_sizes[i]) for i in range(len(exon_starts))]
    if strand == '-':
        exons = exons[::-1]
    exon_seq = []
    for i in range(len(exons)):
        gtf_lines.append([chrom, 'FLAIR', 'exon', exons[i][0] + 1, exons[i][1], score, strand, '.',
                         'gene_id "' + gene + '"; transcript_id "' + iso + '"; exon_number ' + str(i + 1)])
        this_exon_seq = genome.fetch(chrom, exons[i][0], exons[i][1])
        if strand == '-':
            this_exon_seq = get_reverse_complement(this_exon_seq)
        exon_seq.append(this_exon_seq)
    return '\t'.join([str(x) for x in bed_line]) + '\n', gtf_lines, ''.join(exon_seq)


def combine_annot_w_novel_and_write_files(args, gene_to_juncs_to_ends, genome):
    read_to_final_transcript = {}
    with open(args.output + '.isoforms.bed', 'w') as iso_out, \
            open(args.output + '.isoform.read.map.txt', 'w') as map_out, \
            open(args.output + '.isoforms.gtf', 'w') as gtf_out, \
            open(args.output + '.isoforms.fa', 'w') as seq_out, \
            open(args.output + '.isoform.counts.txt', 'w') as counts_out:

        for gene in gene_to_juncs_to_ends:
            for chrom, strand, juncs in gene_to_juncs_to_ends[gene]:
                ends_list = gene_to_juncs_to_ends[gene][(chrom, strand, juncs)]
                ends_list = collapse_end_groups(args.end_window, ends_list, False)
                # FIXME could try accounting for all reads assigned to isoforms - assign them to closest ends
                # not sure how much of an issue this is
                if juncs != ():
                    if args.no_redundant == 'best_only':
                        ends_list.sort(key=lambda x: [len(x[-1]), x[1] - x[0]], reverse=True)
                        ends_list = [ends_list[0]]
                    elif args.no_redundant == 'longest':
                        ends_list.sort(key=lambda x: [x[1] - x[0]], reverse=True)
                        ends_list = [ends_list[0]]
                    else:
                        ends_list.sort(key=lambda x: [len(x[-1]), x[1] - x[0]], reverse=True)
                        ends_list = ends_list[:args.max_ends]
                gene_to_juncs_to_ends[gene][(chrom, strand, juncs)] = ends_list


        for gene in gene_to_juncs_to_ends:
            gene_tot = 0
            for chrom, strand, juncs in gene_to_juncs_to_ends[gene]:
                for iso_info in gene_to_juncs_to_ends[gene][(chrom, strand, juncs)]:
                    gene_tot += len(iso_info[-1])

            gtf_lines, t_starts, t_ends = [], [], []
            for chrom, strand, juncs in gene_to_juncs_to_ends[gene]:
                ends_list = gene_to_juncs_to_ends[gene][(chrom, strand, juncs)]

                name_to_used_counts = {}
                for iso_info in ends_list:
                    if len(iso_info[-1])/gene_tot >= args.frac_support:
                        marker, iso = iso_info[2]
                        iso = iso.split('-endvar')[0]
                        if iso in name_to_used_counts:
                            name_to_used_counts[iso] += 1
                            iso = iso + '-endvar' + str(name_to_used_counts[iso])
                        else:
                            name_to_used_counts[iso] = 1
                        iso_info[2] = (marker, iso)
                        bed_line, gtf_for_transcript, tseq = get_bed_gtf_out_from_info(iso_info, chrom, strand, juncs, gene, genome)
                        iso_out.write(bed_line)
                        t_starts.append(iso_info[0])
                        t_ends.append(iso_info[1])
                        gtf_lines.extend(gtf_for_transcript)
                        map_out.write(iso + '_' + gene + '\t' + ','.join(iso_info[3]) + '\n')
                        for r in iso_info[3]:
                            read_to_final_transcript[r] = (iso + '_' + gene, chrom, strand)
                        counts_out.write(iso + '_' + gene + '\t' + str(len(iso_info[3])) + '\n')
                        seq_out.write('>' + iso + '_' + gene + '\n')
                        seq_out.write(tseq + '\n')
            gtf_lines.insert(0, [chrom, 'FLAIR', 'gene', min(t_starts) + 1, max(t_ends), '.', gtf_lines[0][6], '.',
                                'gene_id "' + gene + '";'])
            for g in gtf_lines:
                gtf_out.write('\t'.join([str(x) for x in g]) + '\n')
    if args.end_norm_dist:
        with open(args.output + '.read_ends.bed', 'w') as ends_out:
            for suffix in ['.matchannot.ends.tsv', '.novelisos.ends.tsv']:
                transcript_to_reads = {}
                for line in open(args.output + suffix):
                    read_name, transcript, start, end = line.rstrip().split('\t')
                    if transcript not in transcript_to_reads:
                        transcript_to_reads[transcript] = []
                    transcript_to_reads[transcript].append((read_name, start, end))
                for t in transcript_to_reads:
                    if len(transcript_to_reads[t]) >= args.sjc_support:  # FIXME this needs to be adjusted to consider single exons vs junction chains, also frac_support
                        for r, start, end in transcript_to_reads[t]:
                            if r in read_to_final_transcript:
                                t_name, chrom, strand = read_to_final_transcript[r]
                                ends_out.write('\t'.join([chrom, start, end, t_name + '|' + r, '.', strand]) + '\n')


def decide_parallel_mode(parallel_option, genome_aligned_bam):
    # parallel_option format already validated in option validation method
    if parallel_option in {'bychrom', 'byregion'}:
        return parallel_option
    else:
        if parallel_option[-2:] == 'GB':
            threshold = float(parallel_option[5:-2])
        else:
            threshold = float(parallel_option[5:])
        file_size_GB = os.path.getsize(genome_aligned_bam) / 1e+9
        if file_size_GB > threshold:
            return 'byregion'
        else:
            return 'bychrom'

def generate_genomic_clipping_reference(temp_prefix, bam_file, region_chrom, region_start, region_end):
    with open(temp_prefix + '.reads.genomicclipping.txt', 'w') as out:
        for read in bam_file.fetch(region_chrom, int(region_start), int(region_end)):
            if not read.is_secondary and not read.is_supplementary:
                name = read.query_name
                cigar = read.cigartuples
                tot_clipped = 0
                if cigar[0][0] in {4, 5}:
                    tot_clipped += cigar[0][1]
                if cigar[-1][0] in {4, 5}:
                    tot_clipped += cigar[-1][1]
                out.write(name + '\t' + str(tot_clipped) + '\n')


def run_collapse_by_chrom(listofargs):
    args, temp_prefix, splice_site_annot_chrom, juncchain_to_transcript, junc_to_gene, all_annot_SE, all_spliced_exons, \
        gene_to_annot_juncs, gene_to_strand, annot_transcript_to_exons, all_annot_transcripts = listofargs
    # first extract reads for chrom as fasta
    temp_split = temp_prefix.split('/')[-1].split('-')
    region_chrom, region_start, region_end = '-'.join(temp_split[:-2]), temp_split[-2], temp_split[-1]
    pipettor.run([('samtools', 'view', '-h', args.genome_aligned_bam, region_chrom + ':' + region_start + '-' + region_end),
                  ('samtools', 'fasta', '-')],
                 stdout=open(temp_prefix + '.reads.fasta', 'w'))
    # then align reads to transcriptome and run count_sam_transcripts
    genome = pysam.FastaFile(args.genome)
    bam_file = pysam.AlignmentFile(args.genome_aligned_bam, 'rb')

    # if args.trimmedreads:
    logging.info('generating genomic clipping reference')
    generate_genomic_clipping_reference(temp_prefix, bam_file, region_chrom, region_start, region_end)
    logging.info('identifying good match to annot')

    good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs = \
        identify_good_match_to_annot(args, temp_prefix, region_chrom, annot_transcript_to_exons, all_annot_transcripts, genome, junc_to_gene)
    # load splice junctions for chrom
    logging.info('correcting splice junctions')
    
    intervalTree, junctionBoundaryDict = buildIntervalTree(splice_site_annot_chrom, args.ss_window, region_chrom, False)

    sj_to_ends = filter_correct_group_reads(args, temp_prefix, region_chrom, region_start, region_end, bam_file, good_align_to_annot, intervalTree,
                                       junctionBoundaryDict)
    bam_file.close()
    logging.info('generating isoforms')

    firstpass_unfiltered, firstpass_junc_to_name, firstpass_SE = process_juncs_to_firstpass_isos(args, temp_prefix,
                                                                                                 region_chrom, sj_to_ends,
                                                                                                 firstpass_SE)
    logging.info('filtering isoforms')

    firstpass, iso_to_unique_bound = filter_firstpass_isos(args, firstpass_unfiltered, firstpass_junc_to_name, firstpass_SE,
                                    junc_to_gene, sup_annot_transcript_to_juncs, annot_transcript_to_exons)
    temp_to_remove = [temp_prefix + '.reads.fasta', temp_prefix + 'reads.notannotmatch.fasta']
    if not args.no_align_to_annot and len(all_annot_transcripts) > 0:
        temp_to_remove.extend([temp_prefix + '.annotated_transcripts.bed', temp_prefix + '.annotated_transcripts.fa'])
    if len(firstpass.keys()) > 0:
        logging.info('getting gene names and writing firstpass')
        
        if args.end_norm_dist:
            get_gene_names_and_write_firstpass(temp_prefix, region_chrom, firstpass, juncchain_to_transcript, junc_to_gene,
                                          all_annot_SE, gene_to_annot_juncs, gene_to_strand, genome, all_spliced_exons,
                                          normalize_ends=True, add_length_at_ends=int(args.end_norm_dist), unique_bound=iso_to_unique_bound)
        else:
            get_gene_names_and_write_firstpass(temp_prefix, region_chrom, firstpass, juncchain_to_transcript, junc_to_gene,
                                          all_annot_SE, gene_to_annot_juncs, gene_to_strand, genome, all_spliced_exons, unique_bound=iso_to_unique_bound)
        clipping_file = temp_prefix + '.reads.genomicclipping.txt'  # if args.trimmedreads else None
        logging.info('identifying good match to firstpass')
        
        transcriptome_align_and_count(args, temp_prefix + 'reads.notannotmatch.fasta',
                                   temp_prefix + '.firstpass.fa',
                                   temp_prefix + '.firstpass.bed',
                                   temp_prefix + '.novelisos.counts.tsv',
                                   temp_prefix + '.novelisos.read.map.txt', False, clipping_file, temp_prefix + '.firstpass.uniquebound.txt')
        temp_to_remove.extend([temp_prefix + '.firstpass.fa'])
    else:
        with open(temp_prefix + '.firstpass.fa', 'w') as _, \
                open(temp_prefix + '.firstpass.bed', 'w') as _, \
                open(temp_prefix + '.novelisos.counts.tsv', 'w') as _, \
                open(temp_prefix + '.novelisos.read.map.txt', 'w') as _:
            pass
        if args.output_endpos:
            with open(temp_prefix + '.ends.tsv', 'w') as _:
                pass
    genome.close()
    if not args.keep_intermediate:
        for f in temp_to_remove:
            os.remove(f)


def run_collapse_from_bam():
    args = get_args()
    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)
    allchrom = genome.references
    logging.info('making temp dir')
    temp_dir = make_correct_temp_dir()
    logging.info('Getting regions')
    all_regions = []
    if decide_parallel_mode(args.parallel_mode, args.genome_aligned_bam) == 'bychrom':
        for chrom in allchrom:
            chromsize = genome.get_reference_length(chrom)
            all_regions.append((chrom, 0, chromsize))
    else:
        pipettor.run(['flair_partition',
                      '--min_partition_items', '1000',
                      '--threads', str(args.threads),
                      '--bam=' + args.genome_aligned_bam,
                      temp_dir + 'regions.bed'])
        for line in open(temp_dir + 'regions.bed'):
            line = line.rstrip().split('\t')
            chrom, start, end = line[0], int(line[1]), int(line[2])
            all_regions.append((chrom, start, end))
    logging.info(f'Number of regions {len(all_regions)}')
    logging.info('Generating splice site database')
    known_chromosomes, annotation_files = generate_known_SS_database(args, temp_dir)

    regions_to_annot_data = {}
    if args.gtf:
        logging.info('Extracting annotation from GTF')
        regions_to_annot_data = get_annot_info(args.gtf, all_regions)
    logging.info('splitting by chunk')
    chunk_cmds = []
    temp_prefixes = []
    for region_chrom, region_start, region_end in all_regions:
        if region_chrom in known_chromosomes:
            juncchain_to_transcript, junc_to_gene, all_annot_SE, all_spliced_exons, gene_to_annot_juncs, gene_to_strand, \
                annot_transcript_to_exons, all_annot_transcripts = {}, {}, [], {}, {}, {}, {}, []
            if args.gtf:
                juncchain_to_transcript, junc_to_gene, all_annot_SE, all_spliced_exons, gene_to_annot_juncs, gene_to_strand, \
                    annot_transcript_to_exons, all_annot_transcripts \
                    = regions_to_annot_data[(region_chrom, region_start, region_end)].return_data()

            splice_site_annot_chrom = annotation_files[region_chrom]
            temp_prefix = temp_dir + '-'.join([region_chrom, str(region_start), str(region_end)])
            chunk_cmds.append([args, temp_prefix, splice_site_annot_chrom, juncchain_to_transcript,
                              junc_to_gene, all_annot_SE, all_spliced_exons, gene_to_annot_juncs, gene_to_strand,
                              annot_transcript_to_exons, all_annot_transcripts])
            temp_prefixes.append(temp_prefix)
    logging.info('running by chunk')
    mp.set_start_method('fork')
    p = mp.Pool(args.threads)
    child_errs = set()
    c = 1
    for i in p.imap(run_collapse_by_chrom, chunk_cmds):
        logging.info(f'\rdone running chunk {c} of {len(chunk_cmds)}')
        child_errs.add(i)
        c += 1
    p.close()
    p.join()
    if len(child_errs) > 1:
        raise ValueError(child_errs)

    files_to_combine = ['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed',
                      '.novelisos.counts.tsv', '.novelisos.read.map.txt']
    if not args.no_align_to_annot:
        files_to_combine.extend(['.matchannot.counts.tsv', '.matchannot.read.map.txt', '.matchannot.bed'])
        if args.end_norm_dist:
            files_to_combine.append('.matchannot.ends.tsv')
    if args.end_norm_dist:
        files_to_combine.append('.novelisos.ends.tsv')
    combine_temp_files_by_suffix(args, temp_prefixes, files_to_combine)

    if not args.keep_intermediate:
        shutil.rmtree(temp_dir)

    if not args.no_align_to_annot:
        gene_to_juncs_to_ends = process_detected_isos(args,
                                                args.output + '.matchannot.read.map.txt',
                                                args.output + '.matchannot.bed',
                                                'a',
                                                {},
                                                args.output + '.matchannot.ends.tsv')
    else:
        gene_to_juncs_to_ends = {}
    gene_to_juncs_to_ends = process_detected_isos(args,
                                            args.output + '.novelisos.read.map.txt',
                                            args.output + '.firstpass.bed',
                                            'n',
                                            gene_to_juncs_to_ends,
                                            args.output + '.novelisos.ends.tsv')

    combine_annot_w_novel_and_write_files(args, gene_to_juncs_to_ends, genome)
    if not args.keep_intermediate:
        files_to_remove = ['.firstpass.reallyunfiltered.bed',
                           '.firstpass.unfiltered.bed',
                           '.firstpass.bed',
                           '.novelisos.counts.tsv',
                           '.novelisos.read.map.txt']
        if not args.no_align_to_annot:
            files_to_remove += ['.matchannot.bed',
                                '.matchannot.counts.tsv',
                                '.matchannot.read.map.txt']
            if args.end_norm_dist:
                files_to_remove.append('.matchannot.ends.tsv')
        if args.end_norm_dist:
            files_to_remove.append('.novelisos.ends.tsv')
        for f in files_to_remove:
            os.remove(args.output + f)

    if args.predict_cds:
        prodcmd = ('predictProductivity',
                   '-i', args.output + '.isoforms.bed',
                   '-o', args.output + '.isoforms.CDS',
                   '--gtf', args.gtf,
                   '--genome_fasta', args.genome,
                   '--longestORF')
        pipettor.run([prodcmd])
        # os.rename(args.output + '.isoforms.CDS.bed', args.output + '.isoforms.bed')
        os.remove(args.output + '.isoforms.CDS.info.tsv')
    genome.close()


if __name__ == "__main__":
    run_collapse_from_bam()
