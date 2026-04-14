#!/usr/bin/env python3
import argparse
import os
import pipettor
import shutil
import pysam
import logging
from statistics import median
from collections import Counter
from flair import FlairError
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.junction_correct import junction_corrector_factory
from flair.partition_runner import parallel_mode_parse, partition_runner_factory
from flair.bed_to_gtf import bed_to_gtf
from flair.isoform_data import (Junc, Exon, ReadRec, Gene, Isoform, exons_to_juncs,
                                get_bed_exons_from_exons, get_sequence_for_exons, binary_search,
                                convert_to_bed)
from flair.read_processing import (should_process_read, add_corrected_read_to_groups,
                                   generate_genomic_alignment_read_to_clipping_file)
from flair.count_sam_transcripts import TRUST_ENDS_WINDOW
from flair.annotation_data import annot_data_from_gtf
from flair.pycbio.hgdata.bed import Bed

MIN_POLYA_FRAC_DIFF_FOR_SE_STRANDING = 0.1

# FIXME: add object for all file names
# FIXME: use real TSVs
# FIXME: need to document all the files
# FIXME: it seems overkill to discard a single read based on one unsupported junction.
#        These can be recovered from other reads. Simple way is to add all reads valid
#        junctions, although maybe faster to track with introns are support and figure
#        this out in correct.  Another might have a polymorphic different with the reference
#        that changes splice junctions.  Should this be discarded if multiple long-reads
#        support it, but it isn't annotated.  Maybe these can be identified.

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
    parser.add_argument('-f', '--gtf', dest="annot_gtf", default=None,
                        help='GTF annotation file, used for identifying annotated isoforms')

    mutexc = parser.add_mutually_exclusive_group(required=False)
    mutexc.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    mutexc.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from long-read alignment with intron-prospector)')
    parser.add_argument('--junction_support', type=int, default=1,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support.')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')
    parser.add_argument('-w', '--end_window', type=int, default=100,
                        help='window size for comparing TSS/TES (100)')

    parser.add_argument('--sjc_support', type=int, default=1,
                        help='''minimum number of supporting reads for a spliced isoform''')
    parser.add_argument('--se_support', type=int, default=3,
                        help='''minimum number of supporting reads for a single exon isoform''')
    parser.add_argument('--frac_support', type=float, default=0.05,
                        help='''minimum fraction of gene locus support for isoform to be called
                        default: 0.05, only isoforms that make up more than 5 percent of the gene
                        locus are reported. Set to 0 for max recall''')
    parser.add_argument('--no_stringent', default=False, action='store_true',
                        help='''specify if all supporting reads don't need to be full-length
                        (aligned to first and last exons of transcript).  Use this for fragmented libraries,
                        with an understanding that it will impact precision.''')
    parser.add_argument('--no_check_splice', default=False, action='store_true',
                        help='''don't enforce accurate alignment around splice site.
                        Specify this for libraries with high error rates, but it will reduce precision''')
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
    parser.add_argument('--end_norm_dist', type=int,
                        help='specify the number of basepairs to extend transcript ends if you want to '
                             'normalize them across transcripts in a gene and extend them')
    parser.add_argument('--output_endpos', default=False, action='store_true',
                        help='specify if you want to output a separate file with corrected read end positions. '
                             'For development purposes')
    parser.add_argument('--output_bam', default=False, action='store_true',
                        help='output intermediate bams aligned to the transcriptome. '
                             'Only works with --keep_intermediate, for debugging')
    parser.add_argument('--fusion_breakpoints',
                        help='''for fusion detection only - bed file containing locations of fusion breakpoints on the synthetic genome''')
    parser.add_argument('--allow_paralogs', default=False, action='store_true',
                        help='specify if want to allow reads to be assigned to multiple paralogs with equivalent alignment')
    parser.add_argument('--generate_map', default=False, action='store_true',
                        help='''specify this argument to generate a txt file of read-isoform assignments''')
    parser.add_argument('--trust_strand', default=False, action='store_true',
                        help='''specify if you want FLAIR to trust the stranding of the input reads and not attempt strand correction''')
    args = parser.parse_args()
    args.parallel_mode = parallel_mode_parse(parser, args.parallel_mode)
    args.trust_ends = False
    args.remove_internal_priming = False

    if not os.path.exists(args.genome_aligned_bam):
        parser.error(f'Aligned reads file path does not exist: {args.genome_aligned_bam}')
    if not os.path.exists(args.genome):
        parser.error(f'Genome file path does not exist: {args.genome}')
    return args


####
# basic types
####
# (Junc, Exon, exons_to_juncs, ISO_SRC_ANNOT, ISO_SRC_NOVEL, IsoIdSrc imported from flair.isoform_data)

# tolerance for terminal exon boundary comparisons
TERMINAL_EXON_BOUNDARY_TOLERANCE = 20

# margin for single-exon isoform overlap comparisons
SINGLE_EXON_OVERLAP_MARGIN = 10

# expression ratio threshold for filtering overlapping single-exon isoforms
SINGLE_EXON_EXPRESSION_RATIO = 1.2

# overlap fraction thresholds for gene assignment
MIN_ISOFORM_OVERLAP_FRAC = 0.5
MIN_ANNOT_OVERLAP_FRAC = 0

# search window for binary search of single-exon annotations
ANNOT_SE_SEARCH_WINDOW = 2


####
# misc
###
def make_temp_dir(out_prefix):
    # FIXME: use TMPDIR unless directory explicitly specified
    temp_dir = out_prefix + ".intermediate"
    try:
        os.makedirs(temp_dir, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Creation of the directory `{temp_dir}' failed") from exc
    return temp_dir + '/'


####
# transcriptome alignment
####
def get_filter_tome_align_cmd(args, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound):  # noqa: C901 - FIXME: reduce complexity
    # FIXME: convert filter_transcriptome_align.py to a library, however
    # minimap output needs to be piped through filter_transcriptome_align
    # without saving the bam file.

    # count sam transcripts ; the dash at the end means STDIN
    # use 1 thread in because this is already multithreaded here
    # count_cmd = ['filter_transcriptome_align.py', '--sam', '-',
    # FIXME why use filter_transcriptome_align if not multithreading? Remove?
    count_cmd = ['count_sam_transcripts.py', '--sam', '-',
                 '-o', output_name,]
    if clipping_file:
        count_cmd.extend(['--trimmedreads', clipping_file])
    if map_file:
        count_cmd.extend(['--generate_map', map_file])
    if args.output_endpos or is_annot:
        count_cmd.extend(['--output_endpos', output_name.split('.counts.txt')[0] + '.ends.tsv'])
    if args.end_norm_dist:
        count_cmd.extend(['--end_norm_dist', args.end_norm_dist])
    if not is_annot and not args.no_stringent:
        count_cmd.extend(['--stringent'])
    if is_annot:
        count_cmd.extend(['--allow_UTR_indels'])
    if args.output_bam:
        count_cmd.extend(['--output_bam', output_name.split('.counts.txt')[0] + '.bam'])
    if not args.no_check_splice:
        count_cmd.append('--check_splice')
    if not args.no_check_splice or not args.no_stringent or is_annot or args.fusion_breakpoints:
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
    if args.fusion_breakpoints:
        count_cmd += ['--fusion_breakpoints', args.fusion_breakpoints]
    if args.allow_paralogs:
        count_cmd += ['--allow_paralogs']
    # print(' '.join([str(x) for x in count_cmd]))
    return count_cmd


def transcriptome_align_and_count(args, input_reads, align_ref_fasta, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound):
    # minimap (results are piped into count_sam_transcripts.py)
    # '--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    if isinstance(input_reads, str):
        input_reads = [input_reads]
    mm2_cmd = ['minimap2', '-a', '-N', '4', '--MD'] + [align_ref_fasta] + input_reads

    # FIXME add in step to filter out chimeric reads here
    # FIXME really need to go in and check on how count_sam_transcripts is working
    count_cmd = get_filter_tome_align_cmd(args, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound)

    pipettor.run([mm2_cmd, count_cmd])


##
# Transcript end assignment
##

def get_best_ends(curr_group, end_window):
    best_ends = []
    if len(curr_group) > int(end_window):
        all_starts = Counter([x.start for x in curr_group])
        all_ends = Counter([x.end for x in curr_group])
        for read_info in curr_group:
            weighted_score = all_starts[read_info.start] + all_ends[read_info.end]
            best_ends.append((weighted_score, read_info.start, read_info.end))
    else:
        for read_info1 in curr_group:
            score, weighted_score = 0, 0
            for read_info2 in curr_group:
                if abs(read_info1.start - read_info2.start) <= end_window and abs(read_info1.end - read_info2.end) <= end_window:
                    score += 2
                    weighted_score += (((end_window - abs(read_info1.start - read_info2.start)) / end_window) +
                                       ((end_window - abs(read_info1.end - read_info2.end)) / end_window))
            best_ends.append((weighted_score, read_info1.start, read_info1.end))
    best_ends.sort(reverse=True)
    # FIXME: DO I WANT TO ADD CORRECTION TO NEARBY ANNOTATED TSS/TTS????
    return best_ends[0]

def group_reads_by_ends(read_info_list, sort_index, end_window):
    sorted_ends = sorted(read_info_list, key=lambda x: x.start if sort_index == 0 else x.end)
    new_groups, group = [], []
    last_edge = 0
    for iso_info in sorted_ends:
        edge = iso_info.start if sort_index == 0 else iso_info.end
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

# MAIN METHOD
# read_ends is a list containing elements with: (read.start, read.end, read.strand, read.name)
# If the reads are spliced, the group will contain only the info for reads with a shared splice junction
# if the reads are unspliced, the group will contain info for all unspliced reads in a given chromosome/region,
# The output is a list of Isoform objects containing:
#    - weighted_score (represents how many reads have ends similar to this exact position)
#    - start, end, strand, read_id (representative read id)
#    - supporting_reads (list of all read names in group)

def collapse_end_groups(end_window, isoform):
    start_groups = group_reads_by_ends(isoform.reads, 0, end_window)
    all_end_groups, iso_end_groups = [], []
    for start_group in start_groups:
        all_end_groups.extend(group_reads_by_ends(start_group, 1, end_window))
    for end_group in all_end_groups:
        weighted_score, start, end = get_best_ends(end_group, end_window)
        read_end_info = Isoform.regroup(isoform, start, end, end_group)
        iso_end_groups.append(read_end_info)
    return iso_end_groups


def get_isos_with_similar_juncs(juncs, junc_to_names, junc_to_gene):
    """Find isoforms sharing junctions with the given junction set."""
    novel_isos = set()
    for j in juncs:
        if junc_to_names and j in junc_to_names:
            novel_isos.update(junc_to_names[j])
        # if j in junc_to_gene:
        #     annot_isos.update(junc_to_gene[j])
    return novel_isos

def _is_junction_subset(juncs, otheriso_juncs):
    """Check if juncs is a proper subset of otheriso_juncs using string matching."""
    if len(juncs) >= len(otheriso_juncs):
        return False
    iso_juncs_str = str(juncs)[1:-1].rstrip(',')
    otheriso_juncs_str = str(otheriso_juncs)[1:-1]
    return iso_juncs_str in otheriso_juncs_str


def _check_terminal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                 terminal_exon_is_subset, superset_support):
    """Check overlap with terminal exon of other transcript (first or last).
    Only requires sharing the same terminal splice site."""
    if first_exon.end == other_exon.end:
        terminal_exon_is_subset[0] = 1
        superset_support.append(otheriso_score)
    elif last_exon.start == other_exon.start:
        terminal_exon_is_subset[1] = 1
        superset_support.append(otheriso_score)


def _check_internal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                 terminal_exon_is_subset, superset_support, unique_seq_bound):
    """Check overlap with internal exon of other transcript.
    Records unique sequence boundaries and checks containment within tolerance."""
    if first_exon.end == other_exon.end:
        unique_seq_bound.append((0, first_exon.end - other_exon.start))
        if first_exon.start >= (other_exon.start - TERMINAL_EXON_BOUNDARY_TOLERANCE):
            terminal_exon_is_subset[0] = 1
            superset_support.append(otheriso_score)
    if last_exon.start == other_exon.start:
        unique_seq_bound.append((1, other_exon.end - last_exon.start))
        if last_exon.end <= (other_exon.end + TERMINAL_EXON_BOUNDARY_TOLERANCE):
            terminal_exon_is_subset[1] = 1
            superset_support.append(otheriso_score)


def _check_junction_subset(juncs, first_exon, last_exon, otheriso_score, otheriso_juncs, otheriso_exons,
                           terminal_exon_is_subset, superset_support, unique_seq_bound):
    """Check if juncs is a subset of otheriso_juncs and update tracking lists."""
    if not _is_junction_subset(juncs, otheriso_juncs):
        return
    for i, other_exon in enumerate(otheriso_exons):
        is_terminal = (i == 0 or i == len(otheriso_exons) - 1)
        if is_terminal:
            _check_terminal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                         terminal_exon_is_subset, superset_support)
        else:
            _check_internal_exon_overlap(first_exon, last_exon, other_exon, otheriso_score,
                                         terminal_exon_is_subset, superset_support, unique_seq_bound)


def _check_novel_iso_subset(novel_iso_id, all_isoforms,
                            juncs, first_exon, last_exon, terminal_exon_is_subset,
                            superset_support, unique_seq_bound):
    """Check if query isoform is subset of a novel (candidate) isoform."""
    otheriso = all_isoforms[novel_iso_id]
    _check_junction_subset(juncs, first_exon, last_exon, otheriso.score, otheriso.juncs, otheriso.exons,
                           terminal_exon_is_subset, superset_support, unique_seq_bound)

def filter_spliced_iso(filter_type, support, juncs, exons, name, score, annots,
                       junc_to_names, all_isoforms,
                       sup_annot_transcript_to_juncs, strand, annot_id):
    assert isinstance(exons[0], Exon)  # FIXME: debugging

    novel_isos = get_isos_with_similar_juncs(juncs, junc_to_names, annots.junc_to_gene)
    terminal_exon_is_subset = [0, 0]  # first exon is a subset, last exon is a subset
    first_exon, last_exon = exons[0], exons[-1]
    superset_support = []
    unique_seq_bound = []
    for novel_iso_id in novel_isos:
        if novel_iso_id != name:
            _check_novel_iso_subset(novel_iso_id, all_isoforms,
                                    juncs, first_exon, last_exon, terminal_exon_is_subset,
                                    superset_support, unique_seq_bound)
    unique_seq_bound = list(set(unique_seq_bound))
    if strand == '-':
        # just invert the indexes
        for i in range(len(unique_seq_bound)):
            unique_seq_bound[i] = f'{abs(unique_seq_bound[i][0] - 1)}_{unique_seq_bound[i][1]}'
    else:
        for i in range(len(unique_seq_bound)):
            unique_seq_bound[i] = f'{unique_seq_bound[i][0]}_{unique_seq_bound[i][1]}'

    if sum(terminal_exon_is_subset) < 2:  # both first and last exon have to overlap
        return True, unique_seq_bound
    elif filter_type != 'nosubset':
        if score >= support and score > max(superset_support) * 1.2:
            return True, unique_seq_bound
    return False, None

####
# terminal exon normalization
####
class GeneMaxTerminalExonsEnds:
    """Class to collect the maximal terminal exons ends for a gene.
    Exons are groups based on the location of the internal splice junction"""
    def __init__(self, gene_id):
        self.gene_id = gene_id
        self.left_ends = {}
        self.right_ends = {}

    def add_left_end(self, exon):
        if exon.end not in self.left_ends:
            self.left_ends[exon.end] = exon.start
        else:
            self.left_ends[exon.end] = min(self.left_ends[exon.end], exon.start)

    def add_right_end(self, exon):
        if exon.start not in self.right_ends:
            self.right_ends[exon.start] = exon.end
        else:
            self.right_ends[exon.start] = max(self.right_ends[exon.start], exon.end)

    def get_left_end(self, exon):
        return self.left_ends[exon.end]

    def get_right_end(self, exon):
        return self.right_ends[exon.start]

class MaxTerminalExonsEnds:
    """Collection of maximal terminal exons ends by gene.
    A genes terminal exons are grouped by the interior exon splice junction
    location.
    """
    def __init__(self):
        # FIXME: this is temporary.  The code groups by (gene_id, strand)
        # for reasons that are suspected to be bugs in stranding.  We keep
        # this but generate a warning until we are sure it is fixed
        self._by_gene_id = {}  # (gene_id, strand) -> GeneMaxTerminalExonsEnds
        self._gene_id_to_strand = {}
        self._genes_warned = set()

    def _obtain(self, gene_id, strand):
        "get current entry or create a new one"
        gene_key = (gene_id, strand)
        gene_entry = self._by_gene_id.get(gene_key)
        if gene_entry is None:
            gene_entry = GeneMaxTerminalExonsEnds(gene_id)
            self._by_gene_id[gene_key] = gene_entry

        # FIXME: tmp generate strand warning or error. it should be impossible to get here
        # with gene broken like this.
        existing_strand = self._gene_id_to_strand.get(gene_id)
        if existing_strand is None:
            self._gene_id_to_strand[gene_id] = strand
        elif strand != existing_strand:
            raise FlairError(f"BUG: gene id '{gene_id}' has transcripts on both strands")
        elif (strand != existing_strand) and (gene_id not in self._genes_warned):
            # FIXME: this is disabled for not to get hard failure
            self._genes_warned.add(gene_id)
            logging.warning("BUG: gene id '%s' has transcripts on both strands", gene_id)

        return gene_entry

    def add_transcript(self, gene_id, strand, transcript_id, exons):
        # don't normalize ends for single exon transcripts, but still record gene
        # FIXME: do we actually want to add single-exon genes?
        gene_entry = self._obtain(gene_id, strand)
        if len(exons) > 1:
            gene_entry.add_left_end(exons[0])
            gene_entry.add_right_end(exons[-1])

    def fetch(self, gene_id, strand) -> GeneMaxTerminalExonsEnds:
        """return entry or error"""
        return self._by_gene_id[(gene_id, strand)]

def max_terminal_exons_ends_from_annots(annots):
    max_terminal_exons_ends = MaxTerminalExonsEnds()
    for transcript_id, gene_id, strand in annots.transcripts:
        exons = annots.transcript_to_exons[(transcript_id, gene_id)]
        max_terminal_exons_ends.add_transcript(gene_id, strand, transcript_id, exons)
    return max_terminal_exons_ends

def max_terminal_exons_ends_from_iso_infos(iso_to_info):
    max_terminal_exons_ends = MaxTerminalExonsEnds()
    for iso_name in iso_to_info:
        isoform = iso_to_info[iso_name]
        max_terminal_exons_ends.add_transcript(isoform.gene_id, isoform.strand, isoform.transcript_id, isoform.exons)
    return max_terminal_exons_ends

####
# transcriptome reference
####
def normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons,
                                  *, add_length_at_ends=0):
    "updates terminal exons ends"
    gene_terminal_exons = max_terminal_exons_ends.fetch(gene_id, strand)
    exons[0] = Exon(gene_terminal_exons.get_left_end(exons[0]) - add_length_at_ends,
                    exons[0].end)
    exons[-1] = Exon(exons[-1].start,
                     gene_terminal_exons.get_right_end(exons[-1]) + add_length_at_ends)

def generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, normalize_ends, max_terminal_exons_ends,
                                                add_length_at_ends, transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand[(transcript_id, gene_id)] = strand
    exons = list(annots.transcript_to_exons[(transcript_id, gene_id)])
    assert isinstance(exons[0], Exon)  # FIXME tmp debugging
    juncs = exons_to_juncs(exons)
    is_not_subset, unique_seq = filter_spliced_iso('nosubset', 0, juncs, exons, (transcript_id, gene_id),
                                                   0, annots, None, None, None, strand, transcript_id)
    if is_not_subset:
        if normalize_ends:
            normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons,
                                          add_length_at_ends=add_length_at_ends)
            transcript_to_new_exons[(transcript_id, gene_id)] = tuple(exons)
        exons = tuple(exons)
        start, end = exons[0].start, exons[-1].end

        # FIXME: duplicated code
        exon_starts, exon_sizes = get_bed_exons_from_exons(exons, start)
        # FIXME: duplicated use BED class,
        bed_line = [chrom, start, end, transcript_id + '_' + gene_id, '.', strand, start, end, '0', len(exons),
                    ','.join([str(x) for x in exon_sizes]), ','.join([str(x) for x in exon_starts])]
        trans_seq = get_sequence_for_exons(genome, chrom, strand, exons)
        annot_bed_fh.write('\t'.join([str(x) for x in bed_line]) + '\n')
        annot_fa_fh.write('>' + transcript_id + '_' + gene_id + '\n')
        annot_fa_fh.write(''.join(trans_seq) + '\n')
        if len(unique_seq) > 0:
            annot_uniqueseq_fh.write(transcript_id + '_' + gene_id + '\t' + ','.join(unique_seq) + '\n')

def generate_transcriptome_reference_guts(normalize_ends, annots, add_length_at_ends, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand = {}
    transcript_to_new_exons = {}
    max_terminal_exons_ends = None
    if normalize_ends:
        max_terminal_exons_ends = max_terminal_exons_ends_from_annots(annots)

    for transcript_id, gene_id, strand in annots.transcripts:
        generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, normalize_ends, max_terminal_exons_ends, add_length_at_ends,
                                                    transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)
    return transcript_to_strand, transcript_to_new_exons

def generate_transcriptome_reference(temp_prefix, annots, chrom, genome,
                                     normalize_ends=False, add_length_at_ends=0):
    with (open(temp_prefix + '.annotated_transcripts.bed', 'w') as annot_bed_fh,
          open(temp_prefix + '.annotated_transcripts.fa', 'w') as annot_fa_fh,
          open(temp_prefix + '.annotated_transcripts_uniquebound.txt', 'w') as annot_uniqueseq_fh):
        return generate_transcriptome_reference_guts(normalize_ends, annots, add_length_at_ends, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)


def identify_good_match_to_annot(args, temp_prefix, chrom, annots, genome):
    # FIXME: refactor
    # good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs = [], set(), {}
    read_to_transcript = {}
    if not args.no_align_to_annot and len(annots.transcripts) > 0:
        # logging.info('generating transcriptome reference')
        # this part generates the fasta file for the annotation
        if args.end_norm_dist is not None:
            transcript_to_strand, transcript_to_new_exons = \
                generate_transcriptome_reference(temp_prefix, annots, chrom, genome,
                                                 normalize_ends=True,
                                                 add_length_at_ends=args.end_norm_dist)
        else:
            transcript_to_strand, transcript_to_new_exons = \
                generate_transcriptome_reference(temp_prefix, annots, chrom, genome)
        # FIXME: make a TSV
        clipping_file = temp_prefix + '.reads.genomicclipping.txt'
        transcriptome_align_and_count(args, temp_prefix + '.reads.fasta',
                                      temp_prefix + '.annotated_transcripts.fa',
                                      temp_prefix + '.annotated_transcripts.bed',
                                      temp_prefix + '.matchannot.counts.txt',
                                      temp_prefix + '.matchannot.read.map.txt', True,
                                      clipping_file,
                                      temp_prefix + '.annotated_transcripts_uniquebound.txt')
        for line in open(temp_prefix + '.matchannot.ends.tsv'):
            line = line.rstrip().split('\t')
            read, transcript = line[:2]
            if line[2] != 'None':
                startindex, startdist, endindex, enddist = [int(x) for x in line[2:]]
                read_to_transcript[read] = (transcript, startindex, startdist, endindex, enddist)
    else:
        # create empty output files
        # FIXME: why doesn't this create all of them?
        # FIXME: change above logic so file open all happens in one place
        with (open(temp_prefix + '.matchannot.counts.txt', 'w'),
              open(temp_prefix + '.matchannot.read.map.txt', 'w')):
            pass
        if args.output_endpos:
            with open(temp_prefix + '.ends.tsv', 'w') as _:
                pass
    # good_align_to_annot = set(good_align_to_annot)
    # return good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs
    return read_to_transcript


def _correct_and_group_read(read, read_to_annot_transcript, annots, junction_corrector, sj_to_ends, genome):
    """Correct a single read's splice junctions and add it to sj_to_ends groups.

    Spliced and single-exon reads are fundamentally different:
    - Spliced: junctions corrected from annotation or intron support, strand from correction
    - Single-exon: no correction, strand resolved later in group_se_by_overlap
    """
    readrec = ReadRec.from_read(read, genome=genome)

    # annotated spliced: correct junctions and strand from annotation
    if read.query_name in read_to_annot_transcript:
        # FIXME more id assumptions
        tid, startindex, startdist, endindex, enddist = read_to_annot_transcript[read.query_name]
        transcript = '_'.join(tid.split('_')[:-1])
        gene = tid.split('_')[-1]
        exons = annots.transcript_to_exons[(transcript, gene)]
        annot_juncs = [(exons[x].end, exons[x + 1].start) for x in range(len(exons) - 1)]
        if len(annot_juncs) > 0:
            newstart = annot_juncs[startindex][0] - startdist
            newend = annot_juncs[endindex][1] + enddist
            juncs = tuple([Junc(x[0], x[1]) for x in annot_juncs[startindex:endindex + 1]])
            readrec.correct_from_annotation(newstart, newend, annots.gene_to_strand[gene], juncs)
            add_corrected_read_to_groups(readrec, sj_to_ends)
            return

    # unannotated spliced: correct junctions and strand from intron support
    if readrec.juncs:
        if junction_corrector.correct_readrec(readrec):
            add_corrected_read_to_groups(readrec, sj_to_ends)
        else:
            logging.debug(f"read dropped: junction correction failed: {readrec.name}")
        return

    # single-exon: no correction, strand resolved later in group_se_by_overlap
    add_corrected_read_to_groups(readrec, sj_to_ends)

def filter_correct_group_reads(args, temp_prefix, region, bam_file, read_to_annot_transcript, annots,
                               junction_corrector, genome, *, sj_to_ends=None,
                               allow_secondary=False):
    """Filter reads, correct splice junctions, and group by junction chain."""
    if sj_to_ends is None:
        sj_to_ends = {}
    for read in bam_file.fetch(region.name, region.start, region.end):
        if should_process_read(read, region, args.quality, args.keep_sup, allow_secondary):
            _correct_and_group_read(read, read_to_annot_transcript, annots, junction_corrector, sj_to_ends, genome)
    return sj_to_ends


def filter_ends_allow_multiple(isoforms, sjc_support, max_ends):
    """Allow multiple ends per junction chain.
    Returns list of Isoform objects that meet support threshold."""
    if isoforms[0].num_reads < sjc_support:
        # If top candidate doesn't meet threshold, merge all reads into it
        best = isoforms[0]
        for iso in isoforms[1:]:
            best.reads.extend(iso.reads)
        return [best]
    else:
        # Filter to those meeting support threshold and limit to max_ends
        filtered = [x for x in isoforms if x.num_reads >= sjc_support]
        filtered = filtered[:max_ends]  # select only top most supported ends
        return filtered

def filter_ends_single_best(isoforms, no_redundant_mode):
    """Pick single best end from junction chain.
    Returns list with single Isoform object."""
    # best_only uses the default sorting, doesn't require additional action
    if no_redundant_mode == 'longest':
        isoforms.sort(reverse=True, key=lambda x: x.length)

    # Pick single best end and merge all reads into it
    best = isoforms[0]
    for iso in isoforms[1:]:
        best.reads.extend(iso.reads)
    return [best]

def filter_ends_by_redundant_and_support(isoforms, sjc_support, se_support, no_redundant, max_ends):
    """Sort ends, then select best ones based on support and value of args.no_redundant."""
    # First by weighted score, then by length
    if isoforms[0].juncs == ():
        support = se_support
        isoforms.sort(key=lambda x: [x.num_reads * x.genomic_length], reverse=True)
    else:
        support = sjc_support
        isoforms.sort(key=lambda x: [x.num_reads, x.genomic_length], reverse=True)

    junc_support = sum([x.num_reads for x in isoforms])
    if junc_support < support:
        logging.debug(f"isoform group dropped: insufficient support ({junc_support} < {support}): {isoforms[0].chrom}:{isoforms[0].start}-{isoforms[0].end}")
        return []

    if no_redundant == 'none':
        # Allow multiple ends per junction chain
        return filter_ends_allow_multiple(isoforms, support, max_ends)
    else:
        # Pick single best end
        return filter_ends_single_best(isoforms, no_redundant)


def _write_unfiltered_ends(isoforms, fh):
    for iso_readrec in isoforms:
        convert_to_bed(iso_readrec).write(fh)

class CandidateIsoforms:
    """Candidate isoforms before filtering, with junction and exon indices.

    isoforms: dict of isoform_name -> Isoform
    junc_to_names: dict of Junc -> set of isoform_names sharing that junction
    exons: set of Exon (named exons for SE, all exons for spliced)
    """
    def __init__(self):
        self.isoforms = {}
        self.junc_to_names = {}
        self.exons = set()

    def add(self, isoform):
        self.isoforms[isoform.name] = isoform
        if isoform.juncs == ():
            self.exons.add(Exon(isoform.start, isoform.end, isoform.name))
        else:
            for j in isoform.juncs:
                if j not in self.junc_to_names:
                    self.junc_to_names[j] = set()
                self.junc_to_names[j].add(isoform.name)
            self.exons.update(set(isoform.exons))


def _process_junc_ends(args, isoforms, candidates, iso_fh):
    # this assumes single exons are pre-grouped by overlap
    # previously treated single exons separately due to them being in larger groups
    filtered_isoforms = filter_ends_by_redundant_and_support(isoforms, args.sjc_support, args.se_support, args.no_redundant, args.max_ends)
    for isoform in filtered_isoforms:
        candidates.add(isoform)
        convert_to_bed(isoform).write(iso_fh)

def _process_junc(args, isoform, candidates, iso_fh, iso_unfilt_fh):
    # NOTE: Harrison's TED code will be slotted in here to replace collapse_end_groups
    these_firstpass = collapse_end_groups(args.end_window, isoform)
    _write_unfiltered_ends(these_firstpass, iso_unfilt_fh)
    _process_junc_ends(args, these_firstpass, candidates, iso_fh)


def correct_se_strand_polyA(read_group, se_support):
    # strand correction for single exon genes based on location of polyA tail sequence

    # FIXME: shouldnt this check for poly(T)
    left_polyA = [read.polyA[0] for read in read_group]
    right_polyA = [read.polyA[1] for read in read_group]
    num_reads = len(read_group)

    left_polyA_count = sum([1 for x in left_polyA if x > 0])
    right_polyA_count = sum([1 for x in right_polyA if x > 0])

    left_polyA_frac = left_polyA_count / num_reads
    right_polyA_frac = right_polyA_count / num_reads

    if abs(left_polyA_frac - right_polyA_frac) > MIN_POLYA_FRAC_DIFF_FOR_SE_STRANDING:
        if max((left_polyA_count, right_polyA_count)) >= se_support:
            if right_polyA_count > left_polyA_count:
                return '+'
            else:
                return '-'
    return None

def _group_se_reads_by_overlap(reads):
    """Group single-exon reads into clusters by coordinate overlap."""
    read_groups = []
    last_end = -1
    read_group = None
    for r in sorted(reads, key=lambda x: (x.start, x.end)):
        if r.start >= last_end:
            if read_group is not None:
                read_groups.append(read_group)
            last_end = r.end
            read_group = []
        if r.end > last_end:
            last_end = r.end
        read_group.append(r)
    if read_group is not None:
        read_groups.append(read_group)
    return read_groups

def group_se_by_overlap(chrom, isoform, se_support, trust_strand):
    for read_group in _group_se_reads_by_overlap(isoform.reads):
        if trust_strand:
            # get most common read strand for group
            read_strands = [x.strand for x in read_group]
            new_strand = max(set(read_strands), key=read_strands.count)
        else:
            # correct based on polyA
            new_strand = correct_se_strand_polyA(read_group, se_support)
        # filter out single exon groups that fail stranding
        if new_strand is None:
            logging.debug(f"single-exon group dropped: strand could not be determined ({len(read_group)} reads): {chrom}:{read_group[0].start}-{read_group[-1].end}")
        else:
            new_key = (chrom, median([x.start for x in read_group]), median([x.end for x in read_group]), ())
            yield new_key, new_strand, read_group

class IsoformOverlapGroups:
    """Isoforms grouped by junction chain with overlap-clustered ends.

    Key: (chrom, median_start, median_end, juncs) where juncs is () for single-exon.
    Value: Isoform.

    Strand is not part of the key. For spliced isoforms, strand is determined
    during junction correction and stored on the Isoform. For single-exon
    reads, strand cannot be determined from junctions, so overlapping reads are
    grouped by coordinate overlap first, then strand is resolved per group by
    majority vote (trust_strand) or polyA consensus.  This means opposite-strand
    single-exon reads at the same locus merge into one group; the minority
    strand is discarded.
    """
    def __init__(self):
        self._groups = {}

    def add_spliced(self, chrom, juncs, isoform):
        """Add a spliced isoform, keyed by chrom, median ends, and junction chain."""
        self._groups[(chrom, median(isoform.starts), median(isoform.ends), juncs)] = isoform

    def add_se_overlap_groups(self, chrom, isoform, se_support, trust_strand):
        """Split single-exon reads into overlap groups and resolve strand."""
        for new_key, new_strand, read_group in group_se_by_overlap(chrom, isoform, se_support, trust_strand):
            self._groups[new_key] = Isoform.regroup(isoform, newreads=read_group, newstrand=new_strand)

    def __iter__(self):
        return iter(self._groups)

    def __getitem__(self, key):
        return self._groups[key]

    def items(self):
        return self._groups.items()


def group_by_overlap(sj_to_ends, se_support, trust_strand):
    groups = IsoformOverlapGroups()
    for (chrom, juncs), isoform in sj_to_ends.items():
        if len(juncs) > 0:
            groups.add_spliced(chrom, juncs, isoform)
        else:
            groups.add_se_overlap_groups(chrom, isoform, se_support, trust_strand)
    return groups


def process_juncs_to_firstpass_isos(args, temp_prefix, sj_to_ends, annots, region_chrom):
    sjc_with_overlap_groups = group_by_overlap(sj_to_ends, args.se_support, args.trust_strand)
    # FIXME everything below here requires confidence in transcript strand
    genes, novel_gene_isos_to_group = build_genes(sjc_with_overlap_groups, annots)
    for strand in novel_gene_isos_to_group:
        generate_non_gene_iso_groups_strand(genes, novel_gene_isos_to_group, strand, region_chrom, sjc_with_overlap_groups)

    candidates = CandidateIsoforms()
    with open(temp_prefix + '.firstpass.unfiltered.bed', 'w') as iso_fh, \
            open(temp_prefix + '.firstpass.reallyunfiltered.bed', 'w') as iso_unfilt_fh:
        for juncs, isoform in sjc_with_overlap_groups.items():
            _process_junc(args, isoform, candidates, iso_fh, iso_unfilt_fh)
    return candidates

####
# single-exon transcript processing
####
def filter_single_exon_iso(args, single_exon, curr_group, all_isoforms):
    """Check if a single-exon isoform passes filtering against its overlap group."""
    isoform = all_isoforms[single_exon.name]
    expression_comp_with_superset = []
    is_contained = False
    for exon in curr_group:
        if exon != single_exon:
            if ((exon.start - SINGLE_EXON_OVERLAP_MARGIN) <= single_exon.start and
                    single_exon.end <= (exon.end + SINGLE_EXON_OVERLAP_MARGIN)):
                if exon.name != '' or args.filter == 'nosubset':  # is exon from spliced transcript
                    is_contained = True
                    break  # filter out
                else:  # is other single exon - check relative expression
                    other_score = all_isoforms[exon.name].score
                    if isoform.score >= args.sjc_support and other_score * SINGLE_EXON_EXPRESSION_RATIO < isoform.score:
                        expression_comp_with_superset.append(True)
                    else:
                        expression_comp_with_superset.append(False)
    return not is_contained and all(expression_comp_with_superset)


def filter_single_exon_group(args, curr_group, all_isoforms, firstpass):
    """Filter single-exon isoforms in an overlap group against spliced exons."""
    for exon in curr_group:
        if exon.name != '':  # is single exon with name
            if filter_single_exon_iso(args, exon, curr_group, all_isoforms):
                firstpass[exon.name] = all_isoforms[exon.name]
            else:
                logging.debug(f"single-exon isoform dropped: contained or low expression: {exon.name} ({all_isoforms[exon.name].num_reads} reads)")
    return firstpass


def filter_all_single_exon(args, sorted_exons, all_isoforms, firstpass):
    """Group exons by overlap and filter single-exon isoforms."""
    last_end = 0
    curr_group = []

    for exon in sorted_exons:
        if exon.start < last_end:
            curr_group.append(exon)
        else:
            if len(curr_group) > 0:
                firstpass = filter_single_exon_group(args, curr_group, all_isoforms, firstpass)
            curr_group = [exon]
        if exon.end > last_end:
            last_end = exon.end
    if len(curr_group) > 0:
        firstpass = filter_single_exon_group(args, curr_group, all_isoforms, firstpass)

    return firstpass


def filter_firstpass_isos(args, candidates, annots, sup_annot_transcript_to_juncs):
    """Filter candidate isoforms by subset/support criteria.
    Returns (firstpass dict, iso_to_unique_bound dict)."""
    iso_to_unique_bound = {}

    if args.filter == 'ginormous':
        firstpass = dict(candidates.isoforms)
    else:
        firstpass = {}
        for iso_name, isoform in candidates.isoforms.items():
            if isoform.juncs != ():
                if args.filter == 'comprehensive':
                    firstpass[iso_name] = isoform
                else:
                    assert isinstance(isoform.exons[0], Exon)  # FIXME tmp debugging
                    is_not_subset, unique_seq = filter_spliced_iso(args.filter, args.sjc_support, isoform.juncs, isoform.exons,
                                                                   iso_name, isoform.num_reads, annots,
                                                                   candidates.junc_to_names, candidates.isoforms,
                                                                   sup_annot_transcript_to_juncs, isoform.strand, isoform.transcript_id)
                    if not is_not_subset:
                        logging.debug(f"isoform dropped: subset of another isoform: {iso_name} ({isoform.num_reads} reads)")
                    else:
                        firstpass[iso_name] = isoform
                        if len(unique_seq) > 0:
                            iso_to_unique_bound[iso_name] = ','.join(unique_seq)
        # HANDLE SINGLE EXONS SEPARATELY - group first - one traversal of list
        firstpass = filter_all_single_exon(args, sorted(candidates.exons), candidates.isoforms, firstpass)

    return firstpass, iso_to_unique_bound


def combine_temp_files_by_suffix(output, temp_prefixes, suffixes):
    for filesuffix in suffixes:
        with open(output + filesuffix, 'wb') as combined_fh:
            for temp_prefix in temp_prefixes:
                with open(temp_prefix + filesuffix, 'rb') as in_fh:
                    shutil.copyfileobj(in_fh, combined_fh, 1024 * 1024 * 10)


def get_genes_with_shared_juncs(juncs, annots):
    # FIXME: what does this actually return?
    gene_hits = {}
    if juncs != ():
        for j in juncs:
            if j in annots.junc_to_gene:
                for transcript_id, gene_id in annots.junc_to_gene[j]:
                    if gene_id not in gene_hits:
                        gene_hits[gene_id] = [0, -1 * len(annots.gene_to_annot_juncs[gene_id])]
                    gene_hits[gene_id][0] += 1
    return gene_hits


def get_single_exon_gene_overlaps(strand, iso_readrec, annots):
    gene_hits = {}
    exon = iso_readrec.exons[0]
    index = binary_search(exon, annots.all_annot_SE[strand])
    # FIXME: how does this ever work? all_annot_SE is [(start, end, strand, gene_id), ...]
    for annot_exon_info in annots.all_annot_SE[strand][index - ANNOT_SE_SEARCH_WINDOW:index + ANNOT_SE_SEARCH_WINDOW]:
        # FIXME: make overlap a function
        overlap = min(exon.end, annot_exon_info.end) - max(exon.start, annot_exon_info.start)
        if overlap > 0:
            # base coverage of long-read isoform by the annotated isoform
            frac_of_iso = float(overlap) / (exon.end - exon.start)
            # base coverage of the annotated isoform by the long-read isoform
            frac_of_annot = float(overlap) / (annot_exon_info.end - annot_exon_info.start)
            if frac_of_iso > MIN_ISOFORM_OVERLAP_FRAC and frac_of_annot > MIN_ANNOT_OVERLAP_FRAC:
                if annot_exon_info.name not in gene_hits or frac_of_iso > gene_hits[annot_exon_info.name][0]:
                    gene_hits[annot_exon_info.name] = [frac_of_iso, frac_of_annot]
    return gene_hits

def get_spliced_exon_overlaps(strand, exons, annots):
    gene_hits = []
    for annot_gene in annots.spliced_exons[strand]:
        annot_exons = sorted(list(annots.spliced_exons[strand][annot_gene]))
        # check if there is overlap in the genes
        # FIXME: not clear how this checks for overlap
        if (min((annot_exons[-1].end, exons[-1].end)) > max((annot_exons[0].start, exons[0].start))):
            covered_pos = set()
            for ex in exons:
                for aex in annot_exons:
                    for p in range(max((aex.start, ex.start)), min((aex.end, ex.end))):
                        covered_pos.add(p)
            if len(covered_pos) > sum([x.end - x.start for x in exons]) * 0.5:
                gene_hits.append([len(covered_pos), annot_gene, strand])
    return gene_hits

def _get_transcript_gene_from_annot(iso_readrec, annots, annot_name_to_used_counts):
    """Return (transcript_id, gene_id) if iso matches an annotated junction chain, else (None, None)."""
    if iso_readrec.juncs != () and iso_readrec.juncs in annots.juncchain_to_transcript:
        transcript_id, gene_id = annots.juncchain_to_transcript[iso_readrec.juncs]
        if transcript_id in annot_name_to_used_counts:
            annot_name_to_used_counts[transcript_id] += 1
            transcript_id = transcript_id + '-endvar' + str(annot_name_to_used_counts[transcript_id])
        else:
            annot_name_to_used_counts[transcript_id] = 1
        return transcript_id, gene_id
    else:
        return None, None


def _find_gene_id_by_overlap(iso_readrec, annots):
    """Find gene_id for an isoform without a matching junction chain, using junction or exon overlap."""
    # FIXME this all requires that we already trust the strand of the transcript
    if iso_readrec.juncs != ():
        gene_hits = get_genes_with_shared_juncs(iso_readrec.juncs, annots)
    else:
        gene_hits = get_single_exon_gene_overlaps(iso_readrec.strand, iso_readrec, annots)
    if gene_hits:
        return sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)[0][0]
    else:
        # look for exon overlap
        if iso_readrec.strand != 'ambig':
            gene_hits = get_spliced_exon_overlaps(iso_readrec.strand, iso_readrec.exons, annots)
        else:
            gene_hits = get_spliced_exon_overlaps(iso_readrec.strand, iso_readrec.exons, annots)
        if gene_hits:
            gene_hits.sort(reverse=True)
            if iso_readrec.strand == 'ambig':
                iso_readrec.strand = gene_hits[0][2]
            return gene_hits[0][1]
        else:
            return None


def get_gene_name_firstpass(isoform, annots, annot_name_to_used_counts):
    transcript_id, gene_id = _get_transcript_gene_from_annot(isoform, annots, annot_name_to_used_counts)
    if transcript_id is None:
        gene_id = _find_gene_id_by_overlap(isoform, annots)
    return gene_id, transcript_id


def build_genes(firstpass, annots):
    """Assign gene names to firstpass isoforms, building Gene objects.

    Returns (genes, novel_gene_isos_to_group) where:
    - genes: dict of gene_id -> Gene for isoforms matched to known genes
    - novel_gene_isos_to_group: isoforms needing novel gene assignment
    """
    annot_name_to_used_counts = {}
    genes = {}
    novel_gene_isos_to_group = {'+': [], '-': []}
    for iso_key in firstpass:
        isoform = firstpass[iso_key]
        gene_id, isoform_id = get_gene_name_firstpass(isoform, annots, annot_name_to_used_counts)
        isoform.transcript_id = isoform_id
        if gene_id is not None:
            # removing this strand correction breaks the unusual junction (due to underlying variant?) test
            isoform.strand = annots.gene_to_strand[gene_id]
            if gene_id not in genes:
                genes[gene_id] = Gene(gene_id, isoform.chrom, annots.gene_to_strand[gene_id])
            genes[gene_id].add_isoform(isoform)
        else:
            novel_gene_isos_to_group[isoform.strand].append((isoform.start, isoform.end, iso_key))

    return genes, novel_gene_isos_to_group


def _assign_novel_gene_group(genes, chrom, strand, group_start, last_end, curr_group, firstpass):
    """Create a Gene for a group of novel overlapping isoforms."""
    gene_id = f'{chrom}:{group_start}-{last_end}:{strand}'
    gene = Gene(gene_id, chrom, strand)
    genes[gene_id] = gene
    for s, e, n in curr_group:
        gene.add_isoform(firstpass[n])

def generate_non_gene_iso_groups_strand(genes, novel_gene_isos_to_group, strand, chrom, firstpass):
    """Group novel isoforms by coordinate overlap and create Gene objects."""
    transcripts_to_group = sorted(novel_gene_isos_to_group[strand])
    last_end = 0
    group_start = 0
    curr_group = []
    for start, end, iso_name in transcripts_to_group:
        if start < last_end:
            curr_group.append((start, end, iso_name))
        else:
            if len(curr_group) > 0:
                _assign_novel_gene_group(genes, chrom, strand, group_start, last_end, curr_group, firstpass)
            curr_group = [(start, end, iso_name)]
            group_start = start
        if end > last_end:
            last_end = end
    if len(curr_group) > 0:
        _assign_novel_gene_group(genes, chrom, strand, group_start, last_end, curr_group, firstpass)

def write_first_pass_isoforms(iso_name, normalize_ends, isoform, max_terminal_exons_ends, add_length_at_ends, unique_bound, unique_fh, iso_fh, seq_fh, genome):
    # FIXME: do normalization outside of write function
    if normalize_ends and len(isoform.exons) > 1:  # don't normalize ends for single exon transcripts
        normalize_gene_terminal_exons(max_terminal_exons_ends, isoform.gene_id, isoform.strand, isoform.exons,
                                      add_length_at_ends=add_length_at_ends)
        isoform.reset_from_exons(isoform.exons)
    if isoform.transcript_id is None:
        isoform.transcript_id = isoform.name
    isoform.name = isoform.transcript_id + '_' + isoform.gene_id

    if unique_bound and iso_name in unique_bound:
        unique_fh.write(isoform.name + '\t' + unique_bound[iso_name] + '\n')

    convert_to_bed(isoform).write(iso_fh)
    seq_fh.write('>' + isoform.name + '\n')
    seq_fh.write(isoform.get_sequence(genome) + '\n')

def write_firstpass(temp_prefix, chrom, firstpass, annots, genome, *,
                    normalize_ends=False, add_length_at_ends=0, unique_bound=None):

    # generating standardized set of ends for gene
    if normalize_ends:
        max_terminal_exons_ends = max_terminal_exons_ends_from_iso_infos(firstpass)
    else:
        # FIXME: passing None is move obvious to flow control,
        # although making write_first_pass_isoforms less monolithic
        # it does more than writing
        max_terminal_exons_ends = {}

    with (open(temp_prefix + '.firstpass.bed', 'w') as iso_fh,
          open(temp_prefix + '.firstpass.fa', 'w') as seq_fh,
          open(temp_prefix + '.firstpass.uniquebound.txt', 'w') as unique_fh):
        for iso_name in firstpass:
            write_first_pass_isoforms(iso_name, normalize_ends, firstpass[iso_name], max_terminal_exons_ends,
                                      add_length_at_ends, unique_bound, unique_fh, iso_fh, seq_fh, genome)


####
# results output
####

def get_transcirpts_to_reads(temp_prefix, suffix):
    transcript_to_reads = {}
    for line in open(temp_prefix + suffix):
        read_name, transcript_id, start, end = line.rstrip().split('\t')
        if transcript_id not in transcript_to_reads:
            transcript_to_reads[transcript_id] = []
        transcript_to_reads[transcript_id].append((read_name, start, end))
    return transcript_to_reads

def write_transcript_ends_bed(args, temp_prefix, suffix, read_to_final_transcript, ends_fh):
    transcript_to_reads = get_transcirpts_to_reads(temp_prefix, suffix)
    for t in transcript_to_reads:
        if len(transcript_to_reads[t]) >= args.sjc_support:  # FIXME this needs to be adjusted to consider single exons vs junction chains, also frac_support
            for r, start, end in transcript_to_reads[t]:
                if r in read_to_final_transcript:
                    t_name, chrom, strand = read_to_final_transcript[r]
                    Bed(chrom, int(start), int(end), name=t_name + '|' + r,
                        score=0, strand=strand).write(ends_fh)

def write_transcript_ends_beds(args, temp_prefix, read_to_final_transcript, ends_fh):
    # FIXME: these are not real TSVs
    for suffix in ['.matchannot.ends.tsv']:
        write_transcript_ends_bed(args, temp_prefix, suffix, read_to_final_transcript, ends_fh)

def predict_productivity(out_prefix, genome_fasta, gtf):
    cmd = ('predictProductivity',
           '-i', out_prefix + '.isoforms.bed',
           '-o', out_prefix + '.isoforms.CDS',
           '--gtf', gtf,
           '--genome_fasta', genome_fasta,
           '--longestORF')
    pipettor.run(cmd)

def _iso_passes_support_filter(args, iso, num_exons, iso_to_counts, gene_to_tot):
    if iso not in iso_to_counts:
        return False
    else:
        count = iso_to_counts[iso][0]
        if num_exons > 1:
            return (count >= args.sjc_support) and (count / gene_to_tot[iso.split('_')[-1]][0]) >= args.frac_support
        else:
            return (count >= args.se_support) and (count / gene_to_tot[iso.split('_')[-1]][1]) >= args.frac_support


def _run_region(*, partition, gtf_data, junction_corrector, args):  # noqa: C901 - FIXME: reduce complexity
    region = partition.region
    # FIXME confusing name if taking gtf_data,
    # FIXME: should only have region, so why take region arg
    annots = annot_data_from_gtf(gtf_data, region)

    # first extract reads for region as fasta
    pipettor.run([('samtools', 'view', '-h', args.genome_aligned_bam, region.name + ':' + str(region.start) + '-' + str(region.end)),
                  ('samtools', 'fasta', '-')],
                 stdout=partition.output_path('reads.fasta'))

    # then align reads to transcriptome and run count_sam_transcripts
    genome = pysam.FastaFile(args.genome)
    bam_file = pysam.AlignmentFile(args.genome_aligned_bam, 'rb')

    # genomic clipping: amount of clipping (from cigar) at ends of reads when aligned to genome
    # generates file with [read{\t}clipping amount] on each line
    # For comparing with amount of clipping after alignment to transcriptome
    # in order to check whether transcriptome alignment is comparable to or better than genomic alignment,
    # which can be considered to support isoform.
    # used in filter_transcriptome_align

    # logging.info('generating genomic clipping reference')
    num_reads, clipping_file = generate_genomic_alignment_read_to_clipping_file(partition.file_prefix, bam_file, region)

    if num_reads == 0:
        suffixes = ['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed',
                    '.isoforms.bed',
                    '.isoform.read.map.txt', '.isoforms.gtf', '.isoforms.fa', '.isoform.counts.txt']
        if args.end_norm_dist:
            suffixes.extend(['.read_ends.bed', '.matchannot.ends.tsv'])
        if not args.no_align_to_annot:
            suffixes.extend(['.matchannot.counts.txt', '.matchannot.read.map.txt'])
        for s in suffixes:
            out = open(partition.file_prefix + s, 'w')
            out.close()
        return

    # aligning to reference transcriptome, then identifying reads that match well to reference transcripts
    # with filter_transcriptome_align
    # logging.info('identifying good match to annot')
    if not args.no_align_to_annot:
        logging.info('aligning to transcriptome reference')
    read_to_annot_transcript = \
        identify_good_match_to_annot(args, partition.file_prefix, region.name, annots, genome)

    logging.info('correcting and grouping reads, filtering isoforms')

    # takes in bam file, for each read attempts to correct splice junctions (removes unsupported ones), then groups reads by junction chains
    # this also handles read strandedness if necessary
    sj_to_ends = filter_correct_group_reads(args, partition.file_prefix, region, bam_file, read_to_annot_transcript, annots,
                                            junction_corrector, genome)
    bam_file.close()

    # for each junction chain, clusters ends - generates junction chain x ends
    # firstpass objects then does initial filtering by read support and
    # redundant ends also separates single exon isoforms from spliced isoforms
    # (because they're handled differently in future step for identifying
    # annotated gene/isoform names)
    candidates = process_juncs_to_firstpass_isos(args, partition.file_prefix, sj_to_ends, annots, region.name)

    # filter isoforms: remove subsets, generate unique boundary sequences
    firstpass, iso_to_unique_bound = filter_firstpass_isos(args, candidates, annots, {})

    if len(firstpass.keys()) > 0:
        # logging.info('getting gene names and writing firstpass')
        # this section identifies annotated gene and isoform names (primarily based on splice junction matching, secondarily by exon overlap)
        # also adjusts isoform strand, determines novel isoform and gene names
        # also normalizes transcript ends (temporarily extends ends so that transcript end alignment does not drive transcript assignment during transcriptome alignment)
        # writes out bed and fa files
        logging.info('realigning to firstpass and getting final isoforms')
        if args.end_norm_dist is not None:
            write_firstpass(partition.file_prefix, region.name, firstpass, annots, genome,
                            normalize_ends=True, add_length_at_ends=args.end_norm_dist, unique_bound=iso_to_unique_bound)
        else:
            write_firstpass(partition.file_prefix, region.name, firstpass, annots, genome, unique_bound=iso_to_unique_bound)
            # logging.info('identifying good match to firstpass')

        # aligns to firstpass transcriptome, identifies best read -> isoform alignment for each read, then gets read counts per isoform
        transcriptome_align_and_count(args, partition.output_path('reads.fasta'),
                                      partition.output_path('firstpass.fa'),
                                      partition.output_path('firstpass.bed'),
                                      partition.output_path('isoform.counts.txt'),
                                      partition.output_path('isoform.read.map.txt'), True,
                                      partition.output_path('reads.genomicclipping.txt'),
                                      partition.output_path('firstpass.uniquebound.txt'))
    else:
        logging.info('no firstpass isoforms found')
        # create empty files
        with open(partition.output_path('firstpass.fa'), 'w') as _, \
             open(partition.output_path('firstpass.bed'), 'w') as _, \
             open(partition.output_path('isoform.counts.txt'), 'w') as _, \
             open(partition.output_path('isoform.read.map.txt'), 'w') as _, \
             open(partition.output_path('isoform.ends.tsv'), 'w') as _:
            pass

    # FIXME this is messy, shouldn't have to reorganize like this
    final_transcript_objs = {}
    for og_key in firstpass:
        final_transcript_objs[firstpass[og_key].name] = firstpass[og_key]

    iso_to_counts = {}
    gene_to_tot = {}
    for line in open(partition.output_path('isoform.ends.tsv')):
        line = line.rstrip().split('\t')
        read, transcript = line[:2]
        gene = transcript.split('_')[-1]
        if gene not in gene_to_tot:
            gene_to_tot[gene] = [0, 0, 0]
        if transcript not in iso_to_counts:
            # total spliced full-length, total full-length spliced + unspliced, total all
            iso_to_counts[transcript] = [0, 0]
        if line[2] == 'None':  # single exon transcript
            startdist, enddist = int(line[3]), int(line[5])
            if args.trust_ends:
                if startdist <= TRUST_ENDS_WINDOW and enddist <= TRUST_ENDS_WINDOW:
                    iso_to_counts[transcript][0] += 1
                    gene_to_tot[gene][1] += 1
            else:
                tlen = final_transcript_objs[transcript].end - final_transcript_objs[transcript].start
                rlen = tlen - (startdist + enddist)
                # FIXME might not work if end_norm_dist is set - add test, or just remove end_norm_dist
                if rlen > (tlen / 2):
                    iso_to_counts[transcript][0] += 1
                    gene_to_tot[gene][1] += 1
        else:
            startindex, startdist, endindex, enddist = [int(x) for x in line[2:]]
            if startindex == 0 and endindex == len(final_transcript_objs[transcript].juncs) - 1:
                iso_to_counts[transcript][0] += 1
                gene_to_tot[gene][0] += 1
                gene_to_tot[gene][1] += 1
        iso_to_counts[transcript][1] += 1
        gene_to_tot[gene][2] += 1

    with open(partition.output_path('isoform.counts.txt'), 'w') as fh:
        for transcript in iso_to_counts:
            fh.write(f'{transcript}\t{iso_to_counts[transcript][0]}\t{iso_to_counts[transcript][1]}\n')

    with open(partition.output_path('isoforms.bed'), 'w') as iso_fh, \
         open(partition.output_path('isoforms.fa'), 'w') as seq_fh:
        for tname in final_transcript_objs:
            # spliced isos checked against spliced total, single exon checked against full-length total
            if _iso_passes_support_filter(args, tname, len(final_transcript_objs[tname].exons), iso_to_counts, gene_to_tot):
                # print(final_transcript_objs[tname].score, iso_to_counts[tname])
                final_transcript_objs[tname].score = iso_to_counts[tname][0]

                convert_to_bed(final_transcript_objs[tname]).write(iso_fh)
                seq_fh.write('>' + final_transcript_objs[tname].name + '\n')
                seq_fh.write(final_transcript_objs[tname].get_sequence(genome) + '\n')

    bed_to_gtf(partition.output_path('isoforms.bed'), partition.output_path('isoforms.gtf'))

    genome.close()

def combine_chunks(args, output, partitions):
    files_to_combine = ['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed',
                        '.isoforms.bed',
                        '.isoform.read.map.txt', '.isoforms.gtf', '.isoforms.fa', '.isoform.counts.txt']
    if args.end_norm_dist:
        files_to_combine.extend(('.read_ends.bed', '.matchannot.ends.tsv'))
    # if args.predict_cds:
    #     files_to_combine.append('.isoforms.CDS.bed')
    if not args.no_align_to_annot:
        files_to_combine.extend(['.matchannot.counts.txt', '.matchannot.read.map.txt'])
    combine_temp_files_by_suffix(output, [p.file_prefix for p in partitions], files_to_combine)

####
# main
####

def flair_transcriptome():
    # FIXME: split options out that are flags to indicate what to do
    # so args doesn't get passes but we don't have to pass so many options

    args = get_args()

    logging.info('loading genome')
    genome = pysam.FastaFile(args.genome)

    temp_dir = make_temp_dir(args.output)

    annot_gtf_data = None
    if args.annot_gtf:
        logging.info('loading annotation GTF')
        annot_gtf_data = gtf_data_parser(args.annot_gtf, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)

    logging.info('building intron support database')
    junction_corrector = junction_corrector_factory(args.ss_window, args.junction_support,
                                                    annot_gtf_data=annot_gtf_data,
                                                    intron_beds=args.junction_bed,
                                                    star_sj_tabs=args.junction_tab)

    logging.info('partitioning genome')
    runner = partition_runner_factory(args.parallel_mode, genome, args.genome_aligned_bam,
                                      temp_dir, args.annot_gtf, args.threads,
                                      gtf_data=annot_gtf_data, junction_corrector=junction_corrector)
    logging.info(f'number of partitions: {len(runner)}')

    logging.info('running partitions')
    runner.run(_run_region, args=args)
    combine_chunks(args, args.output, runner.partitions)

    # this needs to be done here outside of chunking because needs to load whole annotation gtf, which should only be done once
    if args.predict_cds:
        logging.info('predicting CDS')
        # FIXME: why is this passing in annotation GTF?
        predict_productivity(args.output, args.genome, args.annot_gtf)

    if not args.keep_intermediate:
        shutil.rmtree(temp_dir)

    genome.close()


if __name__ == "__main__":
    flair_transcriptome()
