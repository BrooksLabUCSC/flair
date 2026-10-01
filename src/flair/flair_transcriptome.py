#!/usr/bin/env python3
import os
import pipettor
import shutil
import pysam
import hashlib
import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from statistics import median
from collections import Counter
from flair import FlairError, FlairInputDataError, FlairNotImplementedError
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.junction_correct import junction_corrector_factory, UNKNOWN_STRAND
from flair.partition_runner import parallel_mode_parse, PartitionRunner, partition_regions, combine_temp_files_by_suffix
from flair.io_utils import make_run_temp_dir
from flair.bed_to_gtf import bed_to_gtf
from flair.isoform_data import (Exon, Gene, Isoform, ReadRec, get_bed_exons_from_exons,
                                get_sequence_for_exons, binary_search, convert_to_bed12, convert_to_flair_bed, make_big_bed)
from flair.read_processing import generate_genomic_alignment_read_to_clipping_file
from flair.read_correction import filter_correct_group_reads
from flair.count_sam_transcripts import TRUST_ENDS_WINDOW, run_count_sam_transcripts
from flair.annotation_data import annot_data_from_gtf
from flair.pycbio.hgdata.bed import BedReader
from flair.predictProductivity import predict_prod_temp
from flair.flair_bed import FlairBed
from flair.terminal_exon_ends import TerminalExonEnds
from flair import thread_share

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

FILTER_MODES = ('nosubset', 'bysupport', 'comprehensive', 'ginormous')

@dataclass(frozen=True)
class TranscriptomeOpts:
    """Options shared by the functions that make up a transcriptome run.  This is
    internal to this module: the public entry point takes named parameters, and
    this is what they are bundled into for the partition workers, which receive it
    through pickling."""
    genome_aligned_bam: str
    genome: str
    sample_name: str
    output: str
    annot_gtf: str
    junction_tab: str
    junction_bed: str
    junction_support: int
    ss_window: int
    end_window: int
    sjc_support: int
    single_exon_support: int
    frac_support: float
    directRNA: bool
    trust_strand: bool
    trust_junctions: bool
    trust_ends: bool
    no_stringent: bool
    no_check_splice: bool
    no_align_to_annot: bool
    max_ends: int
    filter: str
    keep_supplementary: bool
    quality: int
    threads: int
    parallel_mode: tuple
    fusion_breakpoints: str
    keep_intermediate: bool
    temp_dir: str
    normalize_ends: bool
    generate_map: bool

def add_subparser(subparsers):
    desc = ('generates confident transcript models directly from a bam file '
            'of aligned long rna-seq reads')
    parser = subparsers.add_parser('transcriptome', help="Build a transcriptome from aligned reads",
                                   description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-b', '--genome_aligned_bam', required=True,
                          help='Sorted and indexed bam file aligned to the genome')
    required.add_argument('-g', '--genome', type=str, required=True,
                          help='FastA of reference genome, can be minimap2 indexed')
    required.add_argument('--sample_name', required=True,
                          help='name of sample - this will be added as metadata to your output files')
    parser.add_argument('-o', '--output',
                        help='output file name base for FLAIR isoforms, defaults to --sample_name')

    parser.add_argument('-f', '--gtf', dest="annot_gtf", default=None,
                        help='GTF annotation file, used for identifying annotated isoforms')
    parser.add_argument('--junction_tab', help='short-read junctions in SJ.out.tab format. '
                                               'Use this option if you aligned your short-reads with STAR, '
                                               'STAR will automatically output this file')
    parser.add_argument('--junction_bed', help='short-read junctions in bed format '
                                               '(can be generated from long-read alignment with intron-prospector)')
    parser.add_argument('--junction_support', type=int, default=2,
                        help='if providing short-read junctions, minimum junction support required to keep junction. '
                             'If your junctions file is in bed format, the score field will be used for read support '
                             '(default: %(default)s)')

    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (default: %(default)s)')
    parser.add_argument('--end_window', type=int, default=100,
                        help='window size for comparing TSS/TES (default: %(default)s)')

    parser.add_argument('--sjc_support', type=int, default=1,
                        help='minimum number of supporting reads for a spliced isoform (default: %(default)s)')
    parser.add_argument('--single_exon_support', type=int, default=3,
                        help='minimum number of supporting reads for a single exon isoform (default: %(default)s)')
    parser.add_argument('--frac_support', type=float, default=0.05,
                        help='minimum fraction of gene locus support for isoform to be called; only isoforms '
                             'that make up more than this fraction of the gene locus are reported. Set to 0 for '
                             'max recall (default: %(default)s)')

    parser.add_argument('--directRNA', action='store_true',
                        help="input is directRNA - this sets trust_strand to True, also doesn't allow large deletions in UTRs (an artifact of cDNA amplification)")
    parser.add_argument('--trust_strand', action='store_true',
                        help='trust the stranding of the input reads and do not attempt strand correction: '
                             'reads keep the strand of their alignment, and spliced reads are corrected only '
                             'with splice junctions on that strand, and are only assigned to transcripts '
                             'on that strand')
    parser.add_argument('--trust_junctions', action='store_true',
                        help='trust the strand of every input splice junction for read correction, without '
                             'checking splice motifs.  By default, an unannotated junction without a GT-AG '
                             'motif, such as GC-AG or AT-AC, is only weak evidence of strand: a read stranded '
                             'only by such junctions gets its strand from gene identification, or is dropped')
    parser.add_argument('--trust_ends', action='store_true',
                        help='trust the ends of the input reads: a more stringent way of requiring read ends '
                             'to match the ends of transcript models')

    parser.add_argument('--no_stringent', action='store_true',
                        help="do not require all supporting reads to be full-length, that is aligned to the "
                             "first and last exons of the transcript. Use this for fragmented libraries, "
                             "with an understanding that it will impact precision")
    parser.add_argument('--no_check_splice', action='store_true',
                        help="do not enforce accurate alignment around splice sites. Specify this for "
                             "libraries with high error rates, but it will reduce precision")
    parser.add_argument('--no_align_to_annot', action='store_true',
                        help="skip the initial alignment to the annotated sequences and detect transcripts "
                             "only from the genomic alignment. Slightly faster but less accurate when the "
                             "annotation is good")

    parser.add_argument('--max_ends', type=int, default=1,
                        help='maximum number of TSS/TES picked per isoform; make higher for more precise '
                             'end detection (default: %(default)s)')
    parser.add_argument('--filter', choices=FILTER_MODES, default='nosubset',
                        help='nosubset: any isoforms that are a proper set of another isoform are removed; '
                             'bysupport: subset isoforms are removed based on support; '
                             'comprehensive: default set plus all subset isoforms; '
                             'ginormous: comprehensive set plus single exon subset isoforms '
                             '(default: %(default)s)')

    parser.add_argument('--keep_supplementary', action='store_true',
                        help='keep supplementary alignments when defining isoforms')
    parser.add_argument('--quality', default=1, type=int,
                        help='minimum mapping quality threshold to consider genomic alignments for '
                             'defining transcripts (default: %(default)s)')
    parser.add_argument('--allow_paralogs', action='store_true',
                        help='NOT IMPLEMENTED: assign reads to multiple paralogs with equivalent '
                             'alignment. Specifying this is an error rather than a no-op')

    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='number of threads to run with, related to --parallel_mode (default: %(default)s)')
    parser.add_argument('--parallel_mode', default='auto:1GB',
                        help='parallelization mode. auto:1GB is an automatic threshold where, if the file '
                             'is less than 1GB, parallelization is done by chromosome, and if it is larger, '
                             'parallelization is done by region of non-overlapping reads. Other modes: '
                             'bychrom, byregion, auto:xGB; for the auto threshold, the size must be in units '
                             'of GB (default: %(default)s)')

    parser.add_argument('--fusion_breakpoints',
                        help='for fusion detection only: bed file containing locations of fusion breakpoints '
                             'on the synthetic genome')

    parser.add_argument('--temp_dir',
                        help='directory for temporary files; each run makes its own directory, named for the '
                             'output, in it.  Many small files are written and removed, so a local disk is much '
                             'faster than network storage (default: $TMPDIR or the system temporary directory)')
    parser.add_argument('--keep_intermediate', action='store_true',
                        help='keep intermediate and temporary files for debugging, in the run\'s directory in '
                             '--temp_dir, which must be given. Intermediate files include the promoter-supported '
                             'reads file and read assignments to firstpass isoforms')

    parser.add_argument('--normalize_ends', action='store_true',
                        help='normalize transcript ends with similar terminal splice sites; only recommended '
                             'when --max_ends is 1')
    parser.add_argument('--generate_map', action='store_true',
                        help='generate a txt file of read-isoform assignments')
    parser.set_defaults(entry=transcriptome_cmd)

def transcriptome_cmd(args):
    if args.allow_paralogs:
        raise FlairNotImplementedError("--allow_paralogs is not implemented: a read with an equally "
                                       "good alignment to several paralogs is assigned to one of "
                                       "them, and nothing downstream does otherwise")
    if args.keep_intermediate and args.temp_dir is None:
        raise FlairInputDataError('--keep_intermediate requires --temp_dir, the directory to keep them in')
    for what, path in (('aligned reads bam', args.genome_aligned_bam), ('genome fasta', args.genome)):
        if not os.path.exists(path):
            raise FlairInputDataError(f'{what} file does not exist: {path}')
    flair_transcriptome(genome_aligned_bam=args.genome_aligned_bam, genome=args.genome,
                        sample_name=args.sample_name,
                        output=args.output if args.output is not None else args.sample_name,
                        annot_gtf=args.annot_gtf, junction_tab=args.junction_tab,
                        junction_bed=args.junction_bed, junction_support=args.junction_support,
                        ss_window=args.ss_window, end_window=args.end_window,
                        sjc_support=args.sjc_support,
                        single_exon_support=args.single_exon_support,
                        frac_support=args.frac_support, directRNA=args.directRNA, trust_strand=args.trust_strand,
                        trust_junctions=args.trust_junctions, trust_ends=args.trust_ends, no_stringent=args.no_stringent,
                        no_check_splice=args.no_check_splice,
                        no_align_to_annot=args.no_align_to_annot, max_ends=args.max_ends,
                        filter=args.filter, keep_supplementary=args.keep_supplementary,
                        quality=args.quality, threads=args.threads,
                        parallel_mode=parallel_mode_parse(args.parallel_mode),
                        fusion_breakpoints=args.fusion_breakpoints,
                        keep_intermediate=args.keep_intermediate, temp_dir=args.temp_dir,
                        normalize_ends=args.normalize_ends, generate_map=args.generate_map)


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
# in gene assignment, a gene's matched splice site this close to the isoform's
# 5' splice site counts as matching the isoform's 5' end
FIVE_PRIME_SS_WINDOW = 20
# in the gene assignment fallback with fixed ends, a spliced isoform's terminal exons
# are taken to be this long, instead of reaching its read ends
FALLBACK_TERMINAL_EXON_LEN = 100


####
# transcriptome alignment
####
def transcriptome_align_and_count(args, input_reads, align_ref_fasta, ref_bed, output_name, map_file, is_annot, clipping_file, unique_bound, directRNA):  # noqa: C901 - FIXME: reduce complexity
    # minimap (results are piped into count_sam_transcripts.py)
    # '--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    if isinstance(input_reads, str):
        input_reads = [input_reads]
    mm2_cmd = ['minimap2', '-a', '-N', '4', '--MD'] + [align_ref_fasta] + input_reads

    # FIXME add in step to filter out chimeric reads here
    # FIXME really need to go in and check on how count_sam_transcripts is working
    trimmedreads = clipping_file or None
    generate_map = map_file or None
    output_endpos = output_name.split('.counts.txt')[0] + '.ends.tsv'  # if (args.output_endpos or is_annot) else None)
    stringent = (not is_annot) and (not args.no_stringent)
    check_splice = not args.no_check_splice
    unique_bound_path = unique_bound if unique_bound and (not args.no_stringent or is_annot) else None

    # minimap2 borrows threads left idle in the thread budget while it runs
    with thread_share.program_threads() as mm2_threads:
        run_count_sam_transcripts(
            mm2_cmd=mm2_cmd[:1] + ['-t', str(mm2_threads)] + mm2_cmd[1:],
            output=output_name,
            trimmedreads=trimmedreads,
            generate_map=generate_map,
            output_endpos=output_endpos,
            stringent=stringent,
            allow_UTR_indels=not directRNA,
            check_splice=check_splice,
            isoforms=ref_bed,
            trust_ends=args.trust_ends,
            unique_bound=unique_bound_path,
            fusion_breakpoints=args.fusion_breakpoints,
            # reads come from samtools fasta in their sequenced orientation, and
            # transcript sequences are sense, so an antisense read aligns reversed
            stranded=args.trust_strand)


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
    """Check if juncs is a contiguous run of junctions inside the longer
    otheriso_juncs.  A junction chain is sorted and has no repeats, so the run can
    only start where juncs' first junction is."""
    if len(juncs) >= len(otheriso_juncs):
        return False
    if len(juncs) == 0:
        return True
    try:
        start = otheriso_juncs.index(juncs[0])
    except ValueError:
        return False
    return tuple(otheriso_juncs[start:start + len(juncs)]) == tuple(juncs)


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
    Records unique sequence boundaries and checks containment within tolerance.
    A boundary is only recorded when the terminal exon extends past the other
    exon; when it is inside it there is no unique sequence to require."""
    if first_exon.end == other_exon.end:
        if first_exon.start < other_exon.start:
            unique_seq_bound.append((0, first_exon.end - other_exon.start))
        if first_exon.start >= (other_exon.start - TERMINAL_EXON_BOUNDARY_TOLERANCE):
            terminal_exon_is_subset[0] = 1
            superset_support.append(otheriso_score)
    if last_exon.start == other_exon.start:
        if last_exon.end > other_exon.end:
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
                       sup_annot_transcript_to_juncs, strand):
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
def _exon_bounds(exons):
    return [(e.start, e.end) for e in exons]

def max_terminal_exons_ends_from_annots(annots):
    max_terminal_exons_ends = TerminalExonEnds()
    for transcript_id, gene_id, strand in annots.transcripts:
        exons = annots.transcript_to_exons[(transcript_id, gene_id)]
        max_terminal_exons_ends.add_transcript(gene_id, strand, _exon_bounds(exons))
    return max_terminal_exons_ends

def max_terminal_exons_ends_from_iso_infos(iso_to_info):
    max_terminal_exons_ends = TerminalExonEnds()
    for iso_name in iso_to_info:
        isoform = iso_to_info[iso_name]
        max_terminal_exons_ends.add_transcript(isoform.gene_id, isoform.strand, _exon_bounds(isoform.exons))
    return max_terminal_exons_ends

####
# transcriptome reference
####
def normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons, chrom_len):
    "updates terminal exons ends"
    start, end = max_terminal_exons_ends.normalized_ends(gene_id, strand, _exon_bounds(exons), chrom_len)
    exons[0] = Exon(start, exons[0].end)
    exons[-1] = Exon(exons[-1].start, end)
    return exons

def generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, normalize_ends, max_terminal_exons_ends,
                                                transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand[(transcript_id, gene_id)] = strand
    exons = list(annots.transcript_to_exons[(transcript_id, gene_id)])
    assert isinstance(exons[0], Exon)  # FIXME tmp debugging
    if normalize_ends and len(exons) > 1:
        normalize_gene_terminal_exons(max_terminal_exons_ends, gene_id, strand, exons, genome.get_reference_length(chrom))
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

def generate_transcriptome_reference_guts(normalize_ends, annots, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh):
    transcript_to_strand = {}
    transcript_to_new_exons = {}
    max_terminal_exons_ends = None
    if normalize_ends:
        max_terminal_exons_ends = max_terminal_exons_ends_from_annots(annots)

    for transcript_id, gene_id, strand in annots.transcripts:
        generate_transcriptome_reference_transcript(strand, transcript_to_strand, transcript_id, gene_id, annots, normalize_ends, max_terminal_exons_ends,
                                                    transcript_to_new_exons, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)
    return transcript_to_strand, transcript_to_new_exons

def generate_transcriptome_reference(temp_prefix, annots, chrom, genome, normalize_ends=False):
    with (open(temp_prefix + '.annotated_transcripts.bed', 'w') as annot_bed_fh,
          open(temp_prefix + '.annotated_transcripts.fa', 'w') as annot_fa_fh,
          open(temp_prefix + '.annotated_transcripts_uniquebound.txt', 'w') as annot_uniqueseq_fh):
        return generate_transcriptome_reference_guts(normalize_ends, annots, chrom, genome, annot_bed_fh, annot_fa_fh, annot_uniqueseq_fh)


def identify_good_match_to_annot(args, temp_prefix, chrom, annots, genome):
    # FIXME: refactor
    # good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs = [], set(), {}
    read_to_transcript = {}
    if not args.no_align_to_annot and len(annots.transcripts) > 0:
        # logging.info('generating transcriptome reference')
        # this part generates the fasta file for the annotation
        transcript_to_strand, transcript_to_new_exons = \
            generate_transcriptome_reference(temp_prefix, annots, chrom, genome, normalize_ends=args.normalize_ends)
        # FIXME: make a TSV
        clipping_file = temp_prefix + '.reads.genomicclipping.txt'
        transcriptome_align_and_count(args, temp_prefix + '.reads.fasta',
                                      temp_prefix + '.annotated_transcripts.fa',
                                      temp_prefix + '.annotated_transcripts.bed',
                                      temp_prefix + '.matchannot.counts.txt',
                                      None, True,
                                      clipping_file,
                                      temp_prefix + '.annotated_transcripts_uniquebound.txt',
                                      args.directRNA)
        for line in open(temp_prefix + '.matchannot.ends.tsv'):
            line = line.rstrip().split('\t')
            read, transcript = line[:2]
            start_sj_index, start_sj_dist, start_tend_dist, end_sj_index, end_sj_dist, end_tend_dist = [int(x) if x != 'None' else None for x in line[2:]]
            # is not None: the value was converted on the line above, so the old
            # comparison with the string 'None' was always true
            if start_sj_index is not None and end_sj_index is not None:  # not a single exon transcript
                read_to_transcript[read] = (transcript, start_sj_index, start_sj_dist, end_sj_index, end_sj_dist)
    # good_align_to_annot = set(good_align_to_annot)
    # return good_align_to_annot, firstpass_SE, sup_annot_transcript_to_juncs
    return read_to_transcript


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

def filter_ends_single_best(isoforms):
    """Pick single best end from junction chain.
    Returns list with single Isoform object."""
    # best_only uses the default sorting, doesn't require additional action
    # Pick single best end and merge all reads into it
    best = isoforms[0]
    for iso in isoforms[1:]:
        best.reads.extend(iso.reads)
    return [best]

def filter_ends_by_redundant_and_support(isoforms, sjc_support, se_support, max_ends, normalize_ends):
    """Sort ends, then select best ones based on support and max_ends"""
    if isoforms[0].juncs == ():
        support = se_support
    else:
        support = sjc_support

    if normalize_ends:  # Only by longest length
        isoforms.sort(key=lambda x: x.genomic_length, reverse=True)
    else:  # First by read support, then by length
        # FIXME: the comment above is what was meant; the code multiplies instead, so
        # this is not a two-level sort and a 1-read 10kb candidate outranks a 10-read
        # 500bp one.  With max_ends 1 this single comparison sets the reported TSS and
        # TES for every junction chain.  Sorting by (num_reads, genomic_length), which
        # is what the comment says, moves exactly one locus in the test suite and makes
        # it worse: on a single-exon chain it picks the modal end pair, which for
        # long reads is a 5' truncation cluster, over the pair that reproduces the
        # annotated 3' end.  So the product is compensating for a real bias in the
        # single-exon case while having little to recommend it for spliced chains,
        # where the end window is only the terminal exons.  Probably wants to branch on
        # isoforms[0].juncs == (), which this function already distinguishes above.
        isoforms.sort(key=lambda x: [x.num_reads * x.genomic_length], reverse=True)

    junc_support = sum([x.num_reads for x in isoforms])
    if junc_support < support:
        logging.debug(f"isoform group dropped: insufficient support ({junc_support} < {support}): {isoforms[0].chrom}:{isoforms[0].start}-{isoforms[0].end}")
        return []

    if max_ends > 1:
        # Allow multiple ends per junction chain
        return filter_ends_allow_multiple(isoforms, support, max_ends)
    else:
        # Pick single best end
        return filter_ends_single_best(isoforms)


def _write_unfiltered_ends(isoforms, fh):
    for iso_readrec in isoforms:
        convert_to_bed12(iso_readrec).write(fh)

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


def _filter_isos_by_redundant_and_support(args, isoforms, candidates, iso_fh):
    # this assumes single exons are pre-grouped by overlap
    # previously treated single exons separately due to them being in larger groups
    filtered_isoforms = filter_ends_by_redundant_and_support(isoforms, args.sjc_support, args.single_exon_support, args.max_ends, args.normalize_ends)
    for isoform in filtered_isoforms:
        candidates.add(isoform)
        convert_to_bed12(isoform).write(iso_fh)

def _generate_candidate_isos(args, isoform, candidates, iso_fh, iso_unfilt_fh):
    # NOTE: Harrison's TED code will be slotted in here to replace collapse_end_groups
    these_firstpass = collapse_end_groups(args.end_window, isoform)
    _write_unfiltered_ends(these_firstpass, iso_unfilt_fh)
    _filter_isos_by_redundant_and_support(args, these_firstpass, candidates, iso_fh)


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

    def remove(self, key):
        del self._groups[key]

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
    sjc_with_overlap_groups = group_by_overlap(sj_to_ends, args.single_exon_support, args.trust_strand)
    # FIXME everything below here requires confidence in transcript strand
    build_genes(sjc_with_overlap_groups, annots, region_chrom, sjc_with_overlap_groups)

    candidates = CandidateIsoforms()
    with open(temp_prefix + '.firstpass.unfiltered.bed', 'w') as iso_fh, \
            open(temp_prefix + '.firstpass.reallyunfiltered.bed', 'w') as iso_unfilt_fh:
        for juncs, isoform in sjc_with_overlap_groups.items():
            _generate_candidate_isos(args, isoform, candidates, iso_fh, iso_unfilt_fh)
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
                if exon.name == '' or args.filter == 'nosubset':  # is exon from spliced transcript
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
                                                                   sup_annot_transcript_to_juncs, isoform.strand)
                    if not is_not_subset:
                        logging.debug(f"isoform dropped: subset of another isoform: {iso_name} ({isoform.num_reads} reads)")
                    else:
                        firstpass[iso_name] = isoform
                        if len(unique_seq) > 0:
                            iso_to_unique_bound[iso_name] = ','.join(unique_seq)
        # HANDLE SINGLE EXONS SEPARATELY - group first - one traversal of list
        firstpass = filter_all_single_exon(args, sorted(candidates.exons), candidates.isoforms, firstpass)

    return firstpass, iso_to_unique_bound

def _gene_span(gene_id, annots):
    """span of the gene's annotated junctions, so a long first or last exon does
    not make a gene look longer or overlap more genes"""
    juncs = annots.gene_to_annot_juncs[gene_id]
    return min(j.start for j in juncs), max(j.end for j in juncs)

def _drop_subset_genes(gene_to_matches):
    """keep the genes whose set of matches to the isoform (junctions or splice
    sites) is not a proper subset of another gene's set"""
    return [g for g, matches in gene_to_matches.items()
            if not any(matches < other for other in gene_to_matches.values())]

def _group_overlapping_genes(genes, annots):
    """group genes whose genomic spans overlap, transitively"""
    groups, group_end = [], None
    for (start, end), gene_id in sorted((_gene_span(g, annots), g) for g in genes):
        if groups and start < group_end:
            groups[-1].append(gene_id)
            group_end = max(group_end, end)
        else:
            groups.append([gene_id])
            group_end = end
    return groups

def _pick_group_gene(group, gene_to_juncs, gene_to_sites, five_prime_ss, strand, annots):
    """Pick one gene from a group of overlapping genes, none of whose matches to
    the isoform contain another's, by the highest score: 2 for each of the
    isoform's junctions the gene has annotated, 1 for each other splice site of
    the isoform it has annotated, and 1 if its matched splice site nearest the
    isoform's 5' end is within FIVE_PRIME_SS_WINDOW of the isoform's 5' splice
    site.  Ties go to the gene with the shortest junction span, then the gene id."""
    def score(gene_id):
        juncs = gene_to_juncs.get(gene_id, ())
        junc_sites = {site for j in juncs for site in (j.start, j.end)}
        sites = gene_to_sites.get(gene_id, set())
        matched = junc_sites | sites
        five_prime_match = min(matched) if strand == '+' else max(matched)
        return (2 * len(juncs) + len(sites - junc_sites)
                + (1 if abs(five_prime_match - five_prime_ss) < FIVE_PRIME_SS_WINDOW else 0))

    def length(gene_id):
        start, end = _gene_span(gene_id, annots)
        return end - start
    return min(group, key=lambda g: (-score(g), length(g), g))

def _genes_sharing_splice_sites(sites, strand, annots):
    """map of same-strand gene -> set of the isoform's splice sites it has
    annotated"""
    gene_to_sites = {}
    for site in sites:
        for gene_id in annots.splice_site_to_genes.get(site, ()):
            if annots.gene_to_strand[gene_id] == strand:
                gene_to_sites.setdefault(gene_id, set()).add(site)
    return gene_to_sites

def get_genes_with_shared_juncs(juncs, exons, strand, annots):
    """Assign a spliced isoform to annotated genes by the junctions and individual
    splice sites it shares with them.  A gene is dropped if its matches, junctions
    and splice sites together, are a proper subset of another gene's.  The
    remaining genes are grouped by genomic overlap of their junction spans: each
    group of overlapping genes, such as a gene cluster sharing exons or a
    readthrough gene and its parts, gives one gene, chosen by _pick_group_gene,
    while genes that do not overlap each other, as for a readthrough isoform
    without a readthrough annotation, are all returned.  Returns a tuple of gene
    ids, empty if no gene shares a junction or splice site.  When some gene shares
    a junction, genes sharing only splice sites compete only within a group with
    such a gene; they don't make groups of their own, which would report a
    neighboring gene that shares a splice site as if the isoform read through it."""
    gene_to_juncs = {}
    for j in juncs:
        if j in annots.junc_to_gene:
            for transcript_id, gene_id in annots.junc_to_gene[j]:
                if gene_id not in gene_to_juncs:
                    gene_to_juncs[gene_id] = set()
                gene_to_juncs[gene_id].add(j)
    sites = {site for j in juncs for site in (j.start, j.end)}
    gene_to_sites = _genes_sharing_splice_sites(sites, strand, annots)
    # a gene with a matching junction has its splice sites, strand or not
    for gene_id, gene_juncs in gene_to_juncs.items():
        gene_to_sites.setdefault(gene_id, set()).update(site for j in gene_juncs for site in (j.start, j.end))
    gene_to_matches = {g: {('junc', j) for j in gene_to_juncs.get(g, ())} | {('site', s) for s in gene_sites}
                       for g, gene_sites in gene_to_sites.items()}
    genes = _drop_subset_genes(gene_to_matches)
    groups = _group_overlapping_genes(genes, annots)
    if gene_to_juncs:
        groups = [group for group in groups if any(g in gene_to_juncs for g in group)]
    # the isoform's first splice site, at the inner edge of its 5' exon
    five_prime_ss = exons[-1].start if strand == '-' else exons[0].end
    return tuple(sorted(_pick_group_gene(group, gene_to_juncs, gene_to_sites, five_prime_ss, strand, annots)
                        for group in groups))


def get_single_exon_gene_overlaps(strand, iso_readrec, annots):
    gene_hits = {}
    exon = iso_readrec.exons[0]
    index = binary_search(exon, annots.all_annot_SE[strand])
    # FIXME: how does this ever work? all_annot_SE is [(start, end, strand, gene_id), ...]
    # max(0, ...): a negative slice start would be read as an offset from the end of
    # the list, which for a long list gives an empty window
    for annot_exon_info in annots.all_annot_SE[strand][max(0, index - ANNOT_SE_SEARCH_WINDOW):index + ANNOT_SE_SEARCH_WINDOW]:
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

def _merge_intervals(intervals):
    "merge overlapping (start, end) intervals, sorted by start"
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged

def get_unspliced_exon_overlaps(strand, exons, annots):
    """Genes on the strand whose single-exon (unspliced) annotated transcripts cover
    more than half of the exons' length, as [bases covered, gene_id, strand]"""
    gene_to_intervals = {}
    for annot_exon in annots.all_annot_SE[strand]:
        gene_to_intervals.setdefault(annot_exon.name, []).append((annot_exon.start, annot_exon.end))
    exon_len = sum(e.end - e.start for e in exons)
    gene_hits = []
    for gene_id, intervals in gene_to_intervals.items():
        covered = sum(max(0, min(e.end, a_end) - max(e.start, a_start))
                      for e in exons for a_start, a_end in _merge_intervals(intervals))
        if covered > exon_len * 0.5:
            gene_hits.append([covered, gene_id, strand])
    return gene_hits

def fixed_end_exons(juncs):
    """exons of a spliced isoform from its junctions, with terminal exons of
    FALLBACK_TERMINAL_EXON_LEN rather than reaching its read ends"""
    exons = [Exon(max(0, juncs[0].start - FALLBACK_TERMINAL_EXON_LEN), juncs[0].start)]
    exons.extend(Exon(juncs[i].end, juncs[i + 1].start) for i in range(len(juncs) - 1))
    exons.append(Exon(juncs[-1].end, juncs[-1].end + FALLBACK_TERMINAL_EXON_LEN))
    return exons

def _get_transcript_gene_from_annot(iso_readrec, annots):
    """Return (transcript_id, gene_id) if iso matches an annotated junction chain, else (None, None).
    Each junction chain is named once, before end variants are split off; the variants
    keep this transcript_id."""
    if iso_readrec.juncs != () and iso_readrec.juncs in annots.juncchain_to_transcript:
        transcript_id, gene_id = annots.juncchain_to_transcript[iso_readrec.juncs]
        return transcript_id, (gene_id, )
    else:
        return None, None


def _find_gene_id_by_splicing(iso_readrec, annots):
    """Genes sharing junctions or splice sites with a spliced isoform, or, for a
    single-exon isoform, the best overlapping single-exon gene, as a tuple of gene
    ids, or None"""
    if iso_readrec.juncs != ():
        gene_hits = get_genes_with_shared_juncs(iso_readrec.juncs, iso_readrec.exons, iso_readrec.strand, annots)
        if gene_hits:
            return gene_hits
    else:
        gene_hits = get_single_exon_gene_overlaps(iso_readrec.strand, iso_readrec, annots)
        if gene_hits:
            return (sorted(gene_hits.items(), key=lambda x: x[1], reverse=True)[0][0], )
    return None


def _find_gene_id_by_exon_overlap(iso_readrec, annots):
    """The gene whose exons best overlap the isoform's, as a one-gene tuple, or None"""
    # look for exon overlap.  A spliced isoform's terminal exons
    # are given a fixed length from its terminal junctions (fixed_end_exons), so this
    # doesn't depend on its read ends.  There was an 'ambig' strand branch here; nothing
    # assigns that strand, the two branches were identical, and annots.spliced_exons is
    # keyed by '+' and '-' only, so it would have raised
    if iso_readrec.juncs != ():
        # spliced genes first, then single-exon genes, such as an unspliced lncRNA
        # the isoform is a spliced form of
        exons = fixed_end_exons(iso_readrec.juncs)
        gene_hits = (get_spliced_exon_overlaps(iso_readrec.strand, exons, annots)
                     or get_unspliced_exon_overlaps(iso_readrec.strand, exons, annots))
    else:
        gene_hits = get_spliced_exon_overlaps(iso_readrec.strand, iso_readrec.exons, annots)
    if gene_hits:
        gene_hits.sort(reverse=True)
        return (gene_hits[0][1], )
    else:
        return None


def _find_gene_id_by_overlap(iso_readrec, annots):
    """Find gene_id for an isoform without a matching junction chain, using junction or exon overlap."""
    # this all requires that we already trust the strand of the transcript
    # returns tuple of matching genes, will go into ref_gene_id field
    return (_find_gene_id_by_splicing(iso_readrec, annots)
            or _find_gene_id_by_exon_overlap(iso_readrec, annots))


def get_gene_name_unknown_strand(isoform, annots):
    """Gene identification for a spliced isoform of UNKNOWN_STRAND, whose junctions
    were only weakly stranded, which also gives its strand.  A matching annotated
    junction chain gives the gene.  Otherwise the isoform is tried on each strand,
    first by shared splicing and then by exon overlap; the first of these to find
    genes, all on one strand, gives the gene.  Returns (gene ids, transcript id,
    strand), with gene ids and strand None if no gene was found or genes were found
    on both strands."""
    transcript_id, gene_id = _get_transcript_gene_from_annot(isoform, annots)
    if transcript_id is not None:
        return gene_id, transcript_id, annots.gene_to_strand[gene_id[0]]
    for find_gene_id in (_find_gene_id_by_splicing, _find_gene_id_by_exon_overlap):
        strand_to_genes = {}
        for strand in ('+', '-'):
            isoform.strand = strand
            gene_id = find_gene_id(isoform, annots)
            if gene_id:
                # gene strand, as genes sharing a junction are found on either strand
                strand_to_genes.setdefault(annots.gene_to_strand[gene_id[0]], gene_id)
        isoform.strand = UNKNOWN_STRAND
        if len(strand_to_genes) == 1:
            strand, gene_id = strand_to_genes.popitem()
            return gene_id, None, strand
        elif len(strand_to_genes) > 1:
            break
    return None, None, None


def get_gene_name_firstpass(isoform, annots):
    transcript_id, gene_id = _get_transcript_gene_from_annot(isoform, annots)
    if transcript_id is None:
        gene_id = _find_gene_id_by_overlap(isoform, annots)
    return gene_id, transcript_id

def add_gene_isoform(genes, gene_id, isoform, strand, is_novel):
    # gene_id is a tuple of annotated gene ids, or one novel locus name.  A bare
    # ','.join put a comma between every character of the novel name
    key = gene_id if isinstance(gene_id, str) else ','.join(gene_id)
    hashed_id = int(hashlib.md5(key.encode('utf-8')).hexdigest(), 16)
    if is_novel:
        gene_id = ()
    if hashed_id not in genes:
        genes[hashed_id] = Gene(hashed_id, gene_id, isoform.chrom, strand)
    # this command also sets the gene_id in the isoform object
    genes[hashed_id].add_isoform(isoform)


def build_genes(firstpass, annots, region_chrom, sjc_with_overlap_groups):
    """Assign gene names to firstpass isoforms, building Gene objects.

    Returns (genes, novel_gene_isos_to_group) where:
    - genes: dict of gene_id -> Gene for isoforms matched to known genes
    - novel_gene_isos_to_group: isoforms needing novel gene assignment
    """
    genes = {}
    novel_gene_isos_to_group = {'+': [], '-': []}
    for iso_key in list(firstpass):
        isoform = firstpass[iso_key]
        if isoform.strand == UNKNOWN_STRAND:
            # only weakly stranded by its junctions; the strand has to come from a
            # gene, as there can't be isoforms of unknown strand in the output
            gene_id, isoform_id, strand = get_gene_name_unknown_strand(isoform, annots)
            if gene_id is None:
                logging.debug(f"isoform dropped: strand unknown and not resolved by gene identification: {iso_key}")
                firstpass.remove(iso_key)
                continue
            isoform.strand = strand
        else:
            gene_id, isoform_id = get_gene_name_firstpass(isoform, annots)
        isoform.ref_transcript_id = isoform_id
        if gene_id is not None:
            # removing this strand correction breaks the unusual junction (due to underlying variant?) test
            # currently just using first gene in list, ideally would use the most 5' annotated gene for this strand correction, if we want to do strand correction at all
            isoform.strand = annots.gene_to_strand[gene_id[0]]
            add_gene_isoform(genes, gene_id, isoform, annots.gene_to_strand[gene_id[0]], is_novel=False)
        else:
            novel_gene_isos_to_group[isoform.strand].append((isoform.start, isoform.end, iso_key))

    for strand in novel_gene_isos_to_group:
        generate_non_gene_iso_groups_strand(genes, novel_gene_isos_to_group, strand, region_chrom, sjc_with_overlap_groups)

    return genes


def _sweep_overlap_groups(spans):
    """group (start, end, key) spans whose coordinates overlap, transitively;
    returns lists of keys"""
    groups, group_end = [], None
    for start, end, key in sorted(spans):
        if groups and start < group_end:
            groups[-1].append(key)
            group_end = max(group_end, end)
        else:
            groups.append([key])
            group_end = end
    return groups

def _junction_span(isoform):
    return isoform.juncs[0].start, isoform.juncs[-1].end

def _group_by_read_span(isos, firstpass):
    "group isoforms by overlap of the spans between their read ends"
    return _sweep_overlap_groups([(firstpass[k].start, firstpass[k].end, k) for k in isos])

def _splice_site_components(spliced, firstpass):
    """connected components of spliced isoforms sharing a splice site, a shared
    junction sharing both of its sites"""
    parent = {k: k for k in spliced}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k
    site_owner = {}
    for k in spliced:
        for site in {site for j in firstpass[k].juncs for site in (j.start, j.end)}:
            if site in site_owner:
                parent[find(k)] = find(site_owner[site])
            else:
                site_owner[site] = k
    components = {}
    for k in spliced:
        components.setdefault(find(k), []).append(k)
    return list(components.values())

def _attach_single_exon_isos(single, groups, firstpass):
    """add each single-exon isoform to the group whose exons it overlaps most;
    returns those overlapping no group's exons"""
    group_exons = [[(e.start, e.end) for k in group for e in firstpass[k].exons] for group in groups]
    unattached = []
    for k in single:
        start, end = firstpass[k].start, firstpass[k].end
        overlaps = [sum(max(0, min(end, e_end) - max(start, e_start)) for e_start, e_end in exons)
                    for exons in group_exons]
        best = max(range(len(groups)), key=lambda i: overlaps[i], default=None)
        if best is not None and overlaps[best] > 0:
            groups[best].append(k)
        else:
            unattached.append(k)
    return unattached

def _group_novel_by_splicing(novel, firstpass):
    """Spliced isoforms sharing a junction or splice site are grouped, and groups
    whose junction spans overlap are merged, all without read ends.  A single-exon
    isoform joins the spliced group whose exons it overlaps most; single-exon
    isoforms overlapping no spliced group are grouped by their read-end spans."""
    spliced = [k for k in novel if firstpass[k].juncs]
    single = [k for k in novel if not firstpass[k].juncs]
    comp_spans = []
    for members in _splice_site_components(spliced, firstpass):
        spans = [_junction_span(firstpass[k]) for k in members]
        comp_spans.append((min(s for s, e in spans), max(e for s, e in spans), tuple(members)))
    groups = [[k for members in merged for k in members] for merged in _sweep_overlap_groups(comp_spans)]
    unattached = _attach_single_exon_isos(single, groups, firstpass)
    return groups + _group_by_read_span(unattached, firstpass)

def generate_non_gene_iso_groups_strand(genes, novel_gene_isos_to_group, strand, chrom, firstpass):
    """Group the strand's isoforms that have no annotated gene into novel genes, by
    _group_novel_by_splicing, and create their Gene objects.  A novel gene is named
    for the span of its isoforms."""
    novel = [iso_key for start, end, iso_key in novel_gene_isos_to_group[strand]]
    used_ids = set()
    for group in _group_novel_by_splicing(novel, firstpass):
        group_start = min(firstpass[k].start for k in group)
        group_end = max(firstpass[k].end for k in group)
        gene_id = f'{chrom}:{group_start}-{group_end}:{strand}'
        # groups can now share a span; a shared name would merge them
        suffix = 1
        while gene_id in used_ids:
            suffix += 1
            gene_id = f'{chrom}:{group_start}-{group_end}:{strand}.{suffix}'
        used_ids.add(gene_id)
        for k in group:
            add_gene_isoform(genes, gene_id, firstpass[k], strand, is_novel=True)

def write_first_pass_isoforms(iso_name, normalize_ends, isoform, max_terminal_exons_ends, unique_bound, unique_fh, iso_fh, seq_fh, genome):
    # FIXME: do normalization outside of write function
    if normalize_ends and len(isoform.exons) > 1:  # don't normalize ends for single exon transcripts
        # kept so the padding can be removed from the final isoform, since clamping
        # to the chromosome means it is not always NORM_END_EXTRA_LEN
        isoform.unpadded_ends = max_terminal_exons_ends.furthest_ends(isoform.gene_id, isoform.strand, _exon_bounds(isoform.exons))
        exons = normalize_gene_terminal_exons(max_terminal_exons_ends, isoform.gene_id, isoform.strand, isoform.exons,
                                              genome.get_reference_length(isoform.chrom))
        isoform.reset_from_exons(exons)
    # if isoform.transcript_id is None:
    #     isoform.transcript_id = isoform.name
    # isoform.name = isoform.transcript_id + '_' + isoform.gene_id

    if unique_bound and iso_name in unique_bound:
        unique_fh.write(isoform.name + '\t' + unique_bound[iso_name] + '\n')

    convert_to_bed12(isoform).write(iso_fh)
    seq_fh.write('>' + isoform.name + '\n')
    seq_fh.write(isoform.get_sequence(genome) + '\n')

def write_firstpass(temp_prefix, chrom, firstpass, annots, genome, *,
                    normalize_ends=False, unique_bound=None):

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
                                      unique_bound, unique_fh, iso_fh, seq_fh, genome)


####
# results output
####

def _iso_passes_support_filter(args, iso, gene, num_exons, iso_to_counts, gene_to_tot):
    if iso not in iso_to_counts:
        return False, 0
    else:
        count = iso_to_counts[iso][0]
        if num_exons > 1:
            return (count >= args.sjc_support) and (count / gene_to_tot[gene][0]) >= args.frac_support, (count / gene_to_tot[gene][0])
        else:
            return (count >= args.single_exon_support) and (count / gene_to_tot[gene][1]) >= args.frac_support, (count / gene_to_tot[gene][1])

def generate_empty_intermediate_files(file_prefix, suffixes):
    for s in suffixes:
        out = open(file_prefix + s, 'w')
        out.close()

def calc_final_iso_support(read_ends_file, final_transcript_objs, trust_ends, no_stringent):
    iso_to_counts = {}
    gene_to_tot = {}
    for line in open(read_ends_file):
        line = line.rstrip().split('\t')
        read, transcript = line[:2]
        start_sj_index, start_sj_dist, start_tend_dist, end_sj_index, end_sj_dist, end_tend_dist = [int(x) if x != 'None' else None for x in line[2:]]
        gene = final_transcript_objs[transcript].gene_id
        if gene not in gene_to_tot:
            # total spliced full-length, total full-length spliced + unspliced, total all
            gene_to_tot[gene] = [0, 0, 0]
        if transcript not in iso_to_counts:
            iso_to_counts[transcript] = [0, 0]
        if start_sj_index is None:  # single exon transcript
            if trust_ends:
                if start_tend_dist <= TRUST_ENDS_WINDOW and end_tend_dist <= TRUST_ENDS_WINDOW:
                    iso_to_counts[transcript][0] += 1
                    gene_to_tot[gene][1] += 1
            else:
                tlen = final_transcript_objs[transcript].end - final_transcript_objs[transcript].start
                rlen = tlen - (start_tend_dist + end_tend_dist)
                if rlen > (tlen / 2):
                    iso_to_counts[transcript][0] += 1
                    gene_to_tot[gene][1] += 1
        else:
            # either reads are full-length or the user has explicitly specified no_stringent
            if (start_sj_index == 0 and end_sj_index == len(final_transcript_objs[transcript].juncs) - 1) or no_stringent:
                iso_to_counts[transcript][0] += 1
                gene_to_tot[gene][0] += 1
                gene_to_tot[gene][1] += 1
            else:  # this is only kept to catch bugs in count_sam_transcripts transcript assignment
                # the FIXME above says count_sam_transcripts no longer emits these, so
                # reaching this means the two disagree about what ends.tsv holds
                raise FlairError(f"{read_ends_file}: read '{read}' on transcript '{transcript}' is not "
                                 f"full length: junctions {start_sj_index} to {end_sj_index} of "
                                 f"{len(final_transcript_objs[transcript].juncs)}")
        iso_to_counts[transcript][1] += 1
        gene_to_tot[gene][2] += 1
    return iso_to_counts, gene_to_tot

def write_final_isoform_output(partition, args, final_transcript_objs, iso_to_counts, gene_to_tot, annots, genome, generate_map):
    transcript_to_reads = {}
    if generate_map:
        for line in open(partition.output_path('countsam.read.map.txt')):
            iso, reads = line.split('\t', 1)
            transcript_to_reads[iso] = reads

    with open(partition.output_path('isoforms.bed'), 'w') as iso_fh, \
         open(partition.output_path('isoforms.fa'), 'w') as seq_fh, \
         open(partition.output_path('isoform.counts.txt'), 'w') as counts_fh, \
         open(partition.output_path('isoform.read.map.txt'), 'w') as map_fh:
        for tname in final_transcript_objs:
            # spliced isos checked against spliced total, single exon checked against full-length total
            isoform = final_transcript_objs[tname]
            passes_support, my_frac_support = _iso_passes_support_filter(args, tname, isoform.gene_id, len(isoform.exons), iso_to_counts, gene_to_tot)
            if passes_support:
                if args.normalize_ends and len(isoform.exons) > 1:
                    # removing additional length from ends
                    # not entirely sure why I need to reset the exons outside of the isoform object, but it doesn't work otherwise
                    exons = isoform.exons
                    exons[0] = Exon(isoform.unpadded_ends[0], exons[0].end)
                    exons[-1] = Exon(exons[-1].start, isoform.unpadded_ends[1])
                    isoform.reset_from_exons(exons)

                thickStart, thickEnd, productivity, aaseq = predict_prod_temp(isoform, annots.start_codon_count,
                                                                              annots.gene_to_cds_starts, annots.transcript_to_nmd_except, genome)
                convert_to_flair_bed(isoform, thickStart=thickStart, thickEnd=thickEnd, read_support=iso_to_counts[tname][0],
                                     frac_support=my_frac_support, productivity=productivity, samples=(args.sample_name,), aaseq_id=aaseq).write(iso_fh)
                seq_fh.write('>' + isoform.name + '\n')
                seq_fh.write(isoform.get_sequence(genome) + '\n')
                counts_fh.write(f'{isoform.name}\t{iso_to_counts[tname][0]}\t{iso_to_counts[tname][1]}\n')
                if generate_map:
                    map_fh.write(f'{isoform.name}\t{transcript_to_reads[tname]}')


def _run_region(*, partition, gtf_data, junction_corrector, args):
    region = partition.region
    # junction chains are interned in class-level caches; with one thread every region
    # runs in this process, so they have to be dropped between regions
    ReadRec.clear_juncs_cache()
    Isoform.clear_juncs_cache()

    # first extract reads for region as fasta
    pipettor.run([('samtools', 'view', '-h', args.genome_aligned_bam, region.name + ':' + str(region.start) + '-' + str(region.end)),
                  ('samtools', 'fasta', '-')],
                 stdout=partition.output_path('reads.fasta'))
    if os.path.getsize(partition.output_path('reads.fasta')) > 0:
        _run_region_reads(partition=partition, region=region, gtf_data=gtf_data,
                          junction_corrector=junction_corrector, args=args)
    else:
        generate_empty_intermediate_files(partition.file_prefix, ['.firstpass.bed', '.isoform.counts.txt', '.isoform.read.map.txt', '.isoforms.bed', '.isoforms.fa', '.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed'])


def _run_region_reads(*, partition, region, gtf_data, junction_corrector, args):
    # FIXME confusing name if taking gtf_data,
    # FIXME: should only have region, so why take region arg
    annots = annot_data_from_gtf(gtf_data, region)

    # then align reads to transcriptome and run count_sam_transcripts
    with pysam.FastaFile(args.genome) as genome, \
         pysam.AlignmentFile(args.genome_aligned_bam, 'rb') as bam_file:
        # genomic clipping: amount of clipping (from cigar) at ends of reads when aligned to genome
        # generates file with [read{\t}clipping amount] on each line
        # For comparing with amount of clipping after alignment to transcriptome
        # in order to check whether transcriptome alignment is comparable to or better than genomic alignment,
        # which can be considered to support isoform.

        # logging.info('generating genomic clipping reference')
        # never zero reads here: _run_region only calls this when samtools fasta found
        # reads, and both skip secondary and supplementary alignments
        _, clipping_file = generate_genomic_alignment_read_to_clipping_file(partition.file_prefix, bam_file, region.name, region.start, region.end)

        # aligning to reference transcriptome, then identifying reads that match well to reference transcripts
        # with filter_transcriptome_align
        # logging.info('identifying good match to annot')
        if not args.no_align_to_annot:
            logging.info('aligning to transcriptome reference')
        read_to_annot_transcript = identify_good_match_to_annot(args, partition.file_prefix, region.name, annots, genome)

        logging.info('correcting and grouping reads, filtering isoforms')

        # takes in bam file, for each read attempts to correct splice junctions (removes unsupported ones), then groups reads by junction chains
        # this also handles read strandedness if necessary
        sj_to_ends = {}
        filter_correct_group_reads(bam_file=bam_file, region=region,
                                   read_to_annot_transcript=read_to_annot_transcript,
                                   annots=annots, junction_corrector=junction_corrector,
                                   genome=genome,
                                   quality=args.quality, keep_sup=args.keep_supplementary,
                                   sj_to_ends=sj_to_ends, trust_strand=args.trust_strand,
                                   check_motifs=not args.trust_junctions)
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
            write_firstpass(partition.file_prefix, region.name, firstpass, annots, genome, unique_bound=iso_to_unique_bound, normalize_ends=args.normalize_ends)

            # aligns to firstpass transcriptome, identifies best read -> isoform alignment for each read, then gets read counts per isoform
            read_map_file = partition.output_path('countsam.read.map.txt') if args.generate_map else None
            transcriptome_align_and_count(args, partition.output_path('reads.fasta'),
                                          partition.output_path('firstpass.fa'),
                                          partition.output_path('firstpass.bed'),
                                          partition.output_path('isoform.counts.txt'),
                                          read_map_file, False,  # say is not annot, requires stringent, returns different end values
                                          partition.output_path('reads.genomicclipping.txt'),
                                          partition.output_path('firstpass.uniquebound.txt'),
                                          args.directRNA)
        else:
            logging.info('no firstpass isoforms found')
            generate_empty_intermediate_files(partition.file_prefix, ['.firstpass.fa', '.firstpass.bed', '.isoform.counts.txt', '.countsam.read.map.txt', '.isoform.ends.tsv'])

        # FIXME this is messy, shouldn't have to reorganize like this
        final_transcript_objs = {}
        for og_key in firstpass:
            final_transcript_objs[firstpass[og_key].name] = firstpass[og_key]

        iso_to_counts, gene_to_tot = calc_final_iso_support(partition.output_path('isoform.ends.tsv'), final_transcript_objs, args.trust_ends, args.no_stringent)
        write_final_isoform_output(partition, args, final_transcript_objs, iso_to_counts, gene_to_tot, annots, genome, args.generate_map)

def combine_chunks(args, output, partitions):
    files_to_combine = ['.isoforms.bed', '.isoforms.fa', '.isoform.counts.txt']
    if args.generate_map:
        files_to_combine.append('.isoform.read.map.txt')
    if args.keep_intermediate:
        files_to_combine.extend(['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed'])
    combine_temp_files_by_suffix(output, [p.file_prefix for p in partitions], files_to_combine)

def get_new_ids(output):
    iso_hash_to_ID, gene_hash_to_ID, aaseq_to_id = {}, {}, {}
    iso_count, gene_count, aaseq_count = 1, 1, 1
    with open(output + '.isoforms.newids.bed', 'w') as fh:
        for bed_rec in BedReader(output + '.isoforms.bed', bedClass=FlairBed):
            if bed_rec.name not in iso_hash_to_ID:
                iso_id = f'FLT{iso_count:08d}'
                iso_hash_to_ID[bed_rec.name] = iso_id
                iso_count += 1
            if bed_rec.gene_id not in gene_hash_to_ID:
                gene_id = f'FLG{gene_count:08d}'
                gene_hash_to_ID[bed_rec.gene_id] = gene_id
                gene_count += 1
            if bed_rec.aaseq_id is not None and bed_rec.aaseq_id not in aaseq_to_id:
                aaseq_id = f'FLP{aaseq_count:08d}'
                aaseq_to_id[bed_rec.aaseq_id] = aaseq_id
                aaseq_count += 1
            bed_rec.name = iso_hash_to_ID[bed_rec.name]
            bed_rec.gene_id = gene_hash_to_ID[bed_rec.gene_id]
            if bed_rec.aaseq_id is not None:
                bed_rec.aaseq_id = aaseq_to_id[bed_rec.aaseq_id]
            bed_rec.write(fh)
    with open(output + '.aaseq.tsv', 'w') as fh:
        fh.write('aaseq_id\taaseq\n')
        for aaseq, id in aaseq_to_id.items():
            fh.write(f'{id}\t{aaseq}\n')
    return iso_hash_to_ID

def fix_ids_txt_file(iso_hash_to_ID, oldfile, newfile):
    with open(newfile, 'w') as fh:
        for line in open(oldfile):
            line = line.split('\t', 1)
            if line[0] in iso_hash_to_ID:  # had to add this check due to having isoforms in this file that were not in bed due to not passing final filters
                line[0] = iso_hash_to_ID[line[0]]
                fh.write('\t'.join(line))

def fix_ids_fa_file(iso_hash_to_ID, oldfile, newfile):
    with open(newfile, 'w') as fh:
        last = False
        for line in open(oldfile):
            if line[0] == '>':
                if line[1:].rstrip('\n') in iso_hash_to_ID:
                    line = '>' + iso_hash_to_ID[line[1:].rstrip('\n')] + '\n'
                    last = True
                else:
                    last = False
            if last:
                fh.write(line)

def fix_iso_labels(output, generate_map):
    iso_hash_to_ID = get_new_ids(output)
    fix_ids_txt_file(iso_hash_to_ID, output + '.isoform.counts.txt', output + '.isoform.counts.newids.txt')
    fix_ids_fa_file(iso_hash_to_ID, output + '.isoforms.fa', output + '.isoforms.newids.fa')
    pipettor.run([('mv', output + '.isoforms.newids.bed', output + '.isoforms.bed')])
    pipettor.run([('mv', output + '.isoform.counts.newids.txt', output + '.isoform.counts.txt')])
    pipettor.run([('mv', output + '.isoforms.newids.fa', output + '.isoforms.fa')])
    if generate_map:
        fix_ids_txt_file(iso_hash_to_ID, output + '.isoform.read.map.txt', output + '.isoform.read.map.newids.txt')
        pipettor.run([('mv', output + '.isoform.read.map.newids.txt', output + '.isoform.read.map.txt')])

####
# main
####

def flair_transcriptome(*, genome_aligned_bam, genome, sample_name, output, annot_gtf,
                        junction_tab, junction_bed, junction_support, ss_window, end_window,
                        sjc_support, single_exon_support, frac_support, directRNA, trust_strand,
                        trust_junctions, trust_ends, no_stringent, no_check_splice, no_align_to_annot,
                        max_ends, filter, keep_supplementary, quality, threads, parallel_mode,
                        fusion_breakpoints, keep_intermediate, temp_dir, normalize_ends, generate_map):
    args = TranscriptomeOpts(genome_aligned_bam=genome_aligned_bam, genome=genome,
                             sample_name=sample_name, output=output, annot_gtf=annot_gtf,
                             junction_tab=junction_tab, junction_bed=junction_bed,
                             junction_support=junction_support, ss_window=ss_window,
                             end_window=end_window, sjc_support=sjc_support,
                             single_exon_support=single_exon_support, frac_support=frac_support,
                             directRNA=directRNA, trust_strand=trust_strand or directRNA,
                             trust_junctions=trust_junctions, trust_ends=trust_ends,
                             no_stringent=no_stringent, no_check_splice=no_check_splice,
                             no_align_to_annot=no_align_to_annot, max_ends=max_ends,
                             filter=filter, keep_supplementary=keep_supplementary,
                             quality=quality, threads=threads, parallel_mode=parallel_mode,
                             fusion_breakpoints=fusion_breakpoints,
                             keep_intermediate=keep_intermediate, temp_dir=temp_dir, normalize_ends=normalize_ends,
                             generate_map=generate_map)

    logging.info('loading genome')
    genome_fa = pysam.FastaFile(args.genome)

    if args.keep_intermediate and args.temp_dir is None:
        raise FlairInputDataError('--keep_intermediate requires --temp_dir, the directory to keep them in')
    temp_dir = make_run_temp_dir(args.output, args.temp_dir)
    logging.info(f'temporary files in {temp_dir}')
    try:
        # partitioning, mostly flair_partition, runs while the annotation is parsed; it
        # gets one thread fewer, since the parsing uses one
        logging.info('partitioning genome')
        with ThreadPoolExecutor(max_workers=1) as executor:
            partitioning = executor.submit(partition_regions, args.parallel_mode, genome_fa, args.genome_aligned_bam,
                                           args.annot_gtf, max(1, args.threads - 1))
            annot_gtf_data = None
            if args.annot_gtf:
                logging.info('loading annotation GTF')
                annot_gtf_data = gtf_data_parser(args.annot_gtf, attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES)

            logging.info('building intron support database')
            junction_corrector = junction_corrector_factory(args.ss_window, args.junction_support,
                                                            annot_gtf_data=annot_gtf_data,
                                                            intron_beds=args.junction_bed,
                                                            star_sj_tabs=args.junction_tab)
            regions, weights = partitioning.result()

        runner = PartitionRunner(regions, temp_dir, gtf_data=annot_gtf_data, junction_corrector=junction_corrector,
                                 threads=args.threads, weights=weights)
        logging.info(f'number of partitions: {len(runner)}')

        logging.info('running partitions')
        runner.run(_run_region, args=args)
        combine_chunks(args, args.output, runner.partitions)

        #  simplify isoform and gene ID hashes in bed file, read map, counts, and fa files
        fix_iso_labels(args.output, args.generate_map)

        # index of column with gene id in extracols, then additional column indexes + names
        bed_to_gtf(args.output + '.isoforms.bed', args.output + '.isoforms.gtf', is_flair_bed=True)

        make_big_bed(genome_fa, temp_dir + 'chrom.sizes', args.output + '.isoforms')

    finally:
        # also on failure, so runs don't leave their files in $TMPDIR
        if args.keep_intermediate:
            logging.info(f'intermediate files kept in {temp_dir}')
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)

    genome_fa.close()
