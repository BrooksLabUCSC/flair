#!/usr/bin/env python3
import argparse
import os
import re
import pipettor
import shutil
import pysam
import logging
import math
from bisect import bisect_left, bisect_right
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from statistics import median
from collections import Counter
from flair import FlairError, FlairInputDataError, resolve_deprecated_option
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.junction_correct import junction_corrector_factory
from flair.annotation_alignment import (_exon_bounds, _write_reads_subset,
                                        identify_good_match_to_annot, normalize_gene_terminal_exons,
                                        transcriptome_align_and_count)
from flair.count_sam_transcripts import MINIMAP2_SECONDARY_ALIGNMENTS
from flair.direct_assignment import (covers_unique_bounds, direct_assignments, parse_unique_bounds, read_unique_bounds,
                                     unassignable_reads, write_direct_assignments)
from flair.partition_runner import parallel_mode_parse, PartitionRunner, partition_regions, combine_temp_files_by_suffix
from flair.io_utils import check_genome_fasta, make_run_temp_dir, open_sorted_bed
from flair.bed_to_gtf import bed_to_gtf
from flair.isoform_data import (Exon, Isoform, ReadRec, convert_to_bed12, convert_to_flair_bed, make_big_bed)
from flair.gene_identification import build_genes, _gene_span
from flair.read_processing import generate_genomic_alignment_read_to_clipping_file
from flair.read_correction import filter_correct_group_reads
from flair.annotation_data import annot_data_from_gtf
from flair.pycbio.hgdata.bed import BedReader
from flair.predictProductivity import predict_productivity
from flair.flair_bed import FlairBed
from flair.aaseq_tsv import write_aaseqs
from flair.isoform_counts_tsv import IsoformCountsWriter
from flair.read_ends_tsv import ReadEndsReader, ReadEndsWriter
from flair.read_map_tsv import ReadMapReader, ReadMapWriter
from flair.output_compress import compressed_path, open_fasta_output
from flair.pycbio.sys import fileOps
from flair.terminal_exon_ends import TerminalExonEnds
from flair.interval_index import IntervalIndex

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

# --filter: one of these, or <N>X, keeping a subset isoform with more than N times
# the reads of the isoforms it is a subset of
FILTER_MODES = ('nosubset', 'comprehensive', 'ginormous')
DEFAULT_FILTER = 'nosubset'
# the --filter modes keeping every subset isoform, so a truncated read is aligned alike
# to its own isoform and every isoform containing it, and the secondary alignments to
# them minimap2 reports in the realignment, for the read's own to be among them: with
# the default 4, of a KIF5B--RET fusion's 123 reads, mostly 5' truncated, 74 were
# assigned among its 20 isoforms, and 117 with 50
KEEP_SUBSET_FILTERS = ('comprehensive', 'ginormous')
KEEP_SUBSET_SECONDARY_ALIGNMENTS = 50
_SUPPORT_RATIO_FILTER_RE = re.compile(r'^(\d+(?:\.\d+)?)[xX]$')


def subset_support_ratio(filter_mode):
    """For a --filter of <N>X, N, the ratio of reads by which a subset isoform must
    exceed the isoforms it is a subset of to be kept; None for the named modes"""
    match = _SUPPORT_RATIO_FILTER_RE.match(filter_mode)
    return float(match.group(1)) if match else None


def parse_filter(filter_mode):
    "validate a --filter value"
    if filter_mode in FILTER_MODES:
        return filter_mode
    ratio = subset_support_ratio(filter_mode)
    if ratio is None or ratio <= 0:
        raise FlairInputDataError(f"invalid --filter '{filter_mode}': expected one of {', '.join(FILTER_MODES)}, "
                                  "or a ratio such as 10X")
    return filter_mode


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
    junction_star: str
    junction_bed: str
    junction_support: int
    ss_window: int
    end_window: int
    sjc_support: int
    single_exon_support: int
    subset_backup_support: int
    frac_support: float
    directRNA: bool
    trust_strand: bool
    total_rna: bool
    trust_junctions: bool
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
                          help='FASTA of reference genome, can be minimap2 indexed; a compressed genome '
                               'must be bgzip and indexed with samtools faidx')
    required.add_argument('--sample_name', required=True,
                          help='name of sample - this will be added as metadata to your output files')
    parser.add_argument('-o', '--output',
                        help='prefix for the output file names; defaults to the '
                             'required --sample_name')

    parser.add_argument('-f', '--gtf', dest="annot_gtf", default=None,
                        help='GTF annotation file, used for identifying annotated isoforms')
    parser.add_argument('--junction_star', help='splice junctions in STAR SJ.out.tab format, as STAR writes '
                                                'when it aligns reads')
    parser.add_argument('--junction_tab', help='deprecated name for --junction_star')
    parser.add_argument('--junction_bed', help='splice junctions as BED6 to BED9 with the number of supporting '
                                               'reads in the score column, as intron-prospector writes from a BAM '
                                               'of either short or long reads.  An alternative format to '
                                               '--junction_star, not a different kind of evidence')
    parser.add_argument('--junction_support', type=int, default=2,
                        help='minimum number of supporting reads to keep a junction, from either --junction_star '
                             'or --junction_bed; for bed the score column holds the count (default: %(default)s)')

    parser.add_argument('--ss_window', type=int, default=10,
                        help='window size for correcting splice sites (default: %(default)s)')
    parser.add_argument('--end_window', type=int, default=100,
                        help='window size for comparing TSS/TES (default: %(default)s)')

    parser.add_argument('--sjc_support', type=int, default=1,
                        help='minimum number of supporting reads for a spliced isoform (default: %(default)s)')
    parser.add_argument('--single_exon_support', type=int, default=3,
                        help='minimum number of supporting reads for a single exon isoform (default: %(default)s)')
    parser.add_argument('--subset_backup_support', type=int, default=10,
                        help='a spliced isoform removed as a subset of other isoforms is reported after all, '
                             'with its reads, if every isoform it is a subset of fails support and it has at '
                             'least this many reads; 0 never reports it (default: %(default)s)')
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
    parser.add_argument('--total_rna', action='store_true',
                        help='the library is total RNA rather than poly(A) selected: single-exon reads ending '
                             'in or before genomic A runs are kept, rather than removed as internally primed, '
                             'and no poly(A) tails are needed.  Requires --trust_strand')
    parser.add_argument('--trust_junctions', action='store_true',
                        help='trust the strand of every input splice junction for read correction, without '
                             'checking splice motifs.  By default, an unannotated junction without a GT-AG '
                             'motif, such as GC-AG or AT-AC, is only weak evidence of strand: a read stranded '
                             'only by such junctions gets its strand from gene identification, or is dropped')
    # removed: each isoform's ends come from the reads assigned to it
    parser.add_argument('--trust_ends', action='store_true', help=argparse.SUPPRESS)

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
                        help='maximum number of TSS/TES picked per spliced junction chain; make higher for more '
                             'precise end detection.  Single-exon isoforms are reported at each of their clustered '
                             'ends, however many (default: %(default)s)')
    parser.add_argument('--filter', default=DEFAULT_FILTER,
                        help='which subset isoforms (a contiguous part of another isoform) to keep. '
                             '<N>X, such as 10X: a subset isoform is kept if it has more than N times the reads '
                             'of the isoforms it is a subset of; '
                             'nosubset: all subset isoforms are removed; '
                             'comprehensive: all spliced subset isoforms are kept; '
                             'ginormous: comprehensive, plus single exon subset isoforms '
                             '(default: %(default)s)')

    parser.add_argument('--keep_supplementary', action='store_true',
                        help='keep supplementary alignments when defining isoforms')
    parser.add_argument('--quality', default=1, type=int,
                        help='minimum mapping quality threshold to consider genomic alignments for '
                             'defining transcripts (default: %(default)s)')

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
                        help='report each spliced isoform with the ends it shares with the isoforms of its gene '
                             'whose terminal exons have splice sites near its own: the furthest of their read '
                             'ends, rather than the ends of the reads assigned to it.  Each spliced junction '
                             'chain gives a single isoform, overriding --max_ends')
    parser.add_argument('--generate_map', action='store_true',
                        help='generate a txt file of read-isoform assignments')
    parser.set_defaults(entry=transcriptome_cmd)

def transcriptome_cmd(args):
    resolve_deprecated_option(args, deprecated="junction_tab", current="junction_star")
    parse_filter(args.filter)
    if args.trust_ends:
        raise FlairInputDataError("--trust_ends has been removed: each isoform's ends come from the ends "
                                  "of the reads assigned to it")
    if args.keep_intermediate and args.temp_dir is None:
        raise FlairInputDataError('--keep_intermediate requires --temp_dir, the directory to keep them in')
    if args.total_rna and not (args.trust_strand or args.directRNA):
        raise FlairInputDataError("--total_rna requires --trust_strand: without poly(A) tails, the alignment strand is "
                                  "the only strand evidence for single-exon reads")
    if not os.path.exists(args.genome_aligned_bam):
        raise FlairInputDataError(f'aligned reads bam file does not exist: {args.genome_aligned_bam}')
    check_genome_fasta(args.genome)
    flair_transcriptome(genome_aligned_bam=args.genome_aligned_bam, genome=args.genome,
                        sample_name=args.sample_name,
                        output=args.output if args.output is not None else args.sample_name,
                        annot_gtf=args.annot_gtf, junction_star=args.junction_star,
                        junction_bed=args.junction_bed, junction_support=args.junction_support,
                        ss_window=args.ss_window, end_window=args.end_window,
                        sjc_support=args.sjc_support,
                        single_exon_support=args.single_exon_support,
                        subset_backup_support=args.subset_backup_support,
                        frac_support=args.frac_support, directRNA=args.directRNA, trust_strand=args.trust_strand, total_rna=args.total_rna,
                        trust_junctions=args.trust_junctions, no_stringent=args.no_stringent,
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

# tolerance for terminal exon boundary comparisons
TERMINAL_EXON_BOUNDARY_TOLERANCE = 20

# a spliced junction chain's best supported start and end
# are the read start and end with the most others within this distance
BEST_END_WINDOW = 25
# a single-exon end variant's least share of its cluster's reads: whether a cluster
# of read ends is an end of the transcript or noise among its own reads, unlike
# frac_support, an isoform's share of its gene's reads (_single_exon_end_variants)
SINGLE_EXON_END_MIN_FRAC = 0.05

# margin for single-exon isoform overlap comparisons
SINGLE_EXON_OVERLAP_MARGIN = 10

# expression ratio threshold for filtering overlapping single-exon isoforms with
# --filter comprehensive; with <N>X, N is used
SINGLE_EXON_EXPRESSION_RATIO = 1.2


####
# transcriptome alignment
####
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


def get_isos_with_similar_juncs(juncs, junc_to_names):
    """Find isoforms sharing junctions with the given junction set."""
    novel_isos = set()
    for j in juncs:
        if junc_to_names and j in junc_to_names:
            novel_isos.update(junc_to_names[j])
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
    A boundary is recorded whatever the terminal exon's end, as the transcript
    built may have different ends from those checked here; unique_bounds_past
    keeps the ones a transcript's terminal exons extend past."""
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


def unique_bounds_past(unique_seq_bound, first_exon, last_exon, strand):
    """The unique sequence boundaries, from the subset checks, that a transcript's
    terminal exons extend past, so it has unique sequence there for reads to
    cover, in count_sam_transcripts' '<side>_<distance from splice site>' form
    with side 0 the transcript's 5' end"""
    bounds = []
    for side, dist in sorted(set(unique_seq_bound)):
        exon = first_exon if side == 0 else last_exon
        if exon.end - exon.start > dist:
            if strand == '-':
                side = 1 - side
            bounds.append(f'{side}_{dist}')
    return bounds

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

def _matches_annotated_ends(annots, strand, exons, truncated, end_window, ss_window):
    """Does a subset isoform end, within end_window, where an annotated transcript
    has confirmed ends.  Truncated on one side, the other sharing a superset's
    terminal exon, it can match just that side's end of a basic transcript with the
    same terminal exon, its splice site within ss_window of the subset's: the end of
    a transcript with another terminal exon, or of a non-basic one, is often a
    fragment's, ending inside the full-length isoform's exon.  Otherwise, or truncated on both sides,
    it has to match both ends of one transcript."""
    if truncated[0] != truncated[1]:
        side, splice_site, pos = (0, exons[0].end, exons[0].start) if truncated[0] else (1, exons[-1].start, exons[-1].end)
        if annots.has_transcript_end(strand, side, splice_site, pos, end_window, ss_window):
            return True
    return annots.has_transcript_ends(strand, exons[0].start, exons[-1].end, end_window)

def _juncs_outside_gene(juncs, gene, annots):
    """Do the junctions reach outside the span of the gene's annotated junctions
    (_gene_span), or of its genes' together, for a readthrough.  A novel gene, or
    one without annotated junctions, has no span to be within."""
    gene_ids = [g for g in (gene.gene_desc if gene is not None else ())
                if g in getattr(annots, 'gene_to_annot_juncs', {})]
    if not gene_ids:
        return True
    spans = [_gene_span(g, annots) for g in gene_ids]
    return juncs[0].start < min(s for s, e in spans) or juncs[-1].end > max(e for s, e in spans)

def filter_spliced_iso(filter_type, support, juncs, exons, name, score, annots,
                       junc_to_names, all_isoforms,
                       strand, end_window=0, ss_window=0, supersets=None):
    """Is a spliced candidate kept by the subset filter, and its unique sequence
    boundaries.  It is a subset of isoforms of its own gene, and of other genes'
    only if its junctions reach outside its own gene (_juncs_outside_gene).  An
    isoform of any gene containing it still decides its unique sequence boundaries,
    which come from any isoform with its junctions, as for a subset: none if its
    ends match an annotated transcript's.  The names of the isoforms it is a subset
    of are added to supersets, if given."""
    assert isinstance(exons[0], Exon)  # FIXME: debugging
    isoform = all_isoforms[name]
    any_gene_superset = _juncs_outside_gene(juncs, isoform.gene, annots)

    novel_isos = get_isos_with_similar_juncs(juncs, junc_to_names)
    first_exon, last_exon = exons[0], exons[-1]
    superset_support = []  # of the isoforms this is a subset of
    contained = False  # in an isoform of any gene
    truncated = [False, False]  # its first or last exon is inside a containing isoform's internal exon
    unique_seq_bound = []
    for novel_iso_id in novel_isos:
        # only a superset with enough support to be reported itself can make this
        # a subset.  A full-length isoform still removes its truncated subsets
        # however many reads they have, but one with too few reads to keep can't
        # take a well-supported isoform with it: its reads, which don't cover the
        # superset's extra exons, would be left unassigned
        if novel_iso_id != name and all_isoforms[novel_iso_id].score >= support:
            # this is a subset only of an isoform containing both its first and its
            # last exon: with the two from different isoforms, no one isoform
            # contains it.  Unique sequence boundaries come from any isoform with its
            # junctions, whose truncated reads could otherwise look like this one's
            terminal_exon_is_subset = [0, 0]  # first exon is a subset, last exon is a subset
            _check_novel_iso_subset(novel_iso_id, all_isoforms,
                                    juncs, first_exon, last_exon, terminal_exon_is_subset,
                                    [], unique_seq_bound)
            if sum(terminal_exon_is_subset) == 2:
                otheriso = all_isoforms[novel_iso_id]
                contained = True
                # its junctions are a contiguous run of the other's: a side short
                # of the other's first or last junction has an internal exon there
                first = otheriso.juncs.index(juncs[0])
                truncated[0] |= first > 0
                truncated[1] |= first + len(juncs) < len(otheriso.juncs)
                # an isoform within its gene's junctions isn't taken to be a fragment
                # of another gene's isoform, as of a readthrough or of an overlapping
                # gene containing it, while one reaching past them, as a fragment
                # spanning a readthrough junction, is.  A subset kept or promoted
                # (promote_backup_subsets) keeps the gene it was assigned
                if any_gene_superset or otheriso.gene_id == isoform.gene_id:
                    superset_support.append(otheriso.score)
                    if supersets is not None:
                        supersets.append(novel_iso_id)
    # the transcript's ends aren't known yet, so these are all the boundaries,
    # for unique_bounds_past to choose from once they are
    unique_seq_bound = sorted(set(unique_seq_bound))

    if not contained:
        return True, unique_seq_bound
    elif _matches_annotated_ends(annots, strand, exons, truncated, end_window, ss_window):
        # a contained isoform whose ends match an annotated transcript's confirmed
        # ends is a known isoform, not a fragment of a longer one.  It gets no unique sequence
        # boundaries: its ends are within TERMINAL_EXON_BOUNDARY_TOLERANCE of the
        # containing isoform's exons, so its reads needn't cover sequence past them, which its
        # transcript may only have from end normalization
        return True, []
    elif len(superset_support) == 0:
        # contained only in other genes' isoforms
        return True, unique_seq_bound
    else:
        ratio = subset_support_ratio(filter_type)
        if ratio is not None and score >= support and score > max(superset_support) * ratio:
            return True, unique_seq_bound
    # a removed subset's boundaries are kept for promote_backup_subsets
    return False, unique_seq_bound

def max_terminal_exons_ends_from_iso_infos(iso_to_info):
    max_terminal_exons_ends = TerminalExonEnds()
    for iso_name in iso_to_info:
        isoform = iso_to_info[iso_name]
        max_terminal_exons_ends.add_transcript(isoform.gene_id, isoform.strand, _exon_bounds(isoform.exons))
    return max_terminal_exons_ends

def rank_end_variants(isoforms):
    """Sort a junction chain's or single-exon cluster's end variants, best first;
    with max_ends above 1, the best max_ends of those with the support are kept"""
    # First by read support, then by length
    # FIXME: the comment above is what was meant; the code multiplies instead, so
    # this is not a two-level sort and a 1-read 10kb candidate outranks a 10-read
    # 500bp one.  It once also chose the ends of every single-exon isoform, where
    # sorting by (num_reads, genomic_length) picked the modal end pair, for long
    # reads a 5' truncation cluster, over the pair reproducing the annotated 3' end
    isoforms.sort(key=lambda x: [x.num_reads * x.genomic_length], reverse=True)

def _densest_end_clustered(positions, window, outer):
    """densest_end, and whether it is a cluster of read ends: another read end
    within window of it"""
    positions = sorted(positions)
    exact = Counter(positions)

    def support(pos):
        return bisect_right(positions, pos + window) - bisect_left(positions, pos - window)
    best = max(support(pos) for pos in exact)
    if best == 1:
        return outer(positions), False
    densest = [pos for pos in exact if support(pos) == best]
    most = max(exact[pos] for pos in densest)
    modal = sorted(pos for pos in densest if exact[pos] == most)
    middle = len(modal) // 2
    return (modal[middle] if len(modal) % 2 else outer(modal[middle - 1], modal[middle])), True

def densest_end(positions, window, outer):
    """The read end with the most read ends within window of it.  In a tight
    cluster of read ends, every position in it has about as many, so ties go to the
    position the most reads end at exactly, then to the middle of the tied
    positions, and only between two middle ones to the outer: outer is min for
    starts and max for ends.  Going to the outer one first moved ends past the
    cluster's mode, out of polyA and CAGE peaks.  If no read end has another within
    window, there is no cluster, and the outer end is taken."""
    return _densest_end_clustered(positions, window, outer)[0]

def normalize_chain_ends(isoform):
    """A spliced junction chain is a single first-pass isoform, ending at the
    furthest start and end of its reads.  Its densest read ends are its best
    supported ends, which the subset check uses, and whether its 3' one is a
    cluster of read ends, which the 3' UTR fragment check uses (_ThreePrimeFragments)."""
    starts, ends = isoform.starts, isoform.ends
    best_start, start_clustered = _densest_end_clustered(starts, BEST_END_WINDOW, min)
    best_end, end_clustered = _densest_end_clustered(ends, BEST_END_WINDOW, max)
    isoform.best_ends = (best_start, best_end)
    isoform.three_prime_clustered = end_clustered if isoform.strand == '+' else start_clustered
    isoform.start, isoform.end = min(starts), max(ends)
    return isoform

def _write_unfiltered_ends(isoforms, fh):
    for iso_readrec in isoforms:
        convert_to_bed12(iso_readrec).write(fh)

class CandidateIsoforms:
    """Candidate isoforms before filtering, with junction and exon indices.

    isoforms: dict of isoform_name -> Isoform
    junc_to_names: dict of Junc -> set of isoform_names sharing that junction
    exons: dict of strand -> set of Exon (named exons for SE, all exons for
           spliced, with the ends the subset check uses)
    """
    def __init__(self):
        self.isoforms = {}
        self.junc_to_names = {}
        self.exons = {}

    def add(self, isoform):
        self.isoforms[isoform.name] = isoform
        strand_exons = self.exons.setdefault(isoform.strand, set())
        if isoform.juncs == ():
            strand_exons.add(Exon(isoform.start, isoform.end, isoform.name))
        else:
            for j in isoform.juncs:
                if j not in self.junc_to_names:
                    self.junc_to_names[j] = set()
                self.junc_to_names[j].add(isoform.name)
            # single-exon isoforms on the same strand are checked against these for
            # being subsets
            strand_exons.update(subset_check_exons(isoform))


def _generate_candidate_isos(args, isoform, candidates, iso_fh, iso_unfilt_fh):
    if isoform.juncs != ():
        # no end groups: one isoform, whatever max_ends is; assign_final_ends splits
        # it into end variants
        these_firstpass = [normalize_chain_ends(isoform)]
        support = args.sjc_support
    else:
        # single-exon isoforms aren't realigned to, so a cluster's end variants are
        # final here, each with its own reads; assign_single_exon_isoforms gives a
        # cluster none of whose variants passes support one isoform instead
        # NOTE: Harrison's TED code will be slotted in here to find the end clusters
        these_firstpass = [variant for variant, _ in _single_exon_end_variants(isoform, isoform.reads, support=args.single_exon_support,
                                                                               min_frac=SINGLE_EXON_END_MIN_FRAC)]
        cluster = (isoform.chrom, isoform.strand, min(isoform.starts), max(isoform.ends))
        for variant in these_firstpass:
            variant.end_variant_cluster = cluster
        support = args.single_exon_support
    _write_unfiltered_ends(these_firstpass, iso_unfilt_fh)
    num_reads = sum(variant.num_reads for variant in these_firstpass)
    if num_reads < support:
        logging.debug(f"isoform group dropped: insufficient support ({num_reads} < {support}): {isoform.chrom}:{isoform.start}-{isoform.end}")
        return
    for variant in these_firstpass:
        candidates.add(variant)
        convert_to_bed12(variant).write(iso_fh)


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

def _internally_primed(read, trust_strand):
    """Is a single-exon read internally primed at an end that may be its 3' end
    (ReadRec.intprim, from isoform_data.internally_primed_end).  With trust_strand,
    only the read's 3' end is checked; otherwise either end, its strand not yet
    known: priming at its start is priming on the - strand."""
    if read.intprim is None:
        return False
    left, right = read.intprim
    if trust_strand:
        return right if read.strand == '+' else left
    return left or right

def _se_read_groups(reads, se_support, trust_strand, total_rna=False):
    """(strand, reads) for each overlap group of single-exon reads.  With
    trust_strand, the reads' strands are their transcripts', and each strand's reads
    are grouped separately, so sense and antisense transcripts at a locus are both
    kept.  Otherwise a group is of reads from both strands, and its strand comes
    from the poly(A) tails of its reads, None if they don't give one.  Internally
    primed reads are left out first (_internally_primed): in a poly(A) selected
    library they are reads of no polyadenylated transcript end, which would otherwise
    make single-exon isoforms in introns and UTRs.  With total_rna, the library
    isn't poly(A) selected, and they are kept."""
    if not total_rna:
        reads = [read for read in reads if not _internally_primed(read, trust_strand)]
    if trust_strand:
        by_strand = {}
        for read in reads:
            by_strand.setdefault(read.strand, []).append(read)
        return [(strand, read_group) for strand, strand_reads in sorted(by_strand.items())
                for read_group in _group_se_reads_by_overlap(strand_reads)]
    return [(correct_se_strand_polyA(read_group, se_support), read_group)
            for read_group in _group_se_reads_by_overlap(reads)]

def group_se_by_overlap(chrom, isoform, se_support, trust_strand, total_rna=False):
    for new_strand, read_group in _se_read_groups(isoform.reads, se_support, trust_strand, total_rna):
        # filter out single exon groups that fail stranding
        if new_strand is None:
            logging.debug(f"single-exon group dropped: strand could not be determined ({len(read_group)} reads): {chrom}:{read_group[0].start}-{read_group[-1].end}")
        else:
            new_key = (chrom, median([x.start for x in read_group]), median([x.end for x in read_group]), (), new_strand)
            yield new_key, new_strand, read_group

class IsoformOverlapGroups:
    """Isoforms grouped by junction chain with overlap-clustered ends.

    Key: (chrom, median_start, median_end, juncs) for spliced isoforms, and
    (chrom, median_start, median_end, (), strand) for single-exon ones.
    Value: Isoform.

    For spliced isoforms, strand is determined during junction correction and
    stored on the Isoform.  For single-exon reads, strand cannot be determined
    from junctions: with trust_strand, each strand's reads are grouped by
    coordinate overlap separately; otherwise overlapping reads are grouped
    regardless of strand, and each group's strand is resolved by polyA
    consensus (see _se_read_groups).
    """
    def __init__(self):
        self._groups = {}

    def add_spliced(self, chrom, juncs, isoform):
        """Add a spliced isoform, keyed by chrom, median ends, and junction chain."""
        self._groups[(chrom, median(isoform.starts), median(isoform.ends), juncs)] = isoform

    def add_se_overlap_groups(self, chrom, isoform, se_support, trust_strand, total_rna=False):
        """Split single-exon reads into overlap groups and resolve strand."""
        for new_key, new_strand, read_group in group_se_by_overlap(chrom, isoform, se_support, trust_strand, total_rna):
            self._groups[new_key] = Isoform.regroup(isoform, newreads=read_group, newstrand=new_strand)

    def __iter__(self):
        return iter(self._groups)

    def __getitem__(self, key):
        return self._groups[key]

    def remove(self, key):
        del self._groups[key]

    def items(self):
        return self._groups.items()


def group_by_overlap(sj_to_ends, se_support, trust_strand, total_rna=False):
    groups = IsoformOverlapGroups()
    for (chrom, juncs), isoform in sj_to_ends.items():
        if len(juncs) > 0:
            groups.add_spliced(chrom, juncs, isoform)
        else:
            groups.add_se_overlap_groups(chrom, isoform, se_support, trust_strand, total_rna)
    return groups


def process_juncs_to_firstpass_isos(args, temp_prefix, sj_to_ends, annots, region_chrom):
    sjc_with_overlap_groups = group_by_overlap(sj_to_ends, args.single_exon_support, args.trust_strand, args.total_rna)
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
    """Check if a single-exon isoform passes filtering against its overlap group.
    The end variants of one cluster of reads aren't compared with each other: they
    were chosen as its ends (_end_variants)."""
    isoform = all_isoforms[single_exon.name]
    expression_comp_with_superset = []
    is_contained = False
    for exon in curr_group:
        if isoform.end_variant_cluster is not None and exon.name != '' and \
                all_isoforms[exon.name].end_variant_cluster == isoform.end_variant_cluster:
            continue
        if exon != single_exon:
            if ((exon.start - SINGLE_EXON_OVERLAP_MARGIN) <= single_exon.start and
                    single_exon.end <= (exon.end + SINGLE_EXON_OVERLAP_MARGIN)):
                if exon.name == '' or args.filter == 'nosubset':  # is exon from spliced transcript
                    is_contained = True
                    break  # filter out
                else:  # is other single exon - check relative expression, by the --filter ratio
                    other_score = all_isoforms[exon.name].score
                    ratio = subset_support_ratio(args.filter) or SINGLE_EXON_EXPRESSION_RATIO
                    if isoform.score >= args.sjc_support and other_score * ratio < isoform.score:
                        expression_comp_with_superset.append(True)
                    else:
                        expression_comp_with_superset.append(False)
    return not is_contained and all(expression_comp_with_superset)


def filter_single_exon_group(args, curr_group, all_isoforms, firstpass, utr_fragment):
    """Filter single-exon isoforms in an overlap group against spliced exons."""
    for exon in curr_group:
        if exon.name != '':  # is single exon with name
            if utr_fragment(all_isoforms[exon.name]):
                logging.debug(f"single-exon isoform dropped: starts in a 3' terminal exon: {exon.name} ({all_isoforms[exon.name].num_reads} reads)")
            elif filter_single_exon_iso(args, exon, curr_group, all_isoforms):
                firstpass[exon.name] = all_isoforms[exon.name]
            else:
                logging.debug(f"single-exon isoform dropped: contained or low expression: {exon.name} ({all_isoforms[exon.name].num_reads} reads)")
    return firstpass


class _ThreePrimeFragments:
    """Which of a strand's single-exon isoforms are 3' UTR fragments of a spliced
    transcript: starting inside its 3' terminal exon, within BEST_END_WINDOW, as
    far as a cluster of read ends reaches, a fragment of its 3' UTR or of an
    extension of it, rather than a transcript of their own.

    The 3' terminal exons are those of the annotated spliced transcripts, which
    catch fragments of a long annotated 3' UTR the sample's spliced isoforms don't
    reach, and of the kept spliced isoforms whose 3' end is a cluster of read ends
    (Isoform.three_prime_clustered): a few reads running through a neighbouring
    transcript don't make it a fragment.  An isoform's are taken to its best
    supported ends, as the containment check takes them (subset_check_exons).  A
    single-exon isoform matching an annotated single-exon transcript, overlapping
    more than half of each, is a transcript of its own."""
    def __init__(self, spliced, annots, strand):
        self.strand = strand
        self.terminal = IntervalIndex()
        self.annotated_single = IntervalIndex()
        exon_lists = [subset_check_exons(isoform) for isoform in spliced.values()
                      if isoform.juncs != () and isoform.strand == strand and isoform.three_prime_clustered]
        if annots is not None:
            exon_lists.extend(annots.transcript_to_exons[(transcript_id, gene_id)]
                              for transcript_id, gene_id, tx_strand in annots.transcripts if tx_strand == strand)
            for exon in annots.all_annot_se[strand]:
                self.annotated_single.add(exon.start, exon.end, exon)
        for exons in exon_lists:
            if len(exons) > 1:
                exon = exons[-1] if strand == '+' else exons[0]
                self.terminal.add(exon.start, exon.end, None)

    def _matches_annotated_single_exon(self, isoform):
        for annot in self.annotated_single.overlap(isoform.start, isoform.end):
            overlap = min(isoform.end, annot.end) - max(isoform.start, annot.start)
            if overlap > (isoform.end - isoform.start) / 2 and overlap > (annot.end - annot.start) / 2:
                return True
        return False

    def __call__(self, isoform):
        five_prime = isoform.start if isoform.strand == '+' else isoform.end - 1
        return (len(self.terminal.overlap(five_prime, five_prime + 1, slack=BEST_END_WINDOW)) > 0
                and not self._matches_annotated_single_exon(isoform))

def filter_all_single_exon(args, sorted_exons, all_isoforms, firstpass, utr_fragment):
    """Group exons by overlap and filter single-exon isoforms."""
    last_end = 0
    curr_group = []

    for exon in sorted_exons:
        if exon.start < last_end:
            curr_group.append(exon)
        else:
            if len(curr_group) > 0:
                firstpass = filter_single_exon_group(args, curr_group, all_isoforms, firstpass, utr_fragment)
            curr_group = [exon]
        if exon.end > last_end:
            last_end = exon.end
    if len(curr_group) > 0:
        firstpass = filter_single_exon_group(args, curr_group, all_isoforms, firstpass, utr_fragment)

    return firstpass


def subset_check_exons(isoform):
    """An isoform's exons with its best supported ends, which is what the subset
    check compares; these differ from its furthest ends"""
    exons = list(isoform.exons)
    if isoform.best_ends is not None:
        exons[0] = Exon(isoform.best_ends[0], exons[0].end)
        exons[-1] = Exon(exons[-1].start, isoform.best_ends[1])
    return exons

def filter_firstpass_isos(args, candidates, annots, removed_subsets=None):
    """Filter candidate isoforms by subset/support criteria.
    Returns (firstpass dict, iso_to_unique_bound dict).  Each spliced candidate removed
    as a subset is added to removed_subsets, if given, by name, with the candidates
    it is a subset of and its unique sequence boundaries (promote_backup_subsets)."""
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
                    supersets = []
                    is_not_subset, unique_seq = filter_spliced_iso(args.filter, args.sjc_support, isoform.juncs, subset_check_exons(isoform),
                                                                   iso_name, isoform.num_reads, annots,
                                                                   candidates.junc_to_names, candidates.isoforms,
                                                                   isoform.strand,
                                                                   end_window=args.end_window, ss_window=args.ss_window,
                                                                   supersets=supersets)
                    if not is_not_subset:
                        logging.debug(f"isoform dropped: subset of another isoform: {iso_name} ({isoform.num_reads} reads)")
                        if removed_subsets is not None:
                            removed_subsets[iso_name] = (isoform, [candidates.isoforms[name] for name in supersets],
                                                         unique_seq)
                    else:
                        firstpass[iso_name] = isoform
                        if len(unique_seq) > 0:
                            iso_to_unique_bound[iso_name] = unique_seq
        # HANDLE SINGLE EXONS SEPARATELY - group first - one traversal of each strand's
        # exons, a single-exon isoform only compared with isoforms on its strand.  3' UTR
        # fragments are of the spliced isoforms kept above, not of candidates the
        # filter dropped
        spliced_kept = dict(firstpass)
        for strand in sorted(candidates.exons):
            firstpass = filter_all_single_exon(args, sorted(candidates.exons[strand]), candidates.isoforms, firstpass,
                                               _ThreePrimeFragments(spliced_kept, annots, strand))

    return firstpass, iso_to_unique_bound

def write_first_pass_isoforms(iso_name, isoform, max_terminal_exons_ends, unique_bound, unique_fh, iso_fh, seq_fh, genome):
    # the exons the subset check compared, with the isoform's best supported ends,
    # before normalization moves them to the gene's furthest
    check_exons = subset_check_exons(isoform)
    # FIXME: do normalization outside of write function
    if len(isoform.exons) > 1:  # don't normalize ends for single exon transcripts
        # kept so the padding can be removed from the final isoform, since clamping
        # to the chromosome means it is not always NORM_END_EXTRA_LEN
        isoform.unpadded_ends = max_terminal_exons_ends.furthest_ends(isoform.gene_id, isoform.strand, _exon_bounds(isoform.exons))
        exons = normalize_gene_terminal_exons(max_terminal_exons_ends, isoform.gene_id, isoform.strand, isoform.exons,
                                              genome.get_reference_length(isoform.chrom))
        isoform.reset_from_exons(exons)

    if unique_bound and iso_name in unique_bound:
        # only boundaries the isoform's own best supported ends reach past, as the
        # subset check saw them: one past them can only be covered by reads running
        # on past the isoform's end, as to the gene's furthest end, which isn't its
        # own, so it would only take the isoform's reads from it
        bounds = unique_bounds_past(unique_bound[iso_name], check_exons[0], check_exons[-1], isoform.strand)
        if len(bounds) > 0:
            unique_fh.write(isoform.name + '\t' + ','.join(bounds) + '\n')

    convert_to_bed12(isoform).write(iso_fh)
    if len(isoform.exons) > 1:  # single-exon isoforms get their reads without realignment
        seq_fh.write('>' + isoform.name + '\n')
        seq_fh.write(isoform.get_sequence(genome) + '\n')

def write_firstpass(temp_prefix, chrom, firstpass, annots, genome, *, unique_bound=None):
    # generating standardized set of ends for gene
    max_terminal_exons_ends = max_terminal_exons_ends_from_iso_infos(firstpass)

    with (open(temp_prefix + '.firstpass.bed', 'w') as iso_fh,
          open(temp_prefix + '.firstpass.fa', 'w') as seq_fh,
          open(temp_prefix + '.firstpass.uniquebound.tsv', 'w') as unique_fh):
        for iso_name in firstpass:
            write_first_pass_isoforms(iso_name, firstpass[iso_name], max_terminal_exons_ends,
                                      unique_bound, unique_fh, iso_fh, seq_fh, genome)


####
# results output
####

def _iso_passes_support_filter(args, iso, gene, num_exons, iso_to_counts, gene_to_tot, longest_supported=False):
    """does an isoform pass support, and its fraction of its gene's reads; a
    single-exon cluster's longest supported variant (_single_exon_end_variants)
    needs single_exon_support alone in a gene without spliced reads, where all of
    the gene's reads, frac_support's denominator, are its own and its fragments'"""
    if iso not in iso_to_counts:
        return False, 0
    else:
        count = iso_to_counts[iso][0]
        if num_exons > 1:
            return (count >= args.sjc_support) and (count / gene_to_tot[gene][0]) >= args.frac_support, (count / gene_to_tot[gene][0])
        else:
            frac = count / gene_to_tot[gene][1]
            exempt = longest_supported and gene_to_tot[gene][0] == 0
            return (count >= args.single_exon_support) and (exempt or frac >= args.frac_support), frac


EMPTY_TSV_WRITERS = {'.isoform.counts.tsv': IsoformCountsWriter,
                     '.isoform.read.map.tsv': ReadMapWriter,
                     '.countsam.read.map.tsv': ReadMapWriter,
                     '.isoform.ends.tsv': ReadEndsWriter}

def generate_empty_intermediate_files(file_prefix, suffixes):
    for s in suffixes:
        writer_class = EMPTY_TSV_WRITERS.get(s)
        if writer_class is None:
            open(file_prefix + s, 'w').close()
        else:
            writer_class(file_prefix + s).close()

def single_exon_reads(sj_to_ends):
    "the reads with no junctions after correction"
    return [read for (chrom, juncs), isoform in sj_to_ends.items() if juncs == () for read in isoform.reads]

def _nearest_end_variant(read, variants):
    return min(variants, key=lambda v: abs(read.start - v.start) + abs(read.end - v.end))

def _longest_supported_read_ends(reads):
    """The ends of the longest read whose start and end each have another read's
    within BEST_END_WINDOW, so a read sticking out past the others doesn't set
    them.  If no read's both are, the longest read, less any more than twice the
    length of the next longest."""
    starts, ends = sorted(read.start for read in reads), sorted(read.end for read in reads)

    def near(positions, pos):
        return bisect_right(positions, pos + BEST_END_WINDOW) - bisect_left(positions, pos - BEST_END_WINDOW)
    supported = [read for read in reads if near(starts, read.start) > 1 and near(ends, read.end) > 1]
    by_length = sorted(supported or reads, key=lambda read: (read.end - read.start, -read.start), reverse=True)
    if len(supported) == 0:
        while len(by_length) > 1 and by_length[0].end - by_length[0].start > 2 * (by_length[1].end - by_length[1].start):
            by_length.pop(0)
    return by_length[0].start, by_length[0].end

def _one_end_variant(isoform, reads):
    """The one isoform a junction chain's or single-exon cluster's reads give, with
    all of them: a spliced one at the densest read starts and ends, a single-exon
    one, with no junctions between its ends to place them, at the ends of one real
    read (_longest_supported_read_ends)"""
    if isoform.juncs == ():
        start, end = _longest_supported_read_ends(reads)
    else:
        start = densest_end([read.start for read in reads], BEST_END_WINDOW, min)
        end = densest_end([read.end for read in reads], BEST_END_WINDOW, max)
    return Isoform.regroup(isoform, start, end, reads)

def _end_variants(isoform, reads, *, max_ends, end_window, support):
    """The isoforms, with their reads, that a junction chain's or single-exon
    cluster's reads give: one with all of them (_one_end_variant), or with
    max_ends above 1, the best ranked max_ends of its end groups with at least
    support reads, each read going to the one whose ends are nearest its own, unless
    fewer than two groups have the support"""
    if max_ends > 1:
        groups = collapse_end_groups(end_window, Isoform.regroup(isoform, newreads=reads))
        rank_end_variants(groups)
        variants = [group for group in groups if group.num_reads >= support][:max_ends]
        if len(variants) > 1:
            variant_reads = {id(v): [] for v in variants}
            for read in reads:
                variant_reads[id(_nearest_end_variant(read, variants))].append(read)
            return [(Isoform.regroup(isoform, v.start, v.end, variant_reads[id(v)]), variant_reads[id(v)]) for v in variants]
    return [(_one_end_variant(isoform, reads), reads)]

def _end_clusters(positions, outer, min_reads):
    """The clusters of read ends: the densest of the read ends left (densest_end),
    while at least min_reads read ends are within BEST_END_WINDOW of it, which are
    then left out"""
    remaining, clusters = sorted(positions), []
    while remaining:
        pos, _ = _densest_end_clustered(remaining, BEST_END_WINDOW, outer)
        if sum(1 for p in remaining if abs(p - pos) <= BEST_END_WINDOW) < min_reads:
            break
        clusters.append(pos)
        remaining = [p for p in remaining if abs(p - pos) > BEST_END_WINDOW]
    return clusters

def _nearest_cluster(pos, clusters):
    "the cluster within BEST_END_WINDOW of a read end, if any"
    near = [c for c in clusters if abs(c - pos) <= BEST_END_WINDOW]
    return min(near, key=lambda c: abs(c - pos)) if near else None

def _overlaps_half_of_each(a, b):
    "do two (start, end) intervals overlap by more than half of each"
    overlap = min(a[1], b[1]) - max(a[0], b[0])
    return overlap > (a[1] - a[0]) / 2 and overlap > (b[1] - b[0]) / 2

def _single_exon_end_variants(isoform, reads, *, support, min_frac):
    """The isoforms, with their reads, that a single-exon cluster's reads give, at
    its clustered ends, whatever max_ends is: each pair of a start and an end
    cluster of read ends (_end_clusters) that at least support reads, and min_frac
    of the cluster's reads, have both their ends in, each read going to the one
    whose ends are nearest its own.  Reads ending anywhere, as truncated reads of a
    long transcript do, make no cluster, so a transcript isn't split into pieces;
    min_frac keeps chance clusters of a deeply covered one out.

    A transcript's full-length reads can be too few for a pair, though each of its
    ends is a cluster, where most of its reads are fragments ending in clusters of
    their own, so the longest supported read's ends (_longest_supported_read_ends)
    give a variant too, kept with support reads alone in a gene without spliced
    reads (Isoform.longest_supported, _iso_passes_support_filter).
    If a clustered pair is the same transcript, overlapping more than half of each,
    the nearest such takes, of the longest supported read's ends, those in a
    cluster, at the cluster, instead: a read sticking out past the cluster is
    outlier support, not an end.  Without a clustered pair, one isoform with all
    the reads (_one_end_variant)."""
    min_reads = max(support, math.ceil(min_frac * len(reads)))
    starts = _end_clusters([read.start for read in reads], min, min_reads)
    ends = _end_clusters([read.end for read in reads], max, min_reads)
    pairs = Counter()
    for read in reads:
        start, end = _nearest_cluster(read.start, starts), _nearest_cluster(read.end, ends)
        if start is not None and end is not None and start < end:
            pairs[(start, end)] += 1
    variants = sorted(pair for pair, n in pairs.items() if n >= min_reads)
    if len(variants) == 0:
        return [(_one_end_variant(isoform, reads), reads)]
    longest = _longest_supported_read_ends(reads)
    same = [variant for variant in variants if _overlaps_half_of_each(longest, variant)]
    if same:
        nearest = min(same, key=lambda v: (abs(longest[0] - v[0]) + abs(longest[1] - v[1]), v))
        start, end = _nearest_cluster(longest[0], starts), _nearest_cluster(longest[1], ends)
        longest = (nearest[0] if start is None else start, nearest[1] if end is None else end)
        variants.remove(nearest)
    variants = sorted(set(variants) | {longest})
    variant_reads = {pair: [] for pair in variants}
    for read in reads:
        variant_reads[min(variants, key=lambda v: (abs(read.start - v[0]) + abs(read.end - v[1]), v))].append(read)
    out = []
    for start, end in variants:
        variant = Isoform.regroup(isoform, start, end, variant_reads[(start, end)])
        variant.longest_supported = (start, end) == longest
        out.append((variant, variant_reads[(start, end)]))
    return out

class _RealignedRead:
    "a read's genomic ends from its alignment in the final realignment"
    __slots__ = ('name', 'start', 'end')

    def __init__(self, name, start, end):
        self.name, self.start, self.end = name, start, end

def _realigned_read(isoform, fields):
    """A read's genomic ends from its alignment to a spliced isoform in the final
    realignment, from count_sam_transcripts' ends line for it, which gives, in
    genomic order, how far it aligns before the isoform's first junction and past
    its last.  The padding of the isoform's terminal exons is genome sequence, so
    these are the read's ends on the genome.  None if the line has no distances."""
    left_dist, right_dist = fields[3], fields[6]
    if left_dist == 'None' or right_dist == 'None':
        return None
    return _RealignedRead(fields[0], isoform.juncs[0].start - int(left_dist), isoform.juncs[-1].end + int(right_dist))

def _gene_spliced_reads(firstpass, lines_by_iso):
    """each gene's spliced reads, frac_support's denominator, which splitting a
    junction chain between end variants doesn't change"""
    totals = {}
    for isoform in firstpass.values():
        if isoform.juncs != ():
            totals[isoform.gene_id] = totals.get(isoform.gene_id, 0) + len(lines_by_iso.get(isoform.name, []))
    return totals

def read_ends_fields(ends_path):
    """The rows of an ends file as the string fields this module passes around, with
    a missing distance spelled None, as count_sam_transcripts' own lines had it"""
    return [[row.read, row.transcript] + [str(v) for v in (row.start_sj_index, row.start_sj_dist,
                                                           row.start_tend_dist, row.end_sj_index,
                                                           row.end_sj_dist, row.end_tend_dist)]
            for row in ReadEndsReader(ends_path)]

def _rewrite_ends_and_map(ends_lines, ends_path, map_path):
    "write the ends lines, and from them the read map"
    with ReadEndsWriter(ends_path) as fh:
        for fields in ends_lines:
            fh.writeRow(fields[:2] + [None if v == 'None' else int(v) for v in fields[2:]])
    if map_path is not None:
        iso_reads = {}
        for fields in ends_lines:
            iso_reads.setdefault(fields[1], []).append(fields[0])
        with ReadMapWriter(map_path) as fh:
            for name, names in iso_reads.items():
                fh.writeReads(name, names)

def _covers_half(read, isoform):
    "does a single-exon read cover more than half of a single-exon isoform"
    return min(isoform.end, read.end) - max(isoform.start, read.start) > (isoform.end - isoform.start) / 2

def _variants_all_unsupported(variant_counts, passes_support, isoform):
    """do several end variants, with variant_counts reads each, all fail support,
    though together they pass it"""
    return len(variant_counts) > 1 and passes_support(isoform, sum(variant_counts)) and \
        not any(passes_support(isoform, n) for n in variant_counts)

def assign_final_ends(firstpass, ends_path, map_path, *, normalize_ends, max_ends, end_window, sjc_support, frac_support):
    """Each spliced isoform's reported ends come from the reads assigned to it in the
    final realignment, or directly (direct_assignment), each read's ends taken from
    its alignment, its normalized and padded transcript being only what reads align
    to: one isoform at the densest read ends, or with max_ends above 1, its end
    variants (_end_variants).  With normalize_ends, it keeps instead the ends it
    shares with its gene's isoforms with similar terminal splice sites, set when the
    first-pass isoforms were written (write_first_pass_isoforms), with all its reads.
    A junction chain whose reads pass sjc_support and frac_support, though none of
    its several end variants' do, gives one isoform with all its reads instead, as
    with max_ends 1 (_one_end_variant).

    Single-exon first-pass isoforms are returned as they are, in their places, for
    assign_single_exon_isoforms.  Rewrites the ends and read map files for the
    spliced isoforms this gives, and returns all of them by name."""
    lines_by_iso = {}
    for fields in read_ends_fields(ends_path):
        lines_by_iso.setdefault(fields[1], []).append(fields)
    gene_spliced_reads = _gene_spliced_reads(firstpass, lines_by_iso)

    def passes_support(isoform, num_reads):
        total = gene_spliced_reads.get(isoform.gene_id, 0)
        return num_reads >= sjc_support and total > 0 and num_reads / total >= frac_support

    final, out_lines = {}, []
    for isoform in firstpass.values():
        lines = lines_by_iso.get(isoform.name, [])
        reads = [] if isoform.juncs == () else [r for r in (_realigned_read(isoform, f) for f in lines) if r is not None]
        if len(reads) == 0 or normalize_ends:
            final[isoform.name] = isoform
            out_lines.extend(lines)
            continue
        line_by_read = {f[0]: f for f in lines}
        variants = _end_variants(isoform, reads, max_ends=max_ends, end_window=end_window, support=sjc_support)
        # stringent assignments always have both distances; any line without them
        # goes to the first variant
        others = {f[0] for f in lines} - {read.name for read in reads}
        variant_counts = [len(variant_reads) + (len(others) if i == 0 else 0) for i, (variant, variant_reads) in enumerate(variants)]
        if _variants_all_unsupported(variant_counts, passes_support, isoform):
            variants = [(_one_end_variant(isoform, reads), reads)]
        for i, (variant, variant_reads) in enumerate(variants):
            variant.unpadded_ends = (variant.start, variant.end)
            final[variant.name] = variant
            names = [read.name for read in variant_reads] + (sorted(others) if i == 0 else [])
            out_lines.extend([f[0], variant.name] + f[2:] for f in (line_by_read[n] for n in names))
    _rewrite_ends_and_map(out_lines, ends_path, map_path)
    return final

def assign_single_exon_isoforms(isoforms, iso_to_counts, gene_to_tot, *, se_support, frac_support):
    """The final single-exon isoforms, by the name of the first-pass isoform each
    comes from.  Single-exon isoforms aren't realigned to: each keeps the reads it
    was built from that cover more than half of it.  A cluster whose reads pass
    support, though none of its several end variants' do, gives one isoform with all
    its reads instead, as with max_ends 1 (_one_end_variant).  Support is se_support
    reads and frac_support of the gene's full-length reads: all its spliced ones, as
    calc_final_iso_support counted them, and its single-exon ones covering their
    isoforms.  Adds the isoforms' reads to iso_to_counts and gene_to_tot, as
    calc_final_iso_support counts them."""
    covering = {isoform.name: [read for read in isoform.reads if _covers_half(read, isoform)] for isoform in isoforms}
    gene_total = {}
    for isoform in isoforms:
        gene_total[isoform.gene_id] = (gene_total.get(isoform.gene_id, gene_to_tot.get(isoform.gene_id, [0, 0, 0])[1])
                                       + len(covering[isoform.name]))

    def passes_support(isoform, num_reads):
        total = gene_total[isoform.gene_id]
        return num_reads >= se_support and total > 0 and num_reads / total >= frac_support

    clusters = {}
    for isoform in isoforms:
        clusters.setdefault(isoform.end_variant_cluster or isoform.name, []).append(isoform)
    final = {}
    for variants in clusters.values():
        # if they all fail, the one isoform is at the longest supported ends anyway
        if _variants_all_unsupported([len(covering[v.name]) for v in variants], passes_support, variants[0]):
            merged = _one_end_variant(variants[0], [read for variant in variants for read in variant.reads])
            merged.reads = [read for read in merged.reads if _covers_half(read, merged)]
            final[variants[0].name] = [merged]
        else:
            for variant in variants:
                variant.reads = covering[variant.name]
                final[variant.name] = [variant]
    for isoform in (isoform for variants in final.values() for isoform in variants):
        if isoform.num_reads > 0:
            iso_to_counts[isoform.name] = [isoform.num_reads, isoform.num_reads]
            gene_reads = gene_to_tot.setdefault(isoform.gene_id, [0, 0, 0])
            gene_reads[1] += isoform.num_reads
            gene_reads[2] += isoform.num_reads
    return final

def _with_single_exon_isoforms(final, single_exon):
    "the final isoforms by name, each single-exon one in the place of the first-pass isoform it comes from"
    out = {}
    for name, isoform in final.items():
        for final_isoform in (single_exon.get(name, []) if isoform.juncs == () else [isoform]):
            out[final_isoform.name] = final_isoform
    return out

def calc_final_iso_support(read_ends_file, final_transcript_objs, no_stringent):
    iso_to_counts = {}
    gene_to_tot = {}
    for row in ReadEndsReader(read_ends_file):
        read, transcript = row.read, row.transcript
        start_sj_index, end_sj_index = row.start_sj_index, row.end_sj_index
        gene = final_transcript_objs[transcript].gene_id
        if gene not in gene_to_tot:
            # total spliced full-length, total full-length spliced + unspliced, total all
            gene_to_tot[gene] = [0, 0, 0]
        if transcript not in iso_to_counts:
            iso_to_counts[transcript] = [0, 0]
        # only spliced isoforms' reads have ends lines (assign_single_exon_isoforms);
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

def promote_backup_subsets(args, removed_subsets, final_transcript_objs, iso_to_counts, gene_to_tot, ends_path):
    """Report a spliced isoform the subset filter removed after all, when every
    isoform it is a subset of fails support: a fragment of a well-supported isoform
    is still removed, but one of a read through or extra exon too weak to report
    doesn't take a well-supported isoform with it.  It gets the reads it was built
    from that no other isoform was assigned and that cover its unique sequence
    boundaries, as an assigned read must (covers_unique_bounds), at least
    sjc_support and subset_backup_support of them, and its ends from theirs
    (_one_end_variant).  It
    keeps the gene and reference transcript build_genes gave its junction chain, as
    any isoform does, rather than taking them from an isoform it is a subset of,
    which could put an annotated transcript in a readthrough gene; its reads count
    toward the gene's, as an assigned read's do.  Adds the promoted isoforms to
    final_transcript_objs, iso_to_counts and gene_to_tot."""
    if args.subset_backup_support <= 0 or not removed_subsets:
        return
    passing = {(isoform.strand, isoform.juncs) for name, isoform in final_transcript_objs.items()
               if isoform.juncs != () and _iso_passes_support_filter(args, name, isoform.gene_id, len(isoform.exons),
                                                                     iso_to_counts, gene_to_tot)[0]}
    assigned = {row.read for row in ReadEndsReader(ends_path)} if os.path.exists(ends_path) else set()
    # longest chains first, so an isoform's supersets are decided before it: a fragment
    # of a promoted subset is still a fragment
    for name, (subset, supersets, unique_seq) in sorted(removed_subsets.items(), key=lambda item: -len(item[1][0].juncs)):
        if (subset.strand, subset.juncs) in passing or any((s.strand, s.juncs) in passing for s in supersets):
            continue
        # truncated reads of the isoforms it is a subset of, which don't reach its
        # own sequence, aren't its reads
        check_exons = subset_check_exons(subset)
        bounds = parse_unique_bounds(unique_bounds_past(unique_seq, check_exons[0], check_exons[-1], subset.strand))
        reads = [read for read in subset.reads
                 if read.name not in assigned and covers_unique_bounds(subset, bounds, read.start, read.end)]
        if len(reads) < max(args.sjc_support, args.subset_backup_support):
            continue
        # a copy of the subset, with its gene, strand and reference transcript
        promoted = _one_end_variant(subset, reads)
        promoted.unpadded_ends = (promoted.start, promoted.end)
        promoted.backup_subset = True
        logging.debug(f"subset isoform promoted, its supersets failing support: {name} ({len(reads)} reads)")
        final_transcript_objs[promoted.name] = promoted
        iso_to_counts[promoted.name] = [len(reads), len(reads)]
        gene_reads = gene_to_tot.setdefault(promoted.gene_id, [0, 0, 0])
        for i in range(3):
            gene_reads[i] += len(reads)
        passing.add((promoted.strand, promoted.juncs))

def unpadded_exons(exons, unpadded_ends, name=None):
    """A spliced isoform's exons with the realignment padding taken off its first and
    last exons, back to unpadded_ends.  An end at or past its exon's other edge is
    left where it is, as it would make the exon empty, or reversed, which Exon refuses
    with an assertion that stopped the run."""
    exons = list(exons)
    start, end = unpadded_ends
    if start < exons[0].end:
        exons[0] = Exon(start, exons[0].end)
    else:
        logging.debug(f"unpadded start {start} not before the first exon's end {exons[0].end}, kept at {exons[0].start}: {name}")
    if end > exons[-1].start:
        exons[-1] = Exon(exons[-1].start, end)
    else:
        logging.debug(f"unpadded end {end} not after the last exon's start {exons[-1].start}, kept at {exons[-1].end}: {name}")
    return exons


def write_final_isoform_output(partition, args, final_transcript_objs, iso_to_counts, gene_to_tot, annots, genome, generate_map):
    transcript_to_reads = {}
    if generate_map:
        for row in ReadMapReader(partition.output_path('countsam.read.map.tsv')):
            transcript_to_reads[row.name] = row.reads

    with open(partition.output_path('isoforms.bed'), 'w') as iso_fh, \
         open(partition.output_path('isoforms.fa'), 'w') as seq_fh, \
         IsoformCountsWriter(partition.output_path('isoform.counts.tsv')) as counts_fh, \
         ReadMapWriter(partition.output_path('isoform.read.map.tsv')) as map_fh:
        for tname in final_transcript_objs:
            # spliced isos checked against spliced total, single exon checked against full-length total
            isoform = final_transcript_objs[tname]
            passes_support, my_frac_support = _iso_passes_support_filter(args, tname, isoform.gene_id, len(isoform.exons), iso_to_counts, gene_to_tot,
                                                                         longest_supported=isoform.longest_supported)
            if passes_support:
                if len(isoform.exons) > 1:
                    # removing the padding the realignment needed from the ends
                    isoform.reset_from_exons(unpadded_exons(isoform.exons, isoform.unpadded_ends, isoform.name))

                thickStart, thickEnd, productivity, aaseq = predict_productivity(isoform, annots.start_codon_count,
                                                                                 annots.gene_to_cds_starts, annots.transcript_to_nmd_except, genome)
                convert_to_flair_bed(isoform, thickStart=thickStart, thickEnd=thickEnd, read_support=iso_to_counts[tname][0],
                                     frac_support=my_frac_support, productivity=productivity, samples=(args.sample_name,), aaseq_id=aaseq).write(iso_fh)
                seq_fh.write('>' + isoform.name + '\n')
                seq_fh.write(isoform.get_sequence(genome) + '\n')
                counts_fh.writeRow((isoform.name, iso_to_counts[tname][0], iso_to_counts[tname][1]))
                if generate_map:
                    # single-exon and promoted subset isoforms have their reads; other spliced ones' are in the read map
                    if isoform.juncs == () or isoform.backup_subset:
                        reads = [read.name for read in isoform.reads]
                    else:
                        reads = transcript_to_reads[tname]
                    map_fh.writeReads(isoform.name, reads)


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
        generate_empty_intermediate_files(partition.file_prefix, ['.firstpass.bed', '.isoform.counts.tsv', '.isoform.read.map.tsv', '.isoforms.bed', '.isoforms.fa', '.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed'])


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

        # never zero reads here: _run_region only calls this when samtools fasta found
        # reads, and both skip secondary and supplementary alignments
        _, clipping_file = generate_genomic_alignment_read_to_clipping_file(partition.file_prefix, bam_file, region.name, region.start, region.end)

        # aligning to reference transcriptome, then identifying reads that match well to reference transcripts
        # with filter_transcriptome_align
        if not args.no_align_to_annot:
            logging.info('aligning to transcriptome reference')
        read_to_annot_transcript = identify_good_match_to_annot(
            partition.file_prefix, region, annots, genome, bam_file, no_align_to_annot=args.no_align_to_annot,
            directRNA=args.directRNA, check_splice=not args.no_check_splice, trust_strand=args.trust_strand,
            fusion_breakpoints=args.fusion_breakpoints)

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
        se_reads = single_exon_reads(sj_to_ends)

        # for each junction chain, clusters ends - generates junction chain x ends
        # firstpass objects then does initial filtering by read support and
        # redundant ends also separates single exon isoforms from spliced isoforms
        # (because they're handled differently in future step for identifying
        # annotated gene/isoform names)
        candidates = process_juncs_to_firstpass_isos(args, partition.file_prefix, sj_to_ends, annots, region.name)

        # filter isoforms: remove subsets, generate unique boundary sequences
        removed_subsets = {}
        firstpass, iso_to_unique_bound = filter_firstpass_isos(args, candidates, annots, removed_subsets=removed_subsets)

        if len(firstpass.keys()) > 0:
            # this section identifies annotated gene and isoform names (primarily based on splice junction matching, secondarily by exon overlap)
            # also adjusts isoform strand, determines novel isoform and gene names
            # also normalizes transcript ends (temporarily extends ends so that transcript end alignment does not drive transcript assignment during transcriptome alignment)
            # writes out bed and fa files
            logging.info('realigning to firstpass and getting final isoforms')
            write_firstpass(partition.file_prefix, region.name, firstpass, annots, genome, unique_bound=iso_to_unique_bound)

            # reads matching a firstpass isoform exactly are assigned without the realignment,
            # and reads it can't assign aren't realigned (see direct_assignment); with
            # fusion_breakpoints, count_sam_transcripts checks that reads cover the breakpoint
            direct, unassignable = {}, set()
            if not args.fusion_breakpoints:
                spliced_reads = [read for (chrom, juncs), isoform in sj_to_ends.items() if juncs != () for read in isoform.reads]
                direct = direct_assignments(spliced_reads, firstpass.values(),
                                            read_unique_bounds(partition.output_path('firstpass.uniquebound.tsv')))
                unassignable = unassignable_reads(spliced_reads, firstpass.values(), args.ss_window)
                logging.info(f'{len(direct)} reads assigned to firstpass isoforms without realignment, '
                             f'{len(unassignable)} not realigned as no firstpass isoform can take them')

            # aligns the other spliced reads to the spliced firstpass isoforms, identifies best read -> isoform
            # alignment for each read, then gets read counts per isoform
            read_map_file = partition.output_path('countsam.read.map.tsv') if args.generate_map else None
            num_realigned = _write_reads_subset(partition.output_path('reads.fasta'), partition.output_path('realign.reads.fasta'),
                                                {read.name for read in se_reads} | direct.keys() | unassignable)
            if num_realigned > 0 and any(isoform.juncs != () for isoform in firstpass.values()):
                transcriptome_align_and_count(partition.output_path('realign.reads.fasta'),
                                              partition.output_path('firstpass.fa'),
                                              partition.output_path('firstpass.bed'),
                                              partition.output_path('isoform.counts.tsv'),
                                              read_map_file, False,  # say is not annot, requires stringent, returns different end values
                                              partition.output_path('reads.genomicclipping.tsv'),
                                              partition.output_path('firstpass.uniquebound.tsv'),
                                              directRNA=args.directRNA, no_stringent=args.no_stringent,
                                              check_splice=not args.no_check_splice, trust_strand=args.trust_strand,
                                              fusion_breakpoints=args.fusion_breakpoints,
                                              secondary_alignments=(KEEP_SUBSET_SECONDARY_ALIGNMENTS
                                                                    if args.filter in KEEP_SUBSET_FILTERS
                                                                    else MINIMAP2_SECONDARY_ALIGNMENTS))
            else:
                generate_empty_intermediate_files(partition.file_prefix, ['.isoform.counts.tsv', '.countsam.read.map.tsv', '.isoform.ends.tsv'])
            with open(partition.output_path('isoform.ends.tsv'), 'a') as ends_fh:
                write_direct_assignments(direct, ReadEndsWriter(None, outFh=ends_fh, writeHeader=False))
        else:
            logging.info('no firstpass isoforms found')
            generate_empty_intermediate_files(partition.file_prefix, ['.firstpass.fa', '.firstpass.bed', '.isoform.counts.tsv', '.countsam.read.map.tsv', '.isoform.ends.tsv'])

        # every spliced isoform's ends, from the reads assigned to it; this also writes
        # the read map, for every way spliced reads were assigned, from the ends
        final_transcript_objs = assign_final_ends(firstpass, partition.output_path('isoform.ends.tsv'),
                                                  partition.output_path('countsam.read.map.tsv') if args.generate_map else None,
                                                  normalize_ends=args.normalize_ends,
                                                  max_ends=args.max_ends, end_window=args.end_window,
                                                  sjc_support=args.sjc_support, frac_support=args.frac_support)
        iso_to_counts, gene_to_tot = calc_final_iso_support(partition.output_path('isoform.ends.tsv'), final_transcript_objs, args.no_stringent)
        # single-exon isoforms after the spliced reads are counted, which their support is a fraction of
        single_exon = assign_single_exon_isoforms([isoform for isoform in final_transcript_objs.values() if isoform.juncs == ()],
                                                  iso_to_counts, gene_to_tot,
                                                  se_support=args.single_exon_support, frac_support=args.frac_support)
        final_transcript_objs = _with_single_exon_isoforms(final_transcript_objs, single_exon)
        promote_backup_subsets(args, removed_subsets, final_transcript_objs, iso_to_counts, gene_to_tot,
                               partition.output_path('isoform.ends.tsv'))
        write_final_isoform_output(partition, args, final_transcript_objs, iso_to_counts, gene_to_tot, annots, genome, args.generate_map)


TSV_SUFFIXES = ('.isoform.counts.tsv', '.isoform.read.map.tsv')

def combine_chunks(args, output, partitions):
    files_to_combine = ['.isoforms.bed', '.isoforms.fa', '.isoform.counts.tsv']
    if args.generate_map:
        files_to_combine.append('.isoform.read.map.tsv')
    if args.keep_intermediate:
        files_to_combine.extend(['.firstpass.reallyunfiltered.bed', '.firstpass.unfiltered.bed', '.firstpass.bed'])
    combine_temp_files_by_suffix(output, [p.file_prefix for p in partitions], files_to_combine,
                                 headers=TSV_SUFFIXES)

def get_new_ids(output):
    iso_hash_to_id, gene_hash_to_id, aaseq_to_id = {}, {}, {}
    iso_count, gene_count, aaseq_count = 1, 1, 1
    with open_sorted_bed(compressed_path(output + '.isoforms.newids.bed')) as fh:
        for bed_rec in BedReader(output + '.isoforms.bed', bedClass=FlairBed):
            if bed_rec.name not in iso_hash_to_id:
                iso_id = f'FLT{iso_count:08d}'
                iso_hash_to_id[bed_rec.name] = iso_id
                iso_count += 1
            if bed_rec.gene_id not in gene_hash_to_id:
                gene_id = f'FLG{gene_count:08d}'
                gene_hash_to_id[bed_rec.gene_id] = gene_id
                gene_count += 1
            if bed_rec.aaseq_id is not None and bed_rec.aaseq_id not in aaseq_to_id:
                aaseq_id = f'FLP{aaseq_count:08d}'
                aaseq_to_id[bed_rec.aaseq_id] = aaseq_id
                aaseq_count += 1
            bed_rec.name = iso_hash_to_id[bed_rec.name]
            bed_rec.gene_id = gene_hash_to_id[bed_rec.gene_id]
            if bed_rec.aaseq_id is not None:
                bed_rec.aaseq_id = aaseq_to_id[bed_rec.aaseq_id]
            bed_rec.write(fh)
    write_aaseqs(compressed_path(output + '.aaseq.tsv'), aaseq_to_id)
    return iso_hash_to_id

def fix_ids_txt_file(iso_hash_to_id, oldfile, newfile):
    """Rename the ids in the first column, keeping the header.  This is the last
    thing written to these two files, so it is where they are compressed."""
    with fileOps.opengz(oldfile) as in_fh, fileOps.AtomicFileOpen(newfile, 'w') as fh:
        fh.write(in_fh.readline())  # header, whose first field is a column name
        for line in in_fh:
            line = line.split('\t', 1)
            if line[0] in iso_hash_to_id:  # had to add this check due to having isoforms in this file that were not in bed due to not passing final filters
                line[0] = iso_hash_to_id[line[0]]
                fh.write('\t'.join(line))

def fix_ids_fa_file(iso_hash_to_id, oldfile, newfile):
    with open_fasta_output(newfile) as fh:
        last = False
        for line in open(oldfile):
            if line[0] == '>':
                if line[1:].rstrip('\n') in iso_hash_to_id:
                    line = '>' + iso_hash_to_id[line[1:].rstrip('\n')] + '\n'
                    last = True
                else:
                    last = False
            if last:
                fh.write(line)

def rename_ids_into_place(iso_hash_to_id, joined, final):
    "rewrite the ids of a joined partition file into its final, maybe compressed, name"
    fix_ids_txt_file(iso_hash_to_id, joined, final)
    if final != joined:
        os.unlink(joined)

def move_into_place(renamed, joined):
    """Put the file written with the final ids where the joined partition output was.
    Both names follow the run's compression, so the joined file is only left behind,
    to be removed here, when compression gave the final file a different name."""
    final = compressed_path(joined)
    os.rename(compressed_path(renamed), final)
    if final != joined:
        os.unlink(joined)

def fix_iso_labels(output, generate_map):
    iso_hash_to_id = get_new_ids(output)
    rename_ids_into_place(iso_hash_to_id, output + '.isoform.counts.tsv',
                          compressed_path(output + '.isoform.counts.tsv'))
    fix_ids_fa_file(iso_hash_to_id, output + '.isoforms.fa', output + '.isoforms.newids.fa')
    move_into_place(output + '.isoforms.newids.bed', output + '.isoforms.bed')
    move_into_place(output + '.isoforms.newids.fa', output + '.isoforms.fa')
    if generate_map:
        rename_ids_into_place(iso_hash_to_id, output + '.isoform.read.map.tsv',
                              compressed_path(output + '.isoform.read.map.tsv'))

####
# main
####

def flair_transcriptome(*, genome_aligned_bam, genome, sample_name, output, annot_gtf,
                        junction_star, junction_bed, junction_support, ss_window, end_window,
                        sjc_support, single_exon_support, frac_support, directRNA, trust_strand, total_rna,
                        trust_junctions, no_stringent, no_check_splice, no_align_to_annot,
                        max_ends, filter, keep_supplementary, quality, threads, parallel_mode,
                        fusion_breakpoints, keep_intermediate, temp_dir, normalize_ends, generate_map,
                        subset_backup_support=10):
    args = TranscriptomeOpts(genome_aligned_bam=genome_aligned_bam, genome=genome,
                             sample_name=sample_name, output=output, annot_gtf=annot_gtf,
                             junction_star=junction_star, junction_bed=junction_bed,
                             junction_support=junction_support, ss_window=ss_window,
                             end_window=end_window, sjc_support=sjc_support,
                             single_exon_support=single_exon_support, frac_support=frac_support,
                             subset_backup_support=subset_backup_support,
                             directRNA=directRNA, trust_strand=trust_strand or directRNA, total_rna=total_rna,
                             trust_junctions=trust_junctions,
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
    if args.total_rna and not args.trust_strand:
        raise FlairInputDataError("--total_rna requires --trust_strand: without poly(A) tails, the alignment strand is "
                                  "the only strand evidence for single-exon reads")
    parse_filter(args.filter)
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
                                                            star_sj_tabs=args.junction_star)
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
        bed_to_gtf(compressed_path(args.output + '.isoforms.bed'),
                   compressed_path(args.output + '.isoforms.gtf'), is_flair_bed=True)

        make_big_bed(genome_fa, temp_dir + 'chrom.sizes', args.output + '.isoforms')

    finally:
        # also on failure, so runs don't leave their files in $TMPDIR
        if args.keep_intermediate:
            logging.info(f'intermediate files kept in {temp_dir}')
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)

    genome_fa.close()
