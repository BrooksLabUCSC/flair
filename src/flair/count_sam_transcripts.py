#!/usr/bin/env python3

import argparse
import logging
import re
import os
from dataclasses import dataclass, field
from flair.remove_internal_priming import removeinternalpriming
import pipettor
import pysam
from flair import FlairInputDataError
from flair.pycbio.hgdata.bed import BedReader


_COUNT_SAM_TRANSCRIPTS_SCRIPT = os.path.realpath(__file__)


def parse_args():
    parser = argparse.ArgumentParser(description='''for counting transcript abundances after
            aligning reads to transcripts; for multiple mappers, only the best alignment
            for each read is used''')
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-s', '--sam', type=argparse.FileType('r'), help='sam file or - for STDIN')
    required.add_argument('-o', '--output', default='counts.txt', help='output file name')
    parser.add_argument('-i', '--isoforms',
                        help='specify isoforms.bed file if --stringent and/or --check_splice is specified')
    parser.add_argument('--stringent', action='store_true',
                        help='only count if read alignment passes stringent criteria')
    parser.add_argument('--check_splice', action='store_true',
                        help='''enforce coverage of 4 out of 6 bp around each splice site and no
            insertions greater than 3 bp at the splice site''')
    parser.add_argument('--trust_ends', action='store_true',
                        help='specify if reads are generated from a long read method with minimal fragmentation')
    parser.add_argument('-t', '--threads', default=4, type=int,
                        help='number of threads to use')
    parser.add_argument('--quality', default=0, type=int,
                        help='minimum quality threshold to consider if ends are to be trusted (0)')
    parser.add_argument('--generate_map',
                        help='''specify an output path for a txt file of which isoform each read is assigned to''')
    parser.add_argument('--output_bam',
                        help='''specify an output path for the bam file aligned to the transcriptome if desired''')
    parser.add_argument('--fusion_dist',
                        help='''minimium distance between separate read alignments on the same chromosome to be
            considered a fusion, otherwise no reads will be assumed to be fusions''')
    parser.add_argument('--remove_internal_priming', default=False, action='store_true',
                        help='specify if want to remove reads with internal priming')
    parser.add_argument('--permissive_last_exons', default=False, action='store_true',
                        help='specify if want to allow reads with internal priming in last exon of transcripts (yes for annot, no for firstpass)')
    parser.add_argument('--intprimingthreshold', type=int, default=12,
                        help='number of bases that are at leas 75%% As required to call read as internal priming')
    parser.add_argument('--intprimingfracAs', type=float, default=0.6,
                        help='number of bases that are at leas 75%% As required to call read as internal priming')
    parser.add_argument('--soft_clipping_buffer', type=int, default=50,
                        help='''number of acceptable bases for transcriptome alignment to increase softclipping by''')
    parser.add_argument('--transcriptomefasta',
                        help='provide transcriptome fasta aligned to if --remove_internal_priming is specified')
    parser.add_argument('--unique_bound',
                        help='text file with boundaries of unique sequence in isoforms that are a subset of other isoforms')
    parser.add_argument('--fusion_breakpoints',
                        help='''[OPTIONAL] fusion detection only - bed file containing locations of fusion breakpoints on the synthetic genome''')
    parser.add_argument('--allow_paralogs', default=False, action='store_true',
                        help='specify if want to allow reads to be assigned to multiple paralogs with equivalent alignment')
    parser.add_argument('--allow_UTR_indels', default=False, action='store_true',
                        help='specify if want to allow reads to include indels in UTRs (more permissive for population variation + A to I editing)')
    parser.add_argument('--trimmedreads',
                        help='[requires file path] specify if your reads are properly trimmed and you want to remove alignments '
                        'with too much softclipping at the ends (improves accuracy when possible). Provide a file of read to level of clipping when aligned to the genome.')
    parser.add_argument('--end_norm_dist', type=int, default=0,
                        help='specify the number of basepairs to extend transcript ends if you want to normalize them across transcripts in a gene and extend them')
    parser.add_argument('--output_endpos',
                        help='if desired, specify path to which to output the genomic position of all read ends after transcriptomic alignment')
    args = parser.parse_args()
    return args


def check_args(args):
    if args.stringent or args.fusion_dist or args.check_splice or args.fusion_breakpoints:
        if not os.path.exists(args.isoforms):
            raise FlairInputDataError(f'A valid isoforms bed file needs to be specified: {args.isoforms}')
    if args.fusion_dist:
        args.trust_ends = True
    return args


MIN_INSERTION_LEN = 3
HALF_SS_WINDOW_SIZE = 6
NUM_MISTAKES_IN_SS_WINDOW = 2
TRUST_ENDS_WINDOW = 50
LARGE_INDEL_TOLDERANCE = 25
REQ_BP_ALIGNED_IN_EDGE_EXONS = 10

@dataclass
class IsoformInfo:
    """Per-transcript information parsed out of an isoforms BED (which is
    upstream-derived from a GTF via gtf_to_bed)."""
    transcript_to_exons: dict = field(default_factory=dict)
    transcript_to_bp_ss_index: dict = field(default_factory=dict)
    transcript_to_genomic_ends: dict = field(default_factory=dict)
    transcript_to_unique_bounds: dict = field(default_factory=dict)


def read_isoforms_bed(*, isoforms, stringent=False, check_splice=False,  # noqa: C901 - FIXME: reduce complexity
                      fusion_dist=False, fusion_breakpoints=None,
                      output_endpos=False, unique_bound=None):
    """Read the isoforms BED (and optional unique-bound TSV) and build an
    IsoformInfo with exon block sizes, genomic ends, fusion-breakpoint
    splice-site indices, and unique-sequence boundaries per transcript."""
    info = IsoformInfo()
    chrtobp = {}
    if fusion_breakpoints:
        for bed in BedReader(fusion_breakpoints, numStdCols=3):
            chrtobp[bed.chrom] = bed.chromStart

    if stringent or check_splice or fusion_dist or fusion_breakpoints or output_endpos:
        for bed in BedReader(isoforms, fixScores=True):
            name, left, right, chrom, strand = bed.name, bed.chromStart, bed.chromEnd, bed.chrom, bed.strand
            if name[:10] == 'fusiongene':
                name = '_'.join(name.split('_')[1:])
            blocksizes = [len(blk) for blk in bed.blocks]
            if strand == '+':
                info.transcript_to_exons[name] = blocksizes
            else:
                info.transcript_to_exons[name] = blocksizes[::-1]
            info.transcript_to_genomic_ends[name] = (left, right, strand)
            if fusion_breakpoints:
                blockstarts = [blk.start - left for blk in bed.blocks]
                bpindex = -1
                for i in range(len(blocksizes) - 1):
                    if left + blockstarts[i] + blocksizes[i] <= chrtobp[chrom] <= left + blockstarts[i + 1]:
                        bpindex = i
                if bpindex >= 0 and strand == '-':
                    bpindex = (len(blocksizes) - 2) - bpindex
                info.transcript_to_bp_ss_index[name] = bpindex
        if unique_bound:
            for line in open(unique_bound):
                name, bounds = line.rstrip().split('\t')
                bounds = [x.split('_') for x in bounds.split(',')]
                leftbounds = [int(x[1]) for x in bounds if x[0] == '0']
                rightbounds = [int(x[1]) for x in bounds if x[0] == '1']
                boundsdict = {'left': None, 'right': None}
                if len(leftbounds) > 0:
                    boundsdict['left'] = info.transcript_to_exons[name][0] - max(leftbounds)
                if len(rightbounds) > 0:
                    boundsdict['right'] = sum(info.transcript_to_exons[name][:-1]) + max(rightbounds)
                info.transcript_to_unique_bounds[name] = boundsdict

    return info


def check_singleexon(read_start, read_end, tlen, end_norm_dist):
    """Decide whether a read covers enough of a single-exon transcript to be counted."""
    if read_end - read_start > (tlen / 2) - end_norm_dist:  # must cover at least 50% of single exon transcript
        return True
    else:
        return False


def check_exonenddist(blocksize, read_edge, transcript_edge, trust_ends, disttoblock, unique_bound):
    """Decide whether a read's alignment extends far enough into a terminal
    exon, honoring trust_ends and unique-bound relaxations."""
    if trust_ends:
        return abs(transcript_edge - read_edge) <= TRUST_ENDS_WINDOW
    elif unique_bound:
        if transcript_edge < read_edge:  # left end of transcript
            return unique_bound - read_edge >= min(REQ_BP_ALIGNED_IN_EDGE_EXONS, unique_bound - 5) and disttoblock > abs(transcript_edge - read_edge)
        else:
            return read_edge - unique_bound >= min(REQ_BP_ALIGNED_IN_EDGE_EXONS, (transcript_edge - unique_bound) - 5) and disttoblock > abs(transcript_edge - read_edge)
    else:
        return disttoblock >= min(REQ_BP_ALIGNED_IN_EDGE_EXONS, blocksize - 5)


def check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends, unique_bound_left, unique_bound_right):
    """Require the read to cover enough of both the first and last exons
    of a multi-exon transcript."""
    left_coverage = check_exonenddist(first_blocksize, read_start, 0, trust_ends, first_blocksize - read_start, unique_bound_left)
    right_coverage = check_exonenddist(last_blocksize, read_end, tlen, trust_ends, read_end - (tlen - last_blocksize), unique_bound_right)
    return right_coverage and left_coverage


def check_stringent(coveredpos, exonpos, tlen, blockstarts, blocksizes, trust_ends, tname, end_norm_dist, transcript_to_unique_bounds):
    """Stringent-mode coverage test: single-exon transcripts need 50%
    coverage; multi-exon transcripts need their terminal exons covered."""
    # FIXME - could add back the 80% of the transcript rule - maybe as an option? needs further testing
    read_start, read_end = blockstarts[0], blockstarts[-1] + blocksizes[-1]
    first_blocksize, last_blocksize = exonpos[0], exonpos[-1]
    # covers enough bases into the first and last exons
    if len(exonpos) == 1:  # single exon transcript
        return check_singleexon(read_start, read_end, tlen, end_norm_dist)
    else:
        if tname in transcript_to_unique_bounds:
            unique_bound_left = transcript_to_unique_bounds[tname]['left']  # can also be None
            unique_bound_right = transcript_to_unique_bounds[tname]['right']

        else:
            unique_bound_left, unique_bound_right = None, None

        return check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends, unique_bound_left, unique_bound_right)


def check_splicesites(coveredpos, exonpos, tstart, tend, tname):
    """Require that every splice site the read covers matches (within a
    tolerance window) and that at least one splice site is covered."""
    currpos = 0
    allerrors = []
    all_ss_res = ['notcov' for x in range(len(exonpos) - 1)]
    for i in range(len(exonpos) - 1):
        elen = exonpos[i]
        currpos += elen
        if tstart < currpos < tend:
            ssvals = coveredpos[currpos - 6:currpos + 4]  # total size = 10, need to check indexing, seems off. This worked in a couple cases but is not systematically tested.
            totinsert = sum([x - 1 for x in ssvals if x > 1])  # value is match = 1 + insertsize
            totmatch = sum([1 for x in ssvals if x >= 1])  # insert at pos still counts as match
            if tname == testtname:
                print(i, currpos, ssvals, totmatch, totinsert)
            if totinsert + (len(ssvals) - totmatch) > NUM_MISTAKES_IN_SS_WINDOW:
                # return False
                all_ss_res[i] = 0
            else:
                all_ss_res[i] = 1
            allerrors.append(totinsert + (len(ssvals) - totmatch))
    # Does cover at least one SJ, does not fail to match any junctions it covers
    return 0 not in all_ss_res and 1 in all_ss_res


def check_fusionbp(coveredpos, exonpos, tstart, tend, tname, transcript_to_bp_ss_index):
    """Require that the read's alignment covers the fusion breakpoint splice
    site of a fusion transcript with few mismatches in the window."""
    if tname not in transcript_to_bp_ss_index or transcript_to_bp_ss_index[tname] == -1:
        return False
    else:
        eindex = transcript_to_bp_ss_index[tname]
        currpos = sum(exonpos[:eindex + 1])
        if tstart < currpos < tend:
            ssvals = coveredpos[currpos - HALF_SS_WINDOW_SIZE:currpos + HALF_SS_WINDOW_SIZE]
            totinsert = sum([x - 1 for x in ssvals if x > 1])  # value is match = 1 + insertsize
            totmatch = sum([1 for x in ssvals if x >= 1])  # insert at pos still counts as match
            if tname == testtname:
                print(currpos, ssvals, totmatch, totinsert)
            if totinsert + (len(ssvals) - totmatch) <= NUM_MISTAKES_IN_SS_WINDOW:
                return True
        return False


def get_matchvals(md, *, stringent, check_splice, fusion_breakpoints):
    """Expand the MD tag into a per-reference-position vector of 1 (match)
    and 0 (mismatch), used later to score splice-site windows."""
    matchvals = []
    if stringent or check_splice or fusion_breakpoints:
        mdblocks = re.findall(r'\d+|\D+', md)
        for b in mdblocks:
            if b[0] != '^':
                if b.isnumeric():
                    matchvals.extend([1] * int(b))
                else:
                    matchvals.append(0)
    return matchvals


def process_cigar(matchvals, cigarblocks, startpos, exoninfo, exon_bounds):  # noqa: C901 - FIXME: reduce complexity
    """Walk the CIGAR operations of a transcript alignment to build a
    per-transcript-position coverage vector, collect query-side clipping,
    record alignment block positions, and flag large indels outside the
    terminal exons."""
    matchpos = 0
    coveredpos = [0] * (startpos - 1)
    query_clipping = []
    tendpos = startpos
    blockstarts, blocksizes = [], []
    if exoninfo:
        lb, rb = exoninfo[0], sum(exoninfo) - exoninfo[-1]
        if exon_bounds:
            if exon_bounds['left']:
                lb = 0
            if exon_bounds['right']:
                rb = sum(exoninfo)
    else:
        lb, rb = None, None
    indel_detected = False
    for btype, blen in cigarblocks:
        if btype in (pysam.CSOFT_CLIP, pysam.CHARD_CLIP):
            query_clipping.append(blen)
        elif btype == pysam.CMATCH:
            coveredpos.extend(matchvals[matchpos:matchpos + blen])
            blockstarts.append(tendpos)
            blocksizes.append(blen)
            matchpos += blen
            tendpos += blen
        elif btype in (pysam.CDEL, pysam.CREF_SKIP):
            coveredpos.extend([0] * blen)
            if blen > LARGE_INDEL_TOLDERANCE:
                if exoninfo:
                    if lb + 1 < tendpos and tendpos + blen < rb - 1:  # not in first or last exon
                        indel_detected = True
                else:
                    indel_detected = True
            tendpos += blen
            # if blen > LARGE_INDEL_TOLDERANCE: return True, None, None, None, None, None
        elif btype == pysam.CINS:
            if len(coveredpos) == 0:
                coveredpos.append(blen)
            else:
                coveredpos[-1] += blen
            if blen > LARGE_INDEL_TOLDERANCE:
                if exoninfo:
                    if lb + 1 < tendpos < rb - 1:  # not in first or last exon
                        indel_detected = True
                else:
                    indel_detected = True
            # if blen > LARGE_INDEL_TOLDERANCE: return True, None, None, None, None, None
    return indel_detected, coveredpos, query_clipping, blockstarts, blocksizes, tendpos


def check_transcript_in_annot(exondict, tname):
    """Look up exon structure for a transcript name, raising a descriptive
    error if the annotation and FASTA name sets disagree."""
    try:
        exoninfo = exondict[tname]
    except KeyError:
        raise Exception(
            f"The transcript name ({tname}) in the annotation fasta do not appear to match the ones in the isoforms file."
            " You may be able to fix this by using gtf_to_bed and bedtools getfasta on your annotation gtf"
            " and using the resulting file as your annotation fasta input to this program")
    except Exception as ex:
        raise Exception("** check_splice FAILED for %s" % (tname)) from ex
    return exoninfo


def check_stringent_and_splice(exoninfo, tname, coveredpos, tlen, blockstarts, blocksizes, tstart, tend,
                               transcript_to_bp_ss_index, transcript_to_unique_bounds,
                               *, stringent, check_splice, fusion_breakpoints, trust_ends, end_norm_dist):
    """Combined filter: an alignment must pass the stringent coverage, splice-site,
    and fusion-breakpoint checks that are enabled."""
    passes_stringent, passes_splice, passes_fusion = True, True, True
    if stringent or check_splice or fusion_breakpoints:
        # single exon genes always get checked
        passes_stringent = check_stringent(coveredpos, exoninfo, tlen, blockstarts, blocksizes,
                                           trust_ends, tname, end_norm_dist,
                                           transcript_to_unique_bounds) if stringent or len(exoninfo) == 1 else True
        # only run if spliced transcript
        passes_splice = check_splicesites(coveredpos, exoninfo, tstart, tend, tname) if check_splice and len(exoninfo) > 1 else True
        passes_fusion = check_fusionbp(coveredpos, exoninfo, tstart, tend, tname, transcript_to_bp_ss_index) if fusion_breakpoints else True
        if tname == testtname:
            print(tname, passes_stringent, passes_splice)
    return passes_stringent and passes_splice and passes_fusion


testtname = 'none'


def identify_corrected_ends(exoninfo, startpos, endpos, transcript_to_genomic_ends, tname, output_endpos, tlen):
    """Locate the intron-relative position of a read's start/end on its best
    transcript so callers can translate transcript coords back to genome."""
    left_intron_index, left_dist, right_intron_index, right_dist = None, None, None, None
    # Can only correct read ends if assigned to spliced transcript
    if output_endpos:
        if len(exoninfo) > 1:
            gtstrand = transcript_to_genomic_ends[tname][2]
            currpos = 0
            for i in range(len(exoninfo) - 1):
                elen = exoninfo[i]
                currpos += elen
                if left_intron_index is None and startpos < currpos:
                    left_intron_index = i
                    left_dist = currpos - startpos
                if currpos < endpos:
                    right_intron_index = i
                    right_dist = endpos - currpos
            if gtstrand == '-':
                left_intron_index, right_intron_index = (len(exoninfo) - 2) - right_intron_index, (len(exoninfo) - 2) - left_intron_index
                left_dist, right_dist = right_dist, left_dist
        else:
            left_dist = startpos
            right_dist = tlen - endpos
    return left_intron_index, left_dist, right_intron_index, right_dist


def get_best_transcript(tinfo, info, genomicclipping,
                        *, stringent, check_splice, fusion_breakpoints, allow_UTR_indels,
                        trimmedreads, soft_clipping_buffer, output_endpos,
                        trust_ends, end_norm_dist):
    """Given all transcript alignments for a single read, apply the filtering
    checks and return the single best assignment (or None if none qualify or
    the top two tie)."""
    # parse CIGAR + MD tag to ID transcript pos covered by alignment
    # get start + end of transcript on read, alignment block positions
    # also save soft/hard clipping at ends of read
    # FIXME: MIN_INSERTION_LEN isn't implemented
    # not positions of insertions larger than MIN_INSERTION_LEN, apply those to check_splice
    # filter out reads with long indels
    # generate list of 0s and 1s - transcript pos with match to query, val > 1 = insertion
    passing_transcripts = []
    for tname in tinfo:
        thist = tinfo[tname]
        # process MD tag here to query positions with mismatches
        # for MD tag, keep track of position of mismatch in all match positions
        if stringent or check_splice or fusion_breakpoints:
            exoninfo = check_transcript_in_annot(info.transcript_to_exons, tname)
        else:
            exoninfo = None
        matchvals = get_matchvals(thist.md, stringent=stringent, check_splice=check_splice,
                                  fusion_breakpoints=fusion_breakpoints)
        terminal_exon_info = exoninfo if allow_UTR_indels else None
        terminal_exon_bounds = info.transcript_to_unique_bounds[thist.name] if thist.name in info.transcript_to_unique_bounds else None
        indel_detected, coveredpos, query_clipping, blockstarts, blocksizes, tendpos = process_cigar(matchvals, thist.cigar, thist.startpos, terminal_exon_info, terminal_exon_bounds)
        if indel_detected:
            logging.debug(f"transcript alignment dropped: indel detected: {tname}")
        elif trimmedreads and genomicclipping is not None and sum(query_clipping) > genomicclipping + soft_clipping_buffer:
            logging.debug(f"transcript alignment dropped: excess soft clipping ({sum(query_clipping)} > {genomicclipping} + {soft_clipping_buffer}): {tname}")
        else:
            if check_stringent_and_splice(exoninfo, thist.name, coveredpos, thist.tlen, blockstarts, blocksizes,
                                          thist.startpos, tendpos, info.transcript_to_bp_ss_index, info.transcript_to_unique_bounds,
                                          stringent=stringent, check_splice=check_splice,
                                          fusion_breakpoints=fusion_breakpoints,
                                          trust_ends=trust_ends, end_norm_dist=end_norm_dist):
                left_intron_index, left_dist, right_intron_index, right_dist = identify_corrected_ends(exoninfo, thist.startpos, tendpos, info.transcript_to_genomic_ends, tname, output_endpos, thist.tlen)
                passing_transcripts.append([-1 * thist.alignscore, -1 * sum(matchvals), sum(query_clipping), thist.tlen, tname, (left_intron_index, left_dist), (right_intron_index, right_dist)])
            else:
                logging.debug(f"transcript alignment dropped: failed stringent/splice check: {tname}")

    # order passing transcripts by alignment score
    # then order by amount of query covered
    # then order by amount of transcript covered
    if len(passing_transcripts) > 0:
        passing_transcripts.sort()

        if len(passing_transcripts) == 1 or passing_transcripts[0][:3] != passing_transcripts[1][:3]:
            return [passing_transcripts[0][-3:], ]
        else:
            logging.debug(f"read dropped: ambiguous multi-mapping, top 2 transcripts tied: {passing_transcripts[0][4]}, {passing_transcripts[1][4]}")
            return None

    else:
        return None


class IsoAln(object):
    """Alignment fields needed to score one transcript alignment of one read."""
    def __init__(self, name=None, p=None, cigar=None, tlen=None, als=None, md=None):
        self.name = name
        self.startpos = p
        self.cigar = cigar
        self.tlen = tlen
        self.alignscore = als
        self.md = md


def parse_sam(sam, info, readstoclipping,  # noqa: C901 - FIXME: reduce complexity
              *, quality, remove_internal_priming, transcriptomefasta,
              permissive_last_exons, intprimingthreshold, intprimingfracAs,
              stringent, check_splice, fusion_breakpoints, allow_UTR_indels,
              trimmedreads, soft_clipping_buffer, output_endpos,
              trust_ends, end_norm_dist):
    """Iterate the SAM stream, group alignments per read, call get_best_transcript,
    and accumulate {transcript: [(read, gt_start, gt_end), ...]}."""
    lastread = None
    curr_transcripts = {}
    transcript_to_reads = {}
    samfile = pysam.AlignmentFile(sam, 'r')
    genome = None
    if remove_internal_priming:
        genome = pysam.FastaFile(transcriptomefasta)

    for read in samfile:
        if not read.is_mapped:
            logging.debug(f"read dropped: unmapped: {read.query_name}")
        else:
            readname = read.query_name
            transcript = read.reference_name
            rquality = read.mapping_quality
            if rquality < quality:
                logging.debug(f"read dropped: low quality ({rquality} < {quality}): {readname}")
            elif rquality >= quality:
                # for transcriptome alignment, always take rightmost side on transcript
                if remove_internal_priming:
                    intprim_annot = info.transcript_to_exons if permissive_last_exons else None
                    not_internal_priming = removeinternalpriming(read.reference_name,
                                                                 read.reference_start,
                                                                 read.reference_end, False,
                                                                 genome, None, intprim_annot,
                                                                 intprimingthreshold,
                                                                 intprimingfracAs)
                else:
                    not_internal_priming = True
                if not not_internal_priming:
                    logging.debug(f"read dropped: internal priming on {transcript}: {readname}")
                else:
                    pos = read.reference_start
                    try:
                        alignscore = read.get_tag('AS')
                        mdtag = read.get_tag('MD')
                    except KeyError as ex:
                        raise Exception(f"Missing AS or MD tag in alignment of '{read.query_name}'") from ex
                    cigar = read.cigartuples
                    tlen = samfile.get_reference_length(transcript)
                    if lastread and readname != lastread:
                        if testtname in curr_transcripts:
                            print('\n', lastread, curr_transcripts.keys())
                        clipping = readstoclipping[lastread] if lastread in readstoclipping else None
                        assignedts = get_best_transcript(curr_transcripts, info, clipping,
                                                         stringent=stringent, check_splice=check_splice,
                                                         fusion_breakpoints=fusion_breakpoints,
                                                         allow_UTR_indels=allow_UTR_indels,
                                                         trimmedreads=trimmedreads,
                                                         soft_clipping_buffer=soft_clipping_buffer,
                                                         output_endpos=output_endpos,
                                                         trust_ends=trust_ends, end_norm_dist=end_norm_dist)
                        if not assignedts:
                            logging.debug(f"read dropped: no passing transcript assignment: {lastread}")
                        else:
                            for assignedt, gtstart, gtend in assignedts:
                                if assignedt not in transcript_to_reads:
                                    transcript_to_reads[assignedt] = []
                                transcript_to_reads[assignedt].append((lastread, gtstart, gtend))

                        curr_transcripts = {}
                    curr_transcripts[transcript] = IsoAln(transcript, pos, cigar, tlen, alignscore, mdtag)
                    lastread = readname
    if lastread:
        clipping = readstoclipping[lastread] if lastread in readstoclipping else None
        assignedts = get_best_transcript(curr_transcripts, info, clipping,
                                         stringent=stringent, check_splice=check_splice,
                                         fusion_breakpoints=fusion_breakpoints,
                                         allow_UTR_indels=allow_UTR_indels,
                                         trimmedreads=trimmedreads,
                                         soft_clipping_buffer=soft_clipping_buffer,
                                         output_endpos=output_endpos,
                                         trust_ends=trust_ends, end_norm_dist=end_norm_dist)
        if not assignedts:
            logging.debug(f"read dropped: no passing transcript assignment: {lastread}")
        else:
            for assignedt, gtstart, gtend in assignedts:
                if assignedt not in transcript_to_reads:
                    transcript_to_reads[assignedt] = []
                transcript_to_reads[assignedt].append((lastread, gtstart, gtend))

    return transcript_to_reads


def write_output(args, transcripttoreads):
    if args.output_endpos:
        endout = open(args.output_endpos, 'w')

    if args.generate_map:
        mapout = open(args.generate_map, 'w')
    countout = open(args.output, 'wt')
    for t in transcripttoreads:
        # print(t, transcripttoreads[t])
        if args.generate_map:
            mapout.write(t + '\t' + ','.join([x[0] for x in transcripttoreads[t]]) + '\n')
        countout.write(t + '\t' + str(len(transcripttoreads[t])) + '\n')
        if args.output_endpos:
            for r, s, e in transcripttoreads[t]:
                endout.write('\t'.join([str(x) for x in [r, t, s[0], s[1], e[0], e[1]]]) + '\n')
    if args.output_endpos:
        endout.close()


def build_count_sam_transcripts_cmd(*, output, sam='-', threads=4, quality=0,   # noqa: C901 - linear function okay
                                    isoforms=None, stringent=False, check_splice=False,
                                    trust_ends=False, generate_map=None, output_bam=None,
                                    fusion_dist=None, remove_internal_priming=False,
                                    permissive_last_exons=False, intprimingthreshold=12,
                                    intprimingfracAs=0.6, soft_clipping_buffer=50,
                                    transcriptomefasta=None, unique_bound=None,
                                    fusion_breakpoints=None, allow_paralogs=False,
                                    allow_UTR_indels=False, trimmedreads=None,
                                    end_norm_dist=0, output_endpos=None):
    """Build count_sam_transcripts.py argv."""
    # FIXNE: default values should be centralized
    cmd = ['python3', _COUNT_SAM_TRANSCRIPTS_SCRIPT,
           '--sam', str(sam), '-o', str(output),
           '--quality', str(quality)]
    if threads != 4:
        cmd += ['-t', str(threads)]
    if isoforms:
        cmd += ['-i', str(isoforms)]
    if stringent:
        cmd.append('--stringent')
    if check_splice:
        cmd.append('--check_splice')
    if trust_ends:
        cmd.append('--trust_ends')
    if generate_map:
        cmd += ['--generate_map', str(generate_map)]
    if output_bam:
        cmd += ['--output_bam', str(output_bam)]
    if fusion_dist:
        cmd += ['--fusion_dist', str(fusion_dist)]
    if remove_internal_priming:
        cmd += ['--remove_internal_priming',
                '--intprimingthreshold', str(intprimingthreshold),
                '--intprimingfracAs', str(intprimingfracAs)]
        if transcriptomefasta:
            cmd += ['--transcriptomefasta', str(transcriptomefasta)]
    if permissive_last_exons:
        cmd.append('--permissive_last_exons')
    if soft_clipping_buffer != 50:
        cmd += ['--soft_clipping_buffer', str(soft_clipping_buffer)]
    if unique_bound:
        cmd += ['--unique_bound', str(unique_bound)]
    if fusion_breakpoints:
        cmd += ['--fusion_breakpoints', str(fusion_breakpoints)]
    if allow_paralogs:
        cmd.append('--allow_paralogs')
    if allow_UTR_indels:
        cmd.append('--allow_UTR_indels')
    if trimmedreads:
        cmd += ['--trimmedreads', str(trimmedreads)]
    if end_norm_dist:
        cmd += ['--end_norm_dist', str(end_norm_dist)]
    if output_endpos:
        cmd += ['--output_endpos', str(output_endpos)]
    return cmd


def run_count_sam_transcripts(*, output, mm2_cmd=None, sam='-', threads=4, quality=0,
                              isoforms=None, stringent=False, check_splice=False,
                              trust_ends=False, generate_map=None, output_bam=None,
                              fusion_dist=None, remove_internal_priming=False,
                              permissive_last_exons=False, intprimingthreshold=12,
                              intprimingfracAs=0.6, soft_clipping_buffer=50,
                              transcriptomefasta=None, unique_bound=None,
                              fusion_breakpoints=None, allow_paralogs=False,
                              allow_UTR_indels=False, trimmedreads=None,
                              end_norm_dist=0, output_endpos=None):
    """Run count_sam_transcripts.py; if mm2_cmd given, pipe its stdout in as SAM."""
    cmd = build_count_sam_transcripts_cmd(
        output=output, sam=sam, threads=threads, quality=quality,
        isoforms=isoforms, stringent=stringent, check_splice=check_splice,
        trust_ends=trust_ends, generate_map=generate_map, output_bam=output_bam,
        fusion_dist=fusion_dist, remove_internal_priming=remove_internal_priming,
        permissive_last_exons=permissive_last_exons,
        intprimingthreshold=intprimingthreshold, intprimingfracAs=intprimingfracAs,
        soft_clipping_buffer=soft_clipping_buffer,
        transcriptomefasta=transcriptomefasta, unique_bound=unique_bound,
        fusion_breakpoints=fusion_breakpoints, allow_paralogs=allow_paralogs,
        allow_UTR_indels=allow_UTR_indels, trimmedreads=trimmedreads,
        end_norm_dist=end_norm_dist, output_endpos=output_endpos)
    pipeline = [mm2_cmd, cmd] if mm2_cmd else [cmd]
    pipettor.run(pipeline)


if __name__ == '__main__':
    args = parse_args()
    args = check_args(args)
    info = read_isoforms_bed(
        isoforms=args.isoforms, stringent=args.stringent, check_splice=args.check_splice,
        fusion_dist=args.fusion_dist, fusion_breakpoints=args.fusion_breakpoints,
        output_endpos=args.output_endpos, unique_bound=args.unique_bound)
    readstoclipping = {}
    if args.trimmedreads:
        for line in open(args.trimmedreads):
            rname, clipping = line.rstrip().split('\t')
            readstoclipping[rname] = int(clipping)
    transcript_to_reads = parse_sam(args.sam, info, readstoclipping,
                                    quality=args.quality, remove_internal_priming=args.remove_internal_priming,
                                    transcriptomefasta=args.transcriptomefasta,
                                    permissive_last_exons=args.permissive_last_exons,
                                    intprimingthreshold=args.intprimingthreshold,
                                    intprimingfracAs=args.intprimingfracAs,
                                    stringent=args.stringent, check_splice=args.check_splice,
                                    fusion_breakpoints=args.fusion_breakpoints,
                                    allow_UTR_indels=args.allow_UTR_indels,
                                    trimmedreads=args.trimmedreads,
                                    soft_clipping_buffer=args.soft_clipping_buffer,
                                    output_endpos=args.output_endpos,
                                    trust_ends=args.trust_ends, end_norm_dist=args.end_norm_dist)
    write_output(args, transcript_to_reads)
