#!/usr/bin/env python3

import argparse
import logging
import re
import os
from flair.remove_internal_priming import removeinternalpriming
import pysam
from flair import FlairInputDataError
from flair.pycbio.hgdata.bed import BedReader


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


def get_annot_info(args):  # noqa: C901 - FIXME: reduce complexity
    chrtobp = {}
    if args.fusion_breakpoints:
        for bed in BedReader(args.fusion_breakpoints, numStdCols=3):
            chrtobp[bed.chrom] = bed.chromStart

    transcript_to_bp_ss_index = {}
    transcript_to_exons = {}
    transcript_to_genomic_ends = {}
    transcript_to_unique_bounds = {}
    if args.stringent or args.check_splice or args.fusion_dist or args.fusion_breakpoints or args.output_endpos:
        for bed in BedReader(args.isoforms, fixScores=True):
            name, left, right, chrom, strand = bed.name, bed.chromStart, bed.chromEnd, bed.chrom, bed.strand
            if name[:10] == 'fusiongene':
                name = '_'.join(name.split('_')[1:])
            blocksizes = [len(blk) for blk in bed.blocks]
            if strand == '+':
                transcript_to_exons[name] = blocksizes
            else:
                transcript_to_exons[name] = blocksizes[::-1]
            transcript_to_genomic_ends[name] = (left, right, strand)
            if args.fusion_breakpoints:
                blockstarts = [blk.start - left for blk in bed.blocks]
                bpindex = -1
                for i in range(len(blocksizes) - 1):
                    if left + blockstarts[i] + blocksizes[i] <= chrtobp[chrom] <= left + blockstarts[i + 1]:
                        bpindex = i
                if bpindex >= 0 and strand == '-':
                    bpindex = (len(blocksizes) - 2) - bpindex
                transcript_to_bp_ss_index[name] = bpindex
        if args.unique_bound:
            for line in open(args.unique_bound):
                name, bounds = line.rstrip().split('\t')
                bounds = [x.split('_') for x in bounds.split(',')]
                leftbounds = [int(x[1]) for x in bounds if x[0] == '0']
                rightbounds = [int(x[1]) for x in bounds if x[0] == '1']
                boundsdict = {'left': None, 'right': None}
                if len(leftbounds) > 0:
                    boundsdict['left'] = transcript_to_exons[name][0] - max(leftbounds)
                if len(rightbounds) > 0:
                    boundsdict['right'] = sum(transcript_to_exons[name][:-1]) + max(rightbounds)
                transcript_to_unique_bounds[name] = boundsdict

    return transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends, transcript_to_unique_bounds


def check_singleexon(read_start, read_end, tlen, end_norm_dist):
    if read_end - read_start > (tlen / 2) - end_norm_dist:  # must cover at least 50% of single exon transcript
        return True
    else:
        return False


def check_exonenddist(blocksize, read_edge, transcript_edge, trust_ends, disttoblock, unique_bound):
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
    left_coverage = check_exonenddist(first_blocksize, read_start, 0, trust_ends, first_blocksize - read_start, unique_bound_left)
    right_coverage = check_exonenddist(last_blocksize, read_end, tlen, trust_ends, read_end - (tlen - last_blocksize), unique_bound_right)
    return right_coverage and left_coverage


def check_stringent(coveredpos, exonpos, tlen, blockstarts, blocksizes, trust_ends, tname, end_norm_dist, transcript_to_unique_bounds):
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


def get_matchvals(args, md):
    matchvals = []
    if args.stringent or args.check_splice or args.fusion_breakpoints:
        mdblocks = re.findall(r'\d+|\D+', md)
        for b in mdblocks:
            if b[0] != '^':
                if b.isnumeric():
                    matchvals.extend([1] * int(b))
                else:
                    matchvals.append(0)
    return matchvals


def process_cigar(matchvals, cigarblocks, startpos, exoninfo, exon_bounds):  # noqa: C901 - FIXME: reduce complexity
    matchpos = 0
    coveredpos = [0] * (startpos - 1)
    queryclipping = []
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
        if btype in {4, 5}:  # soft or hard clipping:
            queryclipping.append(blen)
        elif btype == 0:  # match
            coveredpos.extend(matchvals[matchpos:matchpos + blen])
            blockstarts.append(tendpos)
            blocksizes.append(blen)
            matchpos += blen
            tendpos += blen
        elif btype in {2, 3}:  # deletion or intron
            coveredpos.extend([0] * blen)
            if blen > LARGE_INDEL_TOLDERANCE:
                if exoninfo:
                    if lb + 1 < tendpos and tendpos + blen < rb - 1:  # not in first or last exon
                        indel_detected = True
                else:
                    indel_detected = True
            tendpos += blen
            # if blen > LARGE_INDEL_TOLDERANCE: return True, None, None, None, None, None
        elif btype == 1:  # insertion
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
    return indel_detected, coveredpos, queryclipping, blockstarts, blocksizes, tendpos


def check_transcript_in_annot(exondict, tname):
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


def check_stringentandsplice(args, exoninfo, tname, coveredpos, tlen, blockstarts, blocksizes, tstart, tend, transcript_to_bp_ss_index, transcript_to_unique_bounds):
    passesstringent, passessplice, passesfusion = True, True, True
    if args.stringent or args.check_splice or args.fusion_breakpoints:
        # single exon genes always get checked
        passesstringent = check_stringent(coveredpos, exoninfo, tlen, blockstarts, blocksizes,
                                          args.trust_ends, tname, args.end_norm_dist,
                                          transcript_to_unique_bounds) if args.stringent or len(exoninfo) == 1 else True
        # only run if spliced transcript
        passessplice = check_splicesites(coveredpos, exoninfo, tstart, tend, tname) if args.check_splice and len(exoninfo) > 1 else True
        passesfusion = check_fusionbp(coveredpos, exoninfo, tstart, tend, tname, transcript_to_bp_ss_index) if args.fusion_breakpoints else True
        if tname == testtname:
            print(tname, passesstringent, passessplice)
    return passesstringent and passessplice and passesfusion


testtname = 'none'


def identify_corrected_ends(exoninfo, startpos, endpos, transcript_to_genomic_ends, tname, output_endpos, tlen):
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


def get_best_transcript(tinfo, args, transcript_to_exons, transcript_to_bp_ss_index, genomicclipping, transcript_to_genomic_ends, transcript_to_unique_bounds):
    # parse CIGAR + MD tag to ID transcript pos covered by alignment
    # get start + end of transcript on read, alignment block positions
    # also save soft/hard clipping at ends of read
    # FIXME: MIN_INSERTION_LEN isn't implemented
    # not positions of insertions larger than MIN_INSERTION_LEN, apply those to check_splice
    # filter out reads with long indels
    # generate list of 0s and 1s - transcript pos with match to query, val > 1 = insertion
    passingtranscripts = []
    for tname in tinfo:
        thist = tinfo[tname]
        # process MD tag here to query positions with mismatches
        # for MD tag, keep track of position of mismatch in all match positions
        if args.stringent or args.check_splice or args.fusion_breakpoints:
            exoninfo = check_transcript_in_annot(transcript_to_exons, tname)
        else:
            exoninfo = None
        matchvals = get_matchvals(args, thist.md)
        terminal_exon_info = exoninfo if args.allow_UTR_indels else None
        terminal_exon_bounds = transcript_to_unique_bounds[thist.name] if thist.name in transcript_to_unique_bounds else None
        indel_detected, coveredpos, queryclipping, blockstarts, blocksizes, tendpos = process_cigar(matchvals, thist.cigar, thist.startpos, terminal_exon_info, terminal_exon_bounds)
        if tname == testtname:
            print('indel', indel_detected)
        if indel_detected:
            logging.debug(f"transcript alignment dropped: indel detected: {tname}")
        elif args.trimmedreads and genomicclipping is not None and sum(queryclipping) > genomicclipping + args.soft_clipping_buffer:
            logging.debug(f"transcript alignment dropped: excess soft clipping ({sum(queryclipping)} > {genomicclipping} + {args.soft_clipping_buffer}): {tname}")
        else:
            if check_stringentandsplice(args, exoninfo, thist.name, coveredpos, thist.tlen, blockstarts, blocksizes, thist.startpos, tendpos, transcript_to_bp_ss_index, transcript_to_unique_bounds):
                left_intron_index, left_dist, right_intron_index, right_dist = identify_corrected_ends(exoninfo, thist.startpos, tendpos, transcript_to_genomic_ends, tname, args.output_endpos, thist.tlen)
                passingtranscripts.append([-1 * thist.alignscore, -1 * sum(matchvals), sum(queryclipping), thist.tlen, tname, (left_intron_index, left_dist), (right_intron_index, right_dist)])
            else:
                logging.debug(f"transcript alignment dropped: failed stringent/splice check: {tname}")

    # order passing transcripts by alignment score
    # then order by amount of query covered
    # then order by amount of transcript covered
    if len(passingtranscripts) > 0:
        passingtranscripts.sort()

        if len(passingtranscripts) == 1 or passingtranscripts[0][:3] != passingtranscripts[1][:3]:
            return [passingtranscripts[0][-3:], ]
        else:
            logging.debug(f"read dropped: ambiguous multi-mapping, top 2 transcripts tied: {passingtranscripts[0][4]}, {passingtranscripts[1][4]}")
            return None

    else:
        return None


class IsoAln(object):
    def __init__(self, name=None, p=None, cigar=None, tlen=None, als=None, md=None):
        self.name = name
        self.startpos = p
        self.cigar = cigar
        self.tlen = tlen
        self.alignscore = als
        self.md = md


def parse_sam(args, transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends, readstoclipping, transcript_to_unique_bounds):  # noqa: C901 - FIXME: reduce complexity
    lastread = None
    curr_transcripts = {}
    transcripttoreads = {}
    samfile = pysam.AlignmentFile(args.sam, 'r')
    genome = None
    if args.remove_internal_priming:
        genome = pysam.FastaFile(args.transcriptomefasta)

    for read in samfile:
        if not read.is_mapped:
            logging.debug(f"read dropped: unmapped: {read.query_name}")
        else:
            readname = read.query_name
            transcript = read.reference_name
            quality = read.mapping_quality
            if quality < args.quality:
                logging.debug(f"read dropped: low quality ({quality} < {args.quality}): {readname}")
            elif quality >= args.quality:
                # for transcriptome alignment, always take rightmost side on transcript
                if args.remove_internal_priming:
                    intprimannot = transcript_to_exons if args.permissive_last_exons else None
                    notinternalpriming = removeinternalpriming(read.reference_name,
                                                               read.reference_start,
                                                               read.reference_end, False,
                                                               genome, None, intprimannot,
                                                               args.intprimingthreshold,
                                                               args.intprimingfracAs)
                else:
                    notinternalpriming = True
                if not notinternalpriming:
                    logging.debug(f"read dropped: internal priming on {transcript}: {readname}")
                else:
                    pos = read.reference_start
                    try:
                        alignscore = read.get_tag('AS')
                        mdtag = read.get_tag('MD')
                    except KeyError as ex:
                        raise Exception(f"Missing AS or MD tag in alignment of '{read.query_name}' in '{args.sam.name}'") from ex
                    cigar = read.cigartuples
                    tlen = samfile.get_reference_length(transcript)
                    if lastread and readname != lastread:
                        if testtname in curr_transcripts:
                            print('\n', lastread, curr_transcripts.keys())
                        thisclipping = readstoclipping[lastread] if lastread in readstoclipping else None
                        assignedts = get_best_transcript(curr_transcripts, args, transcript_to_exons, transcript_to_bp_ss_index, thisclipping, transcript_to_genomic_ends, transcript_to_unique_bounds)
                        if not assignedts:
                            logging.debug(f"read dropped: no passing transcript assignment: {lastread}")
                        else:
                            for assignedt, gtstart, gtend in assignedts:
                                if assignedt not in transcripttoreads:
                                    transcripttoreads[assignedt] = []
                                transcripttoreads[assignedt].append((lastread, gtstart, gtend))

                        curr_transcripts = {}
                    curr_transcripts[transcript] = IsoAln(transcript, pos, cigar, tlen, alignscore, mdtag)
                    lastread = readname
    if lastread:
        thisclipping = readstoclipping[lastread] if lastread in readstoclipping else None
        assignedts = get_best_transcript(curr_transcripts, args, transcript_to_exons, transcript_to_bp_ss_index, thisclipping, transcript_to_genomic_ends, transcript_to_unique_bounds)
        if not assignedts:
            logging.debug(f"read dropped: no passing transcript assignment: {lastread}")
        else:
            for assignedt, gtstart, gtend in assignedts:
                if assignedt not in transcripttoreads:
                    transcripttoreads[assignedt] = []
                transcripttoreads[assignedt].append((lastread, gtstart, gtend))

    return transcripttoreads


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


MIN_INSERTION_LEN = 3
HALF_SS_WINDOW_SIZE = 6
NUM_MISTAKES_IN_SS_WINDOW = 2
TRUST_ENDS_WINDOW = 50
LARGE_INDEL_TOLDERANCE = 25
REQ_BP_ALIGNED_IN_EDGE_EXONS = 10

if __name__ == '__main__':
    args = parse_args()
    args = check_args(args)
    transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends, transcript_to_unique_bounds = get_annot_info(args)
    readstoclipping = {}
    if args.trimmedreads:
        for line in open(args.trimmedreads):
            rname, clipping = line.rstrip().split('\t')
            readstoclipping[rname] = int(clipping)
    transcripttoreads = parse_sam(args, transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends, readstoclipping, transcript_to_unique_bounds)
    write_output(args, transcripttoreads)
