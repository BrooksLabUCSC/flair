#!/usr/bin/env python3

#import sys
import argparse
import re
#import csv
#import math
import os
#from collections import Counter
#from collections import namedtuple
from flair.remove_internal_priming import removeinternalpriming
import pysam
from flair import FlairInputDataError

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
    parser.add_argument('--transcriptomefasta',
                                            help='provide transcriptome fasta aligned to if --remove_internal_priming is specified')
    parser.add_argument('--fusion_breakpoints',
                                            help='''[OPTIONAL] fusion detection only - bed file containing locations of fusion breakpoints on the synthetic genome''')
    parser.add_argument('--allow_paralogs', default=False, action='store_true',
                                            help='specify if want to allow reads to be assigned to multiple paralogs with equivalent alignment')
    parser.add_argument('--trimmedreads',
                        help='specify if your reads are properly trimmed and you want to remove alignments with too much softclipping at the ends (improves accuracy when possible). Provide a file of read to level of clipping when aligned to the genome.')
    parser.add_argument('--endnormdist', type=int, default=0,
                        help='specify the number of basepairs to extend transcript ends if you want to normalize them across transcripts in a gene and extend them')
    parser.add_argument('--output_endpos',
                        help='whether to output the genomic position of all read ends after transcriptomic alignment')
    args = parser.parse_args()
    return args


def check_args(args):
    if args.stringent or args.fusion_dist or args.check_splice or args.fusion_breakpoints:
        if not os.path.exists(args.isoforms):
            raise FlairInputDataError(f'A valid isoforms bed file needs to be specified: {args.isoforms}')
    if args.fusion_dist:
        args.trust_ends = True
    return args


def get_annot_info(args):
    chrtobp = {}
    if args.fusion_breakpoints:
        for line in open(args.fusion_breakpoints):
            line = line.split('\t')
            chr, pos = line[0], int(line[1])
            chrtobp[chr] = pos

    transcript_to_bp_ss_index = {}
    transcript_to_exons = {}
    transcript_to_genomic_ends = {}
    if args.stringent or args.check_splice or args.fusion_dist or args.fusion_breakpoints:
        for line in open(args.isoforms):
            line = line.rstrip().split('\t')
            name, left, right, chrom, strand = line[3], int(line[1]), int(line[2]), line[0], line[5]
            if name[:10] == 'fusiongene': name = '_'.join(name.split('_')[1:])
            blocksizes = [int(n) for n in line[10].rstrip(',').split(',')]
            if strand == '+': transcript_to_exons[name] = blocksizes
            else: transcript_to_exons[name] = blocksizes[::-1]
            transcript_to_genomic_ends[name] = (left, right, strand)
            if args.fusion_breakpoints:
                blockstarts = [int(n) for n in line[11].rstrip(',').split(',')]
                bpindex = -1
                for i in range(len(blocksizes) - 1):
                    if left + blockstarts[i] + blocksizes[i] <= chrtobp[chrom] <= left + blockstarts[i + 1]:
                        bpindex = i
                if bpindex >= 0 and strand == '-': bpindex = (len(blocksizes)-2) - bpindex
                transcript_to_bp_ss_index[name] = bpindex


    return transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends

def check_singleexon(read_start, read_end, tlen, endnormdist):
    if read_end-read_start > (tlen/2)-endnormdist: ##must cover at least 50% of single exon transcript
        return True
    else:
        return False

def check_exonenddist(blocksize, disttoend, trust_ends, disttoblock):
    if trust_ends:
        return disttoend <= TRUST_ENDS_WINDOW
    else:
        return disttoblock >= min(REQ_BP_ALIGNED_IN_EDGE_EXONS, blocksize-5)

def check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends):
    left_coverage = check_exonenddist(first_blocksize, read_start, trust_ends, first_blocksize-read_start)
    right_coverage = check_exonenddist(last_blocksize, tlen-read_end, trust_ends, read_end - (tlen - last_blocksize))
    return right_coverage and left_coverage

def check_stringent(coveredpos, exonpos, tlen, blockstarts, blocksizes, trust_ends, tname, endnormdist):
    matchpos = len([x for x in coveredpos if x == 1])
    # FIXME - could add back the 80% of the transcript rule - maybe as an option? needs further testing
    read_start, read_end = blockstarts[0], blockstarts[-1] + blocksizes[-1]
    first_blocksize, last_blocksize = exonpos[0], exonpos[-1]
    # covers enough bases into the first and last exons
    if len(exonpos) == 1:  # single exon transcript
        return check_singleexon(read_start, read_end, tlen, endnormdist)
    else:
        return check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends)

def check_splicesites(coveredpos, exonpos, tstart, tend, tname):
    currpos = 0
    allerrors = []
    for i in range(len(exonpos)-1):
        elen = exonpos[i]
        currpos += elen
        if tstart < currpos < tend:
            ssvals = coveredpos[currpos - 6:currpos + 4] ##total size = 10, need to check indexing, seems off. This worked in a couple cases but is not systematically tested.
            totinsert = sum([x-1 for x in ssvals if x > 1]) # value is match = 1 + insertsize
            totmatch = sum([1 for x in ssvals if x >= 1]) # insert at pos still counts as match
            if tname == testtname:
                print(i, currpos, ssvals, totmatch, totinsert)
            if totinsert + (len(ssvals)-totmatch) > NUM_MISTAKES_IN_SS_WINDOW:
                return False
            allerrors.append(totinsert + (len(ssvals)-totmatch))
    return True

def check_fusionbp(coveredpos, exonpos, tstart, tend, tname, transcript_to_bp_ss_index):
    if tname not in transcript_to_bp_ss_index or transcript_to_bp_ss_index[tname] == -1:
        return False
    else:
        eindex = transcript_to_bp_ss_index[tname]
        currpos = sum(exonpos[:eindex+1])
        if tstart < currpos < tend:
            ssvals = coveredpos[currpos - HALF_SS_WINDOW_SIZE:currpos + HALF_SS_WINDOW_SIZE]
            totinsert = sum([x-1 for x in ssvals if x > 1]) # value is match = 1 + insertsize
            totmatch = sum([1 for x in ssvals if x >= 1]) # insert at pos still counts as match
            if tname == testtname:
                print(i, currpos, ssvals, totmatch, totinsert)
            if totinsert + (len(ssvals)-totmatch) > NUM_MISTAKES_IN_SS_WINDOW:
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

def process_cigar(args, matchvals, cigarblocks, startpos):
    matchpos = 0
    coveredpos = [0] * (startpos - 1)
    queryclipping = []
    tendpos = startpos
    blockstarts, blocksizes = [], []
    for btype, blen in cigarblocks:
        if btype in {4,5}:# soft or hard clipping:
            queryclipping.append(blen)
        elif btype == 0 and (args.stringent or args.check_splice or args.fusion_breakpoints): # match
            coveredpos.extend(matchvals[matchpos:matchpos + blen])
            blockstarts.append(tendpos)
            blocksizes.append(blen)
            matchpos += blen
            tendpos += blen
        elif btype in {2,3}:# deletion or intron
            if args.stringent or args.check_splice or args.fusion_breakpoints:
                coveredpos.extend([0] * blen)
                tendpos += blen
            if blen > LARGE_INDEL_TOLDERANCE: return True, None, None, None, None, None
        elif btype == 1: # insertion
            if args.stringent or args.check_splice or args.fusion_breakpoints:
                coveredpos[-1] += blen
            if blen > LARGE_INDEL_TOLDERANCE: return True, None, None, None, None, None
    return False, coveredpos, queryclipping, blockstarts, blocksizes, tendpos

def check_transcript_in_annot(exondict, tname):
    try:
        exoninfo = exondict[tname]
    except KeyError:
        raise Exception(
                "The transcript names in the annotation fasta do not appear to match the ones in the isoforms file. You may be able to fix this by using gtf_to_bed and bed_to_sequence on your annotation gtf and using the resulting file as your annotation fasta input to this program")
    except Exception as ex:
        raise Exception("** check_splice FAILED for %s" % (tname)) from ex
    return exoninfo

def check_stringentandsplice(args, transcript_to_exons, tname, coveredpos, tlen, blockstarts, blocksizes, tstart, tend, transcript_to_bp_ss_index):
    passesstringent, passessplice, passesfusion = True, True, True
    if args.stringent or args.check_splice or args.fusion_breakpoints:
        exoninfo = check_transcript_in_annot(transcript_to_exons, tname)
        passesstringent = check_stringent(coveredpos, exoninfo, tlen, blockstarts, blocksizes,
                                                                          args.trust_ends, tname, args.endnormdist) if args.stringent else True
        passessplice = check_splicesites(coveredpos, exoninfo, tstart, tend, tname) if args.check_splice else True
        passesfusion = check_fusionbp(coveredpos, exoninfo, tstart, tend, tname, transcript_to_bp_ss_index) if args.fusion_breakpoints else True
        if tname == testtname:
            print(tname, passesstringent, passessplice)
    return passesstringent, passessplice, passesfusion

testtname = 'none'


def get_best_transcript(tinfo, args, transcript_to_exons, transcript_to_bp_ss_index, genomicclipping, transcript_to_genomic_ends):
    # parse CIGAR + MD tag to ID transcript pos covered by alignment
    # get start + end of transcript on read, alignment block positions
    # also save soft/hard clipping at ends of read
    # not positions of insertions larger than MIN_INSERTION_LEN, apply those to check_splice
    # filter out reads with long indels
    # generate list of 0s and 1s - transcript pos with match to query, val > 1 = insertion
    passingtranscripts = []
    for tname in tinfo:
        thist = tinfo[tname]
        # process MD tag here to query positions with mismatches
        # for MD tag, keep track of position of mismatch in all match positions
        matchvals = get_matchvals(args, thist.md)
        indel_detected, coveredpos, queryclipping, blockstarts, blocksizes, tendpos = process_cigar(args, matchvals, thist.cigar, thist.startpos)
        if tname == testtname:
            print('indel', indel_detected)
        # print(tname, indel_detected)
        # if queryclipping is not None: print(sum(queryclipping) <= genomicclipping + 50, sum(queryclipping), genomicclipping)
        if not indel_detected and (not args.trimmedreads or genomicclipping == None or sum(queryclipping) <= genomicclipping+SOFT_CLIPPING_BUFFER):
        #     if check_stringentandsplice(args, transcript_to_exons, thist.name, coveredpos, thist.tlen, blockstarts, blocksizes, thist.startpos, tendpos, transcript_to_bp_ss_index):

        # if not indel_detected:
        #     hasworseclipping = args.trimmedreads and genomicclipping != None and sum(queryclipping) > genomicclipping + 25
            passesstringent, passessplice, passesfusion = check_stringentandsplice(args, transcript_to_exons, thist.name, coveredpos, thist.tlen, blockstarts, blocksizes, thist.startpos, tendpos, transcript_to_bp_ss_index)
            # print(passesstringent, passessplice, passesfusion)
            # if passessplice and passesfusion:
            if passessplice and passesfusion and passesstringent:
                if args.stringent:
                    ##THIS ONLY WORKS IF STRINGENT IS ALSO ACTIVATED
                    gtstart, gtend, gtstrand = transcript_to_genomic_ends[tname]
                    if gtstrand == '+':
                        outstart = gtstart + thist.startpos
                        outend = gtend - (sum(transcript_to_exons[tname])-tendpos)
                    else:
                        outend = gtend-thist.startpos
                        outstart = gtstart + (sum(transcript_to_exons[tname])-tendpos)
                else:
                    outstart, outend = 0, 0
                passingtranscripts.append([-1 * thist.alignscore, -1 * sum(matchvals), sum(queryclipping), thist.tlen, tname, outstart, outend])
                # passingtranscripts.append([not passesstringent, hasworseclipping, -1 * thist.alignscore, -1 * sum(matchvals), sum(queryclipping), thist.tlen, tname, outstart, outend])
                # passingtranscripts.append([not passesstringent, -1 * thist.alignscore, -1 * sum(matchvals), sum(queryclipping),
                #                             thist.tlen, tname, outstart, outend])
                # passingtranscripts.append(
                #     [not passesstringent, sum(queryclipping) - genomicclipping, -1 * thist.alignscore, -1 * sum(matchvals), sum(queryclipping),
                #      thist.tlen, tname, outstart, outend])
                # passingtranscripts.append([not passesstringent, tname, outstart, outend])
                # passingtranscripts.append([sum(queryclipping) - genomicclipping, -1 * thist.alignscore,
                #      -1 * sum(matchvals), sum(queryclipping),
                #      thist.tlen, tname, outstart, outend])
    # order passing transcripts by alignment score
    # then order by amount of query covered
    # then order by amount of transcript covered
    if len(passingtranscripts) > 0:
        passingtranscripts.sort()
        if len(passingtranscripts) == 1 or passingtranscripts[0][:3] != passingtranscripts[1][:3]:
        # print(passingtranscripts)
        # if len(passingtranscripts) == 1 and passingtranscripts[0][:2] == [False, False]:
        # if len(passingtranscripts) == 1 and passingtranscripts[0][0] == False:
        # if len(passingtranscripts) == 1:
            return [passingtranscripts[0][-3:], ]
        else:
            return None


        # passingtranscripts.sort()
        # if testtname in tinfo:
        #     print(passingtranscripts)
        # if args.allow_paralogs:
        #     bestmetrics = passingtranscripts[0][:3]
        #     besttranscripts = []
        #     for t in passingtranscripts:
        #         if t[:3] == bestmetrics: besttranscripts.append(t[-3:])
        #     return besttranscripts
        # else: return [passingtranscripts[0][-3:],]
    else: return None

class IsoAln(object):
    def __init__(self, name=None, p=None, cigar=None, tlen=None, als=None, md=None):
        self.name = name
        self.startpos = p
        self.cigar = cigar
        self.tlen = tlen
        self.alignscore = als
        self.md = md

def parse_sam(args, transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends, readstoclipping):
    lastread = None
    curr_transcripts = {}
    transcripttoreads = {}
    samfile = pysam.AlignmentFile(args.sam, 'r')
    genome = None
    if args.remove_internal_priming:
        genome = pysam.FastaFile(args.transcriptomefasta)
    for read in samfile:
        if read.is_mapped:
            readname = read.query_name
            transcript = read.reference_name
            quality = read.mapping_quality
            if quality >= args.quality:
                # for transcriptome alignment, always take rightmost side on transcript
                if args.remove_internal_priming:
                    intprimannot = transcript_to_exons if args.permissive_last_exons else None
                    notinternalpriming = removeinternalpriming(read.reference_name,
                                                                                       read.reference_start,
                                                                                       read.reference_end, False,
                                                                                       genome, None, intprimannot,
                                                                                       args.intprimingthreshold,
                                                                                       args.intprimingfracAs)
                else: notinternalpriming = True
                if notinternalpriming:
                    pos = read.reference_start
                    try:
                        alignscore = read.get_tag('AS')
                        mdtag = read.get_tag('MD')
                    except KeyError as ex:
                        raise Exception(f"Missing AS or MD tag in alignment of '{read.query_name}' in '{args.sam.name}'") from ex
                    cigar = read.cigartuples
                    tlen = samfile.get_reference_length(transcript)
                    if lastread and readname != lastread:
                        if testtname in curr_transcripts: print('\n', lastread, curr_transcripts.keys())
                        thisclipping = readstoclipping[lastread] if lastread in readstoclipping else None
                        assignedts = get_best_transcript(curr_transcripts, args, transcript_to_exons, transcript_to_bp_ss_index, thisclipping, transcript_to_genomic_ends)
                        if assignedts:
                            for assignedt, gtstart, gtend in assignedts:
                                if assignedt not in transcripttoreads: transcripttoreads[assignedt] = []
                                transcripttoreads[assignedt].append(lastread)

                        curr_transcripts = {}
                    curr_transcripts[transcript] = IsoAln(transcript, pos, cigar, tlen, alignscore, mdtag)
                    lastread = readname
    if lastread:
        thisclipping = readstoclipping[lastread] if lastread in readstoclipping else None
        assignedts = get_best_transcript(curr_transcripts, args, transcript_to_exons, transcript_to_bp_ss_index, thisclipping, transcript_to_genomic_ends)
        if assignedts:
            for assignedt, gtstart, gtend in assignedts:
                if assignedt not in transcripttoreads: transcripttoreads[assignedt] = []
                transcripttoreads[assignedt].append(lastread)

    return transcripttoreads

def write_output(args, transcripttoreads):
    if args.generate_map: mapout = open(args.generate_map, 'w')
    countout = open(args.output, 'wt')
    for t in transcripttoreads:
        if args.generate_map:
            mapout.write(t + '\t' + ','.join(transcripttoreads[t]) + '\n')
        countout.write(t + '\t' + str(len(transcripttoreads[t])) + '\n')


MIN_INSERTION_LEN = 3
HALF_SS_WINDOW_SIZE = 6
NUM_MISTAKES_IN_SS_WINDOW = 2
TRUST_ENDS_WINDOW = 50
LARGE_INDEL_TOLDERANCE = 25
REQ_BP_ALIGNED_IN_EDGE_EXONS = 10
SOFT_CLIPPING_BUFFER = 50

if __name__ == '__main__':
    args = parse_args()
    args = check_args(args)
    transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends = get_annot_info(args)
    readstoclipping = {}
    if args.trimmedreads:
        for line in open(args.trimmedreads):
            rname, clipping = line.rstrip().split('\t')
            readstoclipping[rname] = int(clipping)
    transcripttoreads = parse_sam(args, transcript_to_exons, transcript_to_bp_ss_index, transcript_to_genomic_ends, readstoclipping)
    write_output(args, transcripttoreads)
