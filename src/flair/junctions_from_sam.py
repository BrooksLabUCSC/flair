#!/usr/bin/env python3
# preProcess_getASEventReadCounts.py
# Author: Angela Brooks, Alison Tang (modifications for FLAIR)
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""Takes a SAM file as input and will create three files that are used as input
to getASEventReadCounts.py

- *_junctions.bed - Junctions in BED format
- *_genome_reads.txt - Text file containing genome reads.  Reads are assumed to be
                 the same length.  Format is query_name, chr, start, end
Additional files if paired-end reads exist:
- *_paired_end_junctions2qname.txt - File containing the start and end of
  introns that are inferred from junctions.  Only from paired-end reads.
- *_paired_end_as_single_genome.txt - Text file containing genome reads that
                                      are from paired-end reads.
Modified for FLAIR-correct input - Alison Tang 2019
"""
# *_intron_exon_junctions.txt - Text file that is output from
# intron_exon_jctn_counts.py.  Used for quantifying intron retention,
# alternative donor, and alternative acceptor events.

import sys
import optparse
import math
import re
import pysam
import logging
from flair import FlairInputDataError

#############
# CONSTANTS #
#############
#IE_SCRIPT = "/h/angela/bin/intron_exon_jctn_counts.py"

DEF_NAME = "junctions_from_sam"
K_RGB_STR = "0,0,0"
N_RGB_STR = "0,60,120"
# Corresponds to reads equally split between 2 different offsets.
DEF_CONFIDENCE = 1.0

# Some strange junction alignments occur towards the beginning of the chromosome
# which gives errors later. These are likely artifacts, anyway. Junction
# alignments at the beginning of the chromosome must have at least this number
# of bases in the resulting alignments
DEF_BEG_JCN_OVERHANG = 25

SHELL = "/bin/tcsh"

MAX_CHAR = 60

CIGAR_OPER = {0:"M",
              1:"I",
              2:"D",
              3:"N",
              4:"S",
              5:"H",
              6:"P",
              7:"=",
              8:"X"}

#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            self.print_help()
            raise FlairInputDataError(f"{option} option not supplied")


class JcnInfo:
    """
    This object holds all information necessary for junction information
    The junctions will be in BED format
    """
    def __init__(self, name, chr, chromStart, chromEnd, strand, first_block, second_block, intron_start, intron_end, multiJcnBlock=None):

        self.name = name

        self.chr = chr
        if not chr.startswith("chr"):
            self.chr = "chr" + chr

        self.leftmost_start = chromStart
        self.rightmost_end = chromEnd

        total_length = chromEnd - chromStart
        self.second_block_start = total_length - second_block

        self.longest_first_block = first_block
        self.longest_second_block = second_block

        self.strand = strand

        self.block_list = [first_block]

        # Used for calculating entropy when reads span multiple jcns.
        if multiJcnBlock:
            self.multiJcnBlock_list = [multiJcnBlock]
        else:
            self.multiJcnBlock_list = []

        self.intron_start = intron_start
        self.intron_end = intron_end

    def updateJcnInfo(self, name, chr, chromStart, chromEnd, strand, first_block, second_block, intron_start, intron_end, verbosity=False, multiJcnBlock=None):
        # Check that name, chromosome, strand are the same
        if not chr.startswith("chr"):
            chr = "chr" + chr

        if self.name != name:
            raise FlairInputDataError(f"Not the same name of the junctions: {self.name}, {name}")
        if self.chr != chr:
            raise FlairInputDataError(f"Error with chromosome: {self.name}, {name}")

        if self.strand != strand and self.strand != '.':
            if verbosity:
                print("Warning: Junction reads aligning to both strands for %s. Setting strand to '.'" % self.name)
            self.strand = "."

        if self.intron_start != intron_start:
            raise FlairInputDataError(f"Error with intron start: {self.intron_start}, {intron_start}")
        if self.intron_end != intron_end:
            raise FlairInputDataError(f"Error with intron end: {self.intron_end}, {intron_end}")

        # Check blocks start
        if (self.rightmost_end - self.longest_second_block) != (chromEnd - second_block):
            raise FlairInputDataError(f"Error with second block in  {self.name}, {name}")
        if (self.leftmost_start + self.longest_first_block) != (chromStart + first_block):
            raise FlairInputDataError(f"Error with first block {self.name}, {name}")

        if chromStart < self.leftmost_start:
            self.leftmost_start = chromStart
            self.longest_first_block = first_block

            # Change second block start
            total_length = self.rightmost_end - self.leftmost_start
            self.second_block_start = total_length - self.longest_second_block

        if chromEnd > self.rightmost_end:
            self.rightmost_end = chromEnd
            self.longest_second_block = second_block

            # Change second block start
            total_length = self.rightmost_end - self.leftmost_start
            self.second_block_start = total_length - self.longest_second_block

        self.block_list.append(first_block)

        if multiJcnBlock:
            self.multiJcnBlock_list.append(multiJcnBlock)


###############
# END CLASSES #
###############
########
# MAIN #
########
def main():

    opt_parser = OptionParser()

    # Add Options. Required options should have default=None
    opt_parser.add_option("-s",
                          dest="sam_file",
                          type="string",
                          help="""SAM/BAM file of read alignments to junctions and
                                  the genome. More than one file can be listed,
                                  but comma-delimited, e.g file_1.bam,file_2.bam""",
                          default=None)
    opt_parser.add_option("--unique",
                          dest="unique_only",
                          action="store_true",
                          help="""Only keeps uniquely aligned reads. Looks at NH
                                  tag to be 1 for this information.""",
                          default=False)
    opt_parser.add_option("-n",
                          dest="name",
                          type="string",
                          help="""Name prefixed used for output BED file.
                          Default=%s""" % DEF_NAME,
                          default=DEF_NAME)
    opt_parser.add_option("-l",
                          dest="read_length",
                          type="int",
                          help="""Expected read length if all reads should be of
                                  the same length""",
                          default=None)
    opt_parser.add_option("-c",
                          dest="confidence_score",
                          type="float",
                          help="""The mininmum entropy score a junction
                                  has to have in order to be considered
                                  confident. The entropy score =
                                  -Shannon Entropy. Default=%s""" % DEF_CONFIDENCE,
                          default=DEF_CONFIDENCE)
    opt_parser.add_option("-j",
                          dest="forced_junctions",
                          type="string",
                          help="""File containing intron coordinates
                                  that correspond to junctions that will be
                                  kept regardless of the confidence score.""",
                          default=None)
    opt_parser.add_option("-v",
                          dest="verbose",
                          action='store_true',
                          help="""Will run the program with junction strand ambiguity messages""",
                          default=False)
    # opt_parser.add_option("-o",
    #                       dest="output_dir",
    #                       type="string",
    #                       help="Directory to place all output files.",
    #                       default=None)
#   opt_parser.add_option("-p",
#                         dest="paired_end_exists",
#                         action="store_true",
#                         help="Flag to indicate there are paired-end reads.",
#                         default=False)

    (options, args) = opt_parser.parse_args()

    # validate the command line arguments
    opt_parser.check_required("-s")
    # opt_parser.check_required("-o")

    sam_files = []
    sam_file_names = options.sam_file.split(",")
    for sam_file in sam_file_names:
        if sam_file.endswith(".sam"):
            sam_files.append(open(sam_file))
        elif sam_file.endswith(".bam"):
            sam_files.append(pysam.Samfile(sam_file, "rb"))
    #        sam_file = gzip.open(options.sam_file)
        else:
            opt_parser.print_help()
            raise FlairInputDataError("Error in -s: Expecting .sam or .bam file.")

    unique_only = options.unique_only

    name = options.name

    read_len = options.read_length
    confidence_score = options.confidence_score

    # out_dir = options.output_dir

#    paired_end_exists = options.paired_end_exists
    # Still working on the paired end version
    paired_end_exists = False

    # Create ouptut directory if it doesn't exist.
    # if not os.path.exists(out_dir):
    #     os.mkdir(out_dir)
    #     out_dir = os.path.abspath(out_dir)
    #     print "Creating output directory: %s" % out_dir
    # else:
    #     out_dir = os.path.abspath(out_dir)

    # if not out_dir.endswith("/"):
    #     out_dir += "/"

    forced_jcns = None
    if options.forced_junctions:
        forced_jcns = getForcedJunctions(options.forced_junctions)

    # Creating output files.
    junction_bed_file_name = "%s_junctions.bed" % (name)
    junction_bed_file = open(junction_bed_file_name, "w")
    # genome_file_name = "%s%s_genome_reads.txt.gz" % (out_dir, name)
    # genome_file = gzip.open(genome_file_name, "wb")

    # Paired-end reads exist.
    if paired_end_exists:
        pass
        # jcn2qname_file_name = "%s%s_paired_end_junctions2qname.txt" % (out_dir, name)
        # jcn2qname_file = open(jcn2qname_file_name, "w")

        # paired_end_genome_file_name = "%s%s_paired_end_as_single_genome.txt" % (out_dir, name)
        # paired_end_genome_file = open(paired_end_genome_file_name, "w")

    # Creating names of files for intron/exon junction counts.
#   ie_coord_filename = "%s%s_intron_exon_junction_coords.out" % (out_dir, name)
#   ie_read_assoc_filename = "%s%s_intron_exon_junction_coords_w_read.out" % (out_dir, name)
#   ie_read_count_filename = "%s%s_intron_exon_junction_counts.txt" % (out_dir, name)

    # Dictionary to find JcnInfo to update as junction reads are encountered
    jcn2JcnInfo = {}

    jcn2pairedReads = {}

    deletionFlag = False
    insertionFlag = False
    softclipFlag = False
    hardclipFlag = False

    truncation_warn = False

    print("Parsing sam/bam file")
    for sam_file in sam_files:
        for line in sam_file:
            if options.sam_file.endswith(".bam"):
                # I realize this is bad style, but the original code was written
                # for parsing SAM files. So less updating is necessary if I convert
                # the AlignedRead object back to a SAM line
                line = convert2SAMLine(sam_file, line)

            line = formatLine(line)

            # Ignore headers
            if line.startswith("@"):
                continue

            sam_elems = line.split("\t")

            if len(sam_elems) < 11:
                raise FlairInputDataError("Error in SAM file: Expecting at least 11 columns in SAM file.")

            q_name = sam_elems[0]

            flag = int(sam_elems[1])
            chr = sam_elems[2]

            # Ignore unmapped
            if chr == "*":
                continue

            if not chr.startswith("chr"):
                chr = "chr" + chr

            chr_start = int(sam_elems[3])

            cigar = sam_elems[5]

            # Ignore unmapped
            if cigar == "*":
                continue

            tags = sam_elems[11:]
            if unique_only:
                isMultiFlag = isMultiMapped(tags)
                if isMultiFlag is None:
                    print("Could not find tag for %s" % q_name)
                    continue

                if isMultiFlag:
                    continue

            m_count = cigar.count("M")

            i_count = cigar.count("I")
            if i_count > 0:
                if not insertionFlag:
                    print("Not supporting insertions, yet e.g., %s" % cigar)
                insertionFlag = True
                continue

            d_count = cigar.count("D")
            if d_count > 0:
                if not deletionFlag:
                    print("Not supporting deletions, yet e.g., %s" % cigar)
                deletionFlag = True
                continue

            s_count = cigar.count("S")
            if s_count > 0:
                if not softclipFlag:
                    print("Not supporting softclipping, yet e.g., %s" % cigar)
                softclipFlag = True
                continue

            h_count = cigar.count("H")
            if h_count > 0:
                if not hardclipFlag:
                    print("Not supporting hardclipping, yet e.g., %s" % cigar)
                hardclipFlag = True
                continue

            if '=' in cigar:
                cigar = extended_to_simple_cigar(cigar)
                m_count = cigar.count('M')

            # Check if it is a genome read
            if m_count == 1: # A GENOME READ
                pass
                this_read_len = int(cigar.rstrip("M"))

                if read_len:
                    if this_read_len != read_len:
                        raise ValueError(f"Expecting reads of length: {read_len} not {this_read_len}")

                #chr_end = chr_start + this_read_len - 1

                #genome_line = "%s\t%s\t%d\t%d\n" % (q_name, chr, chr_start, chr_end)
                # if paired_end_exists:
                #     if isPairedRead(flag):
                #         paired_end_genome_file.write(genome_line)
                #     else:
                #         genome_file.write(genome_line)
                # else:
                #     genome_file.write(genome_line)

            else: # A JUNCTION READ
                n_count = cigar.count("N")

                if n_count == 0:
                    print("Expecting a junction read: %s" % cigar)
                    continue

                n_split = cigar.split("N")

                # Get the downstrm length
                downstr_len = int(n_split.pop().rstrip("M"))

                # A list to hold the information about each intron. Used in cases
                # where a read aligns to multiple junctions.
                introns_info = []

                first_chr_start = chr_start

                # Get first intron information which also will be used for the
                # upstream length
                upstr_len, intron_len = map(int, n_split.pop(0).split("M"))
                introns_info.append((chr_start, intron_len, None))
                # Updating the chr_start
                chr_start = chr_start + upstr_len + intron_len

                # Get remaining information from additional introns if there are
                # any
                # 3rd element is used in calculating entropy
                for remaining_intron in n_split:
                    exon_len, intron_len = map(int, remaining_intron.split("M"))
                    introns_info.append((chr_start + exon_len - upstr_len,
                                         intron_len,
                                         chr_start - first_chr_start))
                    chr_start = chr_start + exon_len + intron_len

                jcn_tag = None
                jcn_strand = None
                for tag in tags:
                    if tag.startswith("Y0"):
                        jcn_tag = tag[5:]
                        if len(jcn_tag) > MAX_CHAR:
                            jcn_tag = jcn_tag[:MAX_CHAR]
                            if not truncation_warn:
                                truncation_warn = True
                                print("Warning: Y0 tags truncated to %d characters." % MAX_CHAR)

                    if tag.startswith("XS"):
                        tag, almost_strand = tag.split("A")
                        jcn_strand = almost_strand.lstrip(":")

                        if jcn_strand != "+" and jcn_strand != "-" and jcn_strand != ".":
                            raise ValueError(f"Error in strand information for: {line}")

                # Now insert all introns into jcn dictionary.
                for intron_info in introns_info:
                    chr_start = intron_info[0]
                    intron_len = intron_info[1]

                    total_len = upstr_len + intron_len + downstr_len

                    jcn_str = None
                    if not jcn_tag:
                        # Create a junction string based on the 1-based junction
                        # coordinate
                        jcn_str = "%s:%d-%d" % (chr,
                                                chr_start + upstr_len,
                                                chr_start + upstr_len + intron_len - 1)
                    else:
                        # Need to make multiple jcn_str for each intron
                        jcn_str = "%s|%s:%d-%d" % (jcn_tag,
                                                   chr,
                                                  chr_start + upstr_len,
                                                  chr_start + upstr_len + intron_len - 1)

                    # Check for odd junctions that are aligned toward the
                    # beginning of the chromosome, which causes problems later
                    if (chr_start + upstr_len) < DEF_BEG_JCN_OVERHANG:
                        continue

                    if not jcn_strand:
                        jcn_strand = "."

                    # Get BED format information
                    chromStart = chr_start - 1
                    chromEnd = chromStart + total_len

                    # Now add junction to dictionary
                    if jcn_str in jcn2JcnInfo:
                        jcn2JcnInfo[jcn_str].updateJcnInfo(jcn_str,
                                                           chr,
                                                           chromStart, chromEnd,
                                                           jcn_strand, upstr_len,
                                                           downstr_len,
                                                           chr_start + upstr_len,
                                                           chr_start + upstr_len + intron_len - 1,
                                                           options.verbose,
                                                           intron_info[2])
                    else:
                        jcn2JcnInfo[jcn_str] = JcnInfo(jcn_str,
                                                       chr,
                                                       chromStart, chromEnd,
                                                       jcn_strand, upstr_len,
                                                       downstr_len,
                                                       chr_start + upstr_len,
                                                       chr_start + upstr_len + intron_len - 1,
                                                       intron_info[2])

                    # Add read to paired-end dictionary:
                    if paired_end_exists:
                        if isPairedRead(flag):
                            jcn_str = "%s:%d-%d" % (chr,
                                                    chr_start + upstr_len,
                                                    chr_start + upstr_len + intron_len - 1)

                            if jcn_str in jcn2pairedReads:
                                jcn2pairedReads[jcn_str].append(q_name)
                            else:
                                jcn2pairedReads[jcn_str] = [q_name]

    # Close the genome file:
    # genome_file.close()
    # if paired_end_exists:
    #     paired_end_genome_file.close()

    # Done processing all SAM lines, now make junction BED file
    print("Making Junction BED File")
#    bed_header = "track name=\"%s_jcn_counts\" itemRgb=\"On\" useScore=1\n" % name
#    junction_bed_file.write(bed_header)
    confident_jcns = set([])
    strandFlag = False
    for jcn_str in jcn2JcnInfo:
        if isConfidentJunction(jcn2JcnInfo[jcn_str].block_list,
                               jcn2JcnInfo[jcn_str].multiJcnBlock_list,
                               confidence_score,
                               jcn2JcnInfo[jcn_str].chr,
                               jcn2JcnInfo[jcn_str].intron_start,
                               jcn2JcnInfo[jcn_str].intron_end,
                               forced_jcns):

            confident_jcns.add("%s:%d-%d" % (jcn2JcnInfo[jcn_str].chr,
                                             jcn2JcnInfo[jcn_str].intron_start,
                                             jcn2JcnInfo[jcn_str].intron_end))
#            # If junction name has ":::", then the type of junction is known
#            if ":::" in jcn_str:
#                rgb_str = getRGB(jcn_str)
#            else:
#                rgb_str = K_RGB_STR

            intron_left = jcn_str[jcn_str.find(':')+1:jcn_str.find('-')]
            intron_right = jcn_str[jcn_str.find('-')+1:]
            if jcn2JcnInfo[jcn_str].strand in {'+', '-'}:
                strandFlag = True
                jcn_strand = jcn2JcnInfo[jcn_str].strand
            else:
                jcn_strand = '.'
            # bed_line = "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t2\t%s\t%s\n" % (jcn2JcnInfo[jcn_str].chr,
            #                                                                 jcn2JcnInfo[jcn_str].leftmost_start,
            #                                                                 jcn2JcnInfo[jcn_str].rightmost_end,
            #                                                                 jcn_str,
            #                                                                 len(jcn2JcnInfo[jcn_str].block_list),
            #                                                                 jcn2JcnInfo[jcn_str].strand,
            #                                                                 jcn2JcnInfo[jcn_str].leftmost_start,
            #                                                                 jcn2JcnInfo[jcn_str].rightmost_end,
            #                                                                 rgb_str,
            #                                                                 ",".join([repr(jcn2JcnInfo[jcn_str].longest_first_block),
            #                                                                           repr(jcn2JcnInfo[jcn_str].longest_second_block)]),
            #                                                                 ",".join(["0",repr(jcn2JcnInfo[jcn_str].second_block_start)]))
            bed_line = '\t'.join([jcn2JcnInfo[jcn_str].chr,
                                                intron_left,
                                                str(int(intron_right)-1),
                                                '.',
                                                str(min(1000, len(jcn2JcnInfo[jcn_str].block_list))),
                                                 jcn_strand]) + '\n'
            junction_bed_file.write(bed_line)

    junction_bed_file.close()
    # Print out junction to read file for paired end reads
    # if paired_end_exists:
    #     for jcn_str in jcn2pairedReads:
    #         if jcn_str in confident_jcns:
    #             outline = "%s\t%s\n" % (jcn_str, ",".join(jcn2pairedReads[jcn_str]))
    #             jcn2qname_file.write(outline)
    #     jcn2qname_file.close()

    if strandFlag == False:
        logging.info('WARNING, no stranded junctions were found.')


############
# END_MAIN #
############

#############
# FUNCTIONS #
#############


def convert2CIGAR(cigar_tuple):
    cigar_str = ""
    for tup in cigar_tuple:
        cigar_str += "%d%s" % (tup[1],
                               CIGAR_OPER[tup[0]])
    return cigar_str


def convert2SAMLine(bam_file, read_obj):
    # Checks due to depracation of variables in pysam
    try:
        rnext = read_obj.rnext
        pnext = read_obj.pnext
        tlen = read_obj.tlen
    except:
        rnext = read_obj.mrnm
        pnext = read_obj.mpos
        tlen = read_obj.isize

    try:
        rname = bam_file.getrname(read_obj.tid)
    except:
        rname = None

    samline = "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s" % (read_obj.qname,
                                                              read_obj.flag,
                                                              rname if rname else "*",
                                                              read_obj.pos + 1,
                                                              read_obj.mapq,
                                                              convert2CIGAR(read_obj.cigar) if read_obj.cigar else "*",
                                                              rnext if rnext else "*",
                                                              pnext + 1 if pnext else 0,
                                                              tlen if tlen else 0,
                                                              read_obj.seq,
                                                              read_obj.qual)

    if read_obj.tags:
        if read_obj.tags == []:
            return samline + "\n"

        tags = convert2Tags(read_obj.tags)
        if tags != "":
            samline += "\t%s" % tags

    return samline + "\n"


def convert2Tags(tag_tuples):
    tag_list = []
    for tag_tup in tag_tuples:
        if tag_tup[0] == "Y0":
            tag_list.append("Y0:Z:%s" % tag_tup[1])
        if tag_tup[0] == "XS":
            tag_list.append("XS:A:%s" % tag_tup[1])
        if tag_tup[0] == "NH":
            tag_list.append("NH:i:%d" % tag_tup[1])

    return "\t".join(tag_list)


def convertFlag(flag):
    # Bitwise operation to get strand information
    if flag & 16:
        return "-"
    else:
        return "+"


def extended_to_simple_cigar(cigar):
    matches = re.findall(r'([0-9]+)([A-Z]|=)', cigar)
    newcigar = ''
    prevop = '='
    prevnum = 0
    for m in matches:
        num, op = int(m[0]), m[1]
        if op not in ['=', 'X', 'N']:
            raise ValueError(f'Unexpectted operator in extended cigar:{cigar}')
        if op in ['=', 'X'] and prevop in ['=', 'X']:
            prevnum += num
        else:
            if prevop != 'N':
                newcigar += str(prevnum)+'M'
            else:
                newcigar += str(prevnum)+'N'
            prevnum = num
            prevop = op
    if prevop != 'N':
        newcigar += str(prevnum)+'M'
    else:
        newcigar += str(prevnum)+'N'
    return newcigar


def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line


def getForcedJunctions(forced_junction_file):
    """
    File has chr start end for introns that would like to be kept in analysis.
    set will be "chr_start_end"
    """
    jcn_file = open(forced_junction_file)
    forced_junctions = set([])

    for line in jcn_file:
        line = formatLine(line)

        lineList = line.split("\t")

        # Slicing list in case strand is included in the file
        if len(lineList) < 3:
            raise FlairInputDataError("Problem with forced junction file. Needs to be tab-delimited: chr start end")

        if not lineList[0].startswith("chr"):
            lineList[0] = "chr" + lineList[0]

        forced_junctions.add("_".join(lineList[0:3]))

    jcn_file.close()

    return forced_junctions


def getRGB(jcn_str):
    elems = jcn_str.split(":::")

    for elem in elems:
        if elem == "A":
            return K_RGB_STR

    return N_RGB_STR


def getPos2Count(blocklist):
    pos2count = {}

    for block in blocklist:
        if block in pos2count:
            pos2count[block] += 1
        else:
            pos2count[block] = 1

    return pos2count


def getShannonIndex(pos2countDict, totalCount):
    """
    Calculates a Shannon-Wiener Index for the positions and counts.  It treats
    each position like a unique species.
    """

    if totalCount == 0:
        return 0

    summation = 0

    for pos in pos2countDict:
        p = float(pos2countDict[pos]) / totalCount

        if p == 0.0:
            continue

        summation += p * math.log(p, 2)

    if summation == 0:
        return 0

    return -summation


def isConfidentJunction(block_list, multiJcnBlockList, confidence_score,
                        jcn_chr, jcn_start, jcn_end, forced_jcns):
    """
    Confidence is given as the entropy score associated with each junction.
    """
    if forced_jcns:
        jcn_str = "%s:%d-%d" % (jcn_chr,
                                jcn_start,
                                jcn_end)

        if jcn_str in forced_jcns:
            return True

    all_blocks = list(block_list)
    all_blocks.extend(multiJcnBlockList)
    pos2count = getPos2Count(all_blocks)
    totalCount = len(all_blocks)

    entropy = getShannonIndex(pos2count, totalCount)

    if entropy >= confidence_score:
        return True
    else:
        return False


def isMultiMapped(tags):
    NHtag_found = False
    for tag in tags:
        if tag.startswith("NH"):
            NHtag_found = True
            if tag.split(":")[-1] != "1":
                return True

    if not NHtag_found:
        return None

    return False


def isPairedRead(flag):
    # If read is first in a pair
    if flag & 64:
        return True

    # If read is second in a pair
    if flag & 128:
        return True

    return False
#################
# END FUNCTIONS #
#################


if __name__ == "__main__":
    main()
