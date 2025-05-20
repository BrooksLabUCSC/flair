#!/usr/bin/env python3

import sys
import argparse
from bisect import bisect_left
import pysam
import transcriptomic_chimeras




# bamfile = 'smallsimfusionschr22.genomealigned.bam'
#
# transcript_starts = {'chr22':[16592810,]}
# samfile = pysam.AlignmentFile(bamfile, 'rb')
# readtoaligns = {}
# for read in samfile:
#     if read.is_mapped and not read.is_secondary:
#         strand = '-' if read.is_reverse else '+'
#         readname = read.query_name
#         refstart, refend = read.reference_start, read.reference_end
#         qstart, qend = read.query_alignment_start, read.query_alignment_end
#         cigar = read.cigartuples
#         ##hard clipping
#         if cigar[0][0] == 5:
#             qstart += cigar[0][1] ###both start and end are incremented because length of sequence is increased overall
#             qend += cigar[0][1]
#         readlen = read.infer_read_length()
#         if readname not in readtoaligns: readtoaligns[readname] = []
#         if strand == '+':
#             readtoaligns[readname].append([(qstart, refstart), (qend, refend)])
#         else:
#             readtoaligns[readname].append([(readlen - qend, refend), (readlen - qstart, refstart)])
#
#         # print(readname, strand, refstart, refend, qstart, qend)
# for read in readtoaligns:
#     readtoaligns[read] = sorted(readtoaligns[read])
#     print(readtoaligns[read])

###need to flip strand in bam file and generate new fasta of reads
