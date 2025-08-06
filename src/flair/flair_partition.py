#!/usr/bin/env python3

# Program that partitions a BED file into non-overlapping regions.
# This is based on the UCSC Browser utility bedPartition.

import argparse
import os
import logging
import pipettor
from flair.pycbio.sys import fileOps, loggingOps
from flair.pycbio.hgdata.bed import BedReader, Bed

def parse_args():
    parser = argparse.ArgumentParser(
        description="Define non-overlapping regions from a BED file"
    )
    parser.add_argument("--min_partition_items", type=int, default=0,
                        help="Minimum number of input items in a partition")
    parser.add_argument("-part_merge_dist", type=int, default=0,
                        help="Combine adjacent non-overlapping partitions separated by this distance")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of cores for parallel sorting")
    parser.add_argument("bed_file", help="Input BED file, maybe compressed")
    parser.add_argument("ranges_bed", help="Output ranges BED file, will be compressed if it ends in ..gz")
    loggingOps.addCmdOptions(parser, defaultLevel=logging.WARN)
    args = parser.parse_args()
    loggingOps.setupFromCmd(args)
    return args

class PartitionCounts:
    "some statistics"
    def __init__(self):
        self.part_count = 0
        self.item_count = 0
        self.min_part_items = float('inf')
        self.max_part_items = 0

    def count(self, item_count):
        self.part_count += 1
        self.item_count += item_count
        self.min_part_items = min(self.min_part_items, item_count)
        self.max_part_items = max(self.max_part_items, item_count)


def make_sort_cmd(nthreads, bed_file):
    """create command to sort BED by chrom start and reversed end,
    which makes it easy to find overlapping records"""

    # force ASCII sorting
    os.environ["LC_COLLATE"] = "C"

    # decompresses cmd or cat
    return [fileOps.decompressCmd(bed_file) + [bed_file],
            ["sort", "-k1,1", "-k2,2n", "-k3,3nr", f"--parallel={nthreads}"]]

def sorted_bed_reader(nthreads, bed_file):
    """Generator to sort and read a BED file, returning BED3 records"""
    with pipettor.Popen(make_sort_cmd(nthreads, bed_file)) as sort_fh:
        for bed in BedReader(sort_fh, numStdCols=3):
            yield bed

def same_chrom(bed, bed_part):
    return bed.chrom == bed_part.chrom

def is_overlapped(bed, bed_part):
    "determine if a bed is in the partition"
    return (same_chrom(bed, bed_part) and
            (bed.chromStart < bed_part.chromEnd) and
            (bed.chromEnd > bed_part.chromStart))

def should_merge_min_size(bed, bed_part, item_count, min_partition_items):
    return (same_chrom(bed, bed_part) and
            (item_count < min_partition_items))

def should_merge_adjacent(bed, bed_part, part_merge_dist):
    return (same_chrom(bed, bed_part) and
            (bed_part.chromEnd < bed.chromStart) and
            ((bed.chromStart - bed_part.chromEnd) < part_merge_dist))

def incl_in_partition(bed, bed_part, item_count, min_partition_items, part_merge_dist):
    return (is_overlapped(bed, bed_part) or
            should_merge_min_size(bed, bed_part, item_count, min_partition_items) or
            should_merge_adjacent(bed, bed_part, part_merge_dist))

def partition_build(bed_reader, bed_part, min_partition_items, part_merge_dist):
    item_count = 1  # already have one in bed_part
    while (bed := next(bed_reader, None)) is not None:
        if incl_in_partition(bed, bed_part, item_count, min_partition_items, part_merge_dist):
            bed_part.chromStart = min(bed_part.chromStart, bed.chromStart)
            bed_part.chromEnd = max(bed_part.chromEnd, bed.chromEnd)
            item_count += 1
        else:
            break  # have a partition
    return bed_part, bed, item_count

def make_part_bed(bed, part_count):
    "create a BED 4 with name out a BED, or None if bed is None"
    if bed is None:
        return None
    else:
        return Bed(bed.chrom, bed.chromStart, bed.chromEnd,
                   f"P{part_count}")

def partition_reader(bed_reader, min_partition_items, part_merge_dist):
    # Start by taking the first BED item
    part_count = 0
    bed_part = make_part_bed(next(bed_reader, None), part_count)

    while bed_part is not None:
        bed_part, next_bed, item_count = partition_build(bed_reader, bed_part, min_partition_items, part_merge_dist)
        yield bed_part, item_count
        part_count += 1
        bed_part = make_part_bed(next_bed, part_count)

def build_partitions(bed_reader, min_partition_items, part_merge_dist,
                     part_fh, part_counts):
    for bed_part, item_count in partition_reader(bed_reader,
                                                 min_partition_items, part_merge_dist):
        part_counts.count(item_count)
        bed_part.write(part_fh)

def report_stats(part_counts):
    logging.info(f"Number of items: {part_counts.item_count}")
    logging.info(f"Number of partitions: {part_counts.part_count}")
    logging.info(f"Min items per partition: {part_counts.min_part_items}")
    logging.info(f"Max items per partition: {part_counts.max_part_items}")
    if part_counts.part_count > 0:
        logging.info(f"Mean items per partition: {part_counts.item_count / part_counts.part_count:.1f}")

def flair_partition(bed_file, ranges_bed, nthreads, min_partition_items, part_merge_dist):
    part_counts = PartitionCounts()

    bed_reader = sorted_bed_reader(nthreads, bed_file)
    with fileOps.opengz(ranges_bed, 'w') as part_fh:
        build_partitions(bed_reader, min_partition_items, part_merge_dist,
                         part_fh, part_counts)
    report_stats(part_counts)


def main():
    args = parse_args()
    flair_partition(args.bed_file, args.ranges_bed, args.threads,
                    args.min_partition_items, args.part_merge_dist)


main()
