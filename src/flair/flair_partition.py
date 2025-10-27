#!/usr/bin/env python3

# Program that partitions a BED file into non-overlapping regions.
# This is based on the UCSC Browser utility bedPartition.

import argparse
import os
import logging
import shutil
import subprocess
import pipettor
from flair.pycbio.sys import fileOps, loggingOps
from flair.pycbio.hgdata.bed import BedReader, Bed

def check_input_files(bed_files, bam_files):
    for f in bed_files + bam_files:
        open(f).close()

def parse_args():
    parser = argparse.ArgumentParser(
        description=("Define non-overlapping regions from BED, SAM/BAM, or GTF files."
                     "  Partitions are made across all input files")
    )
    parser.add_argument("--min_partition_items", type=int, default=0,
                        help="Minimum number of input items in a partition")
    parser.add_argument("-part_merge_dist", type=int, default=0,
                        help="Combine adjacent non-overlapping partitions separated by this distance")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of cores for parallel sorting")
    parser.add_argument("--bed", dest="bed_files", action="append", default=[],
                        help="Input BED file, maybe compressed.  Maybe repeated")
    parser.add_argument("--bam", dest="bam_files", action="append", default=[],
                        help="Input SAM/BAM file.  Maybe repeated.")
    parser.add_argument("--gtf", dest="gtf_files", action="append", default=[],
                        help="Input GTF file.  Maybe repeated.")
    parser.add_argument("ranges_bed", help="Output ranges BED file, will be compressed if it ends in .gz")
    loggingOps.addCmdOptions(parser, defaultLevel=logging.WARN)
    args = parser.parse_args()
    loggingOps.setupFromCmd(args)
    if (len(args.bed_files) + len(args.bam_files)) == 0:
        parser.error("No input files specified; must have at least one --bam= or --bed= option")
    check_input_files(args.bed_files, args.bam_files)
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


def start_sort_process(nthreads):
    """create process to sort BEDs by chrom start and reversed end,
    which makes it easy to find overlapping records.  Returns
    process object."""

    # force ASCII sorting
    os.environ["LC_COLLATE"] = "C"
    cmd = ["sort", "-k1,1", "-k2,2n", "-k3,3nr", f"--parallel={nthreads}"]
    proc = subprocess.Popen(cmd,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True)
    return proc

def finish_sort_process(sort_proc):
    rc = sort_proc.wait()
    if rc != 0:
        err = sort_proc.stderr.read()
        subprocess.CalledProcessError(rc, sort_proc.args, output=None, stderr=err)

def copy_bed_to_sort(bed_file, to_sort_fh):
    with fileOps.opengz(bed_file) as bed_fh:
        shutil.copyfileobj(bed_fh, to_sort_fh)

def copy_bam_to_sort(bam_file, to_sort_fh):
    cmd = ['bedtools', 'bamtobed', '-i', bam_file]
    with pipettor.Popen(cmd) as bed_fh:
        shutil.copyfileobj(bed_fh, to_sort_fh)

def copy_gtf_to_sort(gtf_file, to_sort_fh):
    cmd = ['gtf_to_bed', '--include_gene', gtf_file, '/dev/stdout']
    with pipettor.Popen(cmd) as bed_fh:
        shutil.copyfileobj(bed_fh, to_sort_fh)

def copy_input_to_sort(bed_files, bam_files, gtf_files, to_sort_fh):
    for bed_file in bed_files:
        copy_bed_to_sort(bed_file, to_sort_fh)
    for bam_file in bam_files:
        copy_bam_to_sort(bam_file, to_sort_fh)
    for gtf_file in gtf_files:
        copy_gtf_to_sort(gtf_file, to_sort_fh)

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

def write_partitions(from_sort_fh, min_partition_items, part_merge_dist,
                     part_fh, part_counts):
    bed_reader = BedReader(from_sort_fh, numStdCols=3)

    for bed_part, item_count in partition_reader(bed_reader,
                                                 min_partition_items, part_merge_dist):
        part_counts.count(item_count)
        bed_part.write(part_fh)

def build_partitions(bed_files, bam_files, gtf_files, nthreads, min_partition_items, part_merge_dist,
                     part_fh, part_counts):
    # NOTE: this takes advantage of sort not writing anything to stdout until
    # stdin is closed.  Don't do this at home.

    sort_proc = start_sort_process(nthreads)
    copy_input_to_sort(bed_files, bam_files, gtf_files, sort_proc.stdin)
    sort_proc.stdin.close()
    write_partitions(sort_proc.stdout, min_partition_items, part_merge_dist,
                     part_fh, part_counts)
    finish_sort_process(sort_proc)

def report_stats(part_counts):
    logging.info(f"Number of items: {part_counts.item_count}")
    logging.info(f"Number of partitions: {part_counts.part_count}")
    logging.info(f"Min items per partition: {part_counts.min_part_items}")
    logging.info(f"Max items per partition: {part_counts.max_part_items}")
    if part_counts.part_count > 0:
        logging.info(f"Mean items per partition: {part_counts.item_count / part_counts.part_count:.1f}")

def flair_partition(bed_files, bam_files, gtf_files, ranges_bed, nthreads, min_partition_items, part_merge_dist):
    part_counts = PartitionCounts()

    with fileOps.opengz(ranges_bed, 'w') as part_fh:
        build_partitions(bed_files, bam_files, gtf_files,
                         nthreads, min_partition_items, part_merge_dist, part_fh, part_counts)
    report_stats(part_counts)


def main():
    args = parse_args()
    flair_partition(args.bed_files, args.bam_files, args.gtf_files, args.ranges_bed, args.threads,
                    args.min_partition_items, args.part_merge_dist)


if __name__ == "__main__":
    main()
