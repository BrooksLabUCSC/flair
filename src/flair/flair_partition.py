#!/usr/bin/env python3

# Program that partitions a BED file into non-overlapping regions.
# This is based on the UCSC Browser utility bedPartition.

import argparse
import functools
import os
import re
import logging
import shutil
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
import pipettor
import pysam
from flair import FlairInputDataError, SeqRange, combine_regions
from flair.interval_index import IntervalIndex
from flair.pycbio.sys import cli, fileOps, loggingOps
from flair.pycbio.hgdata.bed import Bed

def check_input_files(bed_files, bam_files, gtf_files, region_bed):
    "fail here rather than inside the sort pipeline, where the message is not the point"
    for f in bed_files + bam_files + gtf_files + ([region_bed] if region_bed is not None else []):
        open(f).close()

def build_parser():
    parser = argparse.ArgumentParser(
        prog='flair_partition',
        description=("Define non-overlapping regions from BED, SAM/BAM, or GTF files."
                     "  Partitions are made across all input files")
    )
    parser.add_argument("--min_partition_items", type=int, default=0,
                        help="Minimum number of input items in a partition")
    parser.add_argument("--part_merge_dist", type=int, default=0,
                        help="Combine adjacent non-overlapping partitions separated by this distance")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of cores for converting inputs, a BAM by chromosome, and sorting")
    parser.add_argument("--bed", dest="bed_files", action="append", default=[],
                        help="Input BED file, maybe compressed.  Maybe repeated")
    parser.add_argument("--bam", dest="bam_files", action="append", default=[],
                        help="Input SAM/BAM file.  Maybe repeated.")
    parser.add_argument("--gtf", dest="gtf_files", action="append", default=[],
                        help="Input GTF file.  Maybe repeated.")
    parser.add_argument("--region_bed",
                        help="BED of regions to restrict the partitions to, such as the genes of "
                        "interest.  Only input items overlapping a region are read, an indexed BAM "
                        "by random access, and the regions are then grown over the items reaching "
                        "them until they stop growing, so each partition holds every item that a "
                        "partition of the whole input would.  Not used with --min_partition_items")
    parser.add_argument("ranges_bed",
                        help="Output ranges BED file, will be compressed if it ends in .gz.  "
                        "It is a BED4 with the number of input items in each partition as a fifth column")
    loggingOps.addCmdOptions(parser, defaultLevel=logging.WARN)
    return parser

def parse_args():
    parser = build_parser()
    args = parser.parse_args()
    loggingOps.setupFromCmd(args)
    if (len(args.bed_files) + len(args.bam_files) + len(args.gtf_files)) == 0:
        parser.error("No input files specified; must have at least one --bam=, --bed= or --gtf= option")
    if (args.region_bed is not None) and (args.min_partition_items > 0):
        parser.error("--min_partition_items does not work with --region_bed: it would count only the "
                     "items read in the regions, not those of the whole input, and so merge partitions "
                     "across the genome between the regions; drop --min_partition_items")
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


class RegionSet:
    """The regions the input is restricted to, merged, with an index by sequence to
    test input items against"""

    def __init__(self, regions):
        self.regions = combine_regions(regions)
        self._by_seq = {}
        for region in self.regions:
            self._by_seq.setdefault(region.name, IntervalIndex()).add(region.start, region.end, region)
        for index in self._by_seq.values():
            index.build()

    def __eq__(self, other):
        return isinstance(other, RegionSet) and (self.regions == other.regions)

    def __len__(self):
        return len(self.regions)

    def overlaps(self, chrom, start, end):
        index = self._by_seq.get(chrom)
        return (index is not None) and (len(index.overlap(start, end)) > 0)

    def write_bed(self, bed_file):
        "the regions as a BED, for samtools --targets-file"
        with open(bed_file, 'w') as fh:
            for region in self.regions:
                fh.write(f"{region.name}\t{region.start}\t{region.end}\n")


def read_region_bed(region_bed):
    "the regions of a BED, merged"
    with fileOps.opengz(region_bed) as fh:
        return RegionSet([SeqRange(item.chrom, item.chromStart, item.chromEnd)
                          for item in read_ranges(fh)])


def check_region_seqs(regions, bam_files):
    "fail on a region naming a sequence a BAM does not have, which would read nothing"
    region_seqs = {region.name for region in regions.regions}
    for bam_file in bam_files:
        with pysam.AlignmentFile(bam_file) as bam:
            missing = sorted(region_seqs - set(bam.references))
        if missing:
            raise FlairInputDataError(f"--region_bed names {', '.join(missing)}, which {bam_file} does "
                                      f"not have; the regions must name sequences the alignments are on")


def sort_beds(nthreads, in_files, sorted_bed):
    """sort BEDs by chrom start and reversed end, which makes it easy to find
    overlapping records.  This writes a file, rather than being read as it
    outputs, so its final merge doesn't run alongside the partitioning and
    exceed the threads"""
    # force ASCII sorting
    env = dict(os.environ, LC_COLLATE="C")
    cmd = ["sort", "-k1,1", "-k2,2n", "-k3,3nr", f"--parallel={nthreads}", "-o", sorted_bed] + list(in_files)
    proc = subprocess.run(cmd, env=env, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=None, stderr=proc.stderr)

###
# Convert inputs to BED, which is sorted for partitioning.  Only the chrom, start,
# and end are used.  Each conversion writes its own file, so they can run in
# parallel.
###
def bed_to_file(bed_file, out_bed, regions=None):
    with fileOps.opengz(bed_file) as bed_fh, open(out_bed, 'w') as out_fh:
        if regions is None:
            shutil.copyfileobj(bed_fh, out_fh)
        else:
            copy_bed_in_regions(bed_fh, out_fh, regions)

def copy_bed_in_regions(bed_fh, out_fh, regions):
    "the records of a BED overlapping the regions, as a BED3"
    for item in read_ranges(bed_fh):
        if regions.overlaps(item.chrom, item.chromStart, item.chromEnd):
            out_fh.write(f"{item.chrom}\t{item.chromStart}\t{item.chromEnd}\n")

def _bam_region(chrom):
    "region for a whole chrom; braces keep a name containing ':' from being parsed as a range"
    return f"{{{chrom}}}" if ':' in chrom else chrom

def bam_to_file(bam_file, out_bed, chrom=None, regions_bed=None):
    "alignments of a BAM, of one chrom of an indexed BAM, or of the regions of a BED, as BED"
    if regions_bed is not None:
        pipettor.run([['samtools', 'view', '-u', '-M', '-L', regions_bed, bam_file],
                      ['bedtools', 'bamtobed', '-i', 'stdin']], stdout=out_bed)
    elif chrom is None:
        pipettor.run(['bedtools', 'bamtobed', '-i', bam_file], stdout=out_bed)
    else:
        pipettor.run([['samtools', 'view', '-u', bam_file, _bam_region(chrom)],
                      ['bedtools', 'bamtobed', '-i', 'stdin']], stdout=out_bed)

def bam_chroms(bam_file):
    """chroms of an indexed BAM that have mapped reads, largest first, so the BAM
    can be converted by chrom in parallel, or None if it has no index"""
    with pysam.AlignmentFile(bam_file) as bam:
        if not (bam.is_bam or bam.is_cram) or not bam.has_index():
            return None
        stats = [s for s in bam.get_index_statistics() if s.mapped > 0]
    return [s.contig for s in sorted(stats, key=lambda s: s.mapped, reverse=True)]

def check_bam_indexed(bam_file):
    "--region_bed reads only the regions of a BAM, which needs its index"
    if bam_chroms(bam_file) is None:
        raise FlairInputDataError(f"--region_bed needs an indexed BAM or CRAM to read the regions of "
                                  f"{bam_file}; index it with `samtools index', or drop --region_bed")


_TRANSCRIPT_ID_RE = re.compile(r'transcript_id\s+"?([^";]+)"?')

@functools.lru_cache(maxsize=None)
def gtf_spans(gtf_file):
    """The span of each GTF transcript, from its exons, by transcript id.  Only the
    span is needed, and converting to BED12 with gtf_to_bed parses the whole GTF into
    records, which took longer than everything else here combined.  Kept, as
    --region_bed converts the inputs more than once"""
    spans = {}  # transcript_id -> [chrom, start, end], start zero-based
    with fileOps.opengz(gtf_file) as fh:
        for line in fh:
            fields = line.split('\t', 8)
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            match = _TRANSCRIPT_ID_RE.search(fields[8])
            if match is None:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            span = spans.get(match.group(1))
            if span is None:
                spans[match.group(1)] = [fields[0], start, end]
            else:
                span[1] = min(span[1], start)
                span[2] = max(span[2], end)
    return spans

def gtf_to_file(gtf_file, out_bed, regions=None):
    "the span of each GTF transcript, as a BED3, restricted to regions if given"
    with open(out_bed, 'w') as out_fh:
        for chrom, start, end in gtf_spans(gtf_file).values():
            if (regions is None) or regions.overlaps(chrom, start, end):
                out_fh.write(f"{chrom}\t{start}\t{end}\n")

class _ThreadSlots:
    "a budget of threads, from which jobs take as many as they use"
    def __init__(self, nthreads):
        self._free = nthreads
        self._total = nthreads
        self._cond = threading.Condition()

    def run(self, nslots, func, *args):
        # a job wanting more than the whole budget runs alone, rather than waiting
        # for threads that will never be free
        nslots = min(nslots, self._total)
        with self._cond:
            self._cond.wait_for(lambda: self._free >= nslots)
            self._free -= nslots
        try:
            return func(*args)
        finally:
            with self._cond:
                self._free += nslots
                self._cond.notify_all()


# converting a BAM chrom is two processes, samtools view and bedtools bamtobed,
# which together use about one and a half cores
_BAM_CHROM_SLOTS = 2

def _bam_jobs(bam_file, nthreads, regions_bed):
    "the conversions of one BAM: its regions, its chroms in parallel, or the whole file"
    if regions_bed is not None:
        check_bam_indexed(bam_file)
        return [(_BAM_CHROM_SLOTS, bam_to_file, (bam_file, None, regions_bed))]
    # with a single thread, a single bamtobed process for the whole BAM
    chroms = bam_chroms(bam_file) if nthreads >= _BAM_CHROM_SLOTS else None
    if chroms is None:
        return [(1, bam_to_file, (bam_file,))]
    return [(_BAM_CHROM_SLOTS, bam_to_file, (bam_file, chrom)) for chrom in chroms]

def convert_inputs(bed_files, bam_files, gtf_files, nthreads, tmp_dir, regions=None):
    """Convert all inputs to BED files in tmp_dir, with an indexed BAM split by
    chrom, using at most nthreads.  Given regions, a RegionSet, only the items
    overlapping one are converted.  Returns the BED files."""
    nthreads = max(1, nthreads)
    regions_bed = None
    if regions is not None:
        regions_bed = os.path.join(tmp_dir, "regions.bed")
        regions.write_bed(regions_bed)
    jobs = []  # (threads used, function, args), slowest first
    for gtf_file in gtf_files:
        jobs.append((1, gtf_to_file, (gtf_file, regions)))
    for bam_file in bam_files:
        jobs.extend(_bam_jobs(bam_file, nthreads, regions_bed))
    for bed_file in bed_files:
        jobs.append((1, bed_to_file, (bed_file, regions)))

    out_beds = [os.path.join(tmp_dir, f"in{i}.bed") for i in range(len(jobs))]
    # the conversions are mostly separate programs, so threads are enough; each
    # waits for the threads it uses to be free
    slots = _ThreadSlots(nthreads)
    with ThreadPoolExecutor(max_workers=nthreads) as executor:
        futures = [executor.submit(slots.run, nslots, func, args[0], out_bed, *args[1:])
                   for (nslots, func, args), out_bed in zip(jobs, out_beds)]
        for future in futures:
            future.result()
    return out_beds

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

class _Range:
    "location of an input item; reading only the first three columns is much faster than parsing BEDs"
    __slots__ = ("chrom", "chromStart", "chromEnd")

    def __init__(self, chrom, chromStart, chromEnd):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd

def read_ranges(fh):
    "generator of _Range objects for the lines of a BED"
    for line in fh:
        if line.isspace() or line.startswith(('#', 'track ', 'browser ')):
            continue
        chrom, start, end = line.split('\t', 3)[:3]
        yield _Range(chrom, int(start), int(end))

def write_partitions(partitions, part_fh, part_counts):
    for bed_part, item_count in partitions:
        part_counts.count(item_count)
        # the item count lets flair order partitions by size
        bed_part.extraCols = (str(item_count),)
        bed_part.write(part_fh)

def scan_partitions(bed_files, bam_files, gtf_files, nthreads, min_partition_items, part_merge_dist,
                    regions, tmp_dir):
    """Convert, sort and partition the inputs, restricted to regions if given,
    yielding each partition as it is read"""
    in_beds = convert_inputs(bed_files, bam_files, gtf_files, nthreads, tmp_dir, regions)
    sorted_bed = os.path.join(tmp_dir, "sorted.bed")
    sort_beds(nthreads, in_beds, sorted_bed)
    with open(sorted_bed) as sorted_fh:
        yield from partition_reader(read_ranges(sorted_fh), min_partition_items, part_merge_dist)

def covered_regions(partitions):
    "the regions the partitions cover, to read the inputs again with"
    return RegionSet([SeqRange(bed_part.chrom, bed_part.chromStart, bed_part.chromEnd)
                      for bed_part, _ in partitions])

def closed_partitions(bed_files, bam_files, gtf_files, nthreads, part_merge_dist, regions, tmp_dir):
    """The partitions over regions, read again over the regions they cover until they
    stop growing.  Each one then holds every input item a partition of the whole input
    would: an item overlapping a grown region overlaps the chain of items that grew it,
    and the regions only grow over items that touch them, so this terminates.
    min_partition_items is not used, as the counts are of the items read, not of the
    whole input"""
    while True:
        logging.info(f"partitioning over {len(regions)} regions")
        partitions = list(scan_partitions(bed_files, bam_files, gtf_files, nthreads, 0,
                                          part_merge_dist, regions, tmp_dir))
        grown = covered_regions(partitions)
        if (len(grown) == 0) or (grown == regions):
            return partitions
        regions = grown

def report_empty_regions(regions, partitions):
    "warn about a region with no input items in it, which makes no partition"
    covered = covered_regions(partitions)
    for region in regions.regions:
        if not covered.overlaps(region.name, region.start, region.end):
            logging.warning(f"no input items in region {region.name}:{region.start}-{region.end}")

def build_partitions(bed_files, bam_files, gtf_files, nthreads, min_partition_items, part_merge_dist,
                     regions, part_fh, part_counts):
    with tempfile.TemporaryDirectory(prefix="flair_partition.") as tmp_dir:
        if regions is None:
            partitions = scan_partitions(bed_files, bam_files, gtf_files, nthreads,
                                         min_partition_items, part_merge_dist, None, tmp_dir)
        else:
            partitions = closed_partitions(bed_files, bam_files, gtf_files, nthreads,
                                           part_merge_dist, regions, tmp_dir)
            report_empty_regions(regions, partitions)
        write_partitions(partitions, part_fh, part_counts)

def report_stats(part_counts):
    logging.info(f"Number of items: {part_counts.item_count}")
    logging.info(f"Number of partitions: {part_counts.part_count}")
    logging.info(f"Min items per partition: {part_counts.min_part_items}")
    logging.info(f"Max items per partition: {part_counts.max_part_items}")
    if part_counts.part_count > 0:
        logging.info(f"Mean items per partition: {part_counts.item_count / part_counts.part_count:.1f}")

def flair_partition(bed_files, bam_files, gtf_files, ranges_bed, nthreads, min_partition_items,
                    part_merge_dist, region_bed=None):
    part_counts = PartitionCounts()
    regions = None
    if region_bed is not None:
        regions = read_region_bed(region_bed)
        check_region_seqs(regions, bam_files)

    with fileOps.opengz(ranges_bed, 'w') as part_fh:
        build_partitions(bed_files, bam_files, gtf_files, nthreads, min_partition_items,
                         part_merge_dist, regions, part_fh, part_counts)
    if (regions is not None) and (part_counts.part_count == 0):
        raise FlairInputDataError(f"no input items in any of the regions in {region_bed}; "
                                  f"check that they name sequences of the input and are in the "
                                  f"zero-based, half-open coordinates of a BED")
    report_stats(part_counts)


def main():
    args = parse_args()
    with cli.ErrorHandler():
        # inside the handler, so a missing input file is reported as one line
        check_input_files(args.bed_files, args.bam_files, args.gtf_files, args.region_bed)
        flair_partition(args.bed_files, args.bam_files, args.gtf_files, args.ranges_bed, args.threads,
                        args.min_partition_items, args.part_merge_dist, region_bed=args.region_bed)


if __name__ == "__main__":
    main()
