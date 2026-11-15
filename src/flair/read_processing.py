"""Shared read-processing logic for FLAIR modules."""

import logging
import os
from collections import namedtuple
from flair.io_utils import run_program
from flair.isoform_data import Isoform, cigar_clipping
from flair.pycbio.sys import fileOps


# the windows of a region to fetch its reads from: ranges, sorted SeqRanges that don't
# overlap, and starts_only, whether only the reads starting in them or all overlapping them
ReadWindows = namedtuple('ReadWindows', ('ranges', 'starts_only'))


def fetch_reads(bam_file, chrom, start, end, windows=None):
    """The alignments of bam_file overlapping chrom:start-end, or given windows,
    ReadWindows, only those of its ranges, each once, fetched over the ranges alone so
    the reads of the rest of the region are never read"""
    if windows is None:
        yield from bam_file.fetch(chrom, start, end)
        return
    for i, window in enumerate(windows.ranges):
        if windows.starts_only:
            # a partition's window: the reads starting in it, which no other window has
            min_start = window.start
        elif i > 0:
            # report regions' windows: the reads not fetched with the previous window.  The
            # windows are sorted and apart, so a read reaching this one that starts before
            # the previous one's end overlaps it, and was fetched with it, or with an earlier
            # window it also reaches
            min_start = windows.ranges[i - 1].end
        else:
            min_start = None
        for read in bam_file.fetch(window.name, window.start, window.end):
            if min_start is None or read.reference_start >= min_start:
                yield read


def should_process_read(read, region, min_quality, keep_sup, allow_outside_range=False):
    """Check if read passes filtering criteria for processing"""
    if read.mapping_quality < min_quality:
        logging.debug(f"read dropped: low quality ({read.mapping_quality} < {min_quality}): {read.query_name}")
        return False
    if read.is_secondary:
        logging.debug(f"read dropped: secondary alignment: {read.query_name}")
        return False
    if read.is_supplementary and not keep_sup:
        logging.debug(f"read dropped: supplementary alignment: {read.query_name}")
        return False
    if read.reference_name != region.name:
        logging.debug(f"read dropped: wrong reference ({read.reference_name} != {region.name}): {read.query_name}")
        return False
    if not (region.start <= read.reference_start and read.reference_end <= region.end) and not allow_outside_range:
        logging.debug(f"read dropped: outside region range ({read.reference_start}-{read.reference_end} not in {region.start}-{region.end}): {read.query_name}")
        return False
    return True


def get_sequence_from_bed(genome, input_bed, output_fa):
    """Sequence for each BED record, with the (strand) suffix bedtools appends to the
    name removed.  The file is written under a temporary name and renamed when it is
    complete, so a run that stops partway leaves no file a later run would read as a
    finished one."""
    with fileOps.AtomicFileCreate(output_fa) as tmp_fa:
        raw_fa = tmp_fa + '.raw'
        bed_cmd = ('bedtools', 'getfasta', '-nameOnly', '-s', '-split',
                   '-fi', genome,
                   '-bed', input_bed,
                   '-fo', raw_fa)
        run_program([bed_cmd])
        _strip_strand_suffix(raw_fa, tmp_fa)
        os.remove(raw_fa)


def _strip_strand_suffix(in_fa, out_fa):
    "copy a FASTA, dropping the (strand) suffix bedtools appends to each record name"
    with open(in_fa) as in_fh, open(out_fa, 'w') as out_fh:
        for line in in_fh:
            out_fh.write(line.split('(')[0] + '\n' if line[0] == '>' else line)


def add_corrected_read_to_groups(corrected_read, sj_to_ends):
    """Add a corrected read to the junction-to-ends mapping.
    Key is (chrom, juncs) where juncs is () for single-exon reads.
    Single-exon strand is resolved later in group_se_by_overlap."""
    junc_key = (corrected_read.chrom, tuple(sorted(corrected_read.juncs)))
    if junc_key not in sj_to_ends:
        sj_to_ends[junc_key] = Isoform.from_readrec(corrected_read)
    sj_to_ends[junc_key].reads.append(corrected_read)


def generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, chrom, start, end, windows=None):
    """write each primary alignment's clipping at its two ends, given windows, ReadWindows,
    only of the reads in them (fetch_reads)"""
    c = 0
    with open(temp_prefix + '.reads.genomicclipping.tsv', 'w') as clipping_fh:
        for read in fetch_reads(bam_file, chrom, start, end, windows):
            # a record placed in the BAM but not aligned has no cigar, so it has no
            # clipping to report
            if (not read.is_secondary) and (not read.is_supplementary) and (read.cigartuples is not None):
                c += 1
                left_clipping, right_clipping = cigar_clipping(read.cigartuples)
                clipping_fh.write(read.query_name + '\t' + str(left_clipping) + '\t' + str(right_clipping) + '\n')
    return c, temp_prefix + '.reads.genomicclipping.tsv'
