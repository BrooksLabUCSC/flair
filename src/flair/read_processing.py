"""Shared read-processing logic for FLAIR modules."""

import logging
import os
import pysam
import pipettor
from flair.isoform_data import Isoform


def should_process_read(read, region, min_quality, keep_sup, allow_secondary, allow_outside_range=False):
    """Check if read passes filtering criteria for processing"""
    if read.mapping_quality < min_quality:
        logging.debug(f"read dropped: low quality ({read.mapping_quality} < {min_quality}): {read.query_name}")
        return False
    if read.is_secondary and not allow_secondary:
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
    name removed."""
    bed_cmd = ('bedtools', 'getfasta', '-nameOnly', '-s', '-split',
               '-fi', genome,
               '-bed', input_bed,
               '-fo', output_fa)
    pipettor.run([bed_cmd])
    # output_fa + '.fixed', not output_fa.split('.fa')[0]: the split truncates at the
    # first '.fa' anywhere in the path, so a directory named x.fa or an output called
    # out.fasta.fa put the temporary file somewhere else entirely
    fixed_fa = output_fa + '.fixed'
    with open(output_fa) as in_fh, open(fixed_fa, 'w') as out_fh:
        for line in in_fh:
            out_fh.write(line.split('(')[0] + '\n' if line[0] == '>' else line)
    os.replace(fixed_fa, output_fa)


def add_corrected_read_to_groups(corrected_read, sj_to_ends):
    """Add a corrected read to the junction-to-ends mapping.
    Key is (chrom, juncs) where juncs is () for single-exon reads.
    Single-exon strand is resolved later in group_se_by_overlap."""
    junc_key = (corrected_read.chrom, tuple(sorted(corrected_read.juncs)))
    if junc_key not in sj_to_ends:
        sj_to_ends[junc_key] = Isoform.from_readrec(corrected_read)
    sj_to_ends[junc_key].reads.append(corrected_read)


def generate_genomic_alignment_read_to_clipping_file(temp_prefix, bam_file, chrom, start, end):
    c = 0
    # use both soft and hard-clipping because some alignments (secondary, maybe supplementary) can be hard-clipped
    clipping_types = (pysam.CIGAR_OPS.CSOFT_CLIP, pysam.CIGAR_OPS.CHARD_CLIP)
    with open(temp_prefix + '.reads.genomicclipping.txt', 'w') as clipping_fh:
        for read in bam_file.fetch(chrom, start, end):
            if not read.is_secondary and not read.is_supplementary:
                c += 1
                name = read.query_name
                cigar = read.cigartuples
                left_clipping, right_clipping = 0, 0
                if cigar[0][0] in clipping_types:
                    left_clipping = cigar[0][1]
                if cigar[-1][0] in clipping_types:
                    right_clipping = cigar[-1][1]
                clipping_fh.write(name + '\t' + str(left_clipping) + '\t' + str(right_clipping) + '\n')
    return c, temp_prefix + '.reads.genomicclipping.txt'
