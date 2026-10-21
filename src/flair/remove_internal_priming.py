#!/usr/bin/env python3

"""
Detect and filter internal priming artifacts in long-read RNA-seq data.

Internal priming occurs when the oligo-dT primer used during cDNA synthesis
binds to a genomic A-rich region within a transcript rather than at the true
polyA tail, producing a truncated read that falsely appears to end at that
internal site.
"""

from bisect import bisect_left


def check_is_near_annot_end(read_3end_pos, annot_ends):
    """
    Check if read_3end_pos is within 200bp of a known annotated transcript end.
    Uses binary search on the sorted annot_ends list.  Returns True if the read
    likely ends at a real polyA site rather than an internal priming artifact.
    """
    pos1 = bisect_left(annot_ends, read_3end_pos)
    if pos1 == len(annot_ends):
        return abs(annot_ends[pos1 - 1] - read_3end_pos) <= 200
    dist_to_end = min(abs(annot_ends[pos1 - 1] - read_3end_pos), abs(annot_ends[pos1] - read_3end_pos))
    return dist_to_end <= 200


def check_internal_priming(read_3end_pos, chrom, genome, req_freq, threshold):
    """
    Check whether the genomic sequence near the read's 3' end contains an A/T-rich
    stretch that could have caused internal priming.

    Fetches a 60bp window centered on read_3end_pos and scans sub-windows of
    varying sizes.  For each window, counts the dominant base (A or T) frequency.
    Returns True if any window of length >= threshold has A/T frequency >= req_freq,
    indicating a genomic polyA/polyT stretch that could cause internal priming.
    """
    genome_seq_near_end = genome.fetch(chrom, max(read_3end_pos - 30, 0), min(read_3end_pos + 30, genome.get_reference_length(chrom))).upper()
    maxlen = 0
    if len(genome_seq_near_end) > threshold * 2:
        half_seq_len = int(len(genome_seq_near_end) / 2)
        # Scan windows anchored at the center (read_3end_pos), growing outward
        # in both directions.  Negative i values produce windows extending left
        # of center; positive i values extend right.  The min/max swap ensures
        # the slice is always [lower:upper] regardless of sign.
        for i in list(range(-1 * half_seq_len, -1 * threshold)) + list(range(threshold, half_seq_len)):
            seq = genome_seq_near_end[min(i + half_seq_len, half_seq_len): max(i + half_seq_len, half_seq_len)]
            # count whichever of A or T is more frequent (handles both strands)
            count = max(seq.count('A'), seq.count('T'))
            freq = count / len(seq) if len(seq) > 0 else 0
            if freq >= req_freq and len(seq) > maxlen:
                maxlen = len(seq)
    return maxlen >= threshold


def passes_internal_priming_filter(refname, refstart, refend, isrev, genome, annot_transcript_ends, annot_exons, threshold, frac_as):
    """
    Determine if a read should be kept (True) or discarded as internally primed (False).

    Two modes based on alignment type:
    - Transcriptomic alignment (annot_transcript_ends is empty/None): uses annot_exons
      to check if the read's 3' end is near the transcript end, then falls through
      to genomic A-rich check.
    - Genomic alignment (annot_transcript_ends provided): checks if read's 3' end is
      near an annotated transcript end (real polyA site).

    In both cases, if no evidence of a real polyA site is found, the read is
    discarded as internally primed (returns False).
    """
    read_3end_pos = refend if not isrev else refstart

    if not annot_transcript_ends:
        # --- Transcriptomic alignment path ---
        # Check if read's 3' end is near the transcript end based on exon structure
        if annot_exons and refname in annot_exons:
            exon_lens = annot_exons[refname]
            # multi-exon: keep if 3' end is within the last exon
            # (sum of exon lengths minus last exon = start of last exon in transcript coords)
            if len(exon_lens) > 1 and read_3end_pos > sum(exon_lens) - exon_lens[-1]:
                return True
            # single-exon: keep if 3' end is within 200bp of transcript end
            elif len(exon_lens) == 1 and read_3end_pos >= exon_lens[0] - 200:
                return True

        # Fall-through: read didn't pass the near-end check (or had no exon annotation).
        # Check genomic sequence for A-rich stretch.  If no A-rich stretch found, keep the read.
        # FIXME: this uses transcript coordinates (refname/read_3end_pos) as if they were
        # genomic coordinates, which is incorrect for transcriptomic alignments.
        # The refname check (from issue #629) is a workaround: if refname isn't a
        # chromosome in the genome, skip the check and keep the read.
        # FIXME: checking for known chromosome should happen much earlier
        if ((refname not in genome.references) or
                (not check_internal_priming(read_3end_pos, refname, genome, frac_as, threshold))):
            return True

    elif annot_transcript_ends and refname in annot_transcript_ends:
        # --- Genomic alignment path ---
        # If read's 3' end is near a known transcript end, this is likely a real
        # polyA site, so keep the read even if the genomic sequence is A-rich.
        # FIXME: if not near an annotated end, falls through to return False
        # without checking genomic A-richness — asymmetric with the transcriptomic
        # path above, which does check.  Reads on chromosomes with no annotated
        # ends also fall through to False without any A-richness check.
        is_near_annot_end = check_is_near_annot_end(read_3end_pos, annot_transcript_ends[refname])
        if is_near_annot_end:
            return True
    return False
