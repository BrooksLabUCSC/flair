#!/usr/bin/env python3

"""
Detect and filter internal priming artifacts in long-read RNA-seq data.

Internal priming occurs when the oligo-dT primer used during cDNA synthesis
binds to a genomic A-rich region within a transcript rather than at the true
polyA tail, producing a truncated read that falsely appears to end at that
internal site.
"""

from bisect import bisect_left


def checkIsNearAnnotEnd(read3endpos, annotends):
    """
    Check if read3endpos is within 200bp of a known annotated transcript end.
    Uses binary search on the sorted annotends list.  Returns True if the read
    likely ends at a real polyA site rather than an internal priming artifact.
    """
    pos1 = bisect_left(annotends, read3endpos)
    if pos1 == len(annotends):
        return abs(annotends[pos1 - 1] - read3endpos) <= 200
    disttoend = min(abs(annotends[pos1 - 1] - read3endpos), abs(annotends[pos1] - read3endpos))
    return disttoend <= 200


def checkInternalPriming(read3endpos, thischr, genome, reqfreq, threshold):
    """
    Check whether the genomic sequence near the read's 3' end contains an A/T-rich
    stretch that could have caused internal priming.

    Fetches a 60bp window centered on read3endpos and scans sub-windows of
    varying sizes.  For each window, counts the dominant base (A or T) frequency.
    Returns True if any window of length >= threshold has A/T frequency >= reqfreq,
    indicating a genomic polyA/polyT stretch that could cause internal priming.
    """
    genomeseqnearend = genome.fetch(thischr, max(read3endpos - 30, 0), min(read3endpos + 30, genome.get_reference_length(thischr))).upper()
    maxlen = 0
    if len(genomeseqnearend) > threshold * 2:
        halfseqlen = int(len(genomeseqnearend) / 2)
        # Scan windows anchored at the center (read3endpos), growing outward
        # in both directions.  Negative i values produce windows extending left
        # of center; positive i values extend right.  The min/max swap ensures
        # the slice is always [lower:upper] regardless of sign.
        for i in list(range(-1 * halfseqlen, -1 * threshold)) + list(range(threshold, halfseqlen)):
            thisseq = genomeseqnearend[min(i + halfseqlen, halfseqlen): max(i + halfseqlen, halfseqlen)]
            # count whichever of A or T is more frequent (handles both strands)
            thiscount = max(thisseq.count('A'), thisseq.count('T'))
            thisfreq = thiscount / len(thisseq) if len(thisseq) > 0 else 0
            if thisfreq >= reqfreq and len(thisseq) > maxlen:
                maxlen = len(thisseq)
    return maxlen >= threshold

# FIXME: misleading name — returns True when read should be KEPT (no internal
# priming detected), not when internal priming is removed.  A name like
# passesInternalPrimingFilter would be clearer.
def removeinternalpriming(refname, refstart, refend, isrev, genome, annottranscriptends, annotexons, threshold, fracAs):
    """
    Determine if a read should be kept (True) or discarded as internally primed (False).

    Two modes based on alignment type:
    - Transcriptomic alignment (annottranscriptends is empty/None): uses annotexons
      to check if the read's 3' end is near the transcript end, then falls through
      to genomic A-rich check.
    - Genomic alignment (annottranscriptends provided): checks if read's 3' end is
      near an annotated transcript end (real polyA site).

    In both cases, if no evidence of a real polyA site is found, the read is
    discarded as internally primed (returns False).
    """
    read3endpos = refend if not isrev else refstart

    if not annottranscriptends:
        # --- Transcriptomic alignment path ---
        # Check if read's 3' end is near the transcript end based on exon structure
        if annotexons and refname in annotexons:
            theseexons = annotexons[refname]
            # multi-exon: keep if 3' end is within the last exon
            # (sum of exon lengths minus last exon = start of last exon in transcript coords)
            if len(theseexons) > 1 and read3endpos > sum(theseexons) - theseexons[-1]:
                return True
            # single-exon: keep if 3' end is within 200bp of transcript end
            elif len(theseexons) == 1 and read3endpos >= theseexons[0] - 200:
                return True

        # Fall-through: read didn't pass the near-end check (or had no exon annotation).
        # Check genomic sequence for A-rich stretch.  If no A-rich stretch found, keep the read.
        # FIXME: this uses transcript coordinates (refname/read3endpos) as if they were
        # genomic coordinates, which is incorrect for transcriptomic alignments.
        # The refname check (from issue #629) is a workaround: if refname isn't a
        # chromosome in the genome, skip the check and keep the read.
        # FIXME: checking for known chromosome should happen much earlier
        if ((refname not in genome.references) or
                (not checkInternalPriming(read3endpos, refname, genome, fracAs, threshold))):
            return True

    elif annottranscriptends and refname in annottranscriptends:
        # --- Genomic alignment path ---
        # If read's 3' end is near a known transcript end, this is likely a real
        # polyA site, so keep the read even if the genomic sequence is A-rich.
        # FIXME: if not near an annotated end, falls through to return False
        # without checking genomic A-richness — asymmetric with the transcriptomic
        # path above, which does check.  Reads on chromosomes with no annotated
        # ends also fall through to False without any A-richness check.
        isnearannotend = checkIsNearAnnotEnd(read3endpos, annottranscriptends[refname])
        if isnearannotend:
            return True
    return False
