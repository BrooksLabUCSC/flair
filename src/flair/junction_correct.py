"""
Correction of read splice junctions from external evidence.
"""
import logging
from math import inf
from flair import PosRange
from flair.isoform_data import Junc

##
# Notes:
# - We might not want to prefer annotated junction due to possible NAGNAG junctions
#   where annotations often pick just one.
# - A scoring method for junctions based on weighting should be considered
##

##
# somewhat arbitrary sizes to keep from going off ends
# or overlapping other introns.
##
MIN_INTERNAL_EXON_SIZE = 3   # there is one this small!
MIN_TERMINAL_EXON_SIZE = 32

class JunctionCorrector:
    """Correction of read splice sites from orthogonal evidence

       * flank_window - the number of based +/- a read junction to search for
         an supporting intron match.
       * min_read_support - minimum of reads to support a junction.
    """

    def __init__(self, intron_support, flank_window, min_read_support):
        self.intron_support = intron_support
        self.flank_window = flank_window
        self.min_read_support = min_read_support

    @property
    def chroms(self):
        return self.intron_support.chroms

    def overlap_introns(self, chrom, start, end):
        def _filter_intron(intron):
            return intron.annot_supported or (intron.read_support_cnt > self.min_read_support)
        return list(filter(_filter_intron,
                           self.intron_support.overlap_introns(chrom, start, end, self.flank_window)))

    def correct_readrec(self, readrec):
        """Correct a ReadRec's junctions and strand in place from intron support.
        Returns True if corrected, False if there is no support."""
        new_junctions, strand = _correct_junctions(self, readrec)
        if new_junctions is None:
            return False
        readrec.juncs = tuple(Junc(j.start, j.end) for j in new_junctions)
        readrec.strand = strand
        return True


###
# intron support search
###
def _calc_possible_junction_range(readrec, new_junctions):
    """prevent going off ends of read or overlapping small exons"""
    min_start = readrec.start + MIN_TERMINAL_EXON_SIZE
    if len(new_junctions) > 0:
        # adjust for previous intron
        min_start = max(min_start, new_junctions[-1].end + MIN_INTERNAL_EXON_SIZE)
    max_end = readrec.end - MIN_INTERNAL_EXON_SIZE
    return (min_start, max_end)

def _filter_too_close(readrec, new_junctions, intron_hits):
    """drop introns overlapping the previous intron or making a too short
    an exon at ends"""
    min_start, max_end = _calc_possible_junction_range(readrec, new_junctions)
    return list(filter(lambda ih: (ih.start >= min_start) and (ih.end <= max_end),
                       intron_hits))

def _find_best_intron_support(start, end, intron_hits):
    """Pick an intron as `best'  This prefers annotated introns, then
    read-support introns with the larger number of reads.
    """
    closest_introns = _collect_closest_hits(start, end, intron_hits)
    if len(closest_introns) > 1:
        # sort to prefer annotated or then more reads
        closest_introns.sort(key=lambda ci: (ci.annot_supported, ci.read_support_cnt), reverse=True)
    return closest_introns[0]

def _collect_closest_hits(start, end, intron_hits):
    """Return list by introns with closest total distance from ends.
    Multiple introns are return"""
    min_dist = inf
    closest_introns = []
    for intron in intron_hits:
        dist = abs(start - intron.start) + abs(end - intron.end)
        if dist < min_dist:
            # new minimum
            closest_introns.clear()
            min_dist = dist
        if dist <= min_dist:
            closest_introns.append(intron)
    return closest_introns

def _correct_junction(corrector, readrec, start, end, new_junctions):
    """Add a corrected intron junction. Return intron record used or
    None if not supported."""

    intron_hits = corrector.overlap_introns(readrec.chrom, start, end)
    if intron_hits is None:
        logging.debug(f"Read: '{readrec.name}': no intron support for {readrec.chrom}:{start}-{end}")
        return None
    intron_hits = _filter_too_close(readrec, new_junctions, intron_hits)
    if len(intron_hits) == 0:
        logging.debug(f"Read: '{readrec.name}': supporting introns too close to ends or another intron for "
                      f"{readrec.chrom}:{start}-{end}")
        return None
    best_intron = _find_best_intron_support(start, end, intron_hits)
    new_junctions.append(PosRange(best_intron.start, best_intron.end))
    return best_intron

def _determine_strand(readrec, intron_supports):
    """Determine the strand from the IntronSupport objects use for
    splice junction correction, or None if there are conflicts.
    Junctions with strand of '.' don't go into the calculation.
    """
    strand = None
    for intron_support in intron_supports:
        if intron_support.strand != '.':
            if strand is None:
                strand = intron_support.strand
            elif intron_support.strand != strand:
                logging.debug(f"Read: '{readrec.name}': conflicting strands in junction support")
                return None
    if strand is None:
        # all unknown strand
        logging.debug(f"Read: '{readrec.name}': all of the {len(intron_supports)} intron have unknown splice junction so strand can not be determined")
        return None
    return strand

def _correct_junctions(corrector, readrec):
    """Create a list of new junctions for a read and determine
    strand.  Returns (None, None) if can't be corrected or strands are inconsistent"""
    new_junctions = []
    intron_supports = []
    for junc in readrec.juncs:
        intron_support = _correct_junction(corrector, readrec, junc.start, junc.end,
                                           new_junctions)
        if intron_support is None:
            return None, None
        intron_supports.append(intron_support)

    strand = _determine_strand(readrec, intron_supports)
    if strand is None:
        return None, None
    return new_junctions, strand
