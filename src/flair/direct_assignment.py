"""
Assign reads to first-pass isoforms without the final realignment.

The final realignment aligns reads to the first-pass isoforms and assigns each to
the one it fits best.  A read whose corrected junction chain is a kept isoform's
gets that isoform back from the realignment nearly always, so it is assigned
directly when:

  * its chain is the chain of just one kept spliced isoform;
  * junction correction didn't move its junctions from its alignment's introns,
    and its alignment is clean around them (ReadRec.junctions_moved and
    clean_splice_sites, set by read_correction), which also catches a short exon
    the alignment missed; or its junctions came from its annotation match, which
    checked its splice sites against the annotated transcript, and moved none of
    its first and last ones, as for a missed short exon (junctions_from_annotation,
    and not terminal_junctions_moved).  An annotation match that adds or drops a
    junction at an end of the read's alignment can be a read running unspliced
    past the transcript's splice site;
  * its ends cover the isoform's unique sequence, as count_sam_transcripts
    requires of a subset of another isoform (firstpass.uniquebound.txt).

A spliced read whose chain isn't a kept isoform's can't be assigned by the
realignment, which needs a read to span all of an isoform's junctions, unless its
genome alignment left off a short exon at an end: then it is a truncated part of
a kept isoform, ending at that isoform's next splice site with the rest of the
read clipped.  Other such reads aren't realigned (unassignable_reads).

Single-exon reads are assigned without the realignment by
flair_transcriptome.assign_single_exon_reads.
"""
from flair.count_sam_transcripts import SPLICE_SITE_FLANK


def read_unique_bounds(path):
    """The unique sequence boundaries written with the first-pass isoforms: name ->
    {side: distance from the terminal splice site}, side 0 the transcript's 5' end,
    with the largest distance when a side has several"""
    bounds = {}
    for line in open(path):
        name, sides = line.rstrip('\n').split('\t')
        iso_bounds = bounds.setdefault(name, {})
        for side_dist in sides.split(','):
            side, dist = (int(x) for x in side_dist.split('_'))
            iso_bounds[side] = max(iso_bounds.get(side, 0), dist)
    return bounds


def covers_unique_bounds(isoform, iso_bounds, start, end):
    """Does a read from start to end reach SPLICE_SITE_FLANK past each of the
    isoform's unique sequence boundaries, each a distance into a terminal exon from
    its splice site"""
    for side, dist in iso_bounds.items():
        # the transcript's 5' end is the genomic start on +
        if (side == 0) == (isoform.strand == '+'):
            if start > isoform.juncs[0].start - dist - SPLICE_SITE_FLANK:
                return False
        elif end < isoform.juncs[-1].end + dist + SPLICE_SITE_FLANK:
            return False
    return True


def _juncs_key(juncs):
    return tuple((j.start, j.end) for j in sorted(juncs))


def direct_assignments(corrected_reads, isoforms, unique_bounds):
    """The reads that can be assigned without the realignment, as read name ->
    (isoform, read).  corrected_reads are the spliced reads after junction
    correction, and isoforms the first-pass isoforms."""
    by_chain = {}
    for isoform in isoforms:
        if isoform.juncs != ():
            by_chain.setdefault(_juncs_key(isoform.juncs), []).append(isoform)
    assignments = {}
    for read in corrected_reads:
        as_aligned = not read.junctions_moved and read.clean_splice_sites
        if not (as_aligned or (read.junctions_from_annotation and not read.terminal_junctions_moved)):
            continue
        matches = by_chain.get(_juncs_key(read.juncs), [])
        if len(matches) == 1 and covers_unique_bounds(matches[0], unique_bounds.get(matches[0].name, {}), read.start, read.end):
            assignments[read.name] = (matches[0], read)
    return assignments


def _could_extend_into(read, by_junction, window):
    """Could a read be a truncated part of a kept isoform: its chain a contiguous run
    of the isoform's junctions, its alignment clipped at an end within window of the
    isoform's next splice site, where the genome alignment may have left off a
    short exon"""
    juncs = _juncs_key(read.juncs)
    for chain, i in by_junction.get(juncs[0], ()):
        if chain[i:i + len(juncs)] != juncs:
            continue
        if i > 0 and read.clipping[0] > 0 and abs(read.start - chain[i - 1][1]) <= window:
            return True
        after = i + len(juncs)
        if after < len(chain) and read.clipping[1] > 0 and abs(read.end - chain[after][0]) <= window:
            return True
    return False


def unassignable_reads(corrected_reads, isoforms, window):
    """The names of the spliced reads the realignment can't assign: their chain isn't
    a kept isoform's, and they can't be a truncated part of one (_could_extend_into)"""
    kept = {_juncs_key(isoform.juncs) for isoform in isoforms if isoform.juncs != ()}
    by_junction = {}
    for chain in kept:
        for i, junc in enumerate(chain):
            by_junction.setdefault(junc, []).append((chain, i))
    return {read.name for read in corrected_reads
            if _juncs_key(read.juncs) not in kept and not _could_extend_into(read, by_junction, window)}


def write_direct_assignments(assignments, ends_fh):
    """Write the reads assigned without the realignment as count_sam_transcripts
    writes its ends: in genomic order, each end's junction index, distance from
    that junction, and distance from the transcript's end"""
    for name, (isoform, read) in assignments.items():
        last = len(isoform.juncs) - 1
        fields = (name, isoform.name,
                  0, isoform.juncs[0].start - read.start, max(read.start - isoform.start, 0),
                  last, read.end - isoform.juncs[-1].end, max(isoform.end - read.end, 0))
        ends_fh.write('\t'.join(str(x) for x in fields) + '\n')
