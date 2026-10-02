"""
Choose the reads that need aligning to the annotated transcripts.

The alignment to annotated transcripts corrects a read's junctions and ends from
the transcript it matches.  For many reads that can't change anything, and the
read can go straight to junction correction instead:

  * unspliced reads, whose match isn't used, as it gives no junctions;
  * spliced reads whose junctions are a contiguous run of one annotated
    transcript's junctions, and whose genome alignment around each splice site
    has no indels and at most MAX_SPLICE_SITE_MISMATCHES mismatches.  Junction
    correction keeps those junctions, as annotated introns match them exactly.
    A short exon missed by the genome alignment shows up as an indel or
    mismatches next to a junction, so those reads are still aligned.

Either way, a read with an end within END_SPLICE_SITE_DIST of an annotated splice
site is still aligned: its end may continue across that junction, which the
transcript alignment would find.
"""
from bisect import bisect_left
import pysam

# bases on each side of a splice site in which the genome alignment must have no
# indels, as count_sam_transcripts' HALF_SS_WINDOW_SIZE
SPLICE_SITE_WINDOW = 6
# mismatches allowed in the window on each side of a splice site
MAX_SPLICE_SITE_MISMATCHES = 1
# a read end at most this far from an annotated splice site is still aligned
END_SPLICE_SITE_DIST = 10

_REF_AND_QUERY = (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF)


class AnnotJunctionIndex:
    """The junctions of annotated transcripts, for finding reads whose junctions
    are a contiguous run of one transcript's, and read ends near splice sites."""
    def __init__(self, annots):
        self._tx_juncs = {}
        self._junc_index = {}
        for key, exons in annots.transcript_to_exons.items():
            juncs = tuple((exons[i].end, exons[i + 1].start) for i in range(len(exons) - 1))
            self._tx_juncs[key] = juncs
            for i, junc in enumerate(juncs):
                self._junc_index.setdefault(junc, []).append((key, i))
        # intron starts and ends, which a read's right and left ends may be near
        self._intron_starts = sorted({j[0] for j in self._junc_index})
        self._intron_ends = sorted({j[1] for j in self._junc_index})

    def is_annotated_chain(self, juncs):
        "are juncs a contiguous run of one annotated transcript's junctions"
        for key, i in self._junc_index.get(juncs[0], ()):
            if self._tx_juncs[key][i:i + len(juncs)] == juncs:
                return True
        return False

    @staticmethod
    def _near(positions, pos, dist):
        i = bisect_left(positions, pos - dist)
        return i < len(positions) and positions[i] <= pos + dist

    def end_near_splice_site(self, start, end, dist=END_SPLICE_SITE_DIST):
        """is a read's start near an intron end (the start of an exon), or its end
        near an intron start, of any annotated transcript"""
        return self._near(self._intron_ends, start, dist) or self._near(self._intron_starts, end, dist)


def _alignment_blocks(read):
    """CIGAR operations as (op, reference start, reference end, query start), with
    an insertion at the reference position it precedes"""
    blocks = []
    rpos, qpos = read.reference_start, 0
    for op, length in read.cigartuples:
        if op in _REF_AND_QUERY:
            blocks.append((op, rpos, rpos + length, qpos))
            rpos += length
            qpos += length
        elif op in (pysam.CDEL, pysam.CREF_SKIP):
            blocks.append((op, rpos, rpos + length, qpos))
            rpos += length
        elif op == pysam.CINS:
            blocks.append((op, rpos, rpos, qpos))
            qpos += length
        elif op == pysam.CSOFT_CLIP:
            qpos += length
    return blocks


def read_introns(read):
    "the read's introns, from the reference skips in its CIGAR"
    return tuple((b[1], b[2]) for b in _alignment_blocks(read) if b[0] == pysam.CREF_SKIP)


def splice_sites_cleanly_aligned(read, introns, genome,
                                 window=SPLICE_SITE_WINDOW, max_mismatches=MAX_SPLICE_SITE_MISMATCHES):
    """Is the read aligned across the window on each side of each splice site with
    no insertion or deletion, and at most max_mismatches mismatches per side."""
    blocks = _alignment_blocks(read)
    seq = read.query_sequence
    for start, end in introns:
        for lo, hi in ((start - window, start), (end, end + window)):
            covered = mismatches = 0
            for op, bstart, bend, bquery in blocks:
                if op == pysam.CREF_SKIP:
                    continue
                if op == pysam.CINS:
                    if lo <= bstart <= hi:
                        return False
                    continue
                ov_lo, ov_hi = max(bstart, lo), min(bend, hi)
                if ov_lo >= ov_hi:
                    continue
                if op == pysam.CDEL:
                    return False
                ref = genome.fetch(read.reference_name, ov_lo, ov_hi).upper()
                qoff = bquery + ov_lo - bstart
                mismatches += sum(1 for k in range(ov_hi - ov_lo) if seq[qoff + k].upper() != ref[k])
                covered += ov_hi - ov_lo
            if covered < hi - lo or mismatches > max_mismatches:
                return False
    return True


def needs_annotation_alignment(read, junc_index, genome):
    "does aligning this read to the annotated transcripts possibly change its correction"
    if junc_index.end_near_splice_site(read.reference_start, read.reference_end):
        return True
    introns = read_introns(read)
    if len(introns) == 0:
        return False
    return not (junc_index.is_annotated_chain(introns)
                and splice_sites_cleanly_aligned(read, introns, genome))


def reads_skipping_annotation_alignment(bam_file, region, annots, genome):
    """Names of the primary alignments in region that don't need aligning to the
    annotated transcripts"""
    junc_index = AnnotJunctionIndex(annots)
    skip = set()
    for read in bam_file.fetch(region.name, region.start, region.end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not needs_annotation_alignment(read, junc_index, genome):
            skip.add(read.query_name)
    return skip
