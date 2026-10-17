from types import SimpleNamespace
import pysam
from flair.isoform_data import Exon
from flair.annotation_precheck import (AnnotJunctionIndex, read_introns, splice_sites_cleanly_aligned,
                                       needs_annotation_alignment)

# T1: exons 100-200, 300-400, 500-600; T2 skips the middle exon
ANNOTS = SimpleNamespace(transcript_to_exons={
    ('T1', 'G1'): (Exon(100, 200), Exon(300, 400), Exon(500, 600)),
    ('T2', 'G1'): (Exon(100, 200), Exon(500, 600)),
})
GENOME_SEQ = ''.join('ACGT'[(i * 7 + i // 3) % 4] for i in range(1000))


_HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 1000}]})


class _Genome:
    def fetch(self, chrom, start, end):
        return GENOME_SEQ[start:end]


def _read(start, cigar, mutate=()):
    "an alignment matching the genome, with query bases changed at reference positions in mutate"
    seg = pysam.AlignedSegment(_HEADER)
    seg.query_name = 'r'
    seg.reference_id = 0
    seg.reference_start = start
    seg.cigartuples = cigar
    seq, rpos = [], start
    for op, length in cigar:
        if op == pysam.CMATCH:
            for k in range(length):
                base = GENOME_SEQ[rpos + k]
                seq.append('A' if (rpos + k in mutate and base != 'A') else ('C' if rpos + k in mutate else base))
            rpos += length
        elif op in (pysam.CDEL, pysam.CREF_SKIP):
            rpos += length
        elif op in (pysam.CINS, pysam.CSOFT_CLIP):
            seq.append('G' * length)
    seg.query_sequence = ''.join(seq)
    return seg


def test_annotated_chain():
    index = AnnotJunctionIndex(ANNOTS)
    assert index.is_annotated_chain(((200, 300), (400, 500)))
    assert index.is_annotated_chain(((400, 500),))
    assert index.is_annotated_chain(((200, 500),))
    # each junction is annotated, but not as a contiguous run of one transcript's
    assert not index.is_annotated_chain(((200, 500), (400, 500)))
    assert not index.is_annotated_chain(((210, 300),))


def test_end_near_splice_site():
    index = AnnotJunctionIndex(ANNOTS)
    assert index.end_near_splice_site(305, 450)       # start near the 300 exon start
    assert index.end_near_splice_site(250, 395)       # end near the 400 exon end
    assert not index.end_near_splice_site(250, 450)
    assert not index.end_near_splice_site(150, 560)   # transcript ends aren't splice sites


def _spliced(mutate=(), extra=()):
    # 150-200, intron 200-300, 300-350
    return _read(150, [(pysam.CMATCH, 50), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 50)] if not extra else extra, mutate)


def test_clean_splice_sites():
    seg = _spliced()
    assert read_introns(seg) == ((200, 300),)
    assert splice_sites_cleanly_aligned(seg, read_introns(seg), _Genome())


def test_one_mismatch_per_site_allowed():
    assert splice_sites_cleanly_aligned(_spliced(mutate={197}), ((200, 300),), _Genome())
    assert not splice_sites_cleanly_aligned(_spliced(mutate={196, 198}), ((200, 300),), _Genome())
    # a mismatch beyond the window doesn't count
    assert splice_sites_cleanly_aligned(_spliced(mutate={190, 191}), ((200, 300),), _Genome())


def test_indel_near_splice_site():
    ins = _read(150, [(pysam.CMATCH, 48), (pysam.CINS, 3), (pysam.CMATCH, 2), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 50)])
    assert not splice_sites_cleanly_aligned(ins, read_introns(ins), _Genome())
    dele = _read(150, [(pysam.CMATCH, 50), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 2), (pysam.CDEL, 2), (pysam.CMATCH, 46)])
    assert not splice_sites_cleanly_aligned(dele, read_introns(dele), _Genome())
    # an insertion right at the junction, as a missed short exon can be
    atj = _read(150, [(pysam.CMATCH, 50), (pysam.CINS, 9), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 50)])
    assert not splice_sites_cleanly_aligned(atj, read_introns(atj), _Genome())
    # an indel away from the splice sites is fine
    far = _read(150, [(pysam.CMATCH, 20), (pysam.CDEL, 3), (pysam.CMATCH, 27), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 50)])
    assert splice_sites_cleanly_aligned(far, read_introns(far), _Genome())


def test_window_must_be_covered():
    # only 3 bases aligned before the junction
    short = _read(197, [(pysam.CMATCH, 3), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 50)])
    assert not splice_sites_cleanly_aligned(short, read_introns(short), _Genome())


def test_needs_annotation_alignment():
    index = AnnotJunctionIndex(ANNOTS)
    assert not needs_annotation_alignment(_spliced(), index, _Genome())
    # unspliced reads aren't aligned, unless an end is near a splice site
    assert not needs_annotation_alignment(_read(220, [(pysam.CMATCH, 60)]), index, _Genome())
    assert needs_annotation_alignment(_read(150, [(pysam.CMATCH, 48)]), index, _Genome())
    # junction not annotated
    novel = _read(150, [(pysam.CMATCH, 50), (pysam.CREF_SKIP, 120), (pysam.CMATCH, 50)])
    assert needs_annotation_alignment(novel, index, _Genome())
    # annotated junction, but mismatches at a splice site
    assert needs_annotation_alignment(_spliced(mutate={300, 301}), index, _Genome())
