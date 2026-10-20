import pysam
from flair.count_sam_transcripts import deletion_on_splice_motif, process_cigar

FLANK = 'ACCATTCAACCTTCACACTC'


def test_deletion_on_a_canonical_motif():
    for left, right in (('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')):
        tseq = FLANK + left + 'C' * 40 + right + FLANK
        assert deletion_on_splice_motif(tseq, len(FLANK), 44)


def test_deletion_off_a_motif():
    tseq = FLANK + 'CA' + 'C' * 40 + 'TT' + FLANK
    assert not deletion_on_splice_motif(tseq, len(FLANK), 44)


def test_deletion_placed_off_its_motif_in_a_repeat():
    # the intron GT...AG follows a G it repeats at its end, AG|G, so the aligner may
    # equally delete G GT...A, one base to the left, which isn't on the motif
    tseq = FLANK + 'G' + 'GT' + 'C' * 40 + 'AG' + 'C' + FLANK
    start = len(FLANK) + 1
    assert deletion_on_splice_motif(tseq, start, 44)
    assert tseq[start - 1] == tseq[start - 1 + 44]
    assert deletion_on_splice_motif(tseq, start - 1, 44)
    assert not deletion_on_splice_motif(FLANK + 'C' * 44 + FLANK, len(FLANK), 44)


# exons of 300, 100 and 300 bases; tolerated large indels are in the first 250 or
# the last 250 bases
EXONS = [300, 100, 300]


def _indel(tseq, deletion_start, length=60):
    cigar = [(pysam.CMATCH, deletion_start), (pysam.CDEL, length), (pysam.CMATCH, 700 - deletion_start - length)]
    return process_cigar([1] * (700 - length), cigar, 0, EXONS, tseq)[0]


def _tseq(motif_at, left='GT', right='AG', length=60):
    tseq = ['C'] * 700
    tseq[motif_at:motif_at + 2] = left
    tseq[motif_at + length - 2:motif_at + length] = right
    return ''.join(tseq)


def test_terminal_exon_deletion_tolerated_off_a_splice_motif():
    # a skip across double-stranded RNA
    assert not _indel(_tseq(500, 'CA', 'TT'), 500)
    assert not _indel(None, 500)


def test_terminal_exon_deletion_on_a_splice_motif_is_an_intron():
    assert _indel(_tseq(500), 500)
    assert _indel(_tseq(100), 100)


def test_internal_deletion_flagged_either_way():
    assert _indel(_tseq(320, 'CA', 'TT'), 320)
