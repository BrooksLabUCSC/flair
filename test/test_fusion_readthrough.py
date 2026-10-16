"""
check_too_close rejects readthrough: fusion pieces on one chromosome and strand, the
downstream piece starting less than min_dist past the upstream piece in transcript
direction.  Genomic bounds are (chrom, start, end, strand) with start < end on
either strand, in fusion order.
"""
from flair.convert_synthetic_to_genome_bed import check_too_close

MIN_DIST = 100000


def test_plus_strand_readthrough():
    bounds = [('chr14', 37564218, 37579204, '+'), ('chr14', 37595576, 37596008, '+')]
    assert check_too_close(2, bounds, MIN_DIST)


def test_minus_strand_readthrough():
    "TMEM64--LINC01030 in MCF7: 90,444 bp apart, with pieces long enough that their far ends are over min_dist apart"
    bounds = [('chr8', 90645110, 90645718, '-'), ('chr8', 90447230, 90554666, '-')]
    assert check_too_close(2, bounds, MIN_DIST)


def test_far_apart_is_kept():
    assert not check_too_close(2, [('chr1', 1000, 2000, '+'), ('chr1', 202000, 203000, '+')], MIN_DIST)
    assert not check_too_close(2, [('chr1', 202000, 203000, '-'), ('chr1', 1000, 2000, '-')], MIN_DIST)


def test_backward_junction_is_kept():
    "MYO6--SENP6 in MCF7: the downstream gene's first exon joined back into the upstream gene"
    assert not check_too_close(2, [('chr6', 75749188, 75749423, '+'), ('chr6', 75678582, 75715904, '+')], MIN_DIST)
    assert not check_too_close(2, [('chr8', 90447230, 90554666, '-'), ('chr8', 90645110, 90645718, '-')], MIN_DIST)


def test_different_strand_or_chromosome_is_kept():
    assert not check_too_close(2, [('chr1', 1000, 2000, '+'), ('chr1', 5000, 6000, '-')], MIN_DIST)
    assert not check_too_close(2, [('chr1', 1000, 2000, '+'), ('chr2', 5000, 6000, '+')], MIN_DIST)
