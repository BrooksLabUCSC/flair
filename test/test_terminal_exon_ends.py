import pytest
from flair import FlairInputDataError
from flair.terminal_exon_ends import (GeneTerminalExonEnds, TerminalExonEnds,
                                      NORM_END_EXTRA_LEN, NORM_END_SS_WINDOW)

CHROM_LEN = 100000


def _gene():
    gene = GeneTerminalExonEnds('g1')
    # two transcripts sharing splice sites, one reaching further at each end
    gene.add(1000, 1200, 5000, 5500)
    gene.add(900, 1200, 5000, 5400)
    return gene


def test_furthest_ends_same_splice_site():
    assert _gene().furthest_ends(1200, 5000) == (900, 5500)


def test_splice_sites_within_window_are_grouped():
    gene = _gene()
    gene.add(800, 1200 + NORM_END_SS_WINDOW - 1, 5000 - NORM_END_SS_WINDOW, 5600)
    assert gene.furthest_ends(1200, 5000) == (800, 5600)


def test_splice_sites_outside_window_are_not_grouped():
    gene = _gene()
    gene.add(800, 1200 + NORM_END_SS_WINDOW, 5000 - NORM_END_SS_WINDOW - 1, 5600)
    assert gene.furthest_ends(1200, 5000) == (900, 5500)


def test_normalized_ends_are_padded():
    assert _gene().normalized_ends(1200, 5000, CHROM_LEN) == (900 - NORM_END_EXTRA_LEN, 5500 + NORM_END_EXTRA_LEN)


def test_normalized_ends_are_clamped_to_chrom():
    gene = GeneTerminalExonEnds('g1')
    gene.add(10, 200, 900, 990)
    assert gene.normalized_ends(200, 900, 1000) == (0, 1000)


def test_single_exon_transcripts_only_record_gene():
    ends = TerminalExonEnds()
    ends.add_transcript('g1', '+', [(100, 500)])
    assert ends.fetch('g1', '+').first_exon_starts == {}


def test_gene_on_both_strands_is_an_error():
    ends = TerminalExonEnds()
    ends.add_transcript('g1', '+', [(100, 200), (300, 400)])
    with pytest.raises(FlairInputDataError):
        ends.add_transcript('g1', '-', [(100, 200), (300, 400)])
