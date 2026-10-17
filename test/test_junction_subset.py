from flair.isoform_data import Junc
from flair.pycbio.hgdata.bed import BedBlock
from flair.flair_transcriptome import _is_junction_subset

CHAIN = (Junc(100, 200), Junc(300, 400), Junc(500, 600), Junc(700, 800))


def test_contiguous_run_is_subset():
    assert _is_junction_subset(CHAIN[1:3], CHAIN)
    assert _is_junction_subset(CHAIN[:1], CHAIN)
    assert _is_junction_subset(CHAIN[-2:], CHAIN)


def test_non_contiguous_selection_is_not_subset():
    assert not _is_junction_subset((CHAIN[0], CHAIN[2]), CHAIN)


def test_same_or_longer_chain_is_not_subset():
    assert not _is_junction_subset(CHAIN, CHAIN)
    assert not _is_junction_subset(CHAIN + (Junc(900, 1000),), CHAIN)


def test_junction_not_in_chain_is_not_subset():
    assert not _is_junction_subset((Junc(300, 450),), CHAIN)
    assert not _is_junction_subset((Junc(300, 400), Junc(500, 650)), CHAIN)


def test_bed_gaps_as_used_by_quantify():
    gaps = tuple(BedBlock(j.start, j.end) for j in CHAIN)
    assert _is_junction_subset(gaps[1:3], gaps)
    assert not _is_junction_subset((gaps[0], gaps[2]), gaps)
