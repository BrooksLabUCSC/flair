"""
flair fusion --max_dist_to_tss measures a fusion's 5' end to the 5' gene's transcript
starts upstream of it on the gene's strand, allowing starts a little downstream.
"""
from flair.identify_prelim_fusions import dist_to_upstream_tss, TSS_DOWNSTREAM_BUFFER

# ETV6 (+ strand) transcript starts in GENCODE v49
ETV6_STARTS = [11649623, 11649674, 11649747, 11649820, 11649824, 11649854, 11649944, 11650011,
               11752501, 11885937]


def test_downstream_start_does_not_count():
    "SeraCare ETV6--NTRK3 reads begin in ETV6 exon 5, 16 kb upstream of an alternative start"
    assert dist_to_upstream_tss(11869569, '+', ETV6_STARTS) == 11869569 - 11752501


def test_nearest_upstream_start():
    assert dist_to_upstream_tss(11650500, '+', ETV6_STARTS) == 11650500 - 11650011


def test_start_within_the_buffer_downstream_counts():
    assert dist_to_upstream_tss(11649600, '+', ETV6_STARTS) == 0
    assert dist_to_upstream_tss(11649623 - TSS_DOWNSTREAM_BUFFER - 1, '+', ETV6_STARTS) is None


def test_minus_strand_upstream_is_higher():
    starts = [5000, 20000]
    assert dist_to_upstream_tss(12000, '-', starts) == 8000
    assert dist_to_upstream_tss(4800, '-', starts) == 200
    assert dist_to_upstream_tss(20400, '-', starts) == 0
    assert dist_to_upstream_tss(21000, '-', starts) is None
