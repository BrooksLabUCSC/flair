"""
flair fusion --max_dist_to_tss measures a fusion's 5' end to the 5' gene's transcript
starts upstream of it on the gene's strand, allowing starts a little downstream.
"""
from flair.identify_prelim_fusions import dist_to_upstream_tss, fusion_tss_distance, TSS_DOWNSTREAM_BUFFER

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


# ESR1 (+ strand) transcript starts near MCF7's ESR1--CCDC170 reads, GENCODE v49
ESR1_STARTS = [151651283, 151656690, 151690495, 151690529]


def test_median_end_at_a_tss_counts():
    "ESR1--CCDC170: most reads start at the 151690495 TSS, but 3 start 14 kb upstream"
    assert fusion_tss_distance((151676006, 151690531), '+', ESR1_STARTS) == 2
    assert dist_to_upstream_tss(151676006, '+', ESR1_STARTS) == 151676006 - 151656690


def test_distal_end_near_a_tss_counts():
    "5' truncated reads: the most upstream ones reach the promoter, the median doesn't"
    assert fusion_tss_distance((151690600, 151720000), '+', ESR1_STARTS) == 151690600 - 151690529


def test_no_upstream_start_for_either_end():
    assert fusion_tss_distance((151600000, 151610000), '+', ESR1_STARTS) is None
