from flair.isoform_data import Exon
from flair.flair_transcriptome import unique_bounds_past

# a subset's first exon ends at the splice site 1000 shared with a superset's
# internal exon 900-1000, and its last exon starts at 2000, shared with a superset's
# internal exon 2000-2100: each boundary is 100 from the splice site
BOUNDS = [(0, 100), (1, 100), (0, 100)]


def test_only_boundaries_the_terminal_exons_extend_past():
    assert unique_bounds_past(BOUNDS, Exon(850, 1000), Exon(2000, 2050), '+') == ['0_100']
    assert unique_bounds_past(BOUNDS, Exon(950, 1000), Exon(2000, 2150), '+') == ['1_100']
    # ending exactly at the superset exon's edge leaves no unique sequence
    assert unique_bounds_past(BOUNDS, Exon(900, 1000), Exon(2000, 2100), '+') == []


def test_minus_strand_sides_are_transcript_oriented():
    assert unique_bounds_past(BOUNDS, Exon(850, 1000), Exon(2000, 2150), '-') == ['1_100', '0_100']
