from flair.isoform_data import Exon
from flair.predictProductivity import calc_genomic_end_pos

# a two-exon transcript: exons [100, 200) and [300, 350), 150 bases
PLUS_EXONS = [Exon(100, 200), Exon(300, 350)]   # 5' to 3'
MINUS_EXONS = PLUS_EXONS[::-1]                   # 5' to 3' on minus strand
SIZES = [100, 50]


def test_plus_strand_inside_exon():
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 30, '+') == 130
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 120, '+') == 320


def test_minus_strand_inside_exon():
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 30, '-') == 320
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 80, '-') == 170


def test_stop_codon_at_transcript_end():
    # orf_end_pos equal to the transcript length used to give None
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 150, '+') == 350
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 150, '-') == 100
