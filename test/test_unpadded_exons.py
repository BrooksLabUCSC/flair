from flair.flair_transcriptome import unpadded_exons
from flair.isoform_data import Exon

# a three-exon isoform padded 100 bp at each end for the realignment
PADDED = [Exon(900, 1200), Exon(1500, 1600), Exon(2000, 2400)]


def test_padding_removed():
    assert unpadded_exons(PADDED, (1000, 2300)) == [Exon(1000, 1200), Exon(1500, 1600), Exon(2000, 2300)]


def test_input_not_changed():
    exons = list(PADDED)
    unpadded_exons(exons, (1000, 2300))
    assert exons == PADDED


def test_unpadded_start_past_first_exon_kept():
    # an unpadded start at or past the first exon's end would make it empty or
    # reversed: the first exon is left as it is
    assert unpadded_exons(PADDED, (1200, 2300))[0] == Exon(900, 1200)
    assert unpadded_exons(PADDED, (1250, 2300))[0] == Exon(900, 1200)
    assert unpadded_exons(PADDED, (1250, 2300))[-1] == Exon(2000, 2300)


def test_unpadded_end_before_last_exon_kept():
    assert unpadded_exons(PADDED, (1000, 2000))[-1] == Exon(2000, 2400)
    assert unpadded_exons(PADDED, (1000, 1900))[-1] == Exon(2000, 2400)
    assert unpadded_exons(PADDED, (1000, 1900))[0] == Exon(1000, 1200)


def test_two_exons():
    assert unpadded_exons([Exon(0, 300), Exon(500, 900)], (100, 800)) == [Exon(100, 300), Exon(500, 800)]
