from types import SimpleNamespace
from flair.isoform_data import Junc, Exon
from flair.flair_transcriptome import (fixed_end_exons, get_spliced_exon_overlaps,
                                       get_unspliced_exon_overlaps, FALLBACK_TERMINAL_EXON_LEN)

# G1: exons 500-1000, 5000-5100 and a long terminal exon 6000-9000
G1_EXONS = {(500, 1000), (5000, 5100), (6000, 9000)}


def _annots():
    return SimpleNamespace(
        spliced_exons={'+': {'G1': {Exon(s, e) for s, e in G1_EXONS}}, '-': {}},
        # an unspliced lncRNA on each strand
        all_annot_SE={'+': [Exon(20000, 22000, 'LNC')], '-': [Exon(20000, 22000, 'LNC-MINUS')]})


def test_fixed_end_exons():
    exons = fixed_end_exons((Junc(1000, 5000), Junc(5100, 6000)))
    assert exons == [Exon(1000 - FALLBACK_TERMINAL_EXON_LEN, 1000), Exon(5000, 5100),
                     Exon(6000, 6000 + FALLBACK_TERMINAL_EXON_LEN)]
    assert fixed_end_exons((Junc(50, 500),))[0] == Exon(0, 50)


def test_fixed_end_exons_with_spliced_gene():
    # G1's own junctions, so all of the fixed-end exons fall in its exons
    exons = fixed_end_exons((Junc(1000, 5000), Junc(5100, 6000)))
    hits = get_spliced_exon_overlaps('+', exons, _annots())
    assert [(g, covered) for covered, g, _ in hits] == [('G1', 300)]
    assert get_spliced_exon_overlaps('-', exons, _annots()) == []


def test_spliced_within_terminal_exon():
    # spliced inside G1's long terminal exon
    exons = fixed_end_exons((Junc(6500, 7000), Junc(7200, 8000)))
    hits = get_spliced_exon_overlaps('+', exons, _annots())
    assert [g for _, g, _ in hits] == ['G1']


def test_unspliced_gene_overlap():
    # an isoform spliced within the unspliced lncRNA: no spliced gene there, but its
    # fixed-end exons lie inside the lncRNA's exon
    exons = fixed_end_exons((Junc(21000, 21300),))
    assert get_spliced_exon_overlaps('+', exons, _annots()) == []
    assert [(g, covered) for covered, g, _ in get_unspliced_exon_overlaps('+', exons, _annots())] == [('LNC', 200)]
    assert [g for _, g, _ in get_unspliced_exon_overlaps('-', exons, _annots())] == ['LNC-MINUS']


def test_unspliced_gene_less_than_half_covered_is_rejected():
    # only the 5' fixed-end exon (100 of 200 bases) is inside the lncRNA
    exons = fixed_end_exons((Junc(21950, 30000),))
    assert get_unspliced_exon_overlaps('+', exons, _annots()) == []
