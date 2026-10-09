from types import SimpleNamespace
from flair.isoform_data import Junc, Exon, Isoform
from flair.junction_correct import UNKNOWN_STRAND
from flair.flair_transcriptome import (fixed_end_exons, get_spliced_exon_overlaps,
                                       get_unspliced_exon_overlaps, get_gene_name_unknown_strand,
                                       FALLBACK_TERMINAL_EXON_LEN)

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


###
# spliced isoforms of unknown strand get their strand from gene identification
###
def _two_strand_annots():
    # G1 on + as above, and G2 on - with exons 5000-6000 and 8000-8500
    annots = _annots()
    annots.spliced_exons['-'] = {'G2': {Exon(5000, 6000), Exon(8000, 8500)}}
    annots.juncchain_to_transcript = {}
    annots.junc_to_gene = {}
    annots.splice_site_to_genes = {6000: ('G2',)}
    annots.gene_to_strand = {'G1': '+', 'G2': '-', 'LNC': '+', 'LNC-MINUS': '-'}
    annots.gene_to_annot_juncs = {'G1': (Junc(1000, 5000), Junc(5100, 6000)), 'G2': (Junc(6000, 8000),)}
    return annots

def _unknown_iso(juncs, start, end):
    return Isoform('chr1', UNKNOWN_STRAND, juncs, start, end)

def test_unknown_strand_by_exon_overlap():
    # spliced inside G1's long terminal exon, away from G2
    iso = _unknown_iso((Junc(8600, 8700),), 8550, 8800)
    assert get_gene_name_unknown_strand(iso, _two_strand_annots()) == (('G1',), None, '+')
    assert iso.strand == UNKNOWN_STRAND

def test_unknown_strand_by_splice_site_before_exon_overlap():
    # shares G2's splice site at 6000, and its exons also lie in G1's terminal exon
    iso = _unknown_iso((Junc(6000, 7000),), 5900, 7100)
    assert get_gene_name_unknown_strand(iso, _two_strand_annots()) == (('G2',), None, '-')

def test_unknown_strand_genes_on_both_strands():
    # inside the exons of G1 and of G2, sharing no splice site
    iso = _unknown_iso((Junc(8100, 8200),), 8050, 8300)
    assert get_gene_name_unknown_strand(iso, _two_strand_annots()) == (None, None, None)

def test_unknown_strand_no_gene():
    iso = _unknown_iso((Junc(12000, 13000),), 11900, 13100)
    assert get_gene_name_unknown_strand(iso, _two_strand_annots()) == (None, None, None)
