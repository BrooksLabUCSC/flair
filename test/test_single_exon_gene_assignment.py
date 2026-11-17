"""Assignment of a single-exon isoform to an annotated single-exon gene, which reads
every overlapping annotation rather than the few nearest one in a sorted list"""
from types import SimpleNamespace
from flair.annotation_data import AnnotData
from flair.gene_identification import get_single_exon_gene_overlaps
from flair.isoform_data import Exon


def _annots(*annot_exons):
    "annotation of single-exon transcripts on the + strand, given as (start, end, gene)"
    annots = AnnotData()
    annots.all_annot_se['+'] = sorted(Exon(start, end, gene) for start, end, gene in annot_exons)
    annots.index_single_exon_annots()
    return annots


def _isoform(start, end):
    return SimpleNamespace(exons=[Exon(start, end)], strand='+')


def test_the_annotation_containing_the_isoform_is_found_past_the_ones_inside_it():
    # the gene the isoform lies in starts well before five annotations inside it, so a
    # search of the entries nearest the isoform's end never reached it
    annots = _annots((1000, 9000, 'LONG'), (2000, 2100, 'A'), (3000, 3100, 'B'),
                     (4000, 4100, 'C'), (5000, 5100, 'D'), (6000, 6100, 'E'))
    assert list(get_single_exon_gene_overlaps('+', _isoform(1500, 8500), annots)) == ['LONG']


def test_a_gene_covering_less_than_half_of_the_isoform_is_not_a_hit():
    annots = _annots((1000, 2000, 'SHORT'))
    assert get_single_exon_gene_overlaps('+', _isoform(1500, 8500), annots) == {}


def test_an_isoform_in_no_annotation_has_no_gene():
    annots = _annots((1000, 2000, 'A'), (3000, 4000, 'B'))
    assert get_single_exon_gene_overlaps('+', _isoform(2200, 2800), annots) == {}


def test_the_fractions_are_of_the_isoform_and_of_the_annotation():
    annots = _annots((1000, 3000, 'A'))
    hits = get_single_exon_gene_overlaps('+', _isoform(1500, 2500), annots)
    assert hits == {'A': [1.0, 0.5]}


def test_each_overlapping_gene_is_reported():
    annots = _annots((1000, 3000, 'A'), (1100, 2900, 'B'))
    assert sorted(get_single_exon_gene_overlaps('+', _isoform(1200, 2800), annots)) == ['A', 'B']
