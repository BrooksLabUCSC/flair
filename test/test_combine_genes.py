"""
flair combine's genes: isoforms without an annotated gene are grouped into novel
genes as flair transcriptome groups them, from the combined isoforms, so the genes
don't depend on which isoforms a sample happened to put together, nor on the
order of the samples.
"""
import pytest
from flair import FlairError
from flair.flair_bed import FlairBed
from flair.flair_combine import assign_gene_keys, _combined_sort_key


def _bed(name, sample, gene_id, exons, strand='+', ref_genes=()):
    bed = FlairBed('chr1', exons[0][0], exons[-1][1], name, strand=strand, gene_id=gene_id,
                   ref_gene_mappings=ref_genes, samples=(sample,), read_support=10, frac_support=0.5)
    for start, end in exons:
        bed.addBlock(start, end)
    return bed


def _gene_keys(*beds):
    "gene keys of combined isoforms each made of one sample's isoform, in output order"
    combined = sorted([[(bed,)] for bed in beds], key=_combined_sort_key)
    keys, _ = assign_gene_keys(combined)
    return {group[0][0].name: key for group, (key, partners) in zip(combined, keys)}


def test_novel_isoforms_of_different_samples_sharing_a_splice_site_are_one_gene():
    # each sample's isoform is alone in its sample's gene; they share the 3' splice
    # site at 2000
    a = _bed('a', 's1', 'FLG1', [(1000, 1200), (2000, 2300)])
    b = _bed('b', 's2', 'FLG7', [(1500, 1700), (2000, 2300)])
    elsewhere = _bed('c', 's2', 'FLG8', [(9000, 9200), (9500, 9800)])
    keys = _gene_keys(a, b, elsewhere)
    assert keys['a'] == keys['b'] != keys['c']


def test_novel_isoforms_on_other_strands_are_other_genes():
    keys = _gene_keys(_bed('a', 's1', 'FLG1', [(1000, 1200), (2000, 2300)]),
                      _bed('b', 's2', 'FLG1', [(1000, 1200), (2000, 2300)], strand='-'))
    assert keys['a'] != keys['b']


def test_annotated_isoforms_are_in_their_genes():
    keys = _gene_keys(_bed('a', 's1', 'FLG1', [(1000, 1200), (2000, 2300)], ref_genes=('G1',)),
                      _bed('b', 's2', 'FLG3', [(1500, 1700), (2000, 2300)], ref_genes=('G1',)),
                      _bed('c', 's2', 'FLG3', [(1000, 1200), (2500, 2600)], ref_genes=('G2',)))
    assert keys['a'] == keys['b'] == ('ref', ('G1',)) and keys['c'] == ('ref', ('G2',))


def test_different_genes_from_different_samples_are_an_error():
    a = _bed('a', 's1', 'FLG1', [(1000, 1200), (2000, 2300)], ref_genes=('G1',))
    b = _bed('b', 's2', 'FLG1', [(1000, 1200), (2000, 2300)], ref_genes=('G2',))
    with pytest.raises(FlairError, match='s1: G1; s2: G2'):
        assign_gene_keys([[(a,), (b,)]])


def test_different_reference_transcripts_from_different_samples_are_an_error():
    a = _bed('a', 's1', 'FLG1', [(1000, 1200), (2000, 2300)], ref_genes=('G1',))
    b = _bed('b', 's2', 'FLG1', [(1000, 1200), (2000, 2300)], ref_genes=('G1',))
    a.ref_transcript_id, b.ref_transcript_id = 'T1', 'T1-endvar2'
    with pytest.raises(FlairError, match='s1: T1; s2: T1-endvar2'):
        assign_gene_keys([[(a,), (b,)]])
