"""
flair fusion gives each chimeric read piece the shortest gene sharing exons with the
one it was assigned, so a readthrough annotation gives way to the gene it contains.
The layout follows 9p21 in MCF7: a readthrough lncRNA starting with MTAP's exons and
running into CDKN2B-AS1, so reads of the MTAP--CDKN2B-AS1 fusion went to either.
"""
from flair.identify_prelim_fusions import exonic_overlap, exon_sharing_genes, prefer_shorter_gene

GENES = {  # gene -> (chrom, start, end, strand), as genetoinfo without the TSS list
    'MTAP': ('chr9', 1000, 16000, '+'),
    'READTHROUGH': ('chr9', 1000, 60000, '+'),
    'CDKN2B-AS1': ('chr9', 55000, 90000, '+'),
    'CDKN2A': ('chr9', 30000, 40000, '-'),
    'ANTISENSE': ('chr9', 2000, 5000, '-'),
}
EXONS = {
    'MTAP': [(1000, 1200), (5000, 5100), (9000, 9200), (15800, 16000)],
    'READTHROUGH': [(1000, 1200), (5000, 5100), (9000, 9200), (15800, 16000), (30500, 30600), (58000, 60000)],
    'CDKN2B-AS1': [(55000, 55300), (58500, 59000), (89000, 90000)],
    'CDKN2A': [(30000, 30700), (39000, 40000)],
    'ANTISENSE': [(1000, 1200), (4000, 5000)],
}


def _sharing():
    chrom_to_gene_pos = {'chr9': sorted((s, e, strand, g) for g, (_, s, e, strand) in GENES.items())}
    return exon_sharing_genes(chrom_to_gene_pos, EXONS)


def test_exonic_overlap():
    assert exonic_overlap([(1100, 5050)], EXONS['MTAP']) == 100 + 50
    assert exonic_overlap([(1200, 5000)], EXONS['MTAP']) == 0


def test_exon_sharing_is_same_strand_only():
    sharing = _sharing()
    assert sharing['MTAP'] == {'READTHROUGH'}
    assert sharing['READTHROUGH'] == {'MTAP', 'CDKN2B-AS1'}
    assert 'ANTISENSE' not in sharing and 'CDKN2A' not in sharing


def test_readthrough_piece_goes_to_contained_gene():
    "a piece on MTAP's exons 1-4 assigned to the readthrough is MTAP's"
    piece = [(1050, 1200), (5000, 5100), (9000, 9200), (15800, 15950)]
    assert prefer_shorter_gene('READTHROUGH', piece, _sharing(), GENES, EXONS) == 'MTAP'


def test_partner_piece_is_not_pulled_into_the_readthrough():
    "the 3' piece stays CDKN2B-AS1, so the fusion doesn't become one gene"
    piece = [(58500, 59000)]
    assert prefer_shorter_gene('CDKN2B-AS1', piece, _sharing(), GENES, EXONS) == 'CDKN2B-AS1'


def test_readthrough_kept_where_only_it_has_the_exons():
    "a piece on the readthrough's own exon isn't moved to a gene that lacks it"
    piece = [(15800, 16000), (30500, 30600)]
    assert prefer_shorter_gene('READTHROUGH', piece, _sharing(), GENES, EXONS) == 'READTHROUGH'

