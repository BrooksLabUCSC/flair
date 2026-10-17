"""
flair fusion gives each chimeric read piece the shortest gene sharing exons with the
one it was assigned, so a readthrough annotation gives way to the gene it contains.
The layout follows 9p21 in MCF7: a readthrough lncRNA starting with MTAP's exons and
running into CDKN2B-AS1, so reads of the MTAP--CDKN2B-AS1 fusion went to either.
"""
from flair.identify_prelim_fusions import exonic_overlap, exon_sharing_genes, prefer_shorter_gene, is_foldback, getCorrectGene

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


def _piece(qstart, qend, start, end, strand, chrom='chr16'):
    "an id_chimeras piece: read and genome positions in the read's direction"
    if strand == '+':
        return [(qstart, start), (qend, end), 'gene', strand, chrom]
    return [(qstart, end), (qend, start), 'gene', strand, chrom]


def test_foldback_reads_are_detected():
    "SLC7A5 in MCF7: a short piece folded back over the read's main piece"
    overlapping = [_piece(0, 91, 87831293, 87831384, '+'), _piece(85, 1473, 87830016, 87831404, '-')]
    nearby = [_piece(0, 268, 96330704, 96330974, '-'), _piece(244, 1069, 96331069, 96331905, '+')]
    assert is_foldback(overlapping)
    assert is_foldback(nearby)


def test_fusions_are_not_foldbacks():
    same_strand_neighbours = [_piece(0, 500, 1000, 1500, '+'), _piece(500, 1000, 1700, 2200, '+')]
    distant_inversion = [_piece(0, 500, 1000, 1500, '+'), _piece(500, 1000, 50000, 50500, '-')]
    other_chrom = [_piece(0, 500, 1000, 1500, '+'), _piece(500, 1000, 1000, 1500, '-', chrom='chr17')]
    assert not is_foldback(same_strand_neighbours)
    assert not is_foldback(distant_inversion)
    assert not is_foldback(other_chrom)


def _gene_for(blocks, strand='+', juncs=()):
    "getCorrectGene for a piece aligned to blocks, with GENES and EXONS as the annotation"
    chrom_to_gene_pos = {'chr9': sorted((s, e, st, g) for g, (_, s, e, st) in GENES.items())}
    return getCorrectGene(chrom_to_gene_pos, EXONS, {'chr9': {}}, 'chr9', (list(juncs), blocks), strand)


def test_piece_on_exons_gets_the_gene():
    assert _gene_for([(5020, 5100), (9000, 9100)], juncs=[(5100, 9000)]) in ('MTAP', 'READTHROUGH')


def test_intronic_piece_is_not_genic():
    "RPS29 in MCF7: chimeric pieces inside its introns made it a hub with 26 partners"
    assert _gene_for([(2000, 2500)]) == ('chr9', 2000, 2500)


def test_piece_spliced_across_an_exon_without_aligning_to_it_is_not_genic():
    "the piece's span holds MTAP's exon 5000-5100, but none of its aligned bases do"
    assert _gene_for([(3000, 3200), (6000, 6200)], juncs=[(3200, 6000)]) == ('chr9', 3000, 6200)


def test_shorter_gene_not_taken_for_a_piece_off_the_genes_exons():
    "a piece covering none of the assigned gene's exons isn't moved to a nested shorter gene"
    assert prefer_shorter_gene('READTHROUGH', [(20000, 20500)], _sharing(), GENES, EXONS) == 'READTHROUGH'
