"""
flair fusion's preliminary fusions: reads grouped by the genes of their pieces, with
groups whose loci are in the same places merged, so reads of one fusion given
different genes aren't split into several fusions.
"""
from flair.identify_prelim_fusions import consensus_gene, reads_by_genes, _use_transcript_pieces


def _piece(q5, q3, start, end, gene, strand='+', chrom='chr1'):
    "an id_chimeras piece: read and genome positions in the read's direction"
    pos5, pos3 = (start, end) if strand == '+' else (end, start)
    return [(q5, pos5), (q3, pos3), gene, strand, chrom]


def _atp9a(end5, start3, end3, gene3):
    return [_piece(0, 1900, 51639339, end5, 'ATP9A', '-', 'chr20'), _piece(1900, 2200, start3, end3, gene3, '+', 'chr20')]


# ATP9A--ENSG00000234967 in MCF7: reads with one 5' piece, their 3' pieces starting at
# one acceptor, given the gene when they run on into its exons, else a non-genic locus
READS = {'r1': _atp9a(51768693, 47945603, 47945916, ('chr20', 47945603, 47945916)),
         'r2': _atp9a(51768630, 47945603, 47951940, 'G2'),
         'r3': _atp9a(51768599, 47945603, 47949919, ('chr20', 47945603, 47949919)),
         'r4': _atp9a(51768446, 47945603, 47953399, 'G2'),
         'r5': _atp9a(51768437, 47945603, 47951853, 'G2')}
EXONS = {'ATP9A': [(51639339, 51768700)], 'G2': [(47950796, 47953399)], 'SZRD1': [(16350000, 16416000)]}


def test_reads_of_one_fusion_given_different_genes_are_one_fusion():
    assert {genes: sorted(reads) for genes, reads in reads_by_genes(READS, EXONS, 10).items()} == \
        {('ATP9A', 'G2'): ['r1', 'r2', 'r3', 'r4', 'r5']}


def test_fusions_with_a_locus_elsewhere_are_kept_apart():
    # a 3' piece 3.5 kb upstream of the others, further than twice LOCUS_PADDING: another fusion
    reads = dict(READS, r6=_atp9a(51768600, 47941850, 47942094, ('chr20', 47941850, 47942094)),
                 r7=_atp9a(51768610, 47941850, 47942090, ('chr20', 47941850, 47942090)))
    fusions = reads_by_genes(reads, EXONS, 10)
    assert sorted(sorted(reads) for reads in fusions.values()) == [['r1', 'r2', 'r3', 'r4', 'r5'], ['r6', 'r7']]


def test_fusions_with_loci_whose_contigs_would_overlap_are_one_fusion():
    # a 3' piece 1.5 kb upstream: the padded loci overlap, so their contigs would share sequence
    reads = dict(READS, r6=_atp9a(51768600, 47943850, 47944094, ('chr20', 47943850, 47944094)),
                 r7=_atp9a(51768610, 47943850, 47944090, ('chr20', 47943850, 47944090)))
    assert [sorted(reads) for reads in reads_by_genes(reads, EXONS, 10).values()] == [['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7']]


def test_pieces_of_a_read_near_each_other_stay_two_loci():
    # SZRD1 in MCF7: the 3' piece upstream on the same strand, other reads' pieces
    # running between them, isn't made one locus with the 5' piece
    back = [_piece(0, 1000, 16367230, 16395605, 'SZRD1'), _piece(1000, 1500, 16341061, 16341608, ('chr1', 16341061, 16341608))]
    spanning = [_piece(0, 1000, 16341000, 16370000, 'SZRD1'), _piece(1000, 1500, 9000, 9500, 'OTHER', '+', 'chr2')]
    fusions = reads_by_genes({'b1': back, 'b2': back, 's1': spanning}, EXONS, 10)
    assert sorted(fusions[('SZRD1', 'chr1:16341061-16341608')]) == ['b1', 'b2']


def test_consensus_gene_prefers_the_majority_then_an_annotated_gene():
    assert consensus_gene(['A', 'B', 'B']) == 'B'
    assert consensus_gene(['chr1:5-9', 'A']) == 'A'
    assert consensus_gene(['B', 'A']) == 'A'


def test_transcript_pieces_used_for_a_middle_exon_the_genome_alignment_missed():
    genome = [_piece(0, 500, 1000, 1500, 'A'), _piece(620, 1200, 9000, 9580, 'C')]
    transcript = [_piece(0, 500, 1000, 1500, 'A'), _piece(500, 620, 5000, 5120, 'B'), _piece(620, 1200, 9000, 9580, 'C')]
    assert _use_transcript_pieces(transcript, genome)


def test_transcript_pieces_used_for_a_middle_exon_the_genome_alignment_ran_through():
    # ERGIC3--DDX5--KRAS in test-fusion-3gene: the genome's ERGIC3 piece runs on
    # through DDX5's sequence and into KRAS's
    genome = [_piece(0, 1455, 35542020, 35573759, 'ERGIC3'), _piece(919, 5993, 25245395, 25205245, 'KRAS', '-')]
    transcript = [_piece(0, 453, 35542020, 35547418, 'ERGIC3'), _piece(445, 913, 64503572, 64502931, 'DDX5', '-', 'chr17'),
                  _piece(919, 5993, 25245395, 25205245, 'KRAS', '-')]
    assert _use_transcript_pieces(transcript, genome)


def test_transcript_pieces_not_used_for_a_gene_aligned_where_the_genome_alignment_has_another():
    # PKP2--ENSG00000309487 in MCF7: a short KANK1 alignment inside the 3' piece
    genome = [_piece(0, 260, 1000, 1260, 'PKP2'), _piece(256, 2119, 9000, 10863, 'G2')]
    transcript = [_piece(0, 259, 1000, 1259, 'PKP2'), _piece(257, 1068, 9000, 9811, 'G2'),
                  _piece(1207, 1556, 500, 849, 'KANK1', chrom='chr9')]
    assert not _use_transcript_pieces(transcript, genome)
