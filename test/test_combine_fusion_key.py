"""
flair combine's junction chain key for fusion isoforms: each piece's junctions with
its splice sites at the fusion junctions, so isoforms are taken for one only if they
share their fusion junctions too.
"""
from flair.flair_bed import FlairBed
from flair.flair_combine import organize_isos_by_junc_chain


def _piece(name, chrom, strand, exons, pos):
    bed = FlairBed(chrom, exons[0][0], exons[-1][1], name, strand=strand, pos_in_fusion=pos,
                   transcript_class='fusion', samples=('s',), read_support=5, frac_support=1.0)
    for start, end in exons:
        bed.addBlock(start, end)
    return bed


def _chains(*fusions):
    "fusions as name -> pieces, keyed as load_isos_by_id keys them, by junction chain"
    isotoinfo = {}
    for name, pieces in fusions:
        isotoinfo[name] = [((bed.pos_in_fusion, bed.chrom, bed.strand, bed.getGaps()), bed) for bed in pieces]
    chains = {}
    organize_isos_by_junc_chain(isotoinfo, chains)
    return sorted(sorted(beds[0].name for beds in groups) for groups in chains.values())


def _five_prime():
    return [_piece(None, 'chr1', '+', [(100, 200), (300, 400)], 0)]


def test_fusions_spliced_to_different_acceptors_are_different_isoforms():
    # one 5' piece, its single-exon 3' partners on chr2 - starting at 5000 and at 5600
    a = _five_prime() + [_piece(None, 'chr2', '-', [(4000, 5000)], 1)]
    b = _five_prime() + [_piece(None, 'chr2', '-', [(4500, 5600)], 1)]
    for bed in a:
        bed.name = 'a'
    for bed in b:
        bed.name = 'b'
    assert _chains(('a', a), ('b', b)) == [['a'], ['b']]


def test_fusions_sharing_their_fusion_junction_are_one_isoform():
    # the same acceptor, chr2:5000 on -, with different 3' ends
    a = _five_prime() + [_piece(None, 'chr2', '-', [(4000, 5000)], 1)]
    b = _five_prime() + [_piece(None, 'chr2', '-', [(4300, 5000)], 1)]
    for bed in a:
        bed.name = 'a'
    for bed in b:
        bed.name = 'b'
    assert _chains(('a', a), ('b', b)) == [['a', 'b']]
