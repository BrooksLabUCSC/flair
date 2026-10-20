import pysam
from flair.flair_bed import FlairBed
from flair.pycbio.hgdata.bed import BedBlock
from flair.flair_quantify import (direct_assignment_chains, direct_isoform, add_direct_assignments, GeneData,
                                  is_realigned_to, single_exon_index, single_exon_isoform)

GENOME_SEQ = ''.join('ACGT'[(i * 7 + i // 3) % 4] for i in range(2000))
_HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 2000}]})


class _Genome:
    def fetch(self, chrom, start, end):
        return GENOME_SEQ[start:end]


def _iso(name, *exons, pos_in_fusion=None, transcript_class='basic'):
    return FlairBed('chr1', exons[0][0], exons[-1][1], name, strand='+', blocks=[BedBlock(s, e) for s, e in exons],
                    pos_in_fusion=pos_in_fusion, transcript_class=transcript_class)


def _read(start, cigar, mutate=(), sa=False):
    "an alignment matching the genome, with query bases changed at reference positions in mutate"
    seg = pysam.AlignedSegment(_HEADER)
    seg.query_name = 'r'
    seg.reference_id = 0
    seg.reference_start = start
    seg.cigartuples = cigar
    seq, rpos = [], start
    for op, length in cigar:
        if op == pysam.CMATCH:
            seq.extend('A' if rpos + k in mutate and GENOME_SEQ[rpos + k] != 'A' else
                       'C' if rpos + k in mutate else GENOME_SEQ[rpos + k] for k in range(length))
            rpos += length
        elif op in (pysam.CDEL, pysam.CREF_SKIP):
            rpos += length
        elif op in (pysam.CINS, pysam.CSOFT_CLIP):
            seq.append('G' * length)
    seg.query_sequence = ''.join(seq)
    if sa:
        seg.set_tag('SA', 'chr1,1500,+,50M,60,0;')
    return seg


# ASCL1's case: the main isoform, exons 500-1000 and 1200-1400, and a
# low-support isoform with an extra first exon, 100-150, whose chain contains it
MAIN = _iso('main', (500, 1000), (1200, 1400))
UPSTREAM = _iso('upstream', (100, 150), (600, 1000), (1200, 1400))
# a read starting inside the upstream isoform's internal exon, with the main
# isoform's one junction
TSS_READ = [(pysam.CMATCH, 300), (pysam.CREF_SKIP, 200), (pysam.CMATCH, 200)]


def _direct(read, isoforms, trust_ends=False):
    bed = direct_isoform(read, direct_assignment_chains(isoforms), _Genome(), trust_ends)
    return None if bed is None else bed.name


def test_exact_match_is_assigned_though_another_chain_contains_it():
    # realigned, its start short of the main isoform's unique sequence past 600
    # left it to neither isoform
    assert _direct(_read(700, TSS_READ), [MAIN, UPSTREAM]) == 'main'


def test_contained_chain_is_not_assigned_to_the_longer_isoform():
    assert _direct(_read(120, [(pysam.CMATCH, 30), (pysam.CREF_SKIP, 450), (pysam.CMATCH, 400),
                               (pysam.CREF_SKIP, 200), (pysam.CMATCH, 100)]), [MAIN, UPSTREAM]) == 'upstream'
    # a chain no isoform has, an intron off by a base, is realigned
    assert _direct(_read(700, [(pysam.CMATCH, 301), (pysam.CREF_SKIP, 199), (pysam.CMATCH, 200)]), [MAIN, UPSTREAM]) is None


def test_unspliced_read_is_realigned():
    assert _direct(_read(600, [(pysam.CMATCH, 300)]), [MAIN, UPSTREAM]) is None


def test_splice_site_not_cleanly_aligned_is_realigned():
    assert _direct(_read(700, TSS_READ, mutate={997, 998}), [MAIN]) is None
    assert _direct(_read(700, TSS_READ, mutate={998}), [MAIN]) == 'main'
    assert _direct(_read(700, [(pysam.CMATCH, 297), (pysam.CDEL, 2), (pysam.CMATCH, 1), (pysam.CREF_SKIP, 200),
                               (pysam.CMATCH, 200)]), [MAIN]) is None


def test_chain_of_two_isoforms_is_realigned():
    other_ends = _iso('other_ends', (300, 1000), (1200, 1600))
    assert direct_assignment_chains([MAIN, other_ends]) == {}
    assert _direct(_read(700, TSS_READ), [MAIN, other_ends]) is None


def test_fusion_pieces_and_chimeric_reads_are_realigned():
    piece = _iso('fusion', (500, 1000), (1200, 1400), pos_in_fusion=0, transcript_class='fusion')
    assert direct_assignment_chains([piece]) == {}
    assert _direct(_read(700, TSS_READ, sa=True), [MAIN]) is None


def test_trust_ends_needs_ends_near_the_isoforms():
    assert _direct(_read(700, TSS_READ), [MAIN], trust_ends=True) is None
    assert _direct(_read(520, [(pysam.CMATCH, 480), (pysam.CREF_SKIP, 200), (pysam.CMATCH, 180)]), [MAIN],
                   trust_ends=True) == 'main'


def test_direct_assignments_added_to_counts_and_map(tmp_path):
    counts, read_map = tmp_path / 'counts.txt', tmp_path / 'map.txt'
    counts.write_text('main\t2\nupstream\t1\n')
    read_map.write_text('main\tr1,r2\nupstream\tr3\n')
    add_direct_assignments({'main': ['r4'], 'third': ['r5', 'r6']}, counts, read_map)
    assert counts.read_text() == 'main\t3\nupstream\t1\nthird\t2\n'
    assert read_map.read_text() == 'main\tr1,r2,r4\nupstream\tr3\nthird\tr5,r6\n'


def test_direct_assignments_without_realigned_reads_or_map(tmp_path):
    counts = tmp_path / 'counts.txt'
    counts.write_text('')
    add_direct_assignments({'main': ['r1']}, counts, None)
    assert counts.read_text() == 'main\t1\n'


def _gene_data(*genes):
    "gene id -> GeneData, from (gene id, isoforms)"
    data = {}
    for gene_id, isos in genes:
        data[gene_id] = GeneData(gene_id, 'chr1', '+')
        for iso in isos:
            iso.gene_id = gene_id
            data[gene_id].add_isoform_bed(iso)
    return data


def _assigned_se(read, *genes):
    iso = single_exon_isoform(read, single_exon_index(_gene_data(*genes)))
    return None if iso is None else (iso[2], iso[3])


def test_unspliced_read_goes_to_single_exon_isoform_with_nearest_ends():
    near, far = _iso('near', (600, 900)), _iso('far', (500, 1000))
    assert _assigned_se(_read(610, [(pysam.CMATCH, 280)]), ('G1', [near, far])) == ('near', 'G1')
    assert _assigned_se(_read(510, [(pysam.CMATCH, 480)]), ('G1', [near, far])) == ('far', 'G1')


def test_unspliced_read_must_cover_half_of_single_exon_isoform():
    assert _assigned_se(_read(600, [(pysam.CMATCH, 240)]), ('G1', [_iso('se', (500, 1000))])) is None
    assert _assigned_se(_read(600, [(pysam.CMATCH, 260)]), ('G1', [_iso('se', (500, 1000))])) == ('se', 'G1')


def test_unspliced_read_counted_for_one_gene_of_overlapping_genes():
    # the same choice from either gene's reads: G2's isoform, of nearer ends
    genes = (('G1', [_iso('g1', (500, 1000)), MAIN]), ('G2', [_iso('g2', (600, 950))]))
    assert _assigned_se(_read(610, [(pysam.CMATCH, 330)]), *genes) == ('g2', 'G2')


def test_spliced_and_fusion_isoforms_are_realigned_to():
    piece = _iso('piece', (500, 1000), pos_in_fusion=0, transcript_class='fusion')
    assert is_realigned_to(MAIN) and is_realigned_to(piece)
    assert not is_realigned_to(_iso('se', (500, 1000)))
    assert _assigned_se(_read(510, [(pysam.CMATCH, 480)]), ('G1', [MAIN, piece])) is None
