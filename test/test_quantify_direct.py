import pytest
import pysam
from flair import FlairInputDataError
from flair.isoform_data import get_reverse_complement
from flair.flair_bed import FlairBed
from flair.pycbio.hgdata.bed import BedBlock
from flair.flair_quantify import (direct_assignment_chains, direct_isoform, read_count_sam_results, make_batches, GeneData,
                                  batch_reads, is_realigned_to, single_exon_index, single_exon_isoform,
                                  fused_transcript, load_isoform_data)

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


def test_chimeric_reads_are_realigned():
    # a read with a supplementary alignment may be a fusion's
    assert _direct(_read(700, TSS_READ, sa=True), [MAIN]) is None


def test_trust_ends_needs_ends_near_the_isoforms():
    assert _direct(_read(700, TSS_READ), [MAIN], trust_ends=True) is None
    assert _direct(_read(520, [(pysam.CMATCH, 480), (pysam.CREF_SKIP, 200), (pysam.CMATCH, 180)]), [MAIN],
                   trust_ends=True) == 'main'


def test_realigned_reads_added_to_direct_assignments(tmp_path):
    counts, read_map = tmp_path / 'counts.txt', tmp_path / 'map.txt'
    counts.write_text('main\t2\nupstream\t1\n')
    read_map.write_text('main\tr1,r2\nupstream\tr3\n')
    assigned = {'main': ['r4'], 'third': ['r5', 'r6']}
    read_count_sam_results(counts, read_map, assigned)
    assert assigned == {'main': ['r4', 'r1', 'r2'], 'third': ['r5', 'r6'], 'upstream': ['r3']}


def test_realigned_counts_added_without_a_map(tmp_path):
    counts = tmp_path / 'counts.txt'
    counts.write_text('main\t2\nupstream\t1\n')
    assigned = {'main': 1}
    read_count_sam_results(counts, None, assigned)
    assert assigned == {'main': 3, 'upstream': 1}


def _gene(gene_id, start, end, n_isoforms=1, chrom='chr1'):
    gene = GeneData(gene_id, chrom, '+')
    for i in range(n_isoforms):
        gene.add_isoform_bed(FlairBed(chrom, start, end, f'{gene_id}.{i}', strand='+', blocks=[BedBlock(start, end)]))
    return gene


def _batches(genes, num_batches):
    return [(b['chrom'], b['start'], b['end'], b['loci'], b['prev_end'], b['gene_ids'])
            for b in make_batches({g.gene_id: g for g in genes}, num_batches)]


def test_overlapping_genes_are_batched_together():
    genes = [_gene('A', 100, 500), _gene('B', 400, 900), _gene('C', 2000, 2500)]
    assert _batches(genes, 3) == [('chr1', 100, 900, [(100, 900)], None, ['A', 'B']),
                                  ('chr1', 2000, 2500, [(2000, 2500)], 900, ['C'])]


def test_batches_by_chromosome():
    genes = [_gene('A', 100, 500), _gene('B', 1000, 1500), _gene('C', 100, 300, chrom='chr2')]
    assert _batches(genes, 3) == [('chr1', 100, 500, [(100, 500)], None, ['A']),
                                  ('chr1', 1000, 1500, [(1000, 1500)], 500, ['B']),
                                  ('chr2', 100, 300, [(100, 300)], None, ['C'])]


def test_loci_grouped_up_to_a_share_of_the_isoforms():
    genes = [_gene('A', 100, 200, 2), _gene('B', 300, 400, 1), _gene('C', 500, 600, 1), _gene('D', 700, 800, 4)]
    # 8 isoforms over 2 batches: 4 each, a batch closing once it has that many
    assert _batches(genes, 2) == [('chr1', 100, 600, [(100, 200), (300, 400), (500, 600)], None, ['A', 'B', 'C']),
                                  ('chr1', 700, 800, [(700, 800)], 600, ['D'])]


def test_each_read_overlapping_a_locus_taken_by_one_batch(tmp_path):
    genes = {g.gene_id: g for g in (_gene('A', 100, 500), _gene('B', 600, 900), _gene('C', 1200, 1500))}
    reads = [('over_a_and_b', 250, [(pysam.CMATCH, 450)]),          # A's single-exon isoform, taken by A's batch only
             ('spliced_in_b', 610, [(pysam.CMATCH, 40), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 40)]),
             ('in_the_gap', 950, [(pysam.CMATCH, 50), (pysam.CREF_SKIP, 100), (pysam.CMATCH, 50)])]
    bam_path = str(tmp_path / 'reads.bam')
    with pysam.AlignmentFile(bam_path, 'wb', header=_HEADER) as bam:
        for name, start, cigar in reads:
            seg = _read(start, cigar)
            seg.query_name = name
            bam.write(seg)
    pysam.index(bam_path)
    se_index = single_exon_index(genes)
    taken = []
    with pysam.AlignmentFile(bam_path) as bam:
        for batch in make_batches(genes, 3):
            direct, num_realign = batch_reads(bam, batch, {}, _Genome(), False, se_index,
                                              str(tmp_path / 'r.fa'), str(tmp_path / 'clip.txt'))
            realigned = [line[1:].strip() for line in open(tmp_path / 'r.fa') if line.startswith('>')]
            taken.append((batch['gene_ids'], direct, realigned))
            assert num_realign == len(realigned)
    assert taken == [(['A'], {'A.0': ['over_a_and_b']}, []), (['B'], {}, ['spliced_in_b']), (['C'], {}, [])]


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


def test_spliced_isoforms_are_realigned_to():
    assert is_realigned_to(MAIN)
    assert not is_realigned_to(_iso('se', (500, 1000)))
    # a spliced isoform takes no unspliced read without the realignment
    assert _assigned_se(_read(510, [(pysam.CMATCH, 480)]), ('G1', [MAIN])) is None


def _piece(name, chrom, strand, exons, pos, gene_id='FUS'):
    return FlairBed(chrom, exons[0][0], exons[-1][1], name, strand=strand, blocks=[BedBlock(s, e) for s, e in exons],
                    gene_id=gene_id, pos_in_fusion=pos, transcript_class='fusion')


def test_fused_transcript_joins_pieces_in_order_each_in_its_orientation():
    pieces = [_piece('F', 'chr1', '+', [(100, 110), (120, 130)], 0), _piece('F', 'chr1', '-', [(300, 310), (320, 330)], 1)]
    seq, sizes = fused_transcript(pieces, _Genome())
    minus = get_reverse_complement(GENOME_SEQ[300:310] + GENOME_SEQ[320:330])
    assert seq == GENOME_SEQ[100:110] + GENOME_SEQ[120:130] + minus
    assert sizes == [10, 10, 10, 10]


def _write_bed(path, beds):
    with open(path, 'w') as fh:
        for bed in beds:
            bed.write(fh)
    return str(path)


def test_fusion_pieces_kept_apart_from_genes(tmp_path):
    main = _iso('main', (500, 1000), (1200, 1400))
    main.gene_id = 'G1'
    # pieces on two strands, sharing the fusion's gene id, given out of order
    pieces = [_piece('F', 'chr1', '-', [(1500, 1600)], 1), _piece('F', 'chr1', '+', [(100, 200), (300, 400)], 0)]
    gene_data, fusions = load_isoform_data(_write_bed(tmp_path / 'iso.bed', [main] + pieces))
    assert list(gene_data) == ['G1'] and [b.name for b in gene_data['G1'].isoform_beds] == ['main']
    assert [(b.pos_in_fusion, b.strand) for b in fusions['F']] == [(0, '+'), (1, '-')]


def test_fusion_needs_its_pieces_and_its_own_name(tmp_path):
    with pytest.raises(FlairInputDataError, match='pos_in_fusion'):
        load_isoform_data(_write_bed(tmp_path / 'a.bed', [_piece('F', 'chr1', '+', [(100, 200)], 1)]))
    main = _iso('F', (500, 1000), (1200, 1400))
    main.gene_id = 'G1'
    pieces = [_piece('F', 'chr2', '+', [(100, 200)], 0), _piece('F', 'chr2', '-', [(300, 400)], 1)]
    with pytest.raises(FlairInputDataError, match='both a fusion and another isoform'):
        load_isoform_data(_write_bed(tmp_path / 'b.bed', [main] + pieces))
