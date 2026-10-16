from flair.pycbio.hgdata.bed import BedBlock, BedReader
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import AnnotsByChrom
from flair.flair_bed import FlairBed
from flair.convert_synthetic_to_genome_bed import breakpoint_sites, piece_gene, convert_synthetic_isos

# on chr1 +: gene A, exons 1000-1200 and 1500-1700, and gene B, exons 900-1300 and
# 1600-1800, whose first exon holds all of A's, and genes C and D downstream, not
# overlapping (0-based, end exclusive)
GTF = [('A', 'TA', [(1000, 1200), (1500, 1700)]),
       ('B', 'TB', [(900, 1300), (1600, 1800)]),
       ('C', 'TC', [(3000, 3200), (3500, 3700)]),
       ('D', 'TD', [(4000, 4200), (4500, 4600), (4800, 4900)])]


def _annots(tmp_path):
    gtf = tmp_path / 'annotation.gtf'
    with open(gtf, 'w') as fh:
        for gene, transcript, exons in GTF:
            attrs = f'gene_id "{gene}"; transcript_id "{transcript}"; gene_name "{gene}";'
            fh.write(f'chr1\ttest\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t+\t.\t{attrs}\n')
            for start, end in exons:
                fh.write(f'chr1\ttest\texon\t{start + 1}\t{end}\t.\t+\t.\t{attrs}\n')
    return AnnotsByChrom(gtf_data_parser(str(gtf), attrs=GtfAttrsSet.FLAIR, include_features=TRANSCRIPT_EXON_FEATURES))


def _blocks(*exons):
    return [BedBlock(start, end) for start, end in exons]


def test_breakpoint_sites_are_the_edges_at_the_fusion_junctions():
    blocks = _blocks((100, 200), (300, 400))
    assert breakpoint_sites(blocks, '+', 0, 2) == [400]
    assert breakpoint_sites(blocks, '+', 1, 2) == [100]
    assert breakpoint_sites(blocks, '+', 1, 3) == [400, 100]
    assert breakpoint_sites(blocks, '-', 0, 2) == [100]
    assert breakpoint_sites(blocks, '-', 1, 2) == [400]


def test_first_exon_piece_gets_the_gene_of_its_donor(tmp_path):
    # a 5' piece of only A's first exon, which also lies in B's exon, as the
    # preliminary fusion took it from where its reads ended
    assert piece_gene('chr1', '+', _blocks((1050, 1200)), 0, 2, _annots(tmp_path)) == 'A'


def test_last_exon_piece_gets_the_gene_of_its_acceptor(tmp_path):
    assert piece_gene('chr1', '+', _blocks((1500, 1650)), 1, 2, _annots(tmp_path)) == 'A'


def test_piece_without_an_annotated_splice_site_gets_the_gene_of_its_exons(tmp_path):
    # a breakpoint inside B's exon at 1290, past A's: B
    assert piece_gene('chr1', '+', _blocks((1210, 1290)), 0, 2, _annots(tmp_path)) == 'B'


def test_piece_in_two_genes_exons_gets_the_first_by_id_on_a_tie(tmp_path):
    # in both A's exon 1000-1200 and B's 900-1300, each with one transcript
    assert piece_gene('chr1', '+', _blocks((1050, 1180)), 0, 2, _annots(tmp_path)) == 'A'


def test_piece_in_no_genes_exons_is_a_non_genic_locus(tmp_path):
    # left for convert_synthetic_isos to name
    annots = _annots(tmp_path)
    assert piece_gene('chr1', '+', _blocks((2000, 2100)), 1, 2, annots) is None
    # on a chromosome without genes
    assert piece_gene('chr2', '+', _blocks((2000, 2100)), 1, 2, annots) is None


def test_piece_reading_through_two_genes_gets_the_one_at_its_fusion_junction(tmp_path):
    # C's junction, then D's, the fusion junction at D's donor at 4600: D, the gene
    # fused, not C, first by id
    blocks = _blocks((3100, 3200), (3500, 3700), (4000, 4200), (4500, 4600))
    assert piece_gene('chr1', '+', blocks, 0, 2, _annots(tmp_path)) == 'D'


def _convert(tmp_path, *isoforms):
    """convert_synthetic_isos of isoforms, each a list of exons, on the synthetic
    contig of preliminary genes B (chr1:900-1300, +) and C (chr1:1400-2400, +), loci
    0-400 and 400-1400 of it: the bed's (isoform, position, ref_gene_mappings), and
    the isoforms' genes"""
    breakpoints = tmp_path / 'breakpoints.bed'
    breakpoints.write_text('B--C\t400\t400\tbreakpoint-1--B..chr1..900..1300--C..chr1..1400..2400\n')
    synthetic = tmp_path / 'synthetic.bed'
    with open(synthetic, 'w') as fh:
        for n, exons in enumerate(isoforms, start=1):
            FlairBed('B--C', exons[0][0], exons[-1][1], f'FLT{n}', score=5, strand='+', blocks=_blocks(*exons),
                     gene_id='FLG1', ref_gene_mappings=('B--C',), read_support=5, frac_support=1.0, samples=('s',)).write(fh)
    out = tmp_path / 'fusions.bed'
    genes = convert_synthetic_isos(str(synthetic), str(breakpoints), str(out), 0, _annots(tmp_path))
    return sorted((bed.name, bed.pos_in_fusion, bed.ref_gene_mappings) for bed in BedReader(str(out), bedClass=FlairBed)), genes


def test_isoform_whose_pieces_are_of_one_gene_is_left_out(tmp_path):
    # A's first exon to 1200, its donor, then its last exon from 1500, its acceptor:
    # both pieces are A's
    assert _convert(tmp_path, [(150, 300), (500, 650)]) == ([], {})


def test_isoform_of_a_gene_and_a_non_genic_locus_is_kept(tmp_path):
    # the 3' piece, chr1:2300-2380, is in no gene: a non-genic locus, no gene id
    beds, genes = _convert(tmp_path, [(150, 300), (1300, 1380)])
    assert beds == [('FLT1', 0, ('A', None)), ('FLT1', 1, ('A', None))]
    assert genes == {'FLT1': ('A', 'chr1:2300-2380')}


def test_isoforms_with_different_ends_at_a_non_genic_locus_share_its_name(tmp_path):
    # not C, the preliminary gene, nor each piece's own extent
    _, genes = _convert(tmp_path, [(150, 300), (1300, 1380)], [(150, 300), (1300, 1350)], [(150, 300), (1250, 1330)])
    assert set(genes.values()) == {('A', 'chr1:2250-2380')}
