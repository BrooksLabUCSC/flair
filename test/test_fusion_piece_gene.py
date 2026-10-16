from flair.pycbio.hgdata.bed import BedBlock, BedReader
from flair.gtf_io import gtf_data_parser, GtfAttrsSet, TRANSCRIPT_EXON_FEATURES
from flair.annotation_data import AnnotsByChrom
from flair.flair_bed import FlairBed
from flair.convert_synthetic_to_genome_bed import breakpoint_sites, piece_gene, convert_synthetic_isos

# on chr1 +: gene A, exons 1000-1200 and 1500-1700, and gene B, exons 900-1300 and
# 1600-1800, whose first exon holds all of A's (0-based, end exclusive)
GTF = [('A', 'TA', [(1000, 1200), (1500, 1700)]),
       ('B', 'TB', [(900, 1300), (1600, 1800)])]


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
    assert piece_gene('chr1', '+', _blocks((1050, 1200)), 0, 2, 'B', _annots(tmp_path)) == 'A'


def test_last_exon_piece_gets_the_gene_of_its_acceptor(tmp_path):
    assert piece_gene('chr1', '+', _blocks((1500, 1650)), 1, 2, 'B', _annots(tmp_path)) == 'A'


def test_piece_without_an_annotated_splice_site_keeps_the_preliminary_gene(tmp_path):
    # a breakpoint inside the exon, at 1180
    assert piece_gene('chr1', '+', _blocks((1050, 1180)), 0, 2, 'B', _annots(tmp_path)) == 'B'


def test_without_an_annotation_the_preliminary_gene_is_kept():
    assert piece_gene('chr1', '+', _blocks((1050, 1200)), 0, 2, 'B', None) == 'B'


def _convert(tmp_path, exons):
    """convert_synthetic_isos of an isoform with exons on the synthetic contig of
    preliminary genes B (chr1:900-1300, +) and C (chr1:1400-1800, +), loci 0-400 and
    400-800 of it"""
    breakpoints = tmp_path / 'breakpoints.bed'
    breakpoints.write_text('B--C\t400\t400\tbreakpoint-1--B..chr1..900..1300--C..chr1..1400..1800\n')
    synthetic = tmp_path / 'synthetic.bed'
    with open(synthetic, 'w') as fh:
        FlairBed('B--C', exons[0][0], exons[-1][1], 'FLT1', score=5, strand='+', blocks=_blocks(*exons),
                 gene_id='FLG1', ref_gene_mappings=('B--C',), read_support=5, frac_support=1.0, samples=('s',)).write(fh)
    out = tmp_path / 'fusions.bed'
    convert_synthetic_isos(str(synthetic), str(breakpoints), str(out), 0, _annots(tmp_path))
    return [(bed.pos_in_fusion, bed.ref_gene_mappings) for bed in BedReader(str(out), bedClass=FlairBed)]


def test_isoform_whose_pieces_are_of_one_gene_is_left_out(tmp_path):
    # A's first exon to 1200, its donor, then its last exon from 1500, its acceptor:
    # both pieces are A's
    assert _convert(tmp_path, [(150, 300), (500, 650)]) == []


def test_isoform_of_two_genes_is_kept(tmp_path):
    # the 3' piece starts at 1520, no annotated splice site, and keeps C
    assert sorted(_convert(tmp_path, [(150, 300), (520, 650)])) == [(0, ('A', 'C')), (1, ('A', 'C'))]
