"""flair combine of fusion isoforms with an annotation: their genes and proteins"""
import pysam
from flair.pycbio.hgdata.bed import BedBlock
from flair.flair_bed import FlairBed
from flair.flair_combine import fusion_aaseq, set_genes_from_annotation
from test_fusion_piece_gene import _annots


def _genome(tmp_path):
    # chr1 + holds ATG AAA, chr2 holds TTT TAA on its minus strand (TTAAAA)
    fa = tmp_path / 'genome.fa'
    fa.write_text('>chr1\nCCATGAAAGGGG\n>chr2\nCTTAAAACCCC\n')
    pysam.faidx(str(fa))
    return pysam.FastaFile(str(fa))


def _piece(chrom, strand, blocks, thick, pos, productivity='PRO', genes=(), samples=('s',)):
    blocks = [BedBlock(start, end) for start, end in blocks]
    return FlairBed(chrom, blocks[0].start, blocks[-1].end, 'FLT1', score=3, strand=strand, blocks=blocks,
                    thickStart=thick[0], thickEnd=thick[1], gene_id='FLG1', read_support=3, frac_support=1.0,
                    productivity=productivity, transcript_class='fusion', pos_in_fusion=pos, samples=samples,
                    ref_gene_mappings=genes)


def test_fusion_protein_is_its_pieces_coding_parts_joined_in_fusion_order(tmp_path):
    # the 3' piece first, as a sample's bed may list them: still 5' piece's ATG AAA,
    # then the 3' piece's TTT TAA, its stop left off
    beds = [_piece('chr2', '-', [(1, 9)], (1, 7), 1), _piece('chr1', '+', [(0, 10)], (2, 8), 0)]
    assert fusion_aaseq(beds, _genome(tmp_path)) == 'MKF'


def test_fusion_protein_of_a_spliced_piece_skips_its_intron(tmp_path):
    # the 5' piece's second codon split by an intron, ATG A|(intron)|GG
    beds = [_piece('chr1', '+', [(0, 6), (9, 11)], (2, 11), 0), _piece('chr2', '-', [(1, 9)], (1, 7), 1)]
    assert fusion_aaseq(beds, _genome(tmp_path)) == 'MRF'


def test_no_protein_for_an_unproductive_fusion(tmp_path):
    beds = [_piece('chr1', '+', [(0, 10)], (2, 8), 0, 'PTC'), _piece('chr2', '-', [(1, 9)], (1, 7), 1, 'PTC')]
    assert fusion_aaseq(beds, _genome(tmp_path)) is None


def test_fusion_genes_from_the_annotation_whatever_the_samples_gave(tmp_path):
    # on chr1 +, A's first exon to its donor at 1200, then a 3' piece in no gene's
    # exons; one sample named the 5' piece B, whose exon also holds it
    def sample(gene5, name):
        return [_piece('chr1', '+', [(1050, 1200)], (1050, 1050), 0, 'PTC', (gene5, None), (name,)),
                _piece('chr1', '+', [(2000, 2100)], (2000, 2000), 1, 'PTC', (gene5, None), (name,))]
    group = [sample('A', 's1'), sample('B', 's2')]
    assert set_genes_from_annotation(group, _annots(tmp_path), None)
    assert {tuple(bed.ref_gene_mappings) for bed_list in group for bed in bed_list} == {('A', None)}
