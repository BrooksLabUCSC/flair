from flair.pycbio.hgdata.bed import BedBlock, BedReader
from flair.flair_bed import FlairBed
from flair.convert_synthetic_to_genome_bed import convert_synthetic_isos

# a synthetic fusion of GA (chr1:1000-2000, +) and GB (chr2:4000-5000, -), its loci
# at 0-1000 and 1000-2000 of the synthetic chromosome; the isoform has exons
# 100-300 and 500-900 in GA's locus and 1100-1400 and 1600-1900 in GB's, which are
# chr2:4600-4900 and 4100-4400
SYNTH_CHROM = 'GA--GB'
EXONS = [(100, 300), (500, 900), (1100, 1400), (1600, 1900)]


def _convert(tmp_path, thick_start, thick_end, aaseq_id=None):
    breakpoints = tmp_path / 'breakpoints.bed'
    breakpoints.write_text(f'{SYNTH_CHROM}\t1000\t1000\tbreakpoint-1--GA..chr1..1000..2000--GB..chr2..5000..4000\n')
    synthetic = tmp_path / 'synthetic.bed'
    with open(synthetic, 'w') as fh:
        FlairBed(SYNTH_CHROM, 100, 1900, 'FLT1', score=10, strand='+', thickStart=thick_start, thickEnd=thick_end,
                 blocks=[BedBlock(s, e) for s, e in EXONS], gene_id='FLG1', ref_gene_mappings=(SYNTH_CHROM,),
                 read_support=10, frac_support=1.0, productivity='PRO', samples=('s',), aaseq_id=aaseq_id).write(fh)
    out = tmp_path / 'fusions.bed'
    convert_synthetic_isos(str(synthetic), str(breakpoints), str(out), 100)
    return sorted(BedReader(str(out), bedClass=FlairBed), key=lambda bed: bed.pos_in_fusion)


def test_piece_the_coding_region_does_not_reach_has_none(tmp_path):
    # coding 150-700, all in GA: GB's piece used to be drawn all coding
    ga, gb = _convert(tmp_path, 150, 700)
    assert (ga.chrom, ga.chromStart, ga.chromEnd, ga.thickStart, ga.thickEnd) == ('chr1', 1100, 1900, 1150, 1700)
    assert (gb.chrom, gb.chromStart, gb.chromEnd, gb.strand) == ('chr2', 4100, 4900, '-')
    assert gb.thickStart == gb.thickEnd


def test_coding_region_across_the_breakpoint_is_split_between_the_pieces(tmp_path):
    # coding 600-1700: GA's piece from 1600 to its end, GB's from its 5' end (4900,
    # on -) to synthetic 1700, chr2:4300
    ga, gb = _convert(tmp_path, 600, 1700)
    assert (ga.thickStart, ga.thickEnd) == (1600, 1900)
    assert (gb.thickStart, gb.thickEnd) == (4300, 4900)


def test_coding_region_only_in_the_3_prime_piece(tmp_path):
    ga, gb = _convert(tmp_path, 1200, 1800)
    assert ga.thickStart == ga.thickEnd
    assert (gb.thickStart, gb.thickEnd) == (4200, 4800)


def test_noncoding_isoform_has_no_coding_region_in_any_piece(tmp_path):
    for piece in _convert(tmp_path, 100, 100):
        assert piece.thickStart == piece.thickEnd


def test_protein_id_is_carried_to_every_piece(tmp_path):
    assert [piece.aaseq_id for piece in _convert(tmp_path, 150, 700, aaseq_id='FLP00000001')] == ['FLP00000001'] * 2
