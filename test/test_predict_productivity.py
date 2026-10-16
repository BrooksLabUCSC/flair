from types import SimpleNamespace

from flair import predictProductivity
from flair.isoform_data import Exon
from flair.predictProductivity import calc_genomic_end_pos, checked_cds, get_annot_start_codons, predict_productivity

# a two-exon transcript: exons [100, 200) and [300, 350), 150 bases
PLUS_EXONS = [Exon(100, 200), Exon(300, 350)]   # 5' to 3'
MINUS_EXONS = PLUS_EXONS[::-1]                   # 5' to 3' on minus strand
SIZES = [100, 50]


def test_plus_strand_inside_exon():
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 30, '+') == 130
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 120, '+') == 320


def test_minus_strand_inside_exon():
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 30, '-') == 320
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 80, '-') == 170


def test_stop_codon_at_transcript_end():
    # orf_end_pos equal to the transcript length used to give None
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 150, '+') == 350
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 150, '-') == 100


def test_stop_codon_at_exon_end():
    # the stop codon's last base is the last base of the first exon: the ORF ends
    # there, not at the first base of the next exon
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 100, '+') == 200
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 50, '-') == 300


def test_stop_codon_at_exon_start():
    # one base into the next exon
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 101, '+') == 301
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 51, '-') == 199


def test_orf_end_outside_transcript():
    assert calc_genomic_end_pos(PLUS_EXONS, SIZES, 151, '+') is None
    assert calc_genomic_end_pos(MINUS_EXONS, SIZES[::-1], 0, '-') is None


def test_checked_cds():
    assert checked_cds(120, 320, 100, 350) == (120, 320)
    assert checked_cds(100, 100, 100, 350) == (100, 100)      # no CDS
    assert checked_cds(None, 320, 100, 350) is None
    assert checked_cds(120, None, 100, 350) is None
    assert checked_cds(320, 120, 100, 350) is None            # thickStart after thickEnd
    assert checked_cds(90, 320, 100, 350) is None             # outside the transcript
    assert checked_cds(120, 360, 100, 350) is None


def _transcript(strand):
    return SimpleNamespace(gene=SimpleNamespace(gene_desc=['G1']), exons=PLUS_EXONS, strand=strand,
                           ref_transcript_id=None, chrom='chr1', start=100, end=350, name='iso1',
                           get_sequence=lambda genome: 'A' * 150)


def _predict(monkeypatch, strand, orf):
    monkeypatch.setattr(predictProductivity, 'identify_best_orf_from_starts', lambda *args: orf)
    return predict_productivity(_transcript(strand), 1, {'G1': {strand: {120}}}, {}, genome=None)


def test_good_cds_kept(monkeypatch):
    assert _predict(monkeypatch, '+', ['PRO', 120, 320, 150, 0, 'MK']) == (120, 320, 'PRO', 'MK')
    assert _predict(monkeypatch, '-', ['PRO', 320, 120, 150, 0, 'MK']) == (120, 320, 'PRO', 'MK')


def test_missing_orf_end_gives_no_cds(monkeypatch):
    # a missing ORF end used to become thickStart = chromEnd on the minus strand,
    # after thickEnd, which bedToBigBed rejects
    assert _predict(monkeypatch, '-', ['PRO', 320, None, 150, 0, 'MK']) == (100, 100, 'NGO', None)
    assert _predict(monkeypatch, '+', ['PRO', 120, None, 150, 0, 'MK']) == (100, 100, 'NGO', None)


def test_inverted_cds_gives_no_cds(monkeypatch):
    assert _predict(monkeypatch, '+', ['PTC', 320, 120, 150, 0, None]) == (100, 100, 'NGO', None)


def test_start_codons_on_the_transcripts_strand_only():
    # start codons by gene and strand, as annotation_data saves them
    cds_starts = {'G1': {'+': {120}}, 'G2': {'+': {130}, '-': {340}}, 'G3': {'-': {330}}}
    transcript = SimpleNamespace(gene=SimpleNamespace(gene_desc=['G1', 'G2', 'G3', 'NOVEL']), strand='+')
    # other genes' start codons are used when on the same strand
    assert get_annot_start_codons(transcript, cds_starts) == {120, 130}
    transcript.strand = '-'
    assert get_annot_start_codons(transcript, cds_starts) == {340, 330}


def test_no_start_codons_on_the_strand():
    transcript = SimpleNamespace(gene=SimpleNamespace(gene_desc=['G1']), strand='-')
    assert get_annot_start_codons(transcript, {'G1': {'+': {120}}}) == set()
