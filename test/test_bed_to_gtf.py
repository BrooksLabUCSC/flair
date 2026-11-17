"""Conversion of a BED to a GTF, in particular which thick range is a coding region"""
from flair.bed_to_gtf import bed_to_gtf


def _bed_row(chrom_start, chrom_end, thick_start, thick_end):
    "a two-exon BED12 row on the + strand"
    blocks = (chrom_end - chrom_start) - 900
    return '\t'.join(('chr1', str(chrom_start), str(chrom_end), 'tx1', '0', '+',
                      str(thick_start), str(thick_end), '0', '2', f'500,{blocks - 500},',
                      f'0,{(chrom_end - chrom_start) - (blocks - 500)},')) + '\n'


def _features(tmp_path, bed_row):
    "the features of the GTF written from one BED row"
    bed = tmp_path / 'in.bed'
    bed.write_text(bed_row)
    gtf = tmp_path / 'out.gtf'
    bed_to_gtf(query=str(bed), outputfile=str(gtf), force=True)
    return [line.split('\t')[2] for line in gtf.read_text().split('\n') if line]


def test_a_coding_region_covering_the_whole_transcript_is_written(tmp_path):
    # the CDS of a transcript with no UTR has exactly the transcript's bounds, which
    # used to be read as an absent coding region
    features = _features(tmp_path, _bed_row(1000, 3000, 1000, 3000))
    assert features.count('CDS') == 2
    assert features.count('start_codon') == 1
    assert '5UTR' not in features and '3UTR' not in features


def test_a_thick_range_of_no_length_is_not_a_coding_region(tmp_path):
    features = _features(tmp_path, _bed_row(1000, 3000, 1000, 1000))
    assert 'CDS' not in features
    assert features == ['gene', 'transcript', 'exon', 'exon']


def test_a_coding_region_within_the_transcript_has_a_utr_each_side(tmp_path):
    features = _features(tmp_path, _bed_row(1000, 3000, 1100, 2900))
    assert features.count('CDS') == 2
    assert features.count('5UTR') == 1 and features.count('3UTR') == 1
