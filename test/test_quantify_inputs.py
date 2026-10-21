import pytest
import pysam
from flair import FlairInputDataError
from flair.flair_bed import FlairBed
from flair.pycbio.hgdata.bed import BedBlock
from flair.io_utils import read_manifest_rows
from flair.flair_quantify import load_manifest, load_isoform_data
from flair.flair_spliceevents import read_manifest as spliceevents_manifest

_HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 2000}]})


def _bam(path, index=True):
    with pysam.AlignmentFile(str(path), 'wb', header=_HEADER):
        pass
    if index:
        pysam.index(str(path))
    return str(path)


def _manifest(tmp_path, *lines):
    path = tmp_path / 'manifest.tsv'
    path.write_text(''.join(line + '\n' for line in lines))
    return str(path)


def _rows(path):
    return read_manifest_rows(path, (2,), 'a sample name and its bam file', {1: 'bam file'})


def test_manifest_skips_comments_and_blank_lines(tmp_path):
    bam = _bam(tmp_path / 's1.bam')
    manifest = _manifest(tmp_path, '# sample\tbam', '', f's1\t{bam}\r', '   ')
    assert _rows(manifest) == [['s1', bam]]


def test_manifest_line_errors(tmp_path):
    bam = _bam(tmp_path / 's1.bam')
    with pytest.raises(FlairInputDataError, match=r'manifest.tsv:1: a manifest line is a sample name and its bam file, '
                                                  r'tab-separated: '):
        _rows(_manifest(tmp_path, f's1\t{bam}\textra'))
    with pytest.raises(FlairInputDataError, match='tab-separated; fields are separated by tabs, not spaces'):
        _rows(_manifest(tmp_path, f's1 {bam}'))
    with pytest.raises(FlairInputDataError, match='manifest.tsv:1: a manifest line'):
        _rows(_manifest(tmp_path, f'\t{bam}'))
    with pytest.raises(FlairInputDataError, match=r'manifest.tsv:2: bam file .*s2.bam does not exist'):
        _rows(_manifest(tmp_path, f's1\t{bam}', f's2\t{tmp_path}/s2.bam'))


def test_manifest_samples_named_once_and_some(tmp_path):
    bam = _bam(tmp_path / 's1.bam')
    with pytest.raises(FlairInputDataError, match='more than once: s1'):
        _rows(_manifest(tmp_path, f's1\t{bam}', f's1\t{bam}'))
    with pytest.raises(FlairInputDataError, match='no samples'):
        _rows(_manifest(tmp_path, '# nothing'))


def test_quantify_manifest_needs_indexed_bams(tmp_path):
    bam, unindexed = _bam(tmp_path / 's1.bam'), _bam(tmp_path / 's2.bam', index=False)
    assert load_manifest(_manifest(tmp_path, f's1\tA\tb0\t{bam}')) == [('s1', 'A', 'b0', bam)]
    with pytest.raises(FlairInputDataError, match="sample s2's BAM .*s2.bam has no index"):
        load_manifest(_manifest(tmp_path, f's1\tA\tb0\t{bam}', f's2\tB\tb0\t{unindexed}'))
    with pytest.raises(FlairInputDataError, match='a sample id, its condition, its batch'):
        load_manifest(_manifest(tmp_path, f's1\t{bam}'))


def test_spliceevents_manifest(tmp_path):
    bam = _bam(tmp_path / 's1.bam')
    assert spliceevents_manifest(_manifest(tmp_path, '# sample\tbam', f's1\t{bam}')) == [('s1', bam)]


def _iso(name, gene_id, chrom='chr1', strand='+', exons=((100, 200), (300, 400))):
    return FlairBed(chrom, exons[0][0], exons[-1][1], name, strand=strand, blocks=[BedBlock(s, e) for s, e in exons],
                    gene_id=gene_id)


def _bed(tmp_path, *beds):
    path = tmp_path / 'isoforms.bed'
    with open(path, 'w') as fh:
        for bed in beds:
            bed.write(fh)
    return str(path)


def test_isoform_bed_errors(tmp_path):
    with pytest.raises(FlairInputDataError, match="isoform 'T1' has no gene id"):
        load_isoform_data(_bed(tmp_path, _iso('T1', None)))
    with pytest.raises(FlairInputDataError, match="gene id 'G1' has isoforms on chr1 \\+ and on chr2 \\+"):
        load_isoform_data(_bed(tmp_path, _iso('T1', 'G1'), _iso('T2', 'G1', chrom='chr2')))
    with pytest.raises(FlairInputDataError, match='chr1 \\+ and on chr1 -.*current flair combine'):
        load_isoform_data(_bed(tmp_path, _iso('T1', 'G1'), _iso('T2', 'G1', strand='-')))
    with pytest.raises(FlairInputDataError, match="isoform name 'T1' is used more than once on chr1"):
        load_isoform_data(_bed(tmp_path, _iso('T1', 'G1'), _iso('T1', 'G2')))
    with pytest.raises(FlairInputDataError, match='no isoforms'):
        load_isoform_data(_bed(tmp_path))


def test_plain_bed12_is_not_a_flair_isoform_bed(tmp_path):
    path = tmp_path / 'plain.bed'
    path.write_text('chr1\t100\t400\tT1\t0\t+\t100\t100\t0\t2\t100,100,\t0,200,\n')
    with pytest.raises(FlairInputDataError, match='is not a FLAIR isoform bed, of the 25 columns'):
        load_isoform_data(str(path))
