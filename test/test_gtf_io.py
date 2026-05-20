"""
tests of gff_io module
"""
import pytest
from io import StringIO
from flair.gtf_io import (gtf_data_parser, gtf_record_parser, GtfAttrsSet, GtfExon, GtfTranscript, GtfCDS,
                          GtfIdError, GtfParseError, gtf_write_row, FLAIR_ATTRS,
                          EXON_FEATURES, TRANSCRIPT_FEATURES, TRANSCRIPT_EXON_FEATURES, FLAIR_TRANSCRIPT_ATTRS)
from flair import SeqRange

GTF_FILE = "input/basic.annotation.gtf"
KRAS_TRANS_ID = 'ENST00000556131.1'

def test_gtf_data_parser_flair():
    gtf_data = gtf_data_parser(GTF_FILE, attrs=GtfAttrsSet.FLAIR)
    transcript = gtf_data.fetch_transcript(KRAS_TRANS_ID)
    # FIXME: add FLAIR_TRANSCRIPT_ATTRS
    assert frozenset(transcript.attrs.keys()) == FLAIR_TRANSCRIPT_ATTRS

def test_gtf_data_parser_all():
    gtf_data = gtf_data_parser(GTF_FILE, attrs=GtfAttrsSet.ALL)
    transcript = gtf_data.fetch_transcript(KRAS_TRANS_ID)
    for key in FLAIR_ATTRS:
        assert key in transcript.attrs
    assert 'level' in transcript.attrs
    assert 'transcript_support_level' in transcript.attrs

def test_gtf_record_parser_include_exons():
    recs = list(gtf_record_parser(GTF_FILE, include_features=EXON_FEATURES))
    assert len(recs) > 0
    assert all(isinstance(r, GtfExon) for r in recs)

def test_gtf_record_parser_include_transcripts():
    recs = list(gtf_record_parser(GTF_FILE, include_features=TRANSCRIPT_FEATURES))
    assert len(recs) > 0
    assert all(isinstance(r, GtfTranscript) for r in recs)

def test_gtf_record_parser_no_filter():
    recs = list(gtf_record_parser(GTF_FILE))
    types = {type(r) for r in recs}
    assert GtfExon in types
    assert GtfTranscript in types

def test_gtf_data_parser_include_features():
    gtf_data = gtf_data_parser(GTF_FILE, include_features=TRANSCRIPT_EXON_FEATURES)
    transcript = gtf_data.fetch_transcript(KRAS_TRANS_ID)
    assert len(transcript.exons) == 3
    assert len(transcript.cds_recs) == 0

def test_gtf_data_exons():
    gtf_data = gtf_data_parser(GTF_FILE)
    transcript = gtf_data.fetch_transcript(KRAS_TRANS_ID)
    assert len(transcript.exons) == 3
    assert all(isinstance(e, GtfExon) for e in transcript.exons)

def test_gtf_data_cds():
    gtf_data = gtf_data_parser(GTF_FILE)
    transcript = gtf_data.fetch_transcript(KRAS_TRANS_ID)
    assert len(transcript.cds_recs) == 2
    assert all(isinstance(c, GtfCDS) for c in transcript.cds_recs)

def test_parse(basic_gtf_data):
    trans_ids = sorted(basic_gtf_data.iter_transcript_ids())
    assert len(trans_ids) == 59
    assert trans_ids[0:10] == ['ENST00000225792.10', 'ENST00000256078.9', 'ENST00000279052.10', 'ENST00000311936.8', 'ENST00000348547.6',
                               'ENST00000357394.8', 'ENST00000411577.5', 'ENST00000413587.5', 'ENST00000416206.5', 'ENST00000438317.5']

def test_trans_lookup_flair(basic_gtf_data):
    transcript = basic_gtf_data.get_transcript('ENST00000556131.1')
    assert transcript is not None
    assert transcript.transcript_id == 'ENST00000556131.1'
    assert transcript.gene_id == "ENSG00000133703.12"
    assert transcript.gene_name == "KRAS"
    assert transcript.coords == SeqRange(name='chr12', start=25233818, end=25250929, strand='-')
    assert transcript.coords_no_strand == SeqRange(name='chr12', start=25233818, end=25250929)
    assert frozenset(transcript.attrs.keys()) == FLAIR_TRANSCRIPT_ATTRS

def test_trans_lookup_all(basic_gtf_data_all):
    transcript = basic_gtf_data_all.get_transcript('ENST00000556131.1')
    assert transcript is not None
    assert transcript.transcript_id == 'ENST00000556131.1'
    assert transcript.gene_id == "ENSG00000133703.12"
    assert transcript.gene_name == "KRAS"
    assert transcript.coords == SeqRange(name='chr12', start=25233818, end=25250929, strand='-')
    assert transcript.coords_no_strand == SeqRange(name='chr12', start=25233818, end=25250929)
    assert transcript.attrs['level'] == 2
    assert transcript.attrs['transcript_support_level'] == "1"

def test_trans_error(basic_gtf_data):
    transcript = basic_gtf_data.get_transcript('Fred')
    assert transcript is None
    with pytest.raises(GtfIdError, match=r"unknown transcript id `Barney'"):
        basic_gtf_data.fetch_transcript('Barney')

def test_chroms(basic_gtf_data):
    chroms = basic_gtf_data.get_chroms()
    assert chroms == ['chr12', 'chr17', 'chr20']


# overlaps KRAS
KRAS_OVER_RANGE = SeqRange("chr12", 25205000, 25252000)
KRAS_TRANS_IDS = ["ENST00000256078.9", "ENST00000311936.8", "ENST00000556131.1", "ENST00000557334.5"]

def test_iter_overlap(basic_gtf_data):
    trans_ids = sorted([t.transcript_id
                        for t in basic_gtf_data.iter_overlap_transcripts_sr(KRAS_OVER_RANGE)])
    assert trans_ids == KRAS_TRANS_IDS

def test_iter_overlap_strand(basic_gtf_data):
    seq_range = SeqRange(*KRAS_OVER_RANGE[0:3], '-')
    trans_ids = sorted([t.transcript_id
                        for t in basic_gtf_data.iter_overlap_transcripts_sr(seq_range)])
    assert trans_ids == KRAS_TRANS_IDS

def test_iter_overlap_other_strand(basic_gtf_data):
    seq_range = SeqRange(*KRAS_OVER_RANGE[0:3], '+')
    trans_ids = sorted([t.transcript_id
                        for t in basic_gtf_data.iter_overlap_transcripts_sr(seq_range)])
    assert trans_ids == []

def test_rec_str(basic_gtf_data_all):
    trans_str = str(basic_gtf_data_all.fetch_transcript("ENST00000256078.9"))
    assert trans_str == ('chr12\tHAVANA\ttranscript\t25205246\t25250929\t.\t-\t.\t'
                         'gene_id "ENSG00000133703.12"; transcript_id "ENST00000256078.9"; gene_type "protein_coding"; '
                         'gene_name "KRAS"; transcript_type "protein_coding"; transcript_name "KRAS-201"; level 2; '
                         'protein_id "ENSP00000256078.4"; transcript_support_level "1"; hgnc_id "HGNC:6407"; '
                         'tag "RNA_Seq_supported_only"; tag "basic"; tag "appris_principal_4"; tag "CCDS"; ccdsid "CCDS8703.1"; '
                         'havana_gene "OTTHUMG00000171193.4"; havana_transcript "OTTHUMT00000412232.4";')

def test_rec_str1(basic_gtf_data):
    rec_fh = StringIO()

    gtf_write_row(rec_fh, 'chr12', 'HAVANA', 'transcript', 25205246, 25250929, None, '-', None,
                  gene_id="ENSG00000133703.12", transcript_id="ENST00000256078.9", gene_name="KRAS")

    assert rec_fh.getvalue() == ('chr12\tHAVANA\ttranscript\t25205247\t25250929\t.\t-\t.\t'
                                 'gene_id "ENSG00000133703.12"; gene_name "KRAS"; transcript_id "ENST00000256078.9";\n')

def test_rec_str2(basic_gtf_data):
    rec_fh = StringIO()

    gtf_write_row(rec_fh, 'chr12', 'HAVANA', 'transcript', 25205246, 25250929, 10, '-', 1,
                  gene_id="ENSG00000133703.12", transcript_id="ENST00000256078.9", gene_name="KRAS")

    assert rec_fh.getvalue() == ('chr12\tHAVANA\ttranscript\t25205247\t25250929\t10\t-\t1\t'
                                 'gene_id "ENSG00000133703.12"; gene_name "KRAS"; transcript_id "ENST00000256078.9";\n')

def test_rec_str3(basic_gtf_data):
    rec_fh = StringIO()

    gtf_write_row(rec_fh, 'chr12', 'HAVANA', 'transcript', 25205246, 25250929, 100, '-', 2,
                  attrs={"gene_id": "ENSG00000133703.12",
                         "transcript_id": "ENST00000256078.9",
                         "gene_name": "KRAS"})

    assert rec_fh.getvalue() == ('chr12\tHAVANA\ttranscript\t25205247\t25250929\t100\t-\t2\t'
                                 'gene_id "ENSG00000133703.12"; transcript_id "ENST00000256078.9"; gene_name "KRAS";\n')


##
# Whitespace validation tests
##

def _write_gtf(path, gene_id, transcript_id):
    """Write a minimal GTF with one transcript line to path."""
    attrs = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "FOO";'
    path.write_text(f'chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\t{attrs}\n')
    return str(path)


def test_whitespace_space_in_gene_id(tmp_path):
    with pytest.raises(GtfParseError):
        list(gtf_record_parser(_write_gtf(tmp_path / "test.gtf", "ENSG1 BAD", "ENST1")))


def test_whitespace_tab_in_gene_id(tmp_path):
    with pytest.raises(GtfParseError):
        list(gtf_record_parser(_write_gtf(tmp_path / "test.gtf", "ENSG1\tBAD", "ENST1")))


def test_whitespace_space_in_transcript_id(tmp_path):
    with pytest.raises(GtfParseError):
        list(gtf_record_parser(_write_gtf(tmp_path / "test.gtf", "ENSG1", "ENST1 BAD")))


def test_whitespace_tab_in_transcript_id(tmp_path):
    with pytest.raises(GtfParseError):
        list(gtf_record_parser(_write_gtf(tmp_path / "test.gtf", "ENSG1", "ENST1\tBAD")))


def test_no_whitespace_valid(tmp_path):
    recs = list(gtf_record_parser(_write_gtf(tmp_path / "test.gtf", "ENSG1", "ENST1")))
    assert len(recs) == 1
    assert recs[0].gene_id == "ENSG1"
    assert recs[0].transcript_id == "ENST1"
