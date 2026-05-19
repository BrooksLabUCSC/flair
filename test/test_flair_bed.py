"""
tests of flair_bed module
"""
import pytest
from flair.flair_bed import FlairBed, FlairBedReader, parseStrOrNone
from flair.pycbio.hgdata.bed import BedException


def _make_bed(**overrides):
    """Construct a 12-col FlairBed with one block and all extra attrs set."""
    kwargs = dict(chrom="chr1", chromStart=100, chromEnd=200, name="tx1",
                  score=0, strand="+", thickStart=100, thickEnd=200, itemRgb="0",
                  blocks=[(100, 200)],
                  gene_id="geneA", ref_transcript_id="refTx1",
                  ref_gene_mappings=["g1:tx1", "g2:tx2"], read_support=5,
                  frac_support=0.5, productivity="PRO")
    kwargs.update(overrides)
    blocks = kwargs.pop("blocks")
    bed = FlairBed(**{k: v for k, v in kwargs.items() if k != "blocks"})
    if blocks is not None:
        for s, e in blocks:
            bed.addBlock(s, e)
    return bed


def test_parse_str_or_none_empty():
    assert parseStrOrNone("") is None


def test_parse_str_or_none_value():
    assert parseStrOrNone("foo") == "foo"


def test_construct_minimal():
    bed = FlairBed("chr1", 0, 10)
    assert bed.chrom == "chr1"
    assert bed.chromStart == 0
    assert bed.chromEnd == 10
    assert bed.name is None
    assert bed.gene_id is None
    assert bed.ref_transcript_id is None
    assert bed.ref_gene_mappings == []
    assert bed.read_support is None
    assert bed.frac_support is None
    assert bed.productivity is None


def test_construct_full():
    bed = _make_bed()
    assert bed.gene_id == "geneA"
    assert bed.ref_transcript_id == "refTx1"
    assert bed.ref_gene_mappings == ["g1:tx1", "g2:tx2"]
    assert bed.read_support == 5
    assert bed.frac_support == 0.5
    assert bed.productivity == "PRO"


def test_transcript_id_property():
    bed = FlairBed("chr1", 0, 10, name="tx1")
    assert bed.transcript_id == "tx1"
    assert bed.transcript_id == bed.name


def test_transcript_id_setter():
    bed = FlairBed("chr1", 0, 10, name="tx1")
    bed.transcript_id = "tx2"
    assert bed.name == "tx2"
    assert bed.transcript_id == "tx2"


def test_thick_start_passed_to_super():
    """thickStart must reach the underlying Bed."""
    bed = FlairBed("chr1", 0, 100, name="tx", score=0, strand="+",
                   thickStart=10, thickEnd=90)
    assert bed.thickStart == 10
    assert bed.thickEnd == 90


def test_num_columns():
    bed = _make_bed()
    # 12 standard + 6 flair extras
    assert bed.numColumns == 18


def test_to_row_length():
    bed = _make_bed()
    row = bed.toRow()
    assert len(row) == 18


def test_to_row_values():
    bed = _make_bed()
    row = bed.toRow()
    # standard cols 0..11
    assert row[0] == "chr1"
    assert row[1] == "100"
    assert row[2] == "200"
    assert row[3] == "tx1"
    # extras appended in __slots__ order
    assert row[12] == "geneA"
    assert row[13] == "refTx1"
    assert row[14] == "g1:tx1,g2:tx2,"
    assert row[15] == 5
    assert row[16] == 0.5
    assert row[17] == "PRO"


def test_parse_wrong_column_count():
    with pytest.raises(BedException, match="expected at"):
        FlairBed.parse(["chr1", "0", "10"])


def test_parse_round_trip():
    """Build a BED row and parse it back."""
    row = ["chr1", "100", "200", "tx1", "0", "+",
           "100", "200", "0", "1", "100,", "0,",
           "geneA", "refTx1", "g1:tx1,g2:tx2,", "5", "0.5", "PRO"]
    bed = FlairBed.parse(row)
    assert bed.chrom == "chr1"
    assert bed.chromStart == 100
    assert bed.chromEnd == 200
    assert bed.name == "tx1"
    assert bed.transcript_id == "tx1"
    assert bed.gene_id == "geneA"
    assert bed.ref_transcript_id == "refTx1"
    assert bed.ref_gene_mappings == ["g1:tx1", "g2:tx2"]
    assert bed.read_support == 5
    assert bed.frac_support == 0.5
    assert bed.productivity == "PRO"


def test_parse_empty_string_extras_become_none():
    row = ["chr1", "100", "200", "tx1", "0", "+",
           "100", "200", "0", "1", "100,", "0,",
           "", "", "", "0", "0.0", ""]
    bed = FlairBed.parse(row)
    assert bed.gene_id is None
    assert bed.ref_transcript_id is None
    assert bed.ref_gene_mappings == []
    assert bed.read_support == 0
    assert bed.frac_support == 0.0


def test_reader(tmp_path):
    row = ["chr1", "100", "200", "tx1", "0", "+",
           "100", "200", "0", "1", "100,", "0,",
           "geneA", "refTx1", "g1:tx1,g2:tx2,", "5", "0.5", "PRO"]
    p = tmp_path / "test.bed"
    p.write_text("\t".join(row) + "\n")
    beds = list(FlairBedReader(str(p)))
    assert len(beds) == 1
    assert beds[0].chrom == "chr1"
    assert beds[0].gene_id == "geneA"
    assert beds[0].productivity == "PRO"
