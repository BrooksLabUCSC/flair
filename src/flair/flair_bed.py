"""
FLAIR BED record that is used to pass extra fields.
"""
from flair.pycbio.hgdata.bed import Bed, BedException
from flair.pycbio.hgdata.autoSql import strArraySplit, strArrayJoin
from flair.pycbio.tsv.tabFile import TabFileReader

def parseStrOrNone(s):
    return None if len(s) == 0 else s


class FlairBed(Bed):
    """
    BED class that passes along flair-derived attributes.

    This enforces the transcript_id and name columns having the same value
    """
    __slots__ = ("gene_id", "ref_transcript_id", "ref_gene_mappings", "read_support",
                 "frac_support", "productivity")

    def __init__(self, chrom, chromStart, chromEnd, name=None, *, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None, blocks=None,
                 gene_id=None, ref_transcript_id=None, ref_gene_mappings=None,
                 read_support=None, frac_support=None, productivity=None):
        super().__init__(chrom=chrom, chromStart=chromStart, chromEnd=chromEnd,
                         name=name, score=score, strand=strand, thickStart=thickStart, thickEnd=thickEnd,
                         itemRgb=itemRgb, blocks=blocks, numStdCols=12)
        self.gene_id = gene_id
        self.ref_transcript_id = ref_transcript_id
        self.ref_gene_mappings = [] if ref_gene_mappings is None else list(ref_gene_mappings)
        self.read_support = read_support
        self.frac_support = frac_support
        self.productivity = productivity

    @property
    def transcript_id(self):
        return self.name

    @transcript_id.setter
    def transcript_id(self, value):
        self.name = value

    @property
    def numColumns(self):
        """Returns the number of columns in the BED when formatted as a row."""
        return super().numColumns + len(self.__slots__)

    def toRow(self):
        row = super().toRow()
        row.extend([self.gene_id,
                    self.ref_transcript_id,
                    strArrayJoin(self.ref_gene_mappings),
                    self.read_support,
                    self.frac_support,
                    self.productivity])
        return row

    @classmethod
    def _parse(cls, row):
        base = Bed.parse(row[:12], numStdCols=12)
        bed = cls(base.chrom, base.chromStart, base.chromEnd,
                  name=base.name, score=base.score, strand=base.strand,
                  thickStart=base.thickStart, thickEnd=base.thickEnd,
                  itemRgb=base.itemRgb, blocks=base.blocks)
        bed.gene_id = parseStrOrNone(row[12])
        bed.ref_transcript_id = parseStrOrNone(row[13])
        bed.ref_gene_mappings = strArraySplit(row[14])
        bed.read_support = int(row[15])
        bed.frac_support = float(row[16])
        bed.productivity = row[17]
        return bed

    @classmethod
    def parse(cls, row):
        needed_cols = 12 + len(cls.__slots__)
        if len(row) != needed_cols:
            raise BedException("expected at {} columns, found {}: ".format(needed_cols, len(row)))
        try:
            return cls._parse(row)
        except Exception as ex:
            raise BedException(f"parsing of BED row failed: {row}") from ex


def FlairBedReader(fspec):
    """Generator to read BED objects loaded from a tab-file or file-like
    object.  See Bed.parse()."""
    for bed in TabFileReader(fspec, rowClass=lambda r: FlairBed.parse(r),
                             hashAreComments=True, skipBlankLines=True):
        yield bed
