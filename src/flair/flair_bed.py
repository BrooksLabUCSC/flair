"""
FLAIR BED record that is used to pass extra fields.
"""
from flair.pycbio.hgdata.bed import Bed, BedException, defaultIfNone
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
        row.extend([defaultIfNone(self.gene_id, ''),
                    defaultIfNone(self.ref_transcript_id, ''),
                    strArrayJoin(self.ref_gene_mappings),
                    defaultIfNone(self.read_support, ''),
                    defaultIfNone(round(self.frac_support, 4), ''),
                    defaultIfNone(self.productivity, '')])
        return row

    @classmethod
    def _parse(cls, row, fixScores=None):
        base = Bed.parse(row[:12], numStdCols=12, fixScores=fixScores)
        bed = cls(base.chrom, base.chromStart, base.chromEnd,
                  name=base.name, score=base.score, strand=base.strand,
                  thickStart=base.thickStart, thickEnd=base.thickEnd,
                  itemRgb=base.itemRgb, blocks=base.blocks)
        bed.gene_id = parseStrOrNone(row[12])
        bed.ref_transcript_id = parseStrOrNone(row[13])
        bed.ref_gene_mappings = strArraySplit(row[14])
        bed.read_support = int(row[15]) if row[15] != '' else None
        bed.frac_support = float(row[16]) if row[16] != '' else None
        bed.productivity = parseStrOrNone(row[17])
        return bed

    @classmethod
    def parse(cls, row, numStdCols=None, fixScores=None):  # numStdCols is only here for compatibility with BedReader
        needed_cols = 12 + len(cls.__slots__)
        if len(row) != needed_cols:
            raise BedException("expected at {} columns, found {}: ".format(needed_cols, len(row)))
        try:
            return cls._parse(row, fixScores=fixScores)
        except Exception as ex:
            raise BedException(f"parsing of BED row failed: {row}") from ex

    def get_named_extra_attrs(self):
        my_attrs = []
        if self.ref_transcript_id is not None:
            my_attrs.append(('ref_transcript_id', self.ref_transcript_id))
        if self.ref_gene_mappings is not None:
            my_attrs.append(('ref_gene_mappings', self.ref_gene_mappings))
        if self.read_support is not None:
            my_attrs.append(('read_support', self.read_support))
        if self.frac_support is not None:
            my_attrs.append(('frac_support', self.frac_support))
        if self.productivity is not None:
            my_attrs.append(('productivity', self.productivity))
        return my_attrs


def FlairBedReader(fspec):
    """Generator to read BED objects loaded from a tab-file or file-like
    object.  See Bed.parse()."""
    for bed in TabFileReader(fspec, rowClass=lambda r: FlairBed.parse(r),
                             hashAreComments=True, skipBlankLines=True):
        yield bed
