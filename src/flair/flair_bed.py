"""
FLAIR BED record that is used to pass extra fields.
"""
from flair.pycbio.hgdata.bed import Bed, BedException, defaultIfNone
from flair.pycbio.hgdata.autoSql import strArraySplitNone, strArrayJoin
from flair.pycbio.tsv.tabFile import TabFileReader

def parseStrOrNone(s):
    return None if len(s) == 0 else s


def get_strand_rgb(strand, junclen):
    if junclen == 0:
        return "99,99,99"
    elif strand == '+':
        return "27,158,119"
    else:
        return "217,95,2"


class FlairBed(Bed):
    """
    BED class that passes along flair-derived attributes.

    This enforces the transcript_id and name columns having the same value
    """
    __slots__ = ("gene_id", "ref_transcript_id", "ref_gene_mappings", "read_support",
                 "frac_support", "productivity", "transcript_class", "fused_genes", "pos_in_fusion", "samples", "source_isoform", "allele_group", "aaseq_id")

    def __init__(self, chrom, chromStart, chromEnd, name=None, *, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None, blocks=None,
                 gene_id=None, ref_transcript_id=None, ref_gene_mappings=None,
                 read_support=None, frac_support=None, productivity=None, transcript_class='basic',
                 fused_genes=(), pos_in_fusion=None, samples=(), aaseq_id=None, source_isoform=None, allele_group=None):
        super().__init__(chrom=chrom, chromStart=chromStart, chromEnd=chromEnd,
                         name=name, score=score, strand=strand, thickStart=thickStart, thickEnd=thickEnd,
                         itemRgb=itemRgb, blocks=blocks, numStdCols=12)
        self.gene_id = gene_id
        self.ref_transcript_id = ref_transcript_id
        self.ref_gene_mappings = [] if ref_gene_mappings is None else list(ref_gene_mappings)
        self.read_support = read_support
        self.frac_support = frac_support
        self.productivity = productivity
        self.transcript_class = transcript_class
        self.fused_genes = fused_genes
        self.pos_in_fusion = pos_in_fusion
        self.samples = samples
        self.aaseq_id = aaseq_id
        self.source_isoform = source_isoform
        self.allele_group = allele_group
        if itemRgb is None:
            self.itemRgb = self._get_rgb()

    def _get_rgb(self):
        PRODUCTIVITY_COLORS = {"PRO": "103,169,207", "PTC": "239,138,98", "NST": "0,0,0", "NGO": "0,0,0"}
        if self.transcript_class == 'fusion':
            return '179,46,29'
        elif self.transcript_class == 'readthrough':
            return '135,56,181'
        elif self.productivity is not None:
            return PRODUCTIVITY_COLORS[self.productivity]
        else:
            return get_strand_rgb(self.strand, max((self.blockCount, 0)))

    def addBlock(self, start, end):
        super().addBlock(start, end)
        if self.itemRgb == "99,99,99":  # color that is dependent on having few blocks
            self.itemRgb = self._get_rgb()

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
                    defaultIfNone(self.productivity, ''),
                    self.transcript_class,
                    strArrayJoin(self.fused_genes),
                    defaultIfNone(self.pos_in_fusion, ''),
                    strArrayJoin(self.samples),
                    defaultIfNone(self.source_isoform, ''),
                    defaultIfNone(self.allele_group, ''),
                    defaultIfNone(self.aaseq_id, '')
                    ])
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
        bed.ref_gene_mappings = tuple(strArraySplitNone(row[14]))
        bed.read_support = int(row[15]) if row[15] != '' else None
        bed.frac_support = float(row[16]) if row[16] != '' else None
        bed.productivity = parseStrOrNone(row[17])
        bed.transcript_class = row[18]
        bed.fused_genes = tuple(strArraySplitNone(row[19]))
        bed.pos_in_fusion = int(row[20]) if row[20] != '' else None
        bed.samples = tuple(strArraySplitNone(row[21]))
        bed.source_isoform = parseStrOrNone(row[22])
        bed.allele_group = parseStrOrNone(row[23])
        bed.aaseq_id = parseStrOrNone(row[24])
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

    def get_named_extra_attrs(self):  # noqa: C901
        my_attrs = []
        if self.ref_transcript_id is not None:
            my_attrs.append(('ref_transcript_id', self.ref_transcript_id))
        if len(self.ref_gene_mappings) > 0:
            my_attrs.append(('ref_gene_mappings', strArrayJoin(self.ref_gene_mappings)))
        if self.read_support is not None:
            my_attrs.append(('read_support', self.read_support))
        if self.frac_support is not None:
            my_attrs.append(('frac_support', round(self.frac_support, 4)))
        if self.productivity is not None:
            my_attrs.append(('productivity', self.productivity))
        my_attrs.append(('transcript_class', self.transcript_class))
        if len(self.fused_genes) > 0:
            my_attrs.append(('fused_genes', strArrayJoin(self.fused_genes)))
        if self.pos_in_fusion is not None:
            my_attrs.append(('pos_in_fusion', self.pos_in_fusion))
        if len(self.samples) > 0:
            my_attrs.append(('samples', strArrayJoin(self.samples)))
        if self.source_isoform is not None:
            my_attrs.append(('source_isoform', self.source_isoform))
        if self.allele_group is not None:
            my_attrs.append(('allele_group', self.allele_group))
        if self.aaseq_id is not None:
            my_attrs.append(('aaseq_id', self.aaseq_id))
        return my_attrs


def FlairBedReader(fspec):
    """Generator to read BED objects loaded from a tab-file or file-like
    object.  See Bed.parse()."""
    for bed in TabFileReader(fspec, rowClass=lambda r: FlairBed.parse(r),
                             hashAreComments=True, skipBlankLines=True):
        yield bed
