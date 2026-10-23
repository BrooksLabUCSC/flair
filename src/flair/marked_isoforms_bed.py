"""An isoform BED with a column saying whether the isoform retains an intron.

Columns, in order:

  the columns of the input BED, unchanged.  That is twelve for a plain BED and
  twenty five for a FLAIR BED, so the width is taken from the input rather than
  fixed here
  retains_intron   1 when another isoform has an exon where this one has an
                   intron, 0 otherwise

Coordinates are zero based and half open, as they come from a BED.

Written by mark_intron_retention.  The input columns are passed through so the
marked file stays whatever kind of BED it was given.
"""
from flair import FlairInputDataError
from flair.pycbio.tsv import TsvReader, TsvWriter

RETAINS_INTRON_COLUMN = 'retains_intron'

def bed_columns(num_bed_columns):
    "names for the pass-through BED columns, which are addressed by position"
    return [f'col{num}' for num in range(1, num_bed_columns + 1)]

def bed_column_count(bed_rows):
    "the shared width of the input rows, which the marked file keeps"
    widths = {len(row) for row in bed_rows}
    if len(widths) != 1:
        raise FlairInputDataError(
            "the isoform BED rows are not all the same width, found: "
            f"{', '.join(str(width) for width in sorted(widths))}")
    return widths.pop()

class MarkedIsoformsReader(TsvReader):
    """The BED columns are passed through, so they are addressed by the col1, col2
    names the header gives them; only retains_intron is named for what it holds."""
    def __init__(self, marked_bed):
        super().__init__(marked_bed, defaultColType=str,
                         typeMap={RETAINS_INTRON_COLUMN: int})

class MarkedIsoformsWriter(TsvWriter):
    def __init__(self, marked_bed, num_bed_columns):
        super().__init__(marked_bed, columns=bed_columns(num_bed_columns) + [RETAINS_INTRON_COLUMN],
                         defaultColType=str)

    def writeIsoform(self, bed_row, retains_intron):
        self.writeRow(list(bed_row) + [1 if retains_intron else 0])
