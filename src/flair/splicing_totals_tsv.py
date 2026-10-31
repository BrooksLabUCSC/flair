"""The reads each splicing event of flair spliceevents divides its PSIs by,
.spliceevents.totals.tsv: two rows per event.

Columns, in order:

  event                  the event, as in the event files (splicing_events_tsv)
  event_type             the kind of event
  gene                   the gene the event was found in
  gene_name              its name
  additional_genes       other genes sharing the event's junctions
  additional_gene_names  their names
  total                  which total the row is: 'event_reads', the event's reads,
                         which PSIjunc divides by, or 'locus_reads', the reads covering
                         its locus, which PSItot divides by
  <sample>               one integer column per sample

The four gene columns hold several values, comma separated in the file and a list
on the row.  spliceevents_outliers takes its PSIs from the counts file and these.

Each partition writes its part with the header, and only the first part's header is
kept when they are joined.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType, read_header_columns
from flair.splicing_events_tsv import GENE_COLUMNS

ID_COLUMNS = ('event', 'event_type', *GENE_COLUMNS, 'total')

def columns(sample_columns):
    return list(ID_COLUMNS) + list(sample_columns)

def _type_map(sample_columns):
    type_map = {col: commaListType for col in GENE_COLUMNS}
    type_map.update({col: int for col in sample_columns})
    return type_map

class SplicingTotalsReader(TsvReader):
    def __init__(self, splicing_totals_tsv):
        sample_columns = read_header_columns(splicing_totals_tsv)[len(ID_COLUMNS):]
        super().__init__(splicing_totals_tsv, defaultColType=str, typeMap=_type_map(sample_columns))

    @property
    def sample_columns(self):
        return self.columns[len(ID_COLUMNS):]

class SplicingTotalsWriter(TsvWriter):
    def __init__(self, splicing_totals_tsv, sample_columns, *, outFh=None):
        super().__init__(splicing_totals_tsv, columns=columns(sample_columns), defaultColType=str,
                         typeMap=_type_map(sample_columns), outFh=outFh)
