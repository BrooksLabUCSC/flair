"""Splicing events of flair spliceevents, one row per alternative of an event.

Columns, in order:

  event                  the event, <type>|<chrom>:<start>-<end>(<strand>)|<gene>
  alternative            the alternative of the event the row counts, or 'other' for
                         the event's reads of none of its alternatives
  event_type             the kind of event
  gene                   the gene the event was found in
  gene_name              its name
  additional_genes       other genes sharing the event's junctions
  additional_gene_names  their names
  junctions_included     the junctions of the alternative, chrom:start-end
  junctions_excluded     the junctions the alternative skips
  outer_junctions        the junctions flanking the event
  exons                  the exons the alternative includes, chrom:start-end, or the
                         ends of a start or end event, chrom:position
  novel_junctions        the alternative's junctions that aren't annotated
  <sample>               one column per sample

The four gene columns and the five coordinate columns hold several values, comma
separated in the file and a list on the row.  The 'other' row's coordinate
columns are empty.

Coordinates are one based and inclusive, as the event names are.

Three files share these columns and differ only in what the sample columns hold:

  .spliceevents.counts.tsv    reads of the alternative
  .spliceevents.PSIjunc.tsv   its share of the event's reads
  .spliceevents.PSItot.tsv    its share of the reads covering the event's locus

A PSI column is empty when the sample has fewer than --event_support reads to take a
fraction from, rather than 0, which would read as "never included".  The sample
columns stay strings on the row for that reason.

Each partition writes its part of a file with the header, so a part is a valid TSV,
and only the first part's header is kept when they are joined.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType

GENE_COLUMNS = ('gene', 'gene_name', 'additional_genes', 'additional_gene_names')
COORDINATE_COLUMNS = ('junctions_included', 'junctions_excluded', 'outer_junctions', 'exons',
                      'novel_junctions')
EVENT_COLUMNS = ('event', 'alternative', 'event_type', *GENE_COLUMNS, *COORDINATE_COLUMNS)
TYPE_MAP = {col: commaListType for col in GENE_COLUMNS + COORDINATE_COLUMNS}

def columns(sample_columns):
    return list(EVENT_COLUMNS) + list(sample_columns)

class SplicingEventsReader(TsvReader):
    """Reads one of the three files.  The sample columns stay strings, since the
    counts file holds integers and the two PSI files hold fractions or nothing."""
    def __init__(self, splicing_events_tsv):
        super().__init__(splicing_events_tsv, defaultColType=str, typeMap=TYPE_MAP)

    @property
    def sample_columns(self):
        return self.columns[len(EVENT_COLUMNS):]

class SplicingEventsWriter(TsvWriter):
    def __init__(self, splicing_events_tsv, sample_columns, *, outFh=None):
        super().__init__(splicing_events_tsv, columns=columns(sample_columns),
                         defaultColType=str, typeMap=TYPE_MAP, outFh=outFh)
