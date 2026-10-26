"""Alternative splicing events with the reads supporting each side, one row per
side of each event.

Columns, in order:

  feature_id   names one side of one event, inclusion_<event> or exclusion_<event>
  coordinate   the event the two sides belong to; DRIMSeq groups by this the way
               it groups isoforms by gene
  <sample>     one column per sample, the reads supporting this side
  isoform_ids  ids of the isoforms on this side, comma separated in the file
               and a list on the row

Both end in the strand, chr12:25215436-25215560:+, because the same coordinates
on the two strands are two events: without it they share a feature_id, which
DRIMSeq refuses, and a coordinate, which would pool the strands into one test.

Coordinates are zero based and half open, as they come from a BED: a
junction start is the end of the exon before it, a junction end the start of
the exon after.

Written by call_diffsplice_events for alternative 3' and 5' splice sites and
intron retention, and by es_as_inc_excl_to_counts for exon skipping.  DRIMSeq
reads it; flair diffsplice skips an event type whose file holds only the header.
"""
from contextlib import contextmanager
from flair import FlairInputDataError
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType

ID_COLUMNS = ('feature_id', 'coordinate')
ISOFORM_IDS_COLUMN = 'isoform_ids'

def columns(sample_columns):
    return list(ID_COLUMNS) + list(sample_columns) + [ISOFORM_IDS_COLUMN]

class EventQuantReader(TsvReader):
    """Reads an event quant TSV, checking its fixed columns.  Sample columns stay
    strings, since callers decide whether to read them as counts or as fractions."""

    def __init__(self, event_quant_tsv):
        super().__init__(event_quant_tsv, defaultColType=str,
                         typeMap={ISOFORM_IDS_COLUMN: commaListType})
        self._check_columns(event_quant_tsv)

    def _check_columns(self, event_quant_tsv):
        found = tuple(self.columns[:len(ID_COLUMNS)])
        if (found != ID_COLUMNS) or (self.columns[-1] != ISOFORM_IDS_COLUMN):
            raise FlairInputDataError(
                f"{event_quant_tsv}: an event quant TSV must have the columns "
                f"{', '.join(ID_COLUMNS)}, one per sample, then {ISOFORM_IDS_COLUMN}; "
                f"found: {', '.join(self.columns)}")

    @property
    def sample_columns(self):
        "the sample column names, in column order"
        return self.columns[len(ID_COLUMNS):-1]

class EventQuantWriter(TsvWriter):
    def __init__(self, event_quant_tsv, sample_columns, *, outFh=None):
        super().__init__(event_quant_tsv, columns=columns(sample_columns),
                         defaultColType=str, outFh=outFh,
                         typeMap={ISOFORM_IDS_COLUMN: commaListType})
        self.sample_columns = list(sample_columns)

    def writeSide(self, side, feature, event, counts, isoform_ids):
        """One side of an event; side is inclusion or exclusion.  feature names the
        side's own locus, which is the event itself for exon skipping and intron
        retention but the shared anchor for an alternative splice site, where
        several events share one anchor."""
        self.writeRow([f'{side}_{feature}', event] + list(counts) + [isoform_ids])

@contextmanager
def event_quant_writer(event_quant_tsv, sample_columns):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(event_quant_tsv)
    with fileOps.AtomicFileCreate(event_quant_tsv) as tmp_tsv:
        with EventQuantWriter(tmp_tsv, sample_columns) as writer:
            yield writer
