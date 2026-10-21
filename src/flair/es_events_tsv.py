"""Skipped exon events called from an isoform BED, one row per exon.

Columns, in order:

  exon           the skipped exon, chrom:acceptor-donor
  strand         the strand of its gene
  num_inclusion  how many isoforms include the exon
  num_exclusion  how many skip it
  inclusion_isos ids of the including isoforms, comma separated in the file and
                 a list on the row
  exclusion_isos ids of the skipping isoforms, likewise

Coordinates are zero based and half open, as they come from a BED: a
junction start is the end of the exon before it, a junction end the start of
the exon after.

Written by es_as and read by es_as_inc_excl_to_counts, which attaches counts and
writes an event quant TSV; see event_quant_tsv.  A row with no excluding isoform
is not an event and is dropped there.

"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType

COLUMNS = ('exon', 'strand', 'num_inclusion', 'num_exclusion',
           'inclusion_isos', 'exclusion_isos')

TYPE_MAP = {'exon': str, 'strand': str, 'num_inclusion': int, 'num_exclusion': int,
            'inclusion_isos': commaListType, 'exclusion_isos': commaListType}

class EsEventsReader(TsvReader):
    def __init__(self, es_events_tsv):
        super().__init__(es_events_tsv, typeMap=TYPE_MAP)

class EsEventsWriter(TsvWriter):
    def __init__(self, es_events_tsv, *, outFh=None):
        super().__init__(es_events_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)
