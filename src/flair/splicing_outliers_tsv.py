"""Samples whose PSI of an alternative of a splicing event is far from the median of
all samples', from spliceevents_outliers.  One row per outlying sample of an
alternative.

Columns, in order:

  event                  the event, as in the event files (splicing_events_tsv)
  alternative            the alternative of the event
  event_type             the kind of event
  gene                   the gene the event was found in
  gene_name              its name
  additional_genes       other genes sharing the event's junctions
  additional_gene_names  their names
  sample                 the sample that stands out
  medianPSI              the median PSI across samples
  dev(IQR/2)             half the interquartile range, the spread PSI is judged
                         against
  sample_val             this sample's PSI
  tot_not_NA_samples     how many samples had enough reads to take a PSI from
  event_reads;total_locus_reads  this sample's reads of the alternative and the
                         total its PSI divides by, ';' separated in the file and a
                         pair on the row
  delta_PSI_to_med       sample_val minus medianPSI, as an absolute value
  dev_from_med           that difference in units of dev(IQR/2), or 100000 when the
                         other samples agree exactly

The four gene columns hold several values, comma separated in the file and a list
on the row.  The numbers are written as computed, rounded to six places, and read
back as strings.

Written unfiltered as .spliceevents.outliers.tsv and, after the deviation cutoff and
leaving out outliers a simpler event explains, as
.spliceevents.outliers.filtered.tsv.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType, semicolonIntListType
from flair.splicing_events_tsv import GENE_COLUMNS

READ_COUNTS_COLUMN = 'event_reads;total_locus_reads'

COLUMNS = ('event', 'alternative', 'event_type', *GENE_COLUMNS, 'sample', 'medianPSI', 'dev(IQR/2)',
           'sample_val', 'tot_not_NA_samples', READ_COUNTS_COLUMN, 'delta_PSI_to_med', 'dev_from_med')
TYPE_MAP = {READ_COUNTS_COLUMN: semicolonIntListType, **{col: commaListType for col in GENE_COLUMNS}}

class SplicingOutliersReader(TsvReader):
    def __init__(self, splicing_outliers_tsv):
        super().__init__(splicing_outliers_tsv, defaultColType=str, typeMap=TYPE_MAP)

class SplicingOutliersWriter(TsvWriter):
    def __init__(self, splicing_outliers_tsv, *, outFh=None):
        super().__init__(splicing_outliers_tsv, columns=COLUMNS, defaultColType=str, typeMap=TYPE_MAP,
                         outFh=outFh)
