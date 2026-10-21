"""Samples whose PSI for an event is far from the median, from flair spliceevents
--check_outliers.  One row per outlying sample of an event.

Columns, in order:

  eventname                  the event
  eventtype                  the kind of event
  gene                       the gene it was called in
  sample                     the sample that stands out
  medianPSI                  the median PSI across samples
  dev(IQR/2)                 half the interquartile range, the spread PSI is
                             judged against
  sample_val                 this sample's PSI
  tot_not_NA_samples         how many samples had enough reads to take a PSI from
  event_reads;total_locus_reads  this sample's support, the event's reads and
                             the locus total, ';' separated in the file and a
                             pair on the row
  delta_PSI_to_med           sample_val minus medianPSI
  dev_from_med               that difference in units of dev(IQR/2)

Written unfiltered as .diffsplice.outliers.tsv and, after the deviation cutoff,
as .diffsplice.outliers.filtered.tsv.  As with the event files, each partition
writes a headerless part and the header is written once when they are joined.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import semicolonIntListType

READ_COUNTS_COLUMN = 'event_reads;total_locus_reads'

COLUMNS = ('eventname', 'eventtype', 'gene', 'sample', 'medianPSI', 'dev(IQR/2)',
           'sample_val', 'tot_not_NA_samples', READ_COUNTS_COLUMN,
           'delta_PSI_to_med', 'dev_from_med')

class SplicingOutliersReader(TsvReader):
    def __init__(self, splicing_outliers_tsv):
        super().__init__(splicing_outliers_tsv, defaultColType=str,
                         typeMap={READ_COUNTS_COLUMN: semicolonIntListType})

class SplicingOutliersWriter(TsvWriter):
    def __init__(self, splicing_outliers_tsv, *, outFh=None, writeHeader=True):
        super().__init__(splicing_outliers_tsv, columns=COLUMNS, defaultColType=str,
                         typeMap={READ_COUNTS_COLUMN: semicolonIntListType},
                         outFh=outFh, writeHeader=writeHeader)
