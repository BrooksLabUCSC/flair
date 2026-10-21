"""Reads supporting each isoform of a flair transcriptome, one row per isoform.

Columns, in order:

  isoform_id          the isoform
  full_length_reads   reads that run the whole isoform, the ones its support is
                      judged on
  total_reads         every read assigned to it, full length or not

Written per partition as <prefix>.isoform.counts.tsv and joined into one file,
whose isoform ids are then renamed to their FLT ids.

Not the same file as the per-sample isoform.counts.tsv that count_sam_transcripts
writes, which is one count per isoform; see transcript_counts_tsv.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('isoform_id', 'full_length_reads', 'total_reads')
TYPE_MAP = {'isoform_id': str, 'full_length_reads': int, 'total_reads': int}

class IsoformCountsReader(TsvReader):
    def __init__(self, isoform_counts_tsv):
        super().__init__(isoform_counts_tsv, typeMap=TYPE_MAP)

class IsoformCountsWriter(TsvWriter):
    def __init__(self, isoform_counts_tsv, *, outFh=None):
        super().__init__(isoform_counts_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)
