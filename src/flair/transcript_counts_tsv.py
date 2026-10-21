"""Reads assigned to each transcript by count_sam_transcripts, one row per
transcript.

Columns, in order:

  transcript   the transcript
  count        how many reads were assigned to it

Written per sample into the working directory and read back by quantify, which
turns the per-sample files into a counts matrix; see counts_matrix_tsv.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('transcript', 'count')
TYPE_MAP = {'transcript': str, 'count': int}

class TranscriptCountsReader(TsvReader):
    def __init__(self, transcript_counts_tsv):
        super().__init__(transcript_counts_tsv, typeMap=TYPE_MAP)

class TranscriptCountsWriter(TsvWriter):
    def __init__(self, transcript_counts_tsv, *, outFh=None):
        super().__init__(transcript_counts_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)
