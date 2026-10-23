"""Where each read's ends fall on the transcript it was assigned to, written by
count_sam_transcripts --output_endpos.

Columns, in order:

  read             the read name
  transcript       the transcript it was assigned to
  start_sj_index   index of the first splice junction the read covers, empty for
                   a single exon transcript, which has none
  start_sj_dist    distance from the read start to that junction
  start_tend_dist  distance from the read start to the transcript start
  end_sj_index     index of the last splice junction the read covers
  end_sj_dist      distance from the read end to that junction
  end_tend_dist    distance from the read end to the transcript end

One row per read, so a transcript appears once per supporting read.  An empty
start_sj_index is how a single exon transcript shows: there is no junction to
measure from, so only the transcript end distances are meaningful.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter, intOrNoneType

COLUMNS = ('read', 'transcript', 'start_sj_index', 'start_sj_dist', 'start_tend_dist',
           'end_sj_index', 'end_sj_dist', 'end_tend_dist')

TYPE_MAP = {'read': str, 'transcript': str,
            'start_sj_index': intOrNoneType, 'start_sj_dist': intOrNoneType,
            'start_tend_dist': intOrNoneType, 'end_sj_index': intOrNoneType,
            'end_sj_dist': intOrNoneType, 'end_tend_dist': intOrNoneType}

class ReadEndsReader(TsvReader):
    def __init__(self, read_ends_tsv):
        super().__init__(read_ends_tsv, typeMap=TYPE_MAP)

class ReadEndsWriter(TsvWriter):
    def __init__(self, read_ends_tsv, *, outFh=None, writeHeader=True):
        """writeHeader=False appends to a file that already has its header, as the
        reads assigned without the realignment are appended to the ends of a region"""
        super().__init__(read_ends_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh, writeHeader=writeHeader)
