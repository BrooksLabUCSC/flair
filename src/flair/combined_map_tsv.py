"""What each combined id was made from, written by flair combine.

Columns, in order:

  new_id       the FLT or FLG id the combined transcriptome uses
  source_ids   the ids it was built from, each written <samples>:<id> with the
               samples comma separated; '; ' separated in the file, since those
               values contain commas, and a list on the row

flair combine numbers afresh, so these files are the only link between an id in
a combined transcriptome and the per-sample ids behind it; see ids.rst.  One
file for isoforms, one for genes.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import semicolonSpaceListType

COLUMNS = ('new_id', 'source_ids')
TYPE_MAP = {'new_id': str, 'source_ids': semicolonSpaceListType}

class CombinedMapReader(TsvReader):
    def __init__(self, combined_map_txt):
        super().__init__(combined_map_txt, typeMap=TYPE_MAP)

class CombinedMapWriter(TsvWriter):
    def __init__(self, combined_map_txt, *, outFh=None):
        super().__init__(combined_map_txt, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)

    def writeSources(self, new_id, source_ids):
        self.writeRow((new_id, source_ids))
