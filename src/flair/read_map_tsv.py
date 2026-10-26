"""Which reads support a thing, one row per thing.

Columns, in order:

  name         the isoform, allele group or allele the reads support
  reads        the read names, comma separated in the file and a list on the row

Written by flair transcriptome, quantify, alleles and isoalleles, each with its
own file name.
"""
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType

COLUMNS = ('name', 'reads')
TYPE_MAP = {'name': str, 'reads': commaListType}

class ReadMapReader(TsvReader):
    def __init__(self, read_map_tsv):
        super().__init__(read_map_tsv, typeMap=TYPE_MAP)

class ReadMapWriter(TsvWriter):
    def __init__(self, read_map_tsv, *, outFh=None):
        super().__init__(read_map_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)

    def writeReads(self, name, reads):
        "reads are written in the order given, so callers sort when it matters"
        self.writeRow((name, reads))
