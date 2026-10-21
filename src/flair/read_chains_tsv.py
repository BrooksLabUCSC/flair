"""A sample's corrected reads in one partition, written by flair spliceevents' read
correction for its event calling: one row per junction chain.

Columns, in order:

  genes      the chain's annotated genes, 5' to 3', comma separated in the file and a
             list on the row, empty for none
  chrom      the chromosome
  juncs      its junctions, start-end, comma separated in the file and a tuple of
             (start, end) pairs on the row, empty for unspliced reads
  strand     its strand
  ends       the ends of each of its reads, start-end, as juncs
  fragment   for unspliced reads that are fragments of one gene, 'first' for a 5'
             fragment or 'last' for a 3' one, otherwise empty

Coordinates are zero based and half open, as in a BED.

An intermediate, written whole or not at all, since a rerun with --keep_intermediate
skips the samples whose file exists.  Files of an earlier spliceevents have no
header and no fragment column; spliceevents_reads.read_read_chains reads those too.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType, rangeListType

COLUMNS = ('genes', 'chrom', 'juncs', 'strand', 'ends', 'fragment')
TYPE_MAP = {'genes': commaListType, 'juncs': rangeListType, 'ends': rangeListType}

class ReadChainsReader(TsvReader):
    def __init__(self, read_chains_tsv):
        super().__init__(read_chains_tsv, defaultColType=str, typeMap=TYPE_MAP)

class ReadChainsWriter(TsvWriter):
    def __init__(self, read_chains_tsv, *, outFh=None):
        super().__init__(read_chains_tsv, columns=COLUMNS, defaultColType=str, typeMap=TYPE_MAP,
                         outFh=outFh)

@contextmanager
def read_chains_writer(read_chains_tsv):
    "write atomically; the path appears only once complete"
    with fileOps.AtomicFileCreate(read_chains_tsv) as tmp_tsv:
        with ReadChainsWriter(tmp_tsv) as writer:
            yield writer
