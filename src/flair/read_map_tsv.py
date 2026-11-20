"""Which reads support a thing, one row per thing.

Columns, in order:

  isoform_id   the isoform the reads support.  flair isoalleles writes an isoform
               and one of its allele groups, isoform|group, as it does in its
               counts matrix, its BED and its FASTA
  reads        the read names, comma separated in the file and a list on the row

Written by flair transcriptome, quantify and isoalleles, each with its own file name.
"""
import numpy as np
from flair import FlairInputDataError
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import commaListType, read_header_columns

COLUMNS = ('isoform_id', 'reads')
TYPE_MAP = {'isoform_id': str, 'reads': commaListType}

# what the id column was called before FLAIR 3.1, when it was renamed to the name
# every other flair table gives a column of isoform ids
FORMER_ID_COLUMN = 'name'


def check_read_map_columns(read_map_tsv, columns):
    "a read map of a FLAIR before 3.1, whose first column was named name"
    if tuple(columns) == (FORMER_ID_COLUMN, COLUMNS[1]):
        raise FlairInputDataError(
            f'{read_map_tsv}: its columns are {FORMER_ID_COLUMN} and {COLUMNS[1]}, the read map of a '
            f'FLAIR before 3.1, which named the isoform column {FORMER_ID_COLUMN}; make it again with '
            'flair quantify --generate_map, or rename that column to ' + COLUMNS[0])


# repeated ids named in an error, the rest being redundant
MAX_REPORTED_IDS = 3

class ReadMapReader(TsvReader):
    def __init__(self, read_map_tsv):
        check_read_map_columns(read_map_tsv, read_header_columns(read_map_tsv))
        super().__init__(read_map_tsv, typeMap=TYPE_MAP)

class ReadMapWriter(TsvWriter):
    def __init__(self, read_map_tsv, *, outFh=None):
        super().__init__(read_map_tsv, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)

    def writeReads(self, isoform_id, reads):
        "reads are written in the order given, so callers sort when it matters"
        self.writeRow((isoform_id, reads))

def _read_id_hashes(read_map_tsv):
    """The hash of every read id of the map.  Hashes, not the ids: a sample has an id
    per read it counted, and keeping tens of millions of names costs gigabytes."""
    with ReadMapReader(read_map_tsv) as reader:
        return np.fromiter((hash(read) for row in reader for read in row.reads),
                           dtype=np.int64)

def _repeated_hashes(hashes):
    "the hashes that appear more than once, which are candidate repeated ids"
    hashes.sort()
    return frozenset(hashes[:-1][hashes[:-1] == hashes[1:]].tolist())

def _ids_of_hashes(read_map_tsv, hashes):
    "how often each id of the given hashes is used, two hashes being able to collide"
    used = {}
    with ReadMapReader(read_map_tsv) as reader:
        for read in (read for row in reader for read in row.reads if hash(read) in hashes):
            used[read] = used.get(read, 0) + 1
    return used

def check_read_ids_unique(read_map_tsv):
    """Raise if an id names more than one read.  A counted read is written for one
    isoform, so an id in two rows, or twice in one, is two reads sharing a name: the
    reads of a BAM merged from runs that numbered their reads the same way."""
    repeated = _repeated_hashes(_read_id_hashes(read_map_tsv))
    if len(repeated) == 0:
        return
    duplicated = sorted(read for read, used in _ids_of_hashes(read_map_tsv, repeated).items()
                        if used > 1)
    if len(duplicated) > 0:
        named = ', '.join(duplicated[:MAX_REPORTED_IDS])
        raise FlairInputDataError(
            f"{len(duplicated)} read ids name more than one read of {read_map_tsv}: {named}. "
            "An id must name one read of a sample's BAM: label each read with its sample or "
            "run before merging the runs, then align and quantify again")
