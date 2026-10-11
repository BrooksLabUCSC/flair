"""The sample info TSV that flair quantify writes beside a counts matrix.

Columns, in order:

  sample_id    names one counts matrix column
  condition    the condition that sample was assigned in the manifest
  batch        the batch it was assigned

One row per counts column, in the same order as those columns.  The fields are
columns of their own, so they may contain any character.

A counts matrix written before this file existed names each column
sample_condition_batch, the first three manifest fields joined by an underscore.
For such a matrix the fields are recovered from the column name instead, which is
why the manifest documentation forbids underscores in those fields.
"""
from contextlib import contextmanager
from collections import namedtuple
import os
from flair import FlairInputDataError
from flair.pycbio.sys import fileOps
from flair.counts_matrix_tsv import read_sample_columns
from flair.pycbio.tsv import TsvReader, TsvWriter

SAMPLE_INFO_COLUMNS = ('sample_id', 'condition', 'batch')

# fields of a sample_condition_batch column name, for a matrix with no sample info
CONDITION_FIELD = 1
BATCH_FIELD = -1

class SampleInfo(namedtuple('SampleInfo', SAMPLE_INFO_COLUMNS)):
    "one counts matrix column: which sample it holds and how that sample was grouped"
    __slots__ = ()

class SampleInfoReader(TsvReader):
    "Reads a sample info TSV, checking its columns.  Iterating yields SampleInfo."

    def __init__(self, sample_info_tsv):
        super().__init__(sample_info_tsv, defaultColType=str)
        found = tuple(self.columns)
        if found != SAMPLE_INFO_COLUMNS:
            raise FlairInputDataError(
                f"{sample_info_tsv}: expected columns {', '.join(SAMPLE_INFO_COLUMNS)}, "
                f"found: {', '.join(found)}")

    def __iter__(self):
        for row in super().__iter__():
            yield SampleInfo(row.sample_id, row.condition, row.batch)

class SampleInfoWriter(TsvWriter):
    def __init__(self, sample_info_tsv):
        super().__init__(sample_info_tsv, columns=SAMPLE_INFO_COLUMNS,
                         defaultColType=str)

@contextmanager
def sample_info_writer(sample_info_tsv):
    "write sample info atomically; the path appears only once complete"
    fileOps.ensureFileDir(sample_info_tsv)
    with fileOps.AtomicFileCreate(sample_info_tsv) as tmp_tsv:
        with SampleInfoWriter(tmp_tsv) as writer:
            yield writer


# what opengz compresses by, stripped before deriving a sibling name and put back
# after, so a compressed matrix names a compressed sibling
COMPRESS_SUFFIXES = ('.gz', '.bgz', '.bz2', '.Z')

def split_compress_suffix(path):
    "the path without its compression suffix, and that suffix"
    for suffix in COMPRESS_SUFFIXES:
        if path.endswith(suffix):
            return path[:-len(suffix)], suffix
    return path, ''

def sample_info_path(counts_matrix_tsv):
    """flair quantify writes <prefix>.counts.tsv beside <prefix>.sample_info.tsv.  A
    compressed matrix names a compressed sample info file, since quantify writes the
    pair the same way."""
    base, compress = split_compress_suffix(counts_matrix_tsv)
    for suffix in ('.counts.tsv', '.tsv'):
        if base.endswith(suffix):
            return base[:-len(suffix)] + '.sample_info.tsv' + compress
    return base + '.sample_info.tsv' + compress

def write_sample_info(sample_info_tsv, sample_infos):
    with sample_info_writer(sample_info_tsv) as writer:
        for sample_info in sample_infos:
            writer.writeRow(sample_info)

def parse_sample_fields(sample_columns, counts_matrix_tsv):
    "the condition and the batch of each sample column, from the column names"
    try:
        conditions = [col.split('_')[CONDITION_FIELD] for col in sample_columns]
        batches = [col.split('_')[BATCH_FIELD] for col in sample_columns]
    except IndexError as ex:
        raise FlairInputDataError(
            f"{counts_matrix_tsv}: counts columns must be named sample_condition_batch, "
            f"found: {' '.join(sample_columns)}") from ex
    return conditions, batches

def _sample_info_from_columns(sample_columns, counts_matrix_tsv):
    "recover the fields from sample_condition_batch column names"
    conditions, batches = parse_sample_fields(sample_columns, counts_matrix_tsv)
    return [SampleInfo(col, condition, batch)
            for col, condition, batch in zip(sample_columns, conditions, batches)]

def _check_sample_info(sample_infos, sample_columns, sample_info_tsv, counts_matrix_tsv):
    "the sample info rows must describe the counts columns, in the same order"
    named = [si.sample_id for si in sample_infos]
    if named != sample_columns:
        raise FlairInputDataError(
            f"{sample_info_tsv} does not describe the columns of {counts_matrix_tsv}: "
            f"it names {', '.join(named)}, the counts columns are {', '.join(sample_columns)}")

def find_sample_info(counts_matrix_tsv):
    """The sample info file beside a counts matrix, or None.  The two need not agree
    about compression: either may have been written or fetched on its own."""
    derived = sample_info_path(counts_matrix_tsv)
    base, _ = split_compress_suffix(derived)
    for path in (derived, base) + tuple(base + suffix for suffix in COMPRESS_SUFFIXES):
        if os.path.exists(path):
            return path
    return None

def read_sample_info(counts_matrix_tsv):
    """The sample of each counts column and how it was grouped, taken from the sample
    info file that flair quantify writes, or from the column names when a counts
    matrix predates that file."""
    sample_columns = read_sample_columns(counts_matrix_tsv)
    sample_info_tsv = find_sample_info(counts_matrix_tsv)
    if sample_info_tsv is None:
        return _sample_info_from_columns(sample_columns, counts_matrix_tsv)
    with SampleInfoReader(sample_info_tsv) as reader:
        sample_infos = list(reader)
    _check_sample_info(sample_infos, sample_columns, sample_info_tsv, counts_matrix_tsv)
    return sample_infos
