"""The flair quantify counts matrix: the read count of each isoform in each sample.

Columns, in order:

  gene_id      the gene the isoform belongs to
  isoform_id   the isoform
  <sample>     one integer column per sample, named by its sample id

Each row is one isoform of one gene.  Sample column names come from the manifest
and are arbitrary, so they are addressed as a list, CountsRow.counts, in column
order, rather than by name.  flair quantify writes a sample info TSV beside this
file naming the condition and batch of each column; see sample_info_tsv.

flair diffexp and diff_iso_usage also read the counts of flair spliceevents and
flair variantquant as a counts matrix (read_grouped_counts): a group of rows
compared against each other, as an isoform's share of its gene, and a feature
naming each row.  See COUNTS_FORMATS.

A matrix from an earlier FLAIR joined the gene and the isoform into one
isoform_gene column.  That form is refused rather than split, since splitting it
needs a heuristic that is wrong for a gene id containing an underscore.

Reading or writing a path ending .gz or .bz2 compresses it, here as in every other
flair TSV.
"""
from contextlib import contextmanager
from collections import namedtuple
from flair import FlairInputDataError
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvError, TsvReader, TsvRow, TsvWriter
from flair.tsv_column_types import read_header_columns
from flair import splicing_events_tsv, var_counts_tsv

# the id columns every counts matrix starts with, in this order
ID_COLUMNS = ('gene_id', 'isoform_id')

# ids named when an error reports a set of offending ids, the rest being redundant
MAX_REPORTED_IDS = 10

class CountsRow(namedtuple('CountsRow', ('gene_id', 'isoform_id', 'counts'))):
    """One row: one isoform of one gene, and its count in each sample, in column
    order.  Rows are built in memory by quantify, combine and diffexp as well as
    read from a file, so this is a plain record rather than a TsvRow bound to a
    reader."""
    __slots__ = ()

class _CountsMatrixTsvRow(TsvRow):
    "a row as the reader parses it, the sample counts addressed by position"
    @property
    def counts(self):
        return [self[i] for i in range(len(ID_COLUMNS), len(self))]

class CountsMatrixReader(TsvReader):
    """Reads a counts matrix, checking the id columns and converting every sample
    column with countType, int for counts and float for a TPM matrix.  Iterating
    yields CountsRow."""

    def __init__(self, counts_matrix_tsv, *, countType=int):
        super().__init__(counts_matrix_tsv, rowClass=_CountsMatrixTsvRow,
                         typeMap={col: str for col in ID_COLUMNS},
                         defaultColType=countType)
        self._check_id_columns(counts_matrix_tsv)

    def _check_id_columns(self, counts_matrix_tsv):
        found = tuple(self.columns[:len(ID_COLUMNS)])
        if found != ID_COLUMNS:
            raise FlairInputDataError(
                f"{counts_matrix_tsv}: a counts matrix must start with the columns "
                f"{', '.join(ID_COLUMNS)}, found: {', '.join(found)}; a matrix from an "
                "earlier FLAIR, with one isoform_gene column, must be remade with "
                "flair quantify")

    @property
    def sample_columns(self):
        "the sample column names, in column order"
        return self.columns[len(ID_COLUMNS):]

    def __iter__(self):
        for row in super().__iter__():
            yield CountsRow(row.gene_id, row.isoform_id, row.counts)

class CountsMatrixWriter(TsvWriter):
    "Writes a counts matrix for the given samples, in column order."

    def __init__(self, counts_matrix_tsv, sample_columns, *, countType=int):
        super().__init__(counts_matrix_tsv,
                         columns=list(ID_COLUMNS) + list(sample_columns),
                         typeMap={col: str for col in ID_COLUMNS},
                         defaultColType=countType)
        self.sample_columns = list(sample_columns)

    def writeCountsRow(self, row):
        if len(row.counts) != len(self.sample_columns):
            raise FlairInputDataError(
                f"isoform {row.isoform_id} has {len(row.counts)} counts, "
                f"{len(self.sample_columns)} samples are named by the columns")
        self.writeRow([row.gene_id, row.isoform_id] + list(row.counts))

@contextmanager
def counts_matrix_writer(counts_matrix_tsv, sample_columns, *, countType=int):
    "write a counts matrix atomically; the path appears only once complete"
    fileOps.ensureFileDir(counts_matrix_tsv)
    with fileOps.AtomicFileCreate(counts_matrix_tsv) as tmp_tsv:
        with CountsMatrixWriter(tmp_tsv, sample_columns, countType=countType) as writer:
            yield writer

def read_counts_rows(counts_matrix_tsv, *, countType=int):
    "every row of the matrix, in file order"
    with CountsMatrixReader(counts_matrix_tsv, countType=countType) as reader:
        return list(reader)

def read_sample_columns(counts_matrix_tsv):
    "the sample column names, in column order, of any of COUNTS_FORMATS"
    columns = read_header_columns(counts_matrix_tsv)
    return columns[len(counts_format(columns, counts_matrix_tsv).id_columns):]

def read_isoform_ids(counts_matrix_tsv):
    "the isoform id of each row, in file order"
    return [row.isoform_id for row in read_counts_rows(counts_matrix_tsv)]

def write_counts_matrix(counts_matrix_tsv, sample_columns, rows, *, countType=int):
    "write a whole counts matrix at once"
    with counts_matrix_writer(counts_matrix_tsv, sample_columns, countType=countType) as writer:
        for row in rows:
            writer.writeCountsRow(row)

def describe_ids(ids):
    "the first few ids, for naming the offenders in an error without listing them all"
    named = ', '.join(ids[:MAX_REPORTED_IDS])
    return named + (', ...' if len(ids) > MAX_REPORTED_IDS else '')


class CountsFormat(namedtuple('CountsFormat', ('program', 'group_column', 'feature_column', 'id_columns',
                                               'unreported_feature'))):
    """A table of integer counts that the differential usage tests read: its
    id_columns, then a column per sample.  Rows of one group_column value are
    compared against each other, each named by its feature_column value.  Rows of the
    unreported_feature, if any, count toward their group's total but are not tested
    themselves."""
    __slots__ = ()


QUANTIFY_COUNTS = CountsFormat('flair quantify', 'gene_id', 'isoform_id', ID_COLUMNS, None)
# an event's other row holds its reads of none of its alternatives
SPLICEEVENTS_COUNTS = CountsFormat('flair spliceevents', 'event', 'alternative',
                                   splicing_events_tsv.EVENT_COLUMNS, 'other')
# an allele's two rows, the reads showing it and those not
VARIANTQUANT_COUNTS = CountsFormat('flair variantquant', 'var_name', 'is_var',
                                   var_counts_tsv.ID_COLUMNS, None)
COUNTS_FORMATS = (QUANTIFY_COUNTS, SPLICEEVENTS_COUNTS, VARIANTQUANT_COUNTS)

def counts_format(columns, counts_tsv):
    "the format of a table of these columns"
    for fmt in COUNTS_FORMATS:
        if tuple(columns[:len(fmt.id_columns)]) == fmt.id_columns:
            return fmt
    expected = '; '.join(f"{fmt.program}: {', '.join(fmt.id_columns[:2])}" for fmt in COUNTS_FORMATS)
    raise FlairInputDataError(
        f"{counts_tsv} is not a table of counts this reads, which start with the columns of "
        f"{expected}; it starts with: {', '.join(columns[:2])}.  A matrix from an earlier FLAIR, "
        "with one isoform_gene column, must be remade with flair quantify")

class GroupedCounts(namedtuple('GroupedCounts', ('format', 'sample_columns', 'rows', 'unreported'))):
    """The rows of a table of COUNTS_FORMATS as CountsRow: gene_id holds the row's
    group and isoform_id its feature id, unique over the table.  A quantify isoform
    keeps its id, and a spliceevents or variantquant row is named group|feature, as
    spliceevents names its BED records.  unreported is the set of feature ids of rows
    of the format's unreported_feature."""
    __slots__ = ()

def _feature_id(fmt, group, feature):
    return feature if fmt is QUANTIFY_COUNTS else f'{group}|{feature}'

def read_grouped_counts(counts_tsv):
    """Every row of a counts matrix, a flair spliceevents counts table or a flair
    variantquant counts table, in file order, as GroupedCounts"""
    columns = read_header_columns(counts_tsv)
    fmt = counts_format(columns, counts_tsv)
    sample_columns = columns[len(fmt.id_columns):]
    rows, unreported = [], set()
    try:
        with TsvReader(counts_tsv, defaultColType=str, typeMap={col: int for col in sample_columns}) as reader:
            for row in reader:
                group, feature = row[fmt.group_column], row[fmt.feature_column]
                feature_id = _feature_id(fmt, group, feature)
                rows.append(CountsRow(group, feature_id, [row[col] for col in sample_columns]))
                if feature == fmt.unreported_feature:
                    unreported.add(feature_id)
    except TsvError as ex:
        raise FlairInputDataError(f"{counts_tsv}: its sample columns must hold read counts, integers; "
                                  f"a table of PSI values can't be tested ({ex.__cause__})") from ex
    return GroupedCounts(fmt, sample_columns, rows, unreported)
