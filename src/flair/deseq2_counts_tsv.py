"""The counts matrix handed to DESeq2, one row per feature.

Columns, in order:

  <empty>      the feature id, a gene id in the gene table and an isoform id in
               the isoform table.  The column name is empty because R reads the
               file with row.names = 1 and takes this column as the row names
  <sample>     one integer column per sample, named as in the formula matrix

Written by diffexp into its working directory after filtering by expression; not
an output anyone reads directly.
"""
from contextlib import contextmanager
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvReader, TsvWriter
from flair.tsv_column_types import countType

# R takes column one as row names, so it has no name
FEATURE_COLUMN = ''

# the unnamed feature column cannot be a row attribute, so the reader names it
FEATURE_ID = 'feature_id'

def _type_map(feature_column):
    """The column types, taken by both the reader and the writer, so the file is
    parsed as it was formatted.  The feature column has one name in the file and
    another on the row, so the caller names the one it uses."""
    return {feature_column: str}

class Deseq2CountsReader(TsvReader):
    def __init__(self, deseq2_counts_tsv):
        super().__init__(deseq2_counts_tsv, defaultColType=countType,
                         typeMap=_type_map(FEATURE_ID),
                         columnNameMapper=lambda col: FEATURE_ID if col == FEATURE_COLUMN else col)

    @property
    def sample_columns(self):
        return self.columns[1:]

class Deseq2CountsWriter(TsvWriter):
    def __init__(self, deseq2_counts_tsv, sample_columns):
        super().__init__(deseq2_counts_tsv,
                         columns=[FEATURE_COLUMN] + list(sample_columns),
                         typeMap=_type_map(FEATURE_COLUMN),
                         defaultColType=countType)

    def writeFeature(self, feature_id, counts):
        # DESeq2 takes integer counts, and diffexp sums them as floats; countType
        # writes 3.0 as 3, which Deseq2CountsReader reads back as an integer
        self.writeRow([feature_id] + list(counts))

@contextmanager
def deseq2_counts_writer(deseq2_counts_tsv, sample_columns):
    "write atomically; the path appears only once complete"
    fileOps.ensureFileDir(deseq2_counts_tsv)
    with fileOps.AtomicFileCreate(deseq2_counts_tsv) as tmp_tsv:
        with Deseq2CountsWriter(tmp_tsv, sample_columns) as writer:
            yield writer

def write_deseq2_counts(deseq2_counts_tsv, sample_columns, feature_ids, feature_counts):
    with deseq2_counts_writer(deseq2_counts_tsv, sample_columns) as writer:
        for feature_id, counts in zip(feature_ids, feature_counts):
            writer.writeFeature(feature_id, counts)
