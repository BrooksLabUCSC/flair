"""Counts per million, one row per isoform.

Columns, in order:

  gene_id      the gene the isoform belongs to
  isoform_id   the isoform
  <sample>     one column per sample, named by its sample id, containing CPM

The same columns as a counts matrix, with the counts replaced by CPM, so it is read
and written through counts_matrix_tsv with float values.  flair quantify writes it
beside the counts matrix.
"""
from flair.counts_matrix_tsv import (CountsMatrixReader, CountsRow,  # noqa: F401
                                     counts_matrix_writer, read_counts_rows,
                                     read_sample_columns)

class CpmReader(CountsMatrixReader):
    "the counts matrix reader with float values"
    def __init__(self, cpm_tsv):
        super().__init__(cpm_tsv, countType=float)

def cpm_writer(cpm_tsv, sample_columns):
    "write a CPM matrix atomically; values are floats rather than counts"
    return counts_matrix_writer(cpm_tsv, sample_columns, countType=float)

def read_cpm_rows(cpm_tsv):
    "every row of the matrix, in file order"
    with CpmReader(cpm_tsv) as reader:
        return list(reader)
